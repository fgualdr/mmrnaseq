#!/usr/bin/env python3
"""Experimental source-level EM on COMPLETE, quality-scored STAR candidates.

This is a conditional mixture model, not a calibrated posterior over all genomic
origins. Candidate expansion is NOT accepted until independently realigned.
Loci are nonoverlapping BED intervals and are used as a proof-of-concept source
universe. Fragments with unassigned/ambiguous candidates are omitted wholesale.
"""
import argparse
import csv
import gzip
import math
import sys
from bisect import bisect_right
from itertools import zip_longest
from collections import defaultdict

import numpy as np


def opener(path):return gzip.open if path.endswith('.gz') else open

class Loci:
    def __init__(self,bed):
        self.data=defaultdict(list)
        self.names=set()
        with open(bed) as f:
            for no,line in enumerate(f,1):
                if not line.strip() or line.startswith('#'):continue
                p=line.rstrip('\n').split('\t')
                if len(p)<4:raise ValueError('BED4 required')
                start,end=int(p[1]),int(p[2]);name=p[3]
                if start<0 or end<=start or name in self.names:raise ValueError(f'Invalid or repeated BED locus: {no}')
                self.data[p[0]].append((start,end,name));self.names.add(name)
        self.starts={}
        for chrom,items in self.data.items():
            items.sort()
            for a,b in zip(items,items[1:]):
                if a[1]>b[0]:raise ValueError(f'BED loci overlap on {chrom}: {a[2]}, {b[2]}')
            self.starts[chrom]=[r[0] for r in items]
    def lookup(self,chrom,pos):
        items=self.data.get(chrom,[])
        if not items:return None
        i=bisect_right(self.starts[chrom],pos)-1
        return items[i][2] if i>=0 and pos<items[i][1] else None
    def source(self,r):
        src=self.lookup(r['r1_chrom'],int(r['r1_pos1'])-1)
        if src is None:return None
        if int(r['paired']):
            if self.lookup(r['r2_chrom'],int(r['r2_pos1'])-1)!=src:return None
        return src


def read_rows(candidates,scores,loci,max_edges):
    with opener(candidates)(candidates,'rt') as f,opener(scores)(scores,'rt') as g:
        rd=csv.DictReader(f,delimiter='\t');sc=csv.DictReader(g,delimiter='\t')
        if not {'fragment_id','candidate_id','log_quality_score'}<=set(sc.fieldnames or []):
            raise ValueError('Quality scores missing required columns')
        current=None;rows=[];count=0
        def group():
            nonlocal current,rows,count
            if not rows:return None
            nh={int(r[0]['nh']) for r in rows}
            his=[int(r[0]['hi']) for r in rows]
            # Fail closed: a fragment must contain all STAR-declared candidates.
            if len(nh)!=1 or len(rows)!=next(iter(nh)) or set(his)!=set(range(1,len(rows)+1)):
                return (current,None,'INCOMPLETE_STAR_SET')
            if any(source is None for _,source,_ in rows):return (current,None,'UNMODELED_SOURCE')
            return (current,rows,'COMPLETE')
        for r,s in zip_longest(rd,sc):
            if r is None or s is None: raise ValueError('Candidate and score tables have unequal lengths')
            if r['fragment_id']!=s['fragment_id'] or r['candidate_id']!=s['candidate_id']:
                raise ValueError('Candidates and quality scores not identically ordered/complete')
            fragment=r['fragment_id']
            if fragment!=current:
                if current is not None:yield group()
                current=fragment;rows=[]
            rows.append((r,loci.source(r),float(s['log_quality_score'])))
            count+=1
            if count>max_edges:raise MemoryError(f'Edge limit {max_edges} exceeded: run a smaller proof-of-concept region')
        if current is not None:yield group()


def unique_source_counts(path,loci):
    """One count per fragment; never one count per mate or exon block."""
    from itertools import groupby
    counts=defaultdict(float)
    if not path:return counts
    with opener(path)(path,'rt') as f:
        rd=csv.DictReader(f,delimiter='\t')
        for frag,rows in groupby(rd,key=lambda r:r['fragment_id']):
            labels={loci.lookup(r['chrom'],int(r['pos0'])) for r in rows}
            if len(labels)==1 and None not in labels:counts[next(iter(labels))]+=1
    return counts


def em(groups,unique_counts,alpha=1.01,max_iter=100,tol=1e-8):
    # Sparse arrays: each observed alignment is one edge; multiple alignments to
    # the same SOURCE get per-source normalisation, not a multiplicity bonus.
    fragments=[];candidates=[];sources=[];logs=[];frag_ix=[];rejected=defaultdict(int)
    for fragment,rows,status in groups:
        if rows is None:rejected[status]+=1;continue
        g=len(fragments);fragments.append(fragment)
        for r,source,score in rows:
            frag_ix.append(g);sources.append(source);logs.append(score);candidates.append(r['candidate_id'])
    nfrag=len(fragments)
    if not nfrag:raise ValueError('No complete modeled multimapping fragments')
    names=sorted(set(sources)|set(unique_counts))
    ix={n:i for i,n in enumerate(names)}
    f=np.asarray(frag_ix,dtype=np.int64)
    s=np.fromiter((ix[x] for x in sources),dtype=np.int64)
    logs=np.asarray(logs,dtype=np.float64)
    unique=np.asarray([unique_counts.get(n,0.) for n in names],dtype=np.float64)
    # Stable likelihood scale per fragment. n(f,source) removes artificial gain
    # from multiple genomic alignment starts within the same biological source.
    max_logs=np.maximum.reduceat(logs,np.r_[0,np.flatnonzero(np.diff(f))+1])
    pair=f*len(names)+s
    _,inverse,multiplicity=np.unique(pair,return_inverse=True,return_counts=True)
    q=np.exp(logs-max_logs[f])/multiplicity[inverse]
    theta=np.ones(len(names),dtype=np.float64)/len(names)
    converged=False
    for iteration in range(1,max_iter+1):
        weights=q*theta[s]
        denom=np.bincount(f,weights=weights,minlength=nfrag)
        if np.any(denom<=0):raise ValueError('Numerically degenerate fragment likelihood')
        resp=weights/denom[f]
        support=np.bincount(s,weights=resp,minlength=len(names))
        new=unique+support+(alpha-1)
        new/=new.sum()
        delta=float(np.max(np.abs(new-theta)))
        theta=new
        if delta<tol:converged=True;break
    # Final E-step with final theta (avoid one-iteration mismatch).
    weights=q*theta[s];denom=np.bincount(f,weights=weights,minlength=nfrag)
    resp=weights/denom[f]
    if not np.allclose(np.bincount(f,weights=resp,minlength=nfrag),1,atol=1e-8):
        raise AssertionError('Fragment posterior does not sum to one')
    return names,theta,unique,fragments,candidates,s,resp,dict(rejected),iteration,converged


def run(args):
    loci=Loci(args.loci)
    uniques=unique_source_counts(args.unique_anchors,loci)
    result=em(read_rows(args.candidates,args.scores,loci,args.max_edges),uniques,args.alpha,args.max_iter,args.tol)
    names,theta,unique,fragments,candidates,s,resp,rejected,iters,converged=result
    with open(args.output_prefix+'.em_locus_abundance.tsv','w') as out:
        w=csv.writer(out,delimiter='\t',lineterminator='\n')
        w.writerow(['locus','theta','trusted_unique_fragments','expected_multimapping_fragments'])
        mass=np.bincount(s,weights=resp,minlength=len(names))
        for n,t,u,m in zip(names,theta,unique,mass):w.writerow([n,f'{t:.12g}',f'{u:.12g}',f'{m:.12g}'])
    with open(args.output_prefix+'.em_candidate_posterior.tsv','w') as out:
        w=csv.writer(out,delimiter='\t',lineterminator='\n')
        w.writerow(['fragment_id','candidate_id','locus','posterior_conditional'])
        for cand,loc,p in zip(candidates,s,resp):
            w.writerow([cand.split('|HI=')[0],cand,names[loc],f'{p:.12g}'])
    with open(args.output_prefix+'.em_qc.tsv','w') as out:
        w=csv.writer(out,delimiter='\t',lineterminator='\n')
        w.writerow(['metric','value'])
        for k,v in {'mode':'experimental_conditional','fragments':len(fragments),'candidate_edges':len(candidates),
                    'unique_anchors_counted_as_fragments':sum(uniques.values()),'iterations':iters,'converged':int(converged),
                    **{'rejected_'+k:v for k,v in rejected.items()}}.items():w.writerow([k,v])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['candidates','scores','loci','output-prefix']:p.add_argument('--'+name,required=True)
    p.add_argument('--unique-anchors',default=None,help='ABlation ONLY: raw unique anchors as fixed fragment counts (mappability-biased); disabled by default')
    p.add_argument('--alpha',type=float,default=1.01)
    p.add_argument('--max-iter',type=int,default=100)
    p.add_argument('--tol',type=float,default=1e-8)
    p.add_argument('--max-edges',type=int,default=2_000_000)
    a=p.parse_args()
    if a.alpha<1 or a.max_edges<1:p.error('alpha must be >=1, max-edges >=1')
    run(a)
