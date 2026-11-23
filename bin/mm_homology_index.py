#!/usr/bin/env python3
"""Construct coordinate-correspondence blocks from reference self-alignment PAF.

Require a cg:Z: CIGAR; indels split homology blocks. This is proposal evidence
only, not an alignment of a sequenced fragment to a proposed genomic location.
"""
import argparse
import csv
import gzip
import re

CIGAR=re.compile(r'(\d+)([MIDNSHP=X])')

def blocks(x,min_block):
    qs,qe=int(x[2]),int(x[3]);ts,te=int(x[7]),int(x[8]);strand=x[4]
    if strand not in ('+','-'):raise ValueError('Invalid PAF strand')
    cg=[v[5:] for v in x[12:] if v.startswith('cg:Z:')]
    if len(cg)!=1:raise ValueError('PAF must contain exactly one cg:Z: CIGAR')
    cigar=cg[0];tokens=list(CIGAR.finditer(cigar))
    if ''.join(t.group(0) for t in tokens)!=cigar:raise ValueError('Invalid PAF CIGAR')
    q=qs if strand=='+' else qe
    t=ts
    for token in tokens:
        n=int(token.group(1));op=token.group(2)
        if op in 'M=X':
            qlo=q if strand=='+' else q-n
            if n>=min_block:yield [x[0],qlo,qlo+n,x[5],t,t+n,strand]
            q+=(n if strand=='+' else -n);t+=n
        elif op=='I':q+=(n if strand=='+' else -n)
        elif op in 'DN':t+=n
        elif op in 'S':q+=(n if strand=='+' else -n)
        elif op in 'HP':pass
        else:raise ValueError('Unsupported CIGAR operation '+op)
    if t!=te or q!=(qe if strand=='+' else qs):
        raise ValueError('PAF coordinates inconsistent with cg CIGAR')

def build(paf,output,min_block=30,min_identity=.9):
    opener=gzip.open if paf.endswith('.gz') else open
    count=0
    with opener(paf,'rt') as inp,open(output,'w') as out:
        w=csv.writer(out,delimiter='\t',lineterminator='\n')
        w.writerow(['source_chrom','source_start0','source_end0','target_chrom','target_start0','target_end0','strand'])
        for line in inp:
            if not line.strip():continue
            x=line.rstrip('\n').split('\t')
            if len(x)<12:raise ValueError('Malformed PAF')
            if x[0]==x[5] and x[2]==x[7] and x[3]==x[8] and x[4]=='+':continue
            if int(x[10])==0 or int(x[9])/int(x[10])<min_identity:continue
            for b in blocks(x,min_block):
                w.writerow(b)
                w.writerow([b[3],b[4],b[5],b[0],b[1],b[2],b[6]])
                count+=2
    return count

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--paf',required=True);p.add_argument('--output',required=True)
    p.add_argument('--min-block',type=int,default=30);p.add_argument('--min-identity',type=float,default=.9)
    a=p.parse_args()
    if not 0<a.min_identity<=1 or a.min_block<1:p.error('Invalid thresholds')
    print('directed_blocks=',build(a.paf,a.output,a.min_block,a.min_identity))
