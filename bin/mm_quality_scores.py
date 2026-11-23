#!/usr/bin/env python3
"""Quality-aware alignment *pseudo*-log-likelihoods for existing STAR candidates.

Reads SAM SEQ/QUAL (which preserve trimmed FASTQ qualities) and an uncompressed,
faidx-indexed reference. This is NOT a calibrated pair-HMM; indels and clipping
need separate calibration. No graph-proposed locations enter this routine.
"""
import argparse
import csv
import gzip
import math
import re
import sys
from itertools import groupby

OPS = re.compile(r'(\d+)([MIDNSHP=X])')

class Faidx:
    def __init__(self, fasta, index=None):
        self.file=open(fasta,'rb')
        self.idx={}
        with open(index or fasta+'.fai') as f:
            for line in f:
                p=line.rstrip('\n').split('\t')
                self.idx[p[0]]=tuple(map(int,p[1:5]))
    def fetch(self,chrom,start,end):
        if chrom not in self.idx:raise ValueError(f'Unknown chromosome {chrom}')
        length,offset,bases,width=self.idx[chrom]
        if start<0 or end>length or end<start:raise ValueError(f'Coordinates outside {chrom}: {start}:{end}')
        if start==end:return ''
        a=offset + (start//bases)*width + start%bases
        b=offset + ((end-1)//bases)*width + ((end-1)%bases)+1
        self.file.seek(a)
        seq=self.file.read(b-a).replace(b'\n',b'').replace(b'\r',b'').upper().decode('ascii')
        if len(seq)!=end-start:raise ValueError('FASTA FAI indexing mismatch')
        return seq
    def close(self):self.file.close()

def score_record(fields,reference,indel_open=-5.,indel_extend=-1.):
    chrom=fields[2];start=int(fields[3])-1;cigar=fields[5]
    seq=fields[9].upper();quals=fields[10]
    if seq=='*' or quals=='*' or len(seq)!=len(quals):raise ValueError('Missing or inconsistent SEQ/QUAL')
    tokens=OPS.findall(cigar)
    if not tokens or ''.join(n+op for n,op in tokens)!=cigar:raise ValueError('Invalid CIGAR: '+cigar)
    i=0;r=start;logp=0.;usable=0
    for n_str,op in tokens:
        n=int(n_str)
        if op in 'M=X':
            target=reference.fetch(chrom,r,r+n)
            for j,g in enumerate(target):
                b=seq[i+j];q=max(0,ord(quals[i+j])-33)
                eps=max(1e-8,min(.75,10**(-q/10)))
                prob=(1-eps) if b==g and b in 'ACGT' else (eps/3 if b in 'ACGT' and g in 'ACGT' else .25)
                logp+=math.log(prob)
            usable+=n;i+=n;r+=n
        elif op=='I':
            logp+=indel_open+indel_extend*(n-1);i+=n
        elif op=='D':
            logp+=indel_open+indel_extend*(n-1);r+=n
        elif op=='N':r+=n
        elif op=='S':i+=n
        elif op in 'HP':pass
        else:raise ValueError('Unsupported CIGAR operation '+op)
    if i!=len(seq):raise ValueError('CIGAR/SEQ length mismatch')
    return logp,usable

def candidate_groups(path):
    opener=gzip.open if path.endswith('.gz') else open
    with opener(path,'rt') as inp:
        rd=csv.DictReader(inp,delimiter='\t')
        for k,grp in groupby(rd,key=lambda r:r['fragment_id']):
            yield k,{int(r['hi']):r['candidate_id'] for r in grp}

def sam_groups(source):
    inp=sys.stdin if source=='-' else open(source)
    def alignments():
        for line in inp:
            if line.startswith('@'):continue
            p=line.rstrip('\n').split('\t')
            if len(p)<11:raise ValueError('Malformed SAM record')
            flag=int(p[1])
            if flag&4 or flag&2048:continue
            yield p[0],p
    yield from ((key,[row for _,row in rows]) for key,rows in groupby(alignments(),key=lambda t:t[0]))
    if inp is not sys.stdin:inp.close()

def compute(candidates,sam,fasta,out,fai=None):
    reference=Faidx(fasta,fai)
    sam_iter=iter(sam_groups(sam))
    sam_row=next(sam_iter,None)
    with open(out,'w') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['fragment_id','candidate_id','log_quality_score','aligned_bases','status'])
        for fragment,expected in candidate_groups(candidates):
            while sam_row and sam_row[0]<fragment:sam_row=next(sam_iter,None)
            if not sam_row or sam_row[0]!=fragment:
                raise ValueError('Candidate fragment absent from SAM: '+fragment)
            bucket={}
            for fields in sam_row[1]:
                hi=next((int(t[5:]) for t in fields[11:] if t.startswith('HI:i:')),1)
                if hi not in expected:continue
                mate=1 if int(fields[1])&64 else (2 if int(fields[1])&128 else 0)
                if (hi,mate) in bucket:raise ValueError('Duplicated mate/HI in SAM')
                bucket[(hi,mate)]=score_record(fields,reference)
            for hi,candidate_id in expected.items():
                mates=sorted(m for hh,m in bucket if hh==hi)
                if mates not in ([0],[1,2]):raise ValueError(f'Incomplete candidate pair: {fragment} HI={hi}')
                value=sum(bucket[(hi,m)][0] for m in mates)
                bases=sum(bucket[(hi,m)][1] for m in mates)
                writer.writerow([fragment,candidate_id,f'{value:.12g}',bases,'STAR_ALIGNED'])
            sam_row=next(sam_iter,None)
    reference.close()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates',required=True);p.add_argument('--sam',required=True,help='Query-name-sorted SAM; use - for stdin')
    p.add_argument('--fasta',required=True,help='Uncompressed FASTA with .fai index')
    p.add_argument('--output',required=True)
    p.add_argument('--fai',default=None,help='Separate index for staged readonly FASTA')
    a=p.parse_args();compute(a.candidates,a.sam,a.fasta,a.output,a.fai)
