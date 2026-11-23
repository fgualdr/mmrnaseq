#!/usr/bin/env python3
"""Propose genomic origins by projecting contiguous read alignments through homology.

Outputs are UNVALIDATED_PROPOSAL, never candidates for posterior/EM automatically.
This local proof-of-concept loads the PAF correspondence map in RAM.
"""
import argparse
import csv
import gzip
import re
from collections import defaultdict

CIGAR=re.compile(r'(\d+)([MIDNSHP=X])')

def propose(homology,candidates,output,limit=500):
    index=defaultdict(list)
    with open(homology) as inp:
        for h in csv.DictReader(inp,delimiter='\t'):
            lo=int(h['source_start0']);hi=int(h['source_end0'])
            index[h['source_chrom']].append((lo,hi,h))
    for chrom in index:index[chrom].sort(key=lambda x:x[0])
    opener=gzip.open if candidates.endswith('.gz') else open
    with opener(candidates,'rt') as inp,open(output,'w') as out:
        w=csv.writer(out,delimiter='\t',lineterminator='\n')
        w.writerow(['fragment_id','seed_candidate_id','target_chrom','target_start0','homology_strand','status'])
        fragment=None;seen=set()
        for row in csv.DictReader(inp,delimiter='\t'):
            if row['fragment_id']!=fragment:
                fragment=row['fragment_id'];seen=set()
            # No spliced/indel/clipped/partial projections. Reference blocks must
            # encompass the WHOLE contiguous alignment, not just its start.
            cigar=row['r1_cigar'];ops=CIGAR.findall(cigar)
            if not ops or ''.join(n+op for n,op in ops)!=cigar or any(op not in 'M=X' for n,op in ops):continue
            read_len=sum(int(n) for n,op in ops)
            chrom=row['r1_chrom'];start=int(row['r1_pos1'])-1;end=start+read_len
            for lo,hi,h in index.get(chrom,[]):
                if len(seen)>=limit:break
                if lo>start:break
                if not (lo<=start and end<=hi):continue
                offset=start-lo
                if h['strand']=='+':target_start=int(h['target_start0'])+offset
                else:target_start=int(h['target_end0'])-offset-read_len
                key=(h['target_chrom'],target_start,h['strand'])
                if key in seen:continue
                seen.add(key)
                w.writerow([fragment,row['candidate_id'],*key,'UNVALIDATED_PROPOSAL'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--homology',required=True);p.add_argument('--candidates',required=True)
    p.add_argument('--output',required=True);p.add_argument('--max-proposals',type=int,default=500)
    a=p.parse_args()
    if a.max_proposals<1:p.error('--max-proposals must be >=1')
    propose(a.homology,a.candidates,a.output,a.max_proposals)
