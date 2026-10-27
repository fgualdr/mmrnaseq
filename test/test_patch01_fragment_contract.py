#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, subprocess, tempfile
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SAM='''qU\t99\tchr1\t101\t255\t30M100N20M\t=\t501\t450\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:1\tNM:i:1
qU\t147\tchr1\t501\t255\t50M\t=\t101\t-450\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tNH:i:1\tHI:i:1\tAS:i:96\tnM:i:1\tNM:i:0
qM\t99\tchr2\t1001\t3\t50M\t=\t1301\t350\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tNH:i:2\tHI:i:1\tAS:i:98\tnM:i:0\tNM:i:0
qM\t147\tchr2\t1301\t3\t50M\t=\t1001\t-350\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tNH:i:2\tHI:i:1\tAS:i:98\tnM:i:0\tNM:i:0
qM\t355\tchr3\t2001\t3\t48M2S\t=\t2301\t350\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tNH:i:2\tHI:i:2\tAS:i:94\tnM:i:2\tNM:i:2
qM\t403\tchr3\t2301\t3\t50M\t=\t2001\t-350\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tNH:i:2\tHI:i:2\tAS:i:94\tnM:i:2\tNM:i:0
qX\t99\tchr4\t301\t20\t50M\t=\t601\t350\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tNH:i:1\tHI:i:1\tAS:i:100\tnM:i:0\tNM:i:0\tSA:Z:chr5,10,+,20M30S,10,2;
qX\t147\tchr4\t601\t20\t50M\t=\t301\t-350\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tNH:i:1\tHI:i:1\tAS:i:100\tnM:i:0\tNM:i:0
'''

def read_tsv(path):
    with open(path,newline='') as h:
        return list(csv.DictReader(h,delimiter='\t'))

def main():
    with tempfile.TemporaryDirectory() as td:
        d=Path(td); sam=d/'x.sam'; sam.write_text(SAM)
        outs={k:d/f'{k}.tsv' for k in ('cand','uniq','anchor','excluded','qc')}
        cmd=['awk','-v','single_end=0','-v','strandedness=reverse',
             '-v',f'candidate_out={outs["cand"]}','-v',f'unique_out={outs["uniq"]}',
             '-v',f'anchor_out={outs["anchor"]}','-v',f'excluded_out={outs["excluded"]}',
             '-v',f'qc_out={outs["qc"]}','-f',str(ROOT/'bin/mm_build_candidates.awk'),str(sam)]
        subprocess.run(cmd,check=True)
        c=read_tsv(outs['cand']); u=read_tsv(outs['uniq']); a=read_tsv(outs['anchor']); x=read_tsv(outs['excluded'])
        assert len(c)==2 and {r['hi'] for r in c}=={'1','2'}
        assert all(r['fragment_id']=='qM' for r in c)
        assert all(r['r1_chrom']!='.' and r['r2_chrom']!='.' for r in c)
        assert all(r['query_bases']=='100' for r in c)
        assert len(u)==1 and u[0]['fragment_id']=='qU' and u[0]['query_bases']=='100'
        weights=defaultdict(float)
        for r in a: weights[r['fragment_id']]+=float(r['weight'])
        assert abs(weights['qU']-1.0)<1e-9, weights
        assert len([r for r in a if r['mate']=='R1'])==2  # splice-aware blocks
        assert any(r['fragment_id']=='qX' for r in x)
    print('PASS PATCH01 exact fragment contract')
if __name__=='__main__': main()
