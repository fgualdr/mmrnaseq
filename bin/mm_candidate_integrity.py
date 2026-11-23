#!/usr/bin/env python3
"""Audit whether each fragment has all STAR-declared HI candidates.

This audits the existing table; it does not modify the validated resolver.
"""
import argparse
import csv
import gzip
from collections import Counter


def validate(source, destination):
    opener = gzip.open if str(source).endswith('.gz') else open
    summary = Counter()
    with opener(source, 'rt') as inp, open(destination, 'w') as out:
        r = csv.DictReader(inp, delimiter='\t')
        required = {'fragment_id','hi','nh','paired','r1_chrom','r2_chrom'}
        if not required.issubset(r.fieldnames or []):
            raise ValueError(f'Missing columns: {required - set(r.fieldnames or [])}')
        w = csv.writer(out, delimiter='\t', lineterminator='\n')
        w.writerow(['fragment_id','candidate_count','declared_nh','complete'])
        previous = None
        his = set()
        nhs = set()
        def emit():
            if previous is None: return
            complete = len(nhs)==1 and next(iter(nhs)) == len(his) and his == set(range(1, next(iter(nhs))+1))
            w.writerow([previous,len(his),','.join(map(str,sorted(nhs))),int(complete)])
            summary['complete' if complete else 'incomplete'] += 1
        for row in r:
            fragment = row['fragment_id']
            if fragment != previous:
                emit()
                if previous is not None and fragment <= previous:
                    raise ValueError('Candidate table is not strictly sorted by fragment_id')
                previous = fragment
                his=set(); nhs=set()
            hi=int(row['hi']); nh=int(row['nh'])
            if hi in his: raise ValueError(f'Duplicate HI={hi} for {fragment}')
            if nh<1 or hi<1 or hi>nh: raise ValueError(f'Invalid HI/NH for {fragment}')
            if int(row['paired']) and (row['r1_chrom']=='.' or row['r2_chrom']=='.'):
                raise ValueError(f'Incomplete pair {fragment} HI={hi}')
            his.add(hi);nhs.add(nh)
        emit()
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates',required=True);p.add_argument('--output',required=True)
    a=p.parse_args()
    print(dict(validate(a.candidates,a.output)))
