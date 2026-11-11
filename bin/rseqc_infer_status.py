#!/usr/bin/env python3
from __future__ import print_function

import argparse
import math
import re
import sys


def parse_args():
    p = argparse.ArgumentParser(description='Validate RSeQC infer_experiment output and classify strandedness.')
    p.add_argument('--sample', required=True)
    p.add_argument('--provided', required=True)
    p.add_argument('--stdout', required=True, dest='stdout_path')
    p.add_argument('--stderr', required=True, dest='stderr_path')
    p.add_argument('--cutoff', type=float, default=30.0,
                   help='Maximum opposite-strand percentage accepted; default 30 => >=70%% dominant strand.')
    p.add_argument('--output', required=True)
    return p.parse_args()


def read_text(path):
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        return handle.read()


def pct(x):
    return float(x) * 100.0


def fmt(x):
    return 'NA' if x is None else ('%.4f' % x)


def main():
    args = parse_args()
    out_text = read_text(args.stdout_path)
    err_text = read_text(args.stderr_path)

    usable = None
    m = re.search(r'Total\s+(\d+)\s+usable reads were sampled', err_text, flags=re.I)
    if m:
        usable = int(m.group(1))

    und = None
    sense = None
    antisense = None
    protocol = None

    if re.search(r'This is PairEnd Data', out_text, flags=re.I):
        protocol = 'paired_end'
        patterns = [
            ('und', r'Fraction of reads failed to determine:\s*([0-9.eE+-]+)'),
            ('sense', r'Fraction of reads explained by\s*"1\+\+,1--,2\+-,2-\+":\s*([0-9.eE+-]+)'),
            ('antisense', r'Fraction of reads explained by\s*"1\+-,1-\+,2\+\+,2--":\s*([0-9.eE+-]+)'),
        ]
    elif re.search(r'This is SingleEnd Data', out_text, flags=re.I):
        protocol = 'single_end'
        patterns = [
            ('und', r'Fraction of reads failed to determine:\s*([0-9.eE+-]+)'),
            ('sense', r'Fraction of reads explained by\s*"\+\+,--":\s*([0-9.eE+-]+)'),
            ('antisense', r'Fraction of reads explained by\s*"\+-, -\+":\s*([0-9.eE+-]+)'),
        ]
        # Standard RSeQC output has no space after the comma; keep a second regex below.
    else:
        patterns = []

    values = {}
    for key, pattern in patterns:
        mm = re.search(pattern, out_text)
        if mm:
            values[key] = pct(mm.group(1))

    if protocol == 'single_end' and 'antisense' not in values:
        mm = re.search(r'Fraction of reads explained by\s*"\+-,-\+":\s*([0-9.eE+-]+)', out_text)
        if mm:
            values['antisense'] = pct(mm.group(1))

    und = values.get('und')
    sense = values.get('sense')
    antisense = values.get('antisense')

    inference_status = 'VALID'
    reason = 'OK'
    inferred = 'NA'
    check = 'INDETERMINATE'

    if usable == 0:
        inference_status = 'NO_USABLE_READS'
        reason = 'RSeQC sampled zero usable reads; check BED12 content/reference compatibility'
    elif protocol is None:
        inference_status = 'INVALID_OUTPUT'
        reason = 'RSeQC output did not report paired-end or single-end protocol'
    elif any(v is None for v in (und, sense, antisense)):
        inference_status = 'INVALID_OUTPUT'
        reason = 'RSeQC output is missing one or more strandedness fractions'
    elif any((not math.isfinite(v)) or v < -1e-6 or v > 100.0001 for v in (und, sense, antisense)):
        inference_status = 'INVALID_OUTPUT'
        reason = 'RSeQC strandedness fractions are outside [0,100]'
    elif abs((und + sense + antisense) - 100.0) > 1.0:
        inference_status = 'INVALID_OUTPUT'
        reason = 'RSeQC strandedness fractions do not sum to approximately 100%'
    else:
        threshold = 100.0 - args.cutoff
        if sense >= threshold:
            inferred = 'forward'
        elif antisense >= threshold:
            inferred = 'reverse'
        else:
            inferred = 'unstranded'

        provided = args.provided.lower()
        if provided in ('forward', 'reverse', 'unstranded'):
            check = 'PASS' if provided == inferred else 'FAIL'
            if check == 'FAIL':
                reason = 'Provided strandedness disagrees with valid RSeQC inference'
        else:
            check = 'INDETERMINATE'
            reason = 'Provided strandedness is not forward/reverse/unstranded'

    header = [
        'Sample', 'Provided strandedness', 'Inference status', 'Check',
        'Inferred strandedness', 'Sense (%)', 'Antisense (%)',
        'Undetermined (%)', 'Usable reads', 'Reason'
    ]
    row = [
        args.sample, args.provided, inference_status, check, inferred,
        fmt(sense), fmt(antisense), fmt(und),
        'NA' if usable is None else str(usable), reason
    ]

    with open(args.output, 'w', encoding='utf-8') as handle:
        handle.write('\t'.join(header) + '\n')
        handle.write('\t'.join(row) + '\n')

    if inference_status != 'VALID':
        print('WARNING: %s strandedness inference is %s: %s' %
              (args.sample, inference_status, reason), file=sys.stderr)
    elif check == 'FAIL':
        print('WARNING: %s strandedness mismatch: provided=%s inferred=%s' %
              (args.sample, args.provided, inferred), file=sys.stderr)


if __name__ == '__main__':
    main()
