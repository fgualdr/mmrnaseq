#!/usr/bin/env bash
set -euo pipefail
repo=${1:-.}

if [[ ! -f "$repo/nextflow.config" ]]; then
    echo "FAIL: '$repo' is not a repository root" >&2
    exit 1
fi

grep -q "id 'nf-boost@0.7.0'" "$repo/nextflow.config" || {
    echo 'FAIL: nf-boost@0.7.0 is not pinned in top-level nextflow.config' >&2
    exit 1
}
grep -q 'nf_boost_cleanup[[:space:]]*=[[:space:]]*false' "$repo/nextflow.config" || {
    echo 'FAIL: nf_boost_cleanup must default to false' >&2
    exit 1
}
if grep -RqsE "id ['\"]nf-boost|^[[:space:]]*boost[[:space:]]*\{" "$repo/conf/ieo.config" "$repo/conf/ieo.slurm.config"; then
    echo 'FAIL: nf-boost must not be declared in IEO override configs' >&2
    exit 1
fi

python3 - "$repo/nextflow_schema.json" <<'PY'
import json, sys
p=sys.argv[1]
d=json.load(open(p))
x=d['definitions']['generic_options']['properties'].get('nf_boost_cleanup')
assert x and x['type']=='boolean' and x['default'] is False
PY

NXF_SYNTAX_PARSER=v1 nextflow config "$repo" >/dev/null

echo 'PASS nf-boost runtime configuration: 0.7.0 pinned, cleanup=false by default'
