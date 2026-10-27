#!/usr/bin/env bash
set -euo pipefail

required_major=25
required_minor=10
required_patch=4
recommended='25.10.8'

if ! command -v nextflow >/dev/null 2>&1; then
    echo 'FAIL: nextflow is not on PATH' >&2
    exit 1
fi
if ! command -v java >/dev/null 2>&1; then
    echo 'FAIL: java is not on PATH' >&2
    exit 1
fi

java_major=$(java -version 2>&1 | awk -F'[".]' 'NR==1 { if ($2=="1") print $3; else print $2 }')
if [[ ! "$java_major" =~ ^[0-9]+$ ]] || (( java_major < 17 )); then
    echo "FAIL: Java >=17 is required; detected Java ${java_major:-unknown}" >&2
    exit 1
fi

nxf_version=$(nextflow -v 2>&1 | sed -nE 's/.*version[[:space:]]+([0-9]+\.[0-9]+\.[0-9]+).*/\1/p' | head -1)
if [[ -z "$nxf_version" ]]; then
    echo 'FAIL: could not parse Nextflow version' >&2
    nextflow -version >&2 || true
    exit 1
fi
IFS=. read -r major minor patch <<< "$nxf_version"

version_ge_required=false
if (( major > required_major )); then
    version_ge_required=true
elif (( major == required_major && minor > required_minor )); then
    version_ge_required=true
elif (( major == required_major && minor == required_minor && patch >= required_patch )); then
    version_ge_required=true
fi
if [[ "$version_ge_required" != true ]]; then
    echo "FAIL: Nextflow >=25.10.4 required; detected $nxf_version" >&2
    echo "Recommended runtime for this checkpoint: $recommended" >&2
    exit 1
fi

parser=${NXF_SYNTAX_PARSER:-v1}
if [[ "$parser" == 'v2' ]]; then
    echo 'FAIL: NXF_SYNTAX_PARSER=v2 is not supported by this checkpoint.' >&2
    echo 'Unset NXF_SYNTAX_PARSER (Nextflow 25.10.x defaults to the legacy/v1 parser).' >&2
    exit 1
fi
if (( major > 26 || (major == 26 && minor >= 4) )); then
    echo "FAIL: detected Nextflow $nxf_version. Nextflow >=26.04 defaults to the v2 parser." >&2
    echo "Use the validated 25.10.x runtime for this patch series; migrate to v2 in a dedicated later refactor." >&2
    exit 1
fi

repo=${1:-.}
if [[ ! -f "$repo/nextflow.config" || ! -f "$repo/main.nf" ]]; then
    echo "FAIL: '$repo' is not an mmrnaseq repository root" >&2
    exit 1
fi

if grep -RqsE "id ['\"]nf-boost|^[[:space:]]*boost[[:space:]]*\{" "$repo/conf/ieo.config" "$repo/conf/ieo.slurm.config"; then
    echo 'FAIL: nf-boost must not be declared in IEO override configs at this checkpoint' >&2
    exit 1
fi

NXF_SYNTAX_PARSER=v1 nextflow config "$repo" >/dev/null

echo "PASS Nextflow runtime compatibility: Nextflow $nxf_version, Java $java_major, parser v1"
