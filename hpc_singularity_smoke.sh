#!/usr/bin/env bash
set -euo pipefail
DOCKERHUB_USER=fgualdr
RUNTIME_REV="${RUNTIME_REV:-nxf3}"
REPRISE_REF="${REPRISE_REF:-6c1d829af2e88c079215da056dbe1b16a9a96967}"
REPRISE_SHORT="${REPRISE_REF:0:12}"
STD_PATH='/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'

images=(
  "docker://${DOCKERHUB_USER}/mmref-seqkit:2.14.0-${RUNTIME_REV}|seqkit|seqkit version"
  "docker://${DOCKERHUB_USER}/mmref-mummer4:4.0.1-${RUNTIME_REV}|nucmer|nucmer --version"
  "docker://${DOCKERHUB_USER}/mmref-reprise:1.0.1-${REPRISE_SHORT}-${RUNTIME_REV}|REPrise|REPrise 2>&1 | head -20"
  "docker://${DOCKERHUB_USER}/mmref-trf:4.09.1-${RUNTIME_REV}|trf|trf -v"
  "docker://${DOCKERHUB_USER}/mmref-cdhit:4.8.1-${RUNTIME_REV}|cd-hit-est|cd-hit-est -h 2>&1 | head -20"
)

for row in "${images[@]}"; do
    IFS='|' read -r image tool version_cmd <<<"$row"
    echo "==> $image"
    singularity exec --env "PATH=$STD_PATH" "$image" bash -c \
      "nxf-container-check; command -v $tool; $version_cmd"
done

echo 'SINGULARITYCE BASIC ACCEPTANCE: PASS'
