# PATCH 01 — Exact SE/PE fragment-candidate contract

Apply only after **the current uploaded baseline** has passed its validation. This ZIP is a delta against that immediately preceding state. Earlier draft patch archives are obsolete.

## Scientific / computational change

Adds a parallel, non-disruptive fragment-level representation of STAR alternatives. It preserves both mates, splice-aware aligned blocks, observed orientation classes, query length, AS/NM tags, and fragment-normalized unique-anchor evidence. The legacy resolver remains active in this checkpoint.

## Files to add or replace

- `bin/mm_build_candidates.awk`
- `conf/modules.config`
- `modules/local/mm_candidate_prep.nf`
- `subworkflows/local/bam_filter_em.nf`
- `test/test_patch01_fragment_contract.py`

## Files to delete

- None

## Apply

From the repository root, copy the ZIP contents over the checkout, preserving paths. If `DELETE_FILES.txt` is present, remove every listed path. Then inspect before committing:

```bash
git status --short
git diff --stat
```

## Validate

```bash
python3 test/test_patch01_fragment_contract.py
```

## Suggested commit

```bash
git add -A
git commit -m "Add exact fragment-level multimapper contract"
```

## Validation scope

The packaged unit/static test checks the patch contract and synthetic known-answer cases. It does not replace a real Nextflow run on the target HPC environment. Test one representative sample (or a small reference for PATCH 08–09) before proceeding.
