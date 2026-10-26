# PATCH 02 regression checkpoint

After a representative PATCH 02 run:

```bash
python tests/test_patch02_contract.py /path/to/sample.candidate_table.tsv
```

This test checks the schema, required STAR scoring fields, coordinate validity, and NH consistency within fragments. It intentionally does not test the legacy EM mathematics; PATCH 03 replaces that implementation.
