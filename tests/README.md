# Tests

Smoke tests for TAPAS. No pytest, no network, no external data, runs in
under 30 seconds and leaves nothing on disk.

## Run

```bash
bash tests/test_examples.sh
```

Exits non-zero on any failure.

## What is covered

| # | Test | Checks |
|---|---|---|
| T1 | package import | `tapas` imports and exposes `run_general`, `run_refine`, `canonical_label`, `symmetry_permutations`, `write_config`, `parse_config_info`, `validate_and_center`, `shells_for_vacancy`, and a non-empty `__version__` |
| T2 | CLI help | `python tapas.py --help` exits 0 |
| T3 | example 01 dry-run | general mode on the demo slab yields 190 unique configurations (deterministic with `--seed 42`) |
| T4 | example 02 dry-run | refinement with `--balanced --dedup-local` on the shipped `input_cfg/` yields 160 unique configurations |
| T5 | real end-to-end run | repeats example 01 without `--dry-run` and verifies that 190 `cfg_NNNNN/` folders are written, each with a `POSCAR` and a `config_info.txt`, plus a top-level `summary.csv` |
| T6 | config_info.txt format | the written `config_info.txt` is in English (no leftover Spanish strings) and round-trips through `parse_config_info` |

All counts are deterministic.
