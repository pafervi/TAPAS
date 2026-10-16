#!/usr/bin/env bash
# Example 02: refinement for a fixed vacancy set.
# Takes input_cfg/ (one configuration from example 01) and exhaustively
# enumerates reduced-cation distributions with --balanced + --dedup-local.
# Produces out_02/ with ~160 unique configurations (from C(12,4)=495 raw).
set -euo pipefail

python ../../tapas.py POSCAR_CeO2_110_demo \
    --refine input_cfg/ \
    --balanced --dedup-local \
    --symprec 0.01 \
    --out out_02
