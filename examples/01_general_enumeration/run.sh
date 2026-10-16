#!/usr/bin/env bash
# Example 01: general enumeration on a didactic CeO2 (110) slab.
# Produces out_01/ with ~190 unique configurations.
set -euo pipefail

python ../../tapas.py POSCAR_CeO2_110_demo \
    --nvac 2 \
    --mode all \
    --vac-layers surf,subsurf,inner \
    --red-shell both \
    --host-coord 3 4 \
    --r-nn 2.6 --r-nnn 4.8 \
    --min-vac-dist 3.5 \
    --max-bonds-lost 1 \
    --max-red-combos 15 \
    --seed 42 \
    --symprec 0.01 \
    --out out_01
