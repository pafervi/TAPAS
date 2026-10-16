# Example 01: general enumeration

Minimal end-to-end demo of the **general enumeration** mode on a didactic
CeO2 (110) p(1x1) slab (Ce12O24, 3 CeO2 layers, ~20 A vacuum). The slab
is not from any production study; its only purpose is to run in seconds
and exercise every stage of the pipeline.

## Input

- `POSCAR_CeO2_110_demo` — pristine slab, generated with ASE as a
  stoichiometric CeO2 (110) 3-layer slab (`make_example_poscar.py` under
  the parent `examples/` folder documents how it was built).

## Run

```bash
bash run.sh
```

or equivalently:

```bash
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
```

## Expected output

`out_01/` with:

- 190 `cfg_NNNNN/` folders (each with a `POSCAR` and a `config_info.txt`),
- a `summary.csv` with provenance header + 190 rows.

Terminal log should report:

```
Eligible anion sites: 24  |  vacancy combinations (filtered): 128
Vacancy sets after symmetry dedup: 14 (from 128)
Unique configurations: 190
Breakdown by vacancy layer:
  inner_2         56
  inner_3         22
  subsurf_bottom  83
  subsurf_top     50
  surf_bottom     29
  surf_top        140
```

Exact counts are reproducible with `--seed 42`. Changing the seed changes
the sampled reduced-cation subset per vacancy set.
