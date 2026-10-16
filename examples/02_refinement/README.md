# Example 02: refinement for a fixed vacancy set

Demo of the **refinement** mode. Takes one of the configurations produced
by example 01 (`input_cfg/`, 2 oxygen vacancies on the didactic CeO2 (110)
slab) and exhaustively enumerates reduced-cation distributions compatible
with that fixed vacancy set, applying the physical `--balanced` constraint
and the local-symmetry `--dedup-local` approximation.

## Inputs

- `POSCAR_CeO2_110_demo` — same pristine slab as example 01.
- `input_cfg/config_info.txt` — vacancy set and NN/NNN pools for one
  specific configuration (2 vacancies: anion indices 12 and 13 in the
  original centered POSCAR).

The `source_poscar` field in `input_cfg/config_info.txt` has been
rewritten to `../POSCAR_CeO2_110_demo` so the example is self-contained.
The explicit positional POSCAR argument on the command line overrides this
anyway, so either path resolution works.

## Run

```bash
bash run.sh
```

or equivalently:

```bash
python ../../tapas.py POSCAR_CeO2_110_demo \
    --refine input_cfg/ \
    --balanced --dedup-local \
    --symprec 0.01 \
    --out out_02
```

## Expected output

`out_02/` with:

- 160 `cfg_NNNNN/` folders,
- a `summary.csv` with provenance header + 160 rows.

Terminal log should report:

```
Fixed vacancies (2): O indices [12, 13]
Reduced-cation pool: 12 sites  (NN=7, NNN_excl=5)
Exhaustive combinations: C(12,4) = 495
|G| = 16  |  |Stab(vac)| = 1  |  orbit vac = 16
Pool orbits under Stab(vac): 12 (singletons=12, max orbit=1)
[WARNING] Stab(vac_set) trivial: no global exact dedup available. ...
vac   12: pool=6, |Stab(vac_i)|=2, canonical pairs=9/15 (reduction 40%)
vac   13: pool=10, |Stab(vac_i)|=2, canonical pairs=27/45 (reduction 40%)
Local canonical product: 243 (from 675 raw, reduction 64.0%)
--dedup-local + --balanced: 160 candidate configurations
Unique configurations: 160
```

For a smaller output set, add `--pair-samples-per-type 1` to the command.
For the baseline (no reductions), drop both `--balanced` and
`--dedup-local` to see the full C(12, 4) = 495 configurations.
