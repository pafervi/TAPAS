# Examples

End-to-end minimal demos of TAPAS. Both run in seconds and use a
didactic CeO2 (110) p(1x1) slab (Ce12O24) built with ASE just for this
purpose.

| Folder | Mode | What it demonstrates |
|---|---|---|
| `01_general_enumeration/` | general | enumerate symmetry-inequivalent 2-vacancy sets on the demo slab and stratify-sample 15 reduced-cation distributions per set |
| `02_refinement/` | `--refine` | take one configuration from example 01, fix its vacancy set and exhaustively enumerate reduced-cation distributions with `--balanced` + `--dedup-local` |

Each example has its own `README.md` with the exact command and the
expected counts.

## Reproducing the demo slab

The POSCAR used by both examples (`POSCAR_CeO2_110_demo`) is produced
by `make_example_poscar.py` (ASE-based, documented inline). It is a
3-layer CeO2 (110) slab with vacuum; included in the examples as a
pre-built file so the demos run without any additional dependencies.

```bash
python make_example_poscar.py          # regenerates POSCAR_CeO2_110_demo
```
