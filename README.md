# TAPAS: Topology-Aware Polaron Arrangement Sampler

Symmetry-aware enumeration and refinement of anion-vacancy / reduced-cation configurations in oxide slabs.

TAPAS generates the set of symmetry-inequivalent configurations obtained by removing `nvac` anion atoms from an oxide slab and placing the resulting reduced cations on candidate host sites around the vacancies. Each configuration is written as a POSCAR with the reduced cations substituted by a **proxy species** (by default `La` for a `Ce`-based host). The proxy imprints the correct local lattice distortion during a force-field pre-relaxation without requiring an electronic-structure description of the localized state. It is substituted back to the host cation before the final DFT calculation, which is outside the scope of this script.

The defaults target fluorite-type oxides such as CeO2 (one O vacancy releases two electrons localizing on two Ce), but `--host-cation`, `--proxy-cation`, `--anion`, and `--red-per-vac` allow the workflow to be applied to other hosts.

The script operates in two modes:

- **General enumeration (default).** Given a pristine slab POSCAR, enumerate every symmetry-inequivalent vacancy set compatible with the requested physical filters and, for each unique vacancy set, stratify-sample reduced-cation distributions across its NN+NNN pool. Produces a broad, physics-neutral initial exploration of the configurational space.
- **Refinement (`--refine CFG_DIR`).** Given a `config_info.txt` produced by a previous general run (typically the most stable configuration identified by subsequent force-field and DFT screening), fix the vacancy set and exhaustively enumerate all reduced-cation distributions compatible with it. Three optional reductions (physical constraint, local-symmetry deduplication, pair-type sampling) can be combined to keep the enumeration tractable at high vacancy counts.

The two modes share all utilities (slab centering, radii calibration, spglib symmetry, canonical labels, POSCAR writing with verified proxy indices, `config_info.txt` output).

## Requirements

- Python >= 3.8
- `numpy`
- `ase`
- `spglib` (strongly recommended; without it the general-mode deduplication falls back to a distance/depth fingerprint that may merge non-equivalent configurations, and refinement-mode options that rely on local symmetry are disabled)

## Repository layout

The codebase is organised as a small Python package `tapas/` plus a thin CLI entry point `tapas.py` at the repository root.

```
TAPAS/
├── tapas.py                       # CLI entry point: arg parsing + dispatch
├── tapas/                         # package (programmatic API)
│   ├── __init__.py                # public API re-exports
│   ├── common.py                  # die, warn, stable_hash, VALID_LAYER_KEYS
│   ├── geometry.py                # slab centering, layer detection, radii
│   │                              # calibration, coordination check, NN/NNN
│   │                              # shell lookups
│   ├── symmetry.py                # spglib permutations, canonical labels,
│   │                              # distance/depth fingerprint fallback
│   ├── io.py                      # POSCAR writers (verified by composition
│   │                              # and coordinate round-trip),
│   │                              # config_info.txt writer and parser
│   ├── general.py                 # run_general: enumeration mode
│   └── refine.py                  # run_refine: refinement mode
│                                  # (--balanced, --dedup-local,
│                                  #  --pair-samples-per-type, --min-red-dist)
├── requirements.txt
├── LICENSE                        # MIT
├── examples/
│   ├── README.md
│   ├── make_example_poscar.py     # ASE-based builder of the demo slab
│   ├── 01_general_enumeration/
│   │   ├── README.md
│   │   ├── run.sh
│   │   └── POSCAR_CeO2_110_demo
│   └── 02_refinement/
│       ├── README.md
│       ├── run.sh
│       ├── POSCAR_CeO2_110_demo
│       └── input_cfg/             # a cfg from example 01, used as input
│           ├── POSCAR
│           └── config_info.txt
└── tests/
    ├── README.md
    └── test_examples.sh           # smoke tests (dry-run, <10 s)
```

**Module dependency graph** (`A -> B` means A imports from B):

```
tapas.py (root CLI)
    │
    ├─> general ─┬─> geometry ─> common
    │            ├─> symmetry ─┐
    │            └─> io ───────┤
    └─> refine ──┬─> geometry ─┤
                 ├─> symmetry ─┤
                 └─> io ───────┘
                               └─> common
```

The three leaf modules (`geometry`, `symmetry`, `io`) only depend on `common`; `general` and `refine` are the two orchestration layers; the root script owns the CLI (argument parsing and dispatch to `run_general` or `run_refine`).

### Data and control flow

The diagram below shows the two operating modes side by side. Both share the pre-processing stage at the top (geometry and symmetry modules) and the verified writer at the bottom; the red dashed block marks the three combinable approximate reductions that only apply in refinement mode.

<p align="center">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1040 880" style="max-width:1040px;width:100%;height:auto" role="img" aria-label="Flow diagram of tapas. A pristine POSCAR enters a shared pre-processing stage, then splits into general enumeration on the left and refinement on the right, which rejoin at a verified writer that emits one cfg folder per unique configuration plus a summary.csv.">
  <defs>
    <marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <polygon points="0,0 10,5 0,10" fill="currentColor"/>
    </marker>
    <marker id="arrAccent" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <polygon points="0,0 10,5 0,10" fill="#d64545"/>
    </marker>
  </defs>
  <g font-size="12" fill="currentColor" font-family="system-ui, sans-serif">
    <rect x="420" y="20" width="200" height="44" rx="6" fill="none" stroke="currentColor"/>
    <text x="520" y="39" text-anchor="middle" font-weight="600">pristine slab POSCAR</text>
    <text x="520" y="55" text-anchor="middle" opacity="0.7">(host + anion, no vacancies)</text>
    <line x1="520" y1="64" x2="520" y2="92" stroke="currentColor" marker-end="url(#arr)"/>
    <text x="212" y="86" font-size="11" opacity="0.6" font-style="italic">shared pre-processing (geometry.py, symmetry.py)</text>
    <rect x="200" y="94" width="640" height="138" rx="8" fill="none" stroke="currentColor" stroke-dasharray="3 3" opacity="0.5"/>
    <rect x="220" y="110" width="180" height="44" rx="5" fill="none" stroke="currentColor"/>
    <text x="310" y="128" text-anchor="middle">validate_and_center</text>
    <text x="310" y="144" text-anchor="middle" opacity="0.75">vacuum detect, z-shift</text>
    <line x1="400" y1="132" x2="430" y2="132" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="432" y="110" width="180" height="44" rx="5" fill="none" stroke="currentColor"/>
    <text x="522" y="128" text-anchor="middle">classify_anion_layers</text>
    <text x="522" y="144" text-anchor="middle" opacity="0.75">surf, subsurf, inner</text>
    <line x1="612" y1="132" x2="642" y2="132" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="644" y="110" width="180" height="44" rx="5" fill="none" stroke="currentColor"/>
    <text x="734" y="128" text-anchor="middle">suggest_radii, coord check</text>
    <text x="734" y="144" text-anchor="middle" opacity="0.75">r_nn, r_nnn</text>
    <rect x="320" y="178" width="400" height="44" rx="5" fill="none" stroke="currentColor"/>
    <text x="520" y="196" text-anchor="middle">symmetry_permutations (spglib)</text>
    <text x="520" y="212" text-anchor="middle" opacity="0.75">atom-index perms for the pristine slab space group</text>
    <line x1="520" y1="232" x2="520" y2="252" stroke="currentColor"/>
    <line x1="260" y1="252" x2="780" y2="252" stroke="currentColor"/>
    <line x1="260" y1="252" x2="260" y2="276" stroke="currentColor"/>
    <line x1="780" y1="252" x2="780" y2="276" stroke="currentColor"/>
    <rect x="80" y="276" width="360" height="30" rx="4" fill="currentColor" opacity="0.08"/>
    <rect x="80" y="276" width="360" height="30" rx="4" fill="none" stroke="currentColor"/>
    <text x="260" y="296" text-anchor="middle" font-weight="700" font-size="13">GENERAL ENUMERATION (default)</text>
    <line x1="260" y1="306" x2="260" y2="328" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="80" y="330" width="360" height="60" rx="6" fill="none" stroke="currentColor"/>
    <text x="260" y="350" text-anchor="middle" font-weight="600">enumerate candidate vacancy sets</text>
    <text x="260" y="366" text-anchor="middle" opacity="0.75">C(allowed, nvac), filtered by</text>
    <text x="260" y="381" text-anchor="middle" opacity="0.75">--mode, --min-vac-dist, --max-bonds-lost</text>
    <line x1="260" y1="390" x2="260" y2="412" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="80" y="414" width="360" height="46" rx="6" fill="none" stroke="currentColor"/>
    <text x="260" y="433" text-anchor="middle" font-weight="600">dedup vacancy sets under Stab(slab)</text>
    <text x="260" y="449" text-anchor="middle" opacity="0.75">canonical lex-min over symmetry perms</text>
    <line x1="260" y1="460" x2="260" y2="482" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="80" y="484" width="360" height="76" rx="6" fill="none" stroke="currentColor"/>
    <text x="260" y="504" text-anchor="middle" font-weight="600">for each unique vacancy set:</text>
    <text x="260" y="520" text-anchor="middle" opacity="0.75">build NN+NNN pool, split into pair-type strata,</text>
    <text x="260" y="534" text-anchor="middle" opacity="0.75">sample ~--max-red-combos per set uniformly</text>
    <text x="260" y="549" text-anchor="middle" opacity="0.75">(deterministic per-set seed from --seed)</text>
    <line x1="260" y1="560" x2="260" y2="582" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="80" y="584" width="360" height="46" rx="6" fill="none" stroke="currentColor"/>
    <text x="260" y="603" text-anchor="middle" font-weight="600">canonical dedup of (vac, red)</text>
    <text x="260" y="619" text-anchor="middle" opacity="0.75">lex-min over symmetry perms, sha1 hash</text>
    <rect x="600" y="276" width="360" height="30" rx="4" fill="currentColor" opacity="0.08"/>
    <rect x="600" y="276" width="360" height="30" rx="4" fill="none" stroke="currentColor"/>
    <text x="780" y="296" text-anchor="middle" font-weight="700" font-size="13">REFINEMENT (--refine CFG_DIR)</text>
    <line x1="780" y1="306" x2="780" y2="328" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="600" y="330" width="360" height="60" rx="6" fill="none" stroke="currentColor"/>
    <text x="780" y="350" text-anchor="middle" font-weight="600">parse CFG_DIR/config_info.txt</text>
    <text x="780" y="366" text-anchor="middle" opacity="0.75">fixed vacancy indices, per-vacancy</text>
    <text x="780" y="381" text-anchor="middle" opacity="0.75">NN/NNN pools (as recorded)</text>
    <line x1="780" y1="390" x2="780" y2="412" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="600" y="414" width="360" height="46" rx="6" fill="none" stroke="currentColor"/>
    <text x="780" y="433" text-anchor="middle" font-weight="600">build union pool (NN ∪ NNN excl.)</text>
    <text x="780" y="449" text-anchor="middle" opacity="0.75">report |G|, |Stab(vac_set)|, pool orbits</text>
    <line x1="780" y1="460" x2="780" y2="482" stroke="currentColor" marker-end="url(#arr)"/>
    <g fill="#d64545" stroke="#d64545">
      <rect x="600" y="484" width="360" height="130" rx="6" fill="none" stroke-dasharray="5 3"/>
      <text x="780" y="502" text-anchor="middle" font-weight="700" fill="#d64545">optional reductions (combinable)</text>
      <text x="612" y="522" opacity="0.95" fill="#d64545">· --balanced</text>
      <text x="625" y="537" opacity="0.75" fill="#d64545">red_per_vac polarons per vacancy, own pool</text>
      <text x="612" y="557" opacity="0.95" fill="#d64545">· --dedup-local</text>
      <text x="625" y="572" opacity="0.75" fill="#d64545">canonicalize under Stab(vac_i) per vacancy</text>
      <text x="612" y="592" opacity="0.95" fill="#d64545">· --pair-samples-per-type N</text>
      <text x="625" y="607" opacity="0.75" fill="#d64545">N reps per pair-type class (needs both above)</text>
    </g>
    <line x1="780" y1="614" x2="780" y2="640" stroke="#d64545" marker-end="url(#arrAccent)"/>
    <rect x="600" y="642" width="360" height="46" rx="6" fill="none" stroke="currentColor"/>
    <text x="780" y="661" text-anchor="middle" font-weight="600">exhaustive enumeration of red-cations</text>
    <text x="780" y="677" text-anchor="middle" opacity="0.75">canonical dedup, optional --min-red-dist filter</text>
    <line x1="260" y1="630" x2="260" y2="710" stroke="currentColor"/>
    <line x1="780" y1="688" x2="780" y2="710" stroke="currentColor"/>
    <line x1="260" y1="710" x2="780" y2="710" stroke="currentColor"/>
    <line x1="520" y1="710" x2="520" y2="738" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="320" y="740" width="400" height="62" rx="6" fill="none" stroke="currentColor"/>
    <text x="520" y="760" text-anchor="middle" font-weight="600">write_config + write_config_info</text>
    <text x="520" y="776" text-anchor="middle" opacity="0.75">POSCAR (anions removed, host3+ written as proxy),</text>
    <text x="520" y="791" text-anchor="middle" opacity="0.75">config_info.txt (verified by coordinate matching)</text>
    <line x1="520" y1="802" x2="520" y2="828" stroke="currentColor" marker-end="url(#arr)"/>
    <rect x="320" y="830" width="400" height="42" rx="6" fill="none" stroke="currentColor"/>
    <text x="520" y="850" text-anchor="middle" font-weight="600">out_dir/cfg_NNNNN/, out_dir/summary.csv</text>
    <text x="520" y="865" text-anchor="middle" opacity="0.75">one folder per unique configuration</text>
  </g>
</svg>
</p>

### Programmatic use

Every utility is importable from the package root:

```python
from tapas import (
    symmetry_permutations, canonical_label, fingerprint_fallback,
    validate_and_center, classify_anion_layers, suggest_radii,
    shells_for_vacancy,
    write_config, parse_config_info,
    run_general, run_refine,
)
```

Use these if you want to reuse pieces of the workflow from a notebook or your own script without going through the CLI.

## How it works

### General mode

#### 1. Input validation and slab centering

The script checks that the POSCAR contains the expected species, detects the vacuum gap along z, and shifts the slab so the vacuum sits at the cell boundary. This ensures layer detection and the top/bottom face assignment are unambiguous.

#### 2. Anion layer detection

Anion atoms are grouped into atomic layers by detecting gaps in their z coordinates. The layer threshold is chosen adaptively (largest jump in the sorted list of pairwise z separations), avoiding a brittle fixed tolerance. Layers are tagged as `surf_top`, `surf_bottom`, `subsurf_top`, `subsurf_bottom`, or `inner_N` depending on their position relative to both slab surfaces. Vacancies are only placed on anion sites belonging to the layers selected by `--vac-layers`.

#### 3. Coordination shells

For each vacancy site, the script computes two coordination shells of the host cation:

- **NN**: host-cation sites within `--r-nn` of the vacancy.
- **NNN**: host-cation sites between `--r-nn` and `--r-nnn`.

Radii are auto-calibrated from the actual anion-to-host-cation distance distribution of the slab if not provided explicitly. The `--host-coord` check verifies that every eligible anion site has an expected number of host-cation neighbors within `--r-nn`, catching mis-set radii before they corrupt the NN/NNN pools.

#### 4. Vacancy combination filtering

All combinations of `nvac` anion sites from the eligible pool are generated and filtered by:

- `--mode`: restricts to vacancies on the same face, on both faces simultaneously, or unrestricted.
- `--min-vac-dist`: discards vacancy sets where any two vacancies are closer than this distance.
- `--max-bonds-lost`: discards vacancy sets where any single host-cation site would lose more than N bonds to vacancies (i.e., is shared by more than N vacancies). This limits the local perturbation per cation.

#### 5. Symmetry deduplication of vacancy sets

Vacancy sets are deduplicated first using the symmetry operations of the pristine host slab (via spglib). Two vacancy sets related by a symmetry operation are represented by a single canonical set; only that representative is processed for reduced-cation placement. This is done before sampling to avoid over-counting mirror-equivalent configurations.

#### 6. Stratified sampling of reduced-cation placements

For each symmetry-unique vacancy set, `--red-per-vac * nvac` reduced-cation sites are chosen from the NN+NNN pool. Because the number of placements grows combinatorially with vacancy count, a stratified sampling strategy is used:

1. **Stratify by NN/NNN composition**: a stratum is defined by (n_NN, n_NNN) with n_NN + n_NNN = n_red. For 2 reduced cations these are: (2,0), (1,1), (0,2).
2. **Uniform budget**: the `--max-red-combos` quota is split equally across viable strata. No NN/NNN bias is imposed a priori.
3. **Deterministic per-set seed**: within each stratum, placements are drawn at random with a seed derived from `--seed` combined with a hash of the vacancy set, so results are fully reproducible and independent across sets.
4. If a stratum contains fewer placements than its quota, all of them are taken (exhaustive within that stratum).

The `--red-shell` argument can restrict which strata are considered.

#### 7. Exact deduplication of (vacancy, reduced-cation) pairs

Each (vacancy set, reduced-cation placement) pair is mapped to its canonical label: the lexicographic minimum of (sorted vacancy indices, sorted reduced-cation indices) over all host symmetry operations. Two pairs that differ only by a symmetry operation produce the same canonical label and are merged. This is exact, not approximate.

### Refinement mode (`--refine CFG_DIR`)

The refinement path is designed to be run *after* the general enumeration and its downstream force-field and DFT screening have identified a specific vacancy set of interest. It takes the `config_info.txt` of that configuration as input, extracts the vacancy indices and the per-vacancy NN/NNN pools, and enumerates reduced-cation distributions exhaustively (no sampling).

#### 1. Parsing the input configuration

The script reads `CFG_DIR/config_info.txt` and extracts:

- the vacancy indices (anion atoms in the original pristine POSCAR, 0-indexed),
- the per-vacancy NN and NNN host-cation pools, as recorded by the generator,
- the path to the source pristine POSCAR (used by default; can be overridden by passing a POSCAR as positional argument).

Pools are **not recomputed from radii**: they are taken verbatim from the input `config_info.txt`, guaranteeing consistency with the pools that defined the original vacancy set.

#### 2. Baseline enumeration

Starting from the union pool `union(NN_i ∪ NNN_i)` across all vacancies, the script enumerates `C(|pool|, red_per_vac × n_vac)` reduced-cation distributions. Each distribution is then deduplicated under the pristine slab space group with the same canonical-label mechanism as the general mode.

#### 3. Symmetry diagnostic

Before enumeration, the script reports:

- `|G|`: number of symmetry operations found by spglib for the pristine slab.
- `|Stab(vac_set)|`: subset of `|G|` that preserves the vacancy set as a set. Only this subgroup can produce exact global deduplication of reduced-cation distributions for the fixed vacancy set.
- Pool orbits under `Stab(vac_set)`: how the pool decomposes into equivalence classes.

If `|Stab(vac_set)| = 1` (common for canonical vacancy representatives selected by the general enumeration), no further exact dedup is possible on top of what the general run already captured.

#### 4. Optional reductions

Three reductions can be combined on top of the baseline canonical-label dedup.

**`--balanced`** enforces a physical constraint: each vacancy contributes exactly `--red-per-vac` reduced cations drawn from its own NN+NNN pool. Rejects distributions where polarons effectively migrate from one vacancy's local environment to another's. This is a model assumption, not a symmetry statement.

**`--dedup-local`** replaces the global canonical label by a per-vacancy local canonical label computed under `Stab(vac_i)`, the stabiliser of each individual vacancy. Formally:

- with `--balanced`: enumerate canonical 2-subset representatives of each vacancy's pool under `Stab(vac_i)`, then take their cartesian product and keep only distinct flat sets.
- without `--balanced`: enumerate all reduced-cation combinations from the union pool and deduplicate by the tuple of per-vacancy local patterns `(red ∩ pool_v)` canonized under `Stab(vac_i)`.

This is an approximation: two configurations that are locally equivalent at every vacancy are grouped even if they are not globally equivalent under `Stab(vac_set)`. The approximation is tightest when the vacancies are well separated and cross-vacancy interactions are small.

**`--pair-samples-per-type N`** (requires `--dedup-local` + `--balanced`): after the local-canonical pair enumeration, group each vacancy's canonical pairs by pair type (NN-NN, NN-NNN, NNN-NNN) and keep only N representatives per type. For N=1 the final count is bounded by `~3^n_vac`, which is the smallest set compatible with the hypothesis that energetics are dominated by pair-type class rather than by the specific atomic choice within a type. This is a stronger approximation than `--dedup-local` alone.

#### 5. Output

As in general mode, one `cfg_NNNNN/` folder per unique configuration under `--out` (default: `<cfg_dir>_refined/`), each containing a `POSCAR` and a `config_info.txt`. A top-level `summary.csv` lists every refined configuration.

## Output format

### `summary.csv`

One row per configuration. Header rows document the command, timestamp, library versions and (in refinement mode) the source configuration.

**General mode (default)** columns:

| Column | Description |
|---|---|
| `config_id` | folder name, e.g. `cfg_00042` |
| `vac_sites` | anion atom indices in the original POSCAR (pipe-separated) |
| `vac_tags` | layer tags of each vacancy (pipe-separated) |
| `vac_faces` | face assignment of each vacancy: `top`, `bottom`, or `mid` |
| `vac_xyz` | Cartesian coordinates of each vacancy in the centered slab |
| `proxy_sites` | host-cation indices replaced by the proxy (pipe-separated) |
| `proxy_shell_labels` | `NN` or `NNN` for each proxy site |
| `pair_type` | NN/NNN composition string, e.g. `NN-NNN` |
| `dedup_label` | 16-character sha1 of the canonical label |
| `truncated` | always `no` in this version |

**General mode with `--no-proxy`** columns (vacancy-only):

| Column | Description |
|---|---|
| `config_id` | folder name |
| `vac_sites` | anion atom indices |
| `vac_tags` | layer tags |
| `vac_faces` | `top`, `bottom`, or `mid` |
| `vac_xyz` | vacancy Cartesian coordinates |

**Refinement mode** columns:

| Column | Description |
|---|---|
| `config_id` | folder name |
| `red_sites` | reduced-cation host-cation indices in the original POSCAR (pipe-separated) |
| `red_shell_labels` | `NN` or `NNN` for each reduced-cation site |
| `pair_type` | NN/NNN composition string |
| `dedup_label_hash` | 16-character sha1 of the canonical label |

The vacancy information in refinement mode is identical for every row (fixed vacancy set) and is instead recorded once in the provenance header of `summary.csv` and in each `cfg_NNNNN/config_info.txt`.

### `cfg_NNNNN/POSCAR`

VASP POSCAR with anions removed and reduced-cation sites written as the proxy species. Sorted by species (`sort=True`), fractional coordinates.

### `cfg_NNNNN/config_info.txt`

Human-readable record of the configuration, with:

- provenance (command, timestamp, library versions, deduplication method and tolerance),
- `config_id`, source POSCAR path, number of vacancies, mode tag (`all`, `same_face`, `both_faces`, or `refinement_exhaustive`), pair type, and the canonical deduplication hash,
- composition line in the written POSCAR,
- per-vacancy entries with anion index, layer tag (general mode only), face assignment, host coordination count, and Cartesian coordinates,
- per-reduced-cation entries with host-cation index in the original POSCAR, shell label (`NN` or `NNN`), Cartesian coordinates, and distances to every vacancy,
- per-vacancy NN and NNN pools (lists of host-cation indices),
- proxy indices in the **written** POSCAR, verified by coordinate matching, in both base-0 and base-1 (VASP) convention,
- post-DFT control checklist.

This file is the input used by the refinement mode to reconstruct the vacancy set and the per-vacancy pools.

## Usage

### General mode

**Single vacancy, all layers, exhaustive (small combinatorial space):**

```bash
python tapas.py POSCAR_slab \
    --host-cation Ce --proxy-cation La --anion O \
    --nvac 1 \
    --mode all \
    --vac-layers surf,subsurf,inner \
    --red-shell both \
    --host-coord 3 4 \
    --r-nn 2.6 --r-nnn 4.8 \
    --symprec 0.01 \
    --out 1vac_la
```

**Multiple vacancies with physical filters and sampling:**

```bash
python tapas.py POSCAR_slab \
    --host-cation Ce --proxy-cation La --anion O \
    --nvac 3 \
    --mode all \
    --vac-layers surf,subsurf,inner \
    --red-shell both \
    --host-coord 3 4 \
    --r-nn 2.6 --r-nnn 4.8 \
    --min-vac-dist 3.5 \
    --max-bonds-lost 1 \
    --max-red-combos 25 \
    --seed 42 \
    --symprec 0.01 \
    --out 3vac_la
```

**Vacancy-only mode (no proxy substitution):**

```bash
python tapas.py POSCAR_slab \
    --host-cation Ce --anion O \
    --nvac 2 \
    --mode all \
    --vac-layers surf,subsurf,inner \
    --min-vac-dist 3.5 \
    --max-bonds-lost 1 \
    --symprec 0.01 \
    --no-proxy \
    --out 2vac_novac
```

**Dry run (count only, no folders written):**

```bash
python tapas.py POSCAR_slab \
    --nvac 2 \
    --vac-layers surf,subsurf,inner \
    --min-vac-dist 3.5 \
    --max-bonds-lost 1 \
    --max-red-combos 25 \
    --seed 42 \
    --dry-run \
    --out 2vac_la
```

### Refinement mode

**Exhaustive baseline** (no reductions, only global canonical-label dedup):

```bash
python tapas.py POSCAR_slab \
    --refine 3vac_la/cfg_12345/ \
    --symprec 0.01 \
    --out cfg_12345_refined
```

**Physical constraint only:**

```bash
python tapas.py POSCAR_slab \
    --refine 3vac_la/cfg_12345/ \
    --balanced \
    --symprec 0.01 \
    --out cfg_12345_refined
```

**Local-symmetry approximation combined with the physical constraint:**

```bash
python tapas.py POSCAR_slab \
    --refine 3vac_la/cfg_12345/ \
    --balanced --dedup-local \
    --symprec 0.01 \
    --out cfg_12345_refined
```

**Aggressive pair-type reduction (small final count, strong approximation):**

```bash
python tapas.py POSCAR_slab \
    --refine 3vac_la/cfg_12345/ \
    --balanced --dedup-local --pair-samples-per-type 1 \
    --symprec 0.01 \
    --out cfg_12345_refined
```

**POSCAR optional in refinement:** if the `source_poscar` field recorded in `config_info.txt` is accessible from the current working directory, the positional POSCAR argument can be omitted:

```bash
python tapas.py --refine 3vac_la/cfg_12345/ --balanced --dedup-local
```

**Dry run (count only):**

```bash
python tapas.py POSCAR_slab \
    --refine 3vac_la/cfg_12345/ \
    --balanced --dedup-local \
    --dry-run
```

## Arguments

### Positional

| Argument | Description |
|---|---|
| `poscar` | Path to the input POSCAR of the pristine slab. Required in general mode. Optional in `--refine` mode (defaults to the `source_poscar` field recorded in the input `config_info.txt`). |

### Chemistry

| Argument | Default | Description |
|---|---|---|
| `--host-cation` | `Ce` | Symbol of the host cation that becomes reduced. |
| `--proxy-cation` | `La` | Symbol of the proxy species written in the output POSCARs to represent a reduced cation. |
| `--anion` | `O` | Symbol of the anion. Vacancies are placed on anion sites. |
| `--red-per-vac` | `2` | Number of reduced cations produced per vacancy. For CeO2 each O vacancy releases 2 electrons, so the default is 2. |

### Vacancy selection (general mode)

| Argument | Default | Description |
|---|---|---|
| `--nvac` | `1` | Number of vacancies per configuration. |
| `--vac-layers` | `surf,subsurf` | Comma-separated list of anion layer groups eligible for vacancies. Valid keys: `surf` (outermost anion layers on both faces), `subsurf` (second anion layer on both faces), `inner` (all remaining layers). |
| `--mode` | `all` | Face constraint for vacancy sets. `all`: no constraint. `same_face`: all vacancies must be on the same face (top or bottom). `both_faces`: the set must include at least one vacancy on each face. |
| `--min-vac-dist` | `0.0` | Minimum vacancy-to-vacancy distance in Angstrom. Vacancy sets with any pair closer than this are discarded. `0` disables the filter. |
| `--max-bonds-lost` | `None` | Maximum number of anion neighbors a single host-cation site may lose. A value of `1` means no host cation is shared by two or more vacancies. Disabled by default. |

### Coordination shells and radii

| Argument | Default | Description |
|---|---|---|
| `--r-nn` | auto | Cutoff in Angstrom for the first coordination shell (NN). Auto-calibrated from the anion-to-host-cation distance distribution if not set. |
| `--r-nnn` | auto | Cutoff in Angstrom for the second coordination shell (NNN). Auto-calibrated if not set. |
| `--host-coord` | `None` | Expected host-cation coordination numbers for eligible anion sites, e.g. `--host-coord 3 4`. If an anion site has a count outside this set the script aborts with a suggested radius correction (unless `--force`). If not set, the coordination check is skipped. |
| `--red-shell` | `both` | Which coordination shell(s) are candidates for reduced-cation placement. `nn`: NN only. `nnn`: NNN only. `both`: either shell. `nnn_pref`: both shells, but at most `nvac` reduced cations in NN. `nnn_strict`: NNN only (equivalent to `nnn`). |

### Sampling (general mode)

| Argument | Default | Description |
|---|---|---|
| `--max-red-combos` | `25` | Total quota of sampled reduced-cation placements per vacancy set, split uniformly across NN/NNN strata. If a stratum has fewer placements than its share, all are taken. |
| `--seed` | `42` | Base random seed for the stratified sampling. The actual seed per vacancy set is derived from this value combined with a hash of the vacancy set, ensuring reproducibility and independence across sets. |

### Symmetry and layers

| Argument | Default | Description |
|---|---|---|
| `--symprec` | `1e-3` | spglib symmetry tolerance in Angstrom. Raise to `0.01` or `0.05` for pre-relaxed slabs where atomic positions are not exactly on symmetry sites. |
| `--layer-tol` | auto | z threshold in Angstrom to separate anion layers. Adaptive by default. |
| `--face-band` | `0.5` | Half-width in Angstrom of the neutral central band along z. Anion sites within this distance of the slab midplane are assigned face `mid` and are not counted as top or bottom for `--mode same_face` / `both_faces`. |

### Refinement mode

| Argument | Default | Description |
|---|---|---|
| `--refine` | `None` | Path to a folder containing a `config_info.txt` whose vacancy set is to be refined. Setting this activates the refinement mode and makes the positional POSCAR argument optional (the pristine POSCAR path is read from `source_poscar` in `config_info.txt` by default). |
| `--balanced` | off | Physical constraint: each vacancy contributes exactly `--red-per-vac` reduced cations drawn from its own NN+NNN pool. Rejects cross-vacancy migrated distributions. |
| `--dedup-local` | off | Approximate local-symmetry deduplication under each vacancy's individual stabiliser `Stab(vac_i)` instead of (or on top of) the full-slab stabiliser `Stab(vac_set)`. Recovers reductions when `Stab(vac_set)` is trivial, at the cost of grouping configurations that are locally but not globally equivalent. Requires spglib. |
| `--pair-samples-per-type` | `None` | Strong approximation: after `--dedup-local` + `--balanced`, keep only N canonical pair representatives per pair-type class (NN-NN, NN-NNN, NNN-NNN) per vacancy. For N=1 the final count is bounded by `~3^n_vac`. |
| `--min-red-dist` | `0.0` | Minimum reduced-cation / reduced-cation distance in Angstrom. Distributions with any pair closer than this are discarded. `0` disables the filter. |

### Run control

| Argument | Default | Description |
|---|---|---|
| `--no-proxy` | off | General mode only. Vacancy-only mode: write POSCARs with anions removed but no proxy substitution. Reduced-cation arguments (`--red-shell`, `--red-per-vac`, `--max-red-combos`, `--seed`) are ignored. One configuration per symmetry-unique vacancy set. |
| `--dry-run` | off | Count configurations without creating any `cfg_NNNNN/` folders or writing `summary.csv`. |
| `--force` | off | Proceed despite failed coordination checks. |
| `--out` | `configs_out` (general) or `<cfg_dir>_refined` (refinement) | Output directory. Created if it does not exist. |

## Notes

- Flags that are refinement-only (`--balanced`, `--dedup-local`, `--pair-samples-per-type`, `--min-red-dist`) are ignored in general mode and the script prints a warning if they are passed outside `--refine`.
- The `config_info.txt` format is stable across both modes: files produced by general mode are directly consumable by refinement mode, and files produced by refinement mode preserve the fields the downstream pipeline expects.
- When spglib is unavailable, general mode falls back to a distance/depth fingerprint for deduplication, and refinement mode disables `--dedup-local` (it still runs with `--balanced` and the global fingerprint fallback).
- Species defaults (`Ce` / `La` / `O`) are chosen for CeO2 but all configuration output remains in the shared format regardless of species; passing `--host-cation Ti --proxy-cation Sc --anion O`, for example, will produce consistent configuration files for a Ti-based oxide with Sc as the proxy.
>>>>>>> 1b08c99 (Initial commit)
