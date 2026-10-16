"""Refinement mode.

Read a config_info.txt produced by a previous general run, fix the
vacancy set, and exhaustively enumerate all reduced-cation distributions
compatible with it. Three optional reductions (physical constraint,
local-symmetry deduplication, pair-type sampling) can be combined on
top of the global canonical-label dedup.
"""

import datetime
import itertools
import math
import os
import sys
import time

import numpy as np
from ase.io import read
from ase.geometry import get_distances

try:
    import spglib
    HAVE_SPGLIB = True
except ImportError:
    HAVE_SPGLIB = False

from .common import die, warn, stable_hash
from .geometry import validate_and_center
from .symmetry import symmetry_permutations, canonical_label
from .io import write_config, write_config_info, parse_config_info


def run_refine(args):
    """Exhaustively enumerate reduced-cation distributions for a fixed
    vacancy set read from a config_info.txt, with optional physical and
    local-symmetry reductions."""
    cfg_info_path = os.path.join(args.refine, "config_info.txt")
    if not os.path.isfile(cfg_info_path):
        die(f"Not found: {cfg_info_path}")

    source_poscar, n_vac, vac_indices, pools = parse_config_info(cfg_info_path, args.anion)

    poscar_path = args.poscar or source_poscar
    if poscar_path is None:
        die("source_poscar not in config_info.txt; provide POSCAR as positional argument")
    if not os.path.isfile(poscar_path):
        die(f"POSCAR not found: {poscar_path}\n"
            "  If running on a different machine, pass the pristine POSCAR as positional arg")

    atoms = read(poscar_path)
    atoms, vacuum = validate_and_center(atoms, args.host_cation, args.proxy_cation,
                                        args.anion, strict=False)
    print(f"  Slab: {atoms.get_chemical_formula()}  |  vacuum: {vacuum:.2f} A")
    print(f"  Fixed vacancies ({n_vac}): {args.anion} indices {vac_indices}")

    # Union pool across all vacancies, NN taking priority
    union_nn, union_nnn = set(), set()
    for v in vac_indices:
        nn, nnn = pools[v]
        union_nn |= nn
        union_nnn |= nnn
    union_nnn -= union_nn
    pool = sorted(union_nn | union_nnn)

    n_red = args.red_per_vac * n_vac
    total_raw = math.comb(len(pool), n_red)
    print(f"  Reduced-cation pool: {len(pool)} sites  "
          f"(NN={len(union_nn)}, NNN_excl={len(union_nnn)})")
    print(f"  Exhaustive combinations: C({len(pool)},{n_red}) = {total_raw:,}")

    if len(pool) < n_red:
        die(f"Pool size ({len(pool)}) < required reduced cations ({n_red})")

    perms = symmetry_permutations(atoms, args.symprec)
    dedup = "exact_spglib" if perms else "none"
    print(f"  Deduplication: {dedup}")

    # Symmetry diagnostic: |Stab(vac_set)| controls how much global dedup is
    # mathematically possible for this canonical vacancy set.
    if perms:
        vac_set = frozenset(vac_indices)
        stab_perms = [p for p in perms
                      if frozenset(int(p[v]) for v in vac_indices) == vac_set]
        n_stab = len(stab_perms)
        print(f"  |G| = {len(perms)}  |  |Stab(vac)| = {n_stab}  "
              f"|  orbit vac = {len(perms) // n_stab if n_stab else '?'}")

        pool_set = set(pool)
        remaining = set(pool)
        pool_orbits = []
        while remaining:
            seed = min(remaining)
            orb = {seed}
            grew = True
            while grew:
                grew = False
                for x in list(orb):
                    for p in stab_perms:
                        img = int(p[x])
                        if img in pool_set and img not in orb:
                            orb.add(img)
                            grew = True
            pool_orbits.append(orb)
            remaining -= orb
        n_singleton = sum(1 for o in pool_orbits if len(o) == 1)
        max_orb = max(len(o) for o in pool_orbits) if pool_orbits else 0
        print(f"  Pool orbits under Stab(vac): {len(pool_orbits)} "
              f"(singletons={n_singleton}, max orbit={max_orb})")
        if n_stab == 1:
            warn("Stab(vac_set) trivial: no global exact dedup available. "
                 "Options: --balanced (physical constraint) or "
                 "--balanced --dedup-local (local-symmetry approximation).")

    out_dir = args.out or (args.refine.rstrip("/\\") + "_refined")
    if not args.dry_run:
        os.makedirs(out_dir, exist_ok=True)

    vac_combo = tuple(vac_indices)
    pos = atoms.positions

    provenance = [
        f"# generated: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"# command: {' '.join(sys.argv)}",
        f"# ase={__import__('ase').__version__} "
        f"spglib={spglib.__version__ if HAVE_SPGLIB else 'NO'}",
        f"# dedup={dedup} symprec={args.symprec}",
        f"# refinement_source={os.path.abspath(args.refine)}",
    ]

    summary_rows = provenance + [
        "config_id,red_sites,red_shell_labels,pair_type,dedup_label_hash"
    ]

    per_vac_pool = {v: (pools[v][0] | pools[v][1]) for v in vac_indices}

    # Validate flag compatibility
    if args.dedup_local and not perms:
        die("--dedup-local requires spglib symmetry operations")

    canonical_pairs_per_vac = None
    stab_per_vac = None
    if args.dedup_local:
        stab_per_vac = {v: [p for p in perms if int(p[v]) == v] for v in vac_indices}

    # --- Precompute per-mode structures ---
    if args.dedup_local and args.balanced:
        print("\n  [APPROX] --dedup-local + --balanced: canonical pairs of "
              f"{args.red_per_vac} reduced sites per vacancy under Stab(vac_i), "
              "cartesian product.\n")
        canonical_pairs_per_vac = []
        raw_local_product = 1
        total_local_reduction = 1
        for v in vac_indices:
            pool_v = sorted(per_vac_pool[v])
            stab_v = stab_per_vac[v]
            seen_local = {}
            for pair in itertools.combinations(pool_v, args.red_per_vac):
                best = None
                for pm in stab_v:
                    trans = tuple(sorted(int(pm[c]) for c in pair))
                    if best is None or trans < best:
                        best = trans
                if best not in seen_local:
                    seen_local[best] = pair
            pairs_v = list(seen_local.values())
            raw_pairs = math.comb(len(pool_v), args.red_per_vac)
            print(f"  vac {v:4d}: pool={len(pool_v)}, |Stab(vac_i)|={len(stab_v)}, "
                  f"canonical pairs={len(pairs_v)}/{raw_pairs} "
                  f"(reduction {100*(1-len(pairs_v)/raw_pairs):.0f}%)")
            canonical_pairs_per_vac.append(pairs_v)
            raw_local_product *= raw_pairs
            total_local_reduction *= len(pairs_v)
        print(f"  Local canonical product: {total_local_reduction:,} "
              f"(from {raw_local_product:,} raw, "
              f"reduction {100*(1-total_local_reduction/raw_local_product):.1f}%)")

        if args.pair_samples_per_type is not None:
            N = args.pair_samples_per_type
            print(f"\n  [STRONG APPROX] --pair-samples-per-type={N}: keeping only N "
                  "representatives per pair type (NN-NN, NN-NNN, NNN-NNN) per vacancy.\n")
            product_after = 1
            for i, v in enumerate(vac_indices):
                nn_set = pools[v][0]
                by_type = {}
                for pair in canonical_pairs_per_vac[i]:
                    type_key = tuple(sorted(("NN" if c in nn_set else "NNN") for c in pair))
                    by_type.setdefault(type_key, []).append(pair)
                trimmed = []
                summary = []
                for tk in sorted(by_type.keys()):
                    ps = by_type[tk]
                    kept = ps[:N]
                    trimmed.extend(kept)
                    summary.append(f"{'-'.join(tk)}: {len(kept)}/{len(ps)}")
                print(f"  vac {v:4d}: {', '.join(summary)}  ->  "
                      f"{len(trimmed)} pairs (from {len(canonical_pairs_per_vac[i])})")
                canonical_pairs_per_vac[i] = trimmed
                product_after *= len(trimmed)
            print(f"  Product after pair-type sampling: {product_after:,}")

    elif args.dedup_local:
        print("\n  [APPROX] --dedup-local without --balanced: dedup by local pattern "
              "(red ∩ pool_v) canonized under Stab(vac_i), combined per vacancy.\n")
        for v in vac_indices:
            print(f"  vac {v:4d}: pool={len(per_vac_pool[v])}, "
                  f"|Stab(vac_i)|={len(stab_per_vac[v])}")

    def local_config_label(red):
        """Per-vacancy canonical local pattern tuple, used when --dedup-local
        is active without --balanced."""
        red_set = set(red)
        per_vac_labels = []
        for v in vac_indices:
            local = tuple(sorted(red_set & per_vac_pool[v]))
            best = None
            for pm in stab_per_vac[v]:
                trans = tuple(sorted(int(pm[c]) for c in local))
                if best is None or trans < best:
                    best = trans
            per_vac_labels.append(best if best is not None else local)
        return tuple(per_vac_labels)

    def combos_iter():
        """Iterate reduced-cation combos according to the active modes."""
        if args.dedup_local and args.balanced:
            seen_sets = set()
            for chunk in itertools.product(*canonical_pairs_per_vac):
                flat = tuple(sorted(set().union(*[set(c) for c in chunk])))
                if len(flat) != n_red:
                    continue
                if flat in seen_sets:
                    continue
                seen_sets.add(flat)
                yield flat
            return
        if args.dedup_local:
            seen_labels = set()
            for red in itertools.combinations(pool, n_red):
                lab = local_config_label(red)
                if lab in seen_labels:
                    continue
                seen_labels.add(lab)
                yield red
            return
        if args.balanced:
            seen_sets = set()
            pools_sorted = [sorted(per_vac_pool[v]) for v in vac_indices]
            for chunk in itertools.product(
                *[itertools.combinations(p, args.red_per_vac) for p in pools_sorted]
            ):
                flat = tuple(sorted(set().union(*[set(c) for c in chunk])))
                if len(flat) != n_red:
                    continue
                if flat in seen_sets:
                    continue
                seen_sets.add(flat)
                yield flat
            return
        yield from itertools.combinations(pool, n_red)

    # Count total for progress reporting
    if args.balanced:
        total_raw = sum(1 for _ in combos_iter())
        mode_name = "--dedup-local + --balanced" if args.dedup_local else "--balanced"
        print(f"  {mode_name}: {total_raw:,} candidate configurations")
    elif args.dedup_local:
        print(f"  --dedup-local: iterating C({len(pool)},{n_red})={total_raw:,} combos "
              f"and deduplicating by local label on the fly")

    # Optional reduced-cation / reduced-cation distance filter
    if args.min_red_dist > 0:
        pool_arr = np.array(pool)
        _, D_pool = get_distances(pos[pool_arr], pos[pool_arr],
                                  cell=atoms.cell, pbc=atoms.pbc)
        pool_idx_of = {p_atom: i for i, p_atom in enumerate(pool)}

    def passes_dist_filter(red):
        if args.min_red_dist <= 0:
            return True
        idxs = [pool_idx_of[c] for c in red]
        for i in range(len(idxs)):
            for j in range(i + 1, len(idxs)):
                if D_pool[idxs[i], idxs[j]] < args.min_red_dist:
                    return False
        return True

    # Build a mock shell_pools dict for config_info output (uses NN/NNN as read)
    shell_pools = {v: (pools[v][0], pools[v][1]) for v in vac_indices}

    seen = set()
    cid = 0
    n_filtered_dist = 0
    t0 = time.time()
    report_every = max(1, total_raw // 400)

    for i_raw, red in enumerate(combos_iter()):
        if i_raw % report_every == 0 or i_raw == total_raw - 1:
            elapsed = time.time() - t0
            rate = (i_raw + 1) / elapsed if elapsed > 0 else 0
            eta = (total_raw - i_raw - 1) / rate if rate > 0 else 0
            print(f"\r  {i_raw+1:,}/{total_raw:,} "
                  f"({100*(i_raw+1)/total_raw:.1f}%)  unique={cid}  "
                  f"{rate:.0f}/s  ETA {eta:.0f}s   ",
                  end="", flush=True)

        if not passes_dist_filter(red):
            n_filtered_dist += 1
            continue

        if perms:
            lab = canonical_label(vac_combo, red, perms)
        else:
            lab = (tuple(sorted(vac_combo)), tuple(sorted(red)))
        h = stable_hash(lab)
        if h in seen:
            continue
        seen.add(h)

        labels = ["NN" if c in union_nn else "NNN" for c in red]
        pair = "-".join(sorted(labels))
        name = f"cfg_{cid:05d}"

        if not args.dry_run:
            cfg_dir = os.path.join(out_dir, name)
            proxy_written, counts, cons = write_config(
                atoms, list(vac_combo), list(red), cfg_dir,
                args.host_cation, args.proxy_cation, args.anion,
            )
            write_config_info(
                os.path.join(cfg_dir, "config_info.txt"),
                atoms, provenance, name,
                os.path.abspath(poscar_path), n_vac,
                mode_tag="refinement_exhaustive", pair_type=pair, dedup_hash=h,
                counts=counts, host=args.host_cation,
                proxy=args.proxy_cation, anion=args.anion,
                vac_combo=vac_combo, info_layers=None, face_fn=None,
                red_sites=red, red_labels=labels,
                shell_pools=shell_pools, proxy_written=proxy_written,
                constraints=cons,
            )

        summary_rows.append(
            f"{name},{'|'.join(map(str, red))},{'|'.join(labels)},{pair},{h}"
        )
        cid += 1

    if not args.dry_run:
        with open(os.path.join(out_dir, "summary.csv"), "w") as f:
            f.write("\n".join(summary_rows))

    removed = total_raw - cid - n_filtered_dist
    print(f"\n\nUnique configurations: {cid}"
          + ("  [dry-run, no folders written]" if args.dry_run
             else f"\n-> {out_dir}/"))
    print(f"From {total_raw:,} raw combinations: "
          f"{n_filtered_dist:,} filtered by --min-red-dist, "
          f"{removed:,} removed by symmetry/local dedup "
          f"({100*(removed+n_filtered_dist)/max(total_raw,1):.1f}% total reduction)")
