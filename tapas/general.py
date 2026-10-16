"""General enumeration mode.

Enumerate every symmetry-inequivalent vacancy set compatible with the
user-selected anion layers and physical filters, and for each unique
vacancy set stratify-sample reduced-cation (polaron) distributions
across the NN+NNN pool. Vacancy sets are deduplicated exactly under
the pristine slab space group; reduced-cation distributions are
deduplicated with the same canonical label applied to the full
(vacancy, reduced) config.
"""

import datetime
import hashlib
import itertools
import math
import os
import random
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

from .common import die, stable_hash, VALID_LAYER_KEYS
from .geometry import (
    validate_and_center,
    suggest_radii,
    classify_anion_layers,
    coordination_check,
    shells_for_vacancy,
)
from .symmetry import symmetry_permutations, canonical_label, fingerprint_fallback
from .io import write_config, write_config_vacancy_only, write_config_info


def run_general(args):
    """Enumerate symmetry-inequivalent vacancy sets and stratify-sample
    reduced-cation distributions per set."""
    keys = [k.strip() for k in args.vac_layers.split(",") if k.strip()]
    bad_keys = [k for k in keys if k not in VALID_LAYER_KEYS]
    if bad_keys:
        die(f"Unrecognized --vac-layers keys: {bad_keys}. Valid: {sorted(VALID_LAYER_KEYS)}")

    atoms = read(args.poscar)
    atoms, vacuum = validate_and_center(atoms, args.host_cation, args.proxy_cation,
                                        args.anion, strict=True)

    n_host = sum(1 for a in atoms if a.symbol == args.host_cation)
    n_anion = sum(1 for a in atoms if a.symbol == args.anion)
    print(f" Slab: {args.host_cation}{n_host}{args.anion}{n_anion}  |  "
          f"detected vacuum: {vacuum:.2f} A")

    r1, r2 = suggest_radii(atoms, args.host_cation, args.anion)
    r_nn = args.r_nn if args.r_nn is not None else r1
    r_nnn = args.r_nnn if args.r_nnn is not None else r2
    if r_nn >= r_nnn:
        die(f"r_nn ({r_nn}) must be smaller than r_nnn ({r_nnn})")
    print(f" Radii: r_nn={r_nn} A, r_nnn={r_nnn} A "
          f"(auto; suggested from the slab: {r1}/{r2})")

    info, layers = classify_anion_layers(atoms, args.anion, args.layer_tol)

    def layer_allowed(tag):
        for k in keys:
            if k == "surf" and tag.startswith("surf"):
                return True
            if k == "subsurf" and tag.startswith("subsurf"):
                return True
            if k == "inner" and tag.startswith("inner"):
                return True
        return False

    allowed = sorted(idx for idx, (li, tag) in info.items() if layer_allowed(tag))
    if not allowed:
        die("No eligible anion site with those --vac-layers")

    if args.host_coord is not None:
        coordination_check(atoms, allowed, r_nn, args.host_cation,
                           allowed_coord=set(args.host_coord), force=args.force)
    else:
        print("  Coordination check skipped (--host-coord not set).")

    z = atoms.positions[:, 2]
    zmid = (z.max() + z.min()) / 2.0
    band = args.face_band

    def face(idx):
        dz = atoms.positions[idx, 2] - zmid
        return "top" if dz > band else ("bottom" if dz < -band else "mid")

    perms = symmetry_permutations(atoms, args.symprec)
    dedup_method = "exact_spglib_symmetry" if perms else "distance_depth_fingerprint"
    print(f"  Deduplication method: {dedup_method}")

    # 1. Candidate vacancy sets with physical filters
    vac_combos = []
    for combo in itertools.combinations(allowed, args.nvac):
        faces = {face(i) for i in combo}
        if args.mode == "same_face" and len(faces - {"mid"}) > 1:
            continue
        if args.mode == "both_faces" and not {"top", "bottom"} <= faces:
            continue
        if args.min_vac_dist > 0 and len(combo) > 1:
            _, D = get_distances(atoms.positions[list(combo)],
                                 atoms.positions[list(combo)],
                                 cell=atoms.cell, pbc=atoms.pbc)
            iu = np.triu_indices(len(combo), 1)
            if D[iu].min() < args.min_vac_dist:
                continue
        if args.max_bonds_lost is not None and len(combo) > 1:
            lost_count = {}
            bad = False
            for v in combo:
                nn_v, _ = shells_for_vacancy(atoms, v, r_nn, r_nnn, args.host_cation)
                for ce in nn_v:
                    lost_count[ce] = lost_count.get(ce, 0) + 1
                    if lost_count[ce] > args.max_bonds_lost:
                        bad = True
                        break
                if bad:
                    break
            if bad:
                continue
        vac_combos.append(combo)

    print(f" Eligible anion sites: {len(allowed)}  |  "
          f"vacancy combinations (filtered): {len(vac_combos)}")

    # 2. Dedup vacancy sets by host symmetry FIRST (before reduced-cation sampling)
    if perms:
        seen_vac, vac_reps = set(), []
        for combo in vac_combos:
            best = None
            for pm in perms:
                lab = tuple(sorted(int(pm[v]) for v in combo))
                if best is None or lab < best:
                    best = lab
            if best not in seen_vac:
                seen_vac.add(best)
                vac_reps.append(combo)
        n_before = len(vac_combos)
        vac_combos = vac_reps
        print(f"  Vacancy sets after symmetry dedup: {len(vac_combos)} "
              f"(from {n_before})")

    # 3. Provenance and output setup
    provenance = [
        f"# generated: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"# command: {' '.join(sys.argv)}",
        f"# ase={__import__('ase').__version__} "
        f"spglib={spglib.__version__ if HAVE_SPGLIB else 'NO'}",
        f"# dedup={dedup_method} symprec={args.symprec}",
        f"# r_nn={r_nn} r_nnn={r_nnn} layer_tol={'auto' if args.layer_tol is None else args.layer_tol}",
        f"# slab_centered=yes vacuum={vacuum:.2f}A",
    ]

    out_dir = args.out or "configs_out"
    os.makedirs(out_dir, exist_ok=True)
    if args.no_proxy:
        rows = provenance + ["config_id,vac_sites,vac_tags,vac_faces,vac_xyz"]
    else:
        rows = provenance + [
            "config_id,vac_sites,vac_tags,vac_faces,vac_xyz,proxy_sites,"
            "proxy_shell_labels,pair_type,dedup_label,truncated"
        ]
    seen, cid = set(), 0
    breakdown_layers = {}

    t_start = time.time()
    total_vc = len(vac_combos)
    report_every = max(1, total_vc // 200)

    for idx_vc, combo in enumerate(vac_combos):
        if idx_vc % report_every == 0 or idx_vc == total_vc - 1:
            elapsed = time.time() - t_start
            rate = (idx_vc + 1) / elapsed if elapsed > 0 else 0
            eta = int((total_vc - (idx_vc + 1)) / rate) if rate > 0 else 0
            pct = ((idx_vc + 1) / total_vc) * 100
            print(f"\r  combos {idx_vc + 1}/{total_vc} ({pct:.1f}%)  "
                  f"unique={cid}  {int(rate)} combo/s  ETA {eta}s   ",
                  end="", flush=True)

        tags = "|".join(info[v][1] for v in combo)
        faces_s = "|".join(face(v) for v in combo)
        vac_xyz = ";".join(
            f"{atoms.positions[v][0]:.3f}/{atoms.positions[v][1]:.3f}/"
            f"{atoms.positions[v][2]:.3f}" for v in combo
        )

        # No-proxy mode: one entry per unique vacancy set
        if args.no_proxy:
            if perms:
                lab = min(tuple(sorted(int(p[v]) for v in combo)) for p in perms)
            else:
                lab = fingerprint_fallback(atoms, list(combo), [])
            h = stable_hash(lab)
            if h in seen:
                continue
            seen.add(h)
            name = f"cfg_{cid:05d}"
            for v in combo:
                breakdown_layers[info[v][1]] = breakdown_layers.get(info[v][1], 0) + 1
            if not args.dry_run:
                cfg_dir = os.path.join(out_dir, name)
                write_config_vacancy_only(atoms, combo, cfg_dir, args.anion)
            rows.append(f"{name},{'|'.join(map(str, combo))},{tags},{faces_s},{vac_xyz}")
            cid += 1
            continue

        # Normal path: build NN/NNN pools per vacancy and sample reduced cations
        shell_map, union_nn, union_nnn = {}, set(), set()
        for v in combo:
            nn, nnn = shells_for_vacancy(atoms, v, r_nn, r_nnn, args.host_cation)
            shell_map[v] = (set(nn), set(nnn))
            union_nn |= set(nn)
            union_nnn |= set(nnn)
        union_nnn -= union_nn

        if args.red_shell == "nn":
            pool = sorted(union_nn)
        elif args.red_shell == "nnn":
            pool = sorted(union_nnn)
        else:
            pool = sorted(union_nn | union_nnn)

        n_red = args.red_per_vac * args.nvac
        if len(pool) < n_red:
            continue

        # Deterministic per-combo RNG for reproducible stratified sampling
        _combo_seed = int(hashlib.sha1(repr(tuple(sorted(combo))).encode())
                          .hexdigest()[:8], 16)
        _rng = random.Random(args.seed + _combo_seed)
        nn_pool = sorted(union_nn)
        nnn_pool = sorted(union_nnn)

        strata = []
        for n_nn in range(0, n_red + 1):
            n_nnn = n_red - n_nn
            if n_nn > len(nn_pool) or n_nnn > len(nnn_pool):
                continue
            if args.red_shell == "nn" and n_nnn > 0:
                continue
            if args.red_shell == "nnn" and n_nn > 0:
                continue
            if args.red_shell == "nnn_pref" and n_nn > args.nvac:
                continue
            if args.red_shell == "nnn_strict" and n_nn > 0:
                continue
            strata.append((n_nn, n_nnn))
        if not strata:
            continue

        per_stratum = max(1, args.max_red_combos // len(strata))

        def sample_stratum(n_nn, n_nnn, k):
            total = math.comb(len(nn_pool), n_nn) * math.comb(len(nnn_pool), n_nnn)
            out = set()
            if total <= k:
                for a in itertools.combinations(nn_pool, n_nn):
                    for b in itertools.combinations(nnn_pool, n_nnn):
                        out.add(tuple(sorted(a + b)))
                return out
            tries = 0
            while len(out) < k and tries < k * 20:
                a = tuple(_rng.sample(nn_pool, n_nn)) if n_nn else ()
                b = tuple(_rng.sample(nnn_pool, n_nnn)) if n_nnn else ()
                out.add(tuple(sorted(a + b)))
                tries += 1
            return out

        candidate_red = []
        for n_nn, n_nnn in strata:
            candidate_red.extend(sample_stratum(n_nn, n_nnn, per_stratum))

        for red in candidate_red:
            labels = ["NN" if c in union_nn else "NNN" for c in red]
            if perms:
                lab = canonical_label(combo, red, perms)
            else:
                lab = fingerprint_fallback(atoms, list(combo), list(red))
            h = stable_hash(lab)
            if h in seen:
                continue
            seen.add(h)

            name = f"cfg_{cid:05d}"
            pair = "-".join(sorted(labels))
            for v in combo:
                breakdown_layers[info[v][1]] = breakdown_layers.get(info[v][1], 0) + 1

            if not args.dry_run:
                cfg_dir = os.path.join(out_dir, name)
                proxy_written, counts, cons = write_config(
                    atoms, combo, red, cfg_dir,
                    args.host_cation, args.proxy_cation, args.anion
                )
                write_config_info(
                    os.path.join(cfg_dir, "config_info.txt"),
                    atoms, provenance, name,
                    os.path.abspath(args.poscar), args.nvac,
                    mode_tag=args.mode, pair_type=pair, dedup_hash=h,
                    counts=counts, host=args.host_cation,
                    proxy=args.proxy_cation, anion=args.anion,
                    vac_combo=combo, info_layers=info, face_fn=face,
                    red_sites=red, red_labels=labels,
                    shell_pools=shell_map, proxy_written=proxy_written,
                    constraints=cons,
                )

            rows.append(
                f"{name},{'|'.join(map(str, combo))},{tags},{faces_s},{vac_xyz},"
                f"{'|'.join(map(str, red))},{'|'.join(labels)},{pair},{h},no"
            )
            cid += 1

    with open(os.path.join(out_dir, "summary.csv"), "w") as f:
        f.write("\n".join(rows))

    print(f"\n\nUnique configurations: {cid}"
          + ("  [dry-run: no folders written]" if args.dry_run else ""))
    print(f"-> {out_dir}/summary.csv")
    if args.max_red_combos < 999999:
        print(f"  (Reduced sites SAMPLED: ~{args.max_red_combos} distributions "
              f"per vacancy set, stratified by pair_type, seed {args.seed}. "
              "Per-set minimum is from the sampled subset, not global.)")

    print("Breakdown by vacancy layer:")
    for l_name in sorted(breakdown_layers.keys()):
        print(f"  {l_name:<30} {breakdown_layers[l_name]}")
