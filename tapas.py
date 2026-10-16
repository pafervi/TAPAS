#!/usr/bin/env python3
"""Command-line entry point for TAPAS (Topology-Aware Polaron Arrangement Sampler).

Symmetry-aware enumeration of oxygen-vacancy / reduced-cation configurations in slabs.

Two operating modes:

(A) GENERAL ENUMERATION (default). Enumerate every symmetry-inequivalent vacancy
    set compatible with the user-selected anion layers and physical filters
    (--nvac, --mode, --vac-layers, --min-vac-dist, --max-bonds-lost), and for
    each unique vacancy set stratify-sample reduced-cation (polaron) distributions
    across the NN+NNN pool, optionally restricted to a specific shell choice
    (--red-shell). Vacancy sets are deduplicated exactly under the pristine slab
    space group (spglib); reduced-cation distributions are deduplicated with the
    same canonical label applied to the full (vac, red) config.

(B) REFINEMENT (--refine CFG_DIR). Read a config_info.txt produced by a previous
    general run, fix the vacancy set, and exhaustively enumerate all reduced-cation
    distributions compatible with that vacancy set (no sampling, no cap). Three
    optional reductions can be combined on top of the global canonical-label dedup:
        --balanced                  enforce the physical constraint that each
                                    vacancy donates its 2 polarons to its own
                                    NN+NNN pool (rejects cross-vacancy migrated
                                    configs).
        --dedup-local               deduplicate under each vacancy's own
                                    stabiliser Stab(vac_i) instead of (or on top
                                    of) the full-slab stabiliser Stab(vac_set);
                                    local-symmetry approximation that recovers
                                    reductions when Stab(vac_set) is trivial.
        --pair-samples-per-type N   requires --dedup-local + --balanced; keep
                                    only N canonical pair representatives per
                                    pair-type class (NN-NN, NN-NNN, NNN-NNN) per
                                    vacancy. Strong approximation that assumes
                                    energetics are dominated by pair type.

Output in both modes: one subfolder per unique configuration containing POSCAR
(vacancies removed, reduced sites substituted as proxy cation) and config_info.txt
(provenance, vacancy indices with xyz, reduced-cation indices with shell labels
and vacancy distances, per-vacancy NN/NNN pools, proxy indices verified by
coordinate matching, post-DFT control notes). A top-level summary.csv with
provenance header lists every written configuration.
"""

import argparse

from tapas import run_general, run_refine
from tapas.common import die, warn


def parse_args():
    """Parse the command-line arguments for both operating modes."""
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )
    p.add_argument("poscar", nargs="?",
                   help="pristine slab POSCAR. Required in general mode; "
                        "in --refine mode optional (default: source_poscar from "
                        "config_info.txt).")

    # Species
    p.add_argument("--host-cation", default="Ce")
    p.add_argument("--proxy-cation", default="La")
    p.add_argument("--anion", default="O")

    # General-mode arguments
    p.add_argument("--nvac", type=int, default=1,
                   help="[general] number of vacancies per configuration")
    p.add_argument("--mode", choices=["same_face", "both_faces", "all"], default="all",
                   help="[general] vacancy face constraint")
    p.add_argument("--vac-layers", default="surf,subsurf",
                   help="[general] comma-separated anion layer groups: "
                        "surf, subsurf, inner")
    p.add_argument("--red-shell", choices=["nn", "nnn", "both", "nnn_pref", "nnn_strict"],
                   default="both",
                   help="[general] reduced-cation shell(s): nn only, nnn only, both, "
                        "nnn_pref (NN <= nvac), nnn_strict (0 NN)")
    p.add_argument("--red-per-vac", type=int, default=2,
                   help="[general] reduced cations produced per vacancy (default: 2)")
    p.add_argument("--max-red-combos", type=int, default=25,
                   help="[general] stratified sample budget per vacancy set")
    p.add_argument("--seed", type=int, default=42,
                   help="[general] RNG seed for stratified sampling")
    p.add_argument("--min-vac-dist", type=float, default=0.0,
                   help="[general] minimum vacancy-vacancy distance (A); 0 disables")
    p.add_argument("--max-bonds-lost", type=int, default=None,
                   help="[general] cap on vacancies sharing a NN host cation")
    p.add_argument("--face-band", type=float, default=0.5,
                   help="[general] half-width of the neutral mid-plane band (A)")
    p.add_argument("--no-proxy", action="store_true",
                   help="[general] vacancy-only mode: write POSCARs with anions removed, "
                        "no proxy substitution; reduced-cation arguments are ignored")

    # Shared geometry / symmetry
    p.add_argument("--r-nn", type=float, default=None,
                   help="first host-anion shell radius (A); auto if omitted")
    p.add_argument("--r-nnn", type=float, default=None,
                   help="second host-anion shell radius (A); auto if omitted")
    p.add_argument("--layer-tol", type=float, default=None,
                   help="z threshold to split anion layers; adaptive if omitted")
    p.add_argument("--host-coord", type=int, nargs="+", default=None,
                   help="expected host-cation coordination numbers for anion sites "
                        "(e.g. 3 4); skips coord check if omitted")
    p.add_argument("--symprec", type=float, default=1e-3,
                   help="spglib symmetry tolerance (raise to 0.01-0.05 for relaxed slabs)")

    # Refinement-mode arguments
    p.add_argument("--refine", default=None, metavar="CFG_DIR",
                   help="REFINEMENT mode: path to a folder containing a config_info.txt "
                        "whose vacancy set is to be refined. Enables --balanced, "
                        "--dedup-local, --pair-samples-per-type, --min-red-dist.")
    p.add_argument("--balanced", action="store_true",
                   help="[refine] each vacancy contributes exactly red-per-vac polarons "
                        "from its own NN+NNN pool (polaron-vacancy locality constraint)")
    p.add_argument("--dedup-local", action="store_true",
                   help="[refine, APPROXIMATION] deduplicate under Stab(vac_i) per vacancy "
                        "instead of (or on top of) Stab(vac_set). Valid when energetics "
                        "are dominated by local environments.")
    p.add_argument("--pair-samples-per-type", type=int, default=None,
                   help="[refine, STRONG APPROXIMATION] with --dedup-local + --balanced, "
                        "keep only N canonical pair representatives per pair type "
                        "(NN-NN, NN-NNN, NNN-NNN) per vacancy. N=1 gives ~3^n_vac configs.")
    p.add_argument("--min-red-dist", type=float, default=0.0,
                   help="[refine] minimum reduced-cation / reduced-cation distance (A); "
                        "0 disables")

    # Output / behaviour
    p.add_argument("--dry-run", action="store_true",
                   help="count unique configs only; do not write any files")
    p.add_argument("--force", action="store_true",
                   help="continue despite failed coordination checks")
    p.add_argument("--out", default=None,
                   help="output directory (general default: configs_out; "
                        "refine default: <cfg_dir>_refined)")

    args = p.parse_args()

    # Mode-specific validation
    if args.refine:
        if args.pair_samples_per_type is not None and not (args.dedup_local and args.balanced):
            die("--pair-samples-per-type requires --dedup-local AND --balanced")
    else:
        if args.poscar is None:
            die("POSCAR is required in general mode (not using --refine)")
        if args.nvac < 1:
            die("--nvac must be >= 1")
        if (args.balanced or args.dedup_local or args.pair_samples_per_type is not None
                or args.min_red_dist > 0):
            warn("--balanced / --dedup-local / --pair-samples-per-type / --min-red-dist "
                 "are refinement-only and will be ignored in general mode.")

    return args


def main():
    """Dispatch to the general or refinement entry point based on --refine."""
    args = parse_args()
    if args.refine:
        run_refine(args)
    else:
        run_general(args)


if __name__ == "__main__":
    main()
