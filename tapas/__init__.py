"""TAPAS: Topology-Aware Polaron Arrangement Sampler.

Symmetry-aware enumeration and refinement of anion-vacancy / reduced-cation
configurations in oxide slabs.

The command-line interface lives in `tapas.py` at the repository root. This
package exposes the programmatic API used by that script and reusable by
notebooks or external code:

    from tapas import run_general, run_refine
    from tapas import canonical_label, symmetry_permutations
    from tapas import write_config, parse_config_info
"""

from .common import VALID_LAYER_KEYS, die, warn, stable_hash
from .geometry import (
    validate_and_center,
    classify_anion_layers,
    suggest_radii,
    coordination_check,
    shells_for_vacancy,
)
from .symmetry import (
    symmetry_permutations,
    canonical_label,
    fingerprint_fallback,
    HAVE_SPGLIB,
)
from .io import (
    write_config,
    write_config_vacancy_only,
    write_config_info,
    parse_config_info,
)
from .general import run_general
from .refine import run_refine

__all__ = [
    "VALID_LAYER_KEYS",
    "die",
    "warn",
    "stable_hash",
    "validate_and_center",
    "classify_anion_layers",
    "suggest_radii",
    "coordination_check",
    "shells_for_vacancy",
    "symmetry_permutations",
    "canonical_label",
    "fingerprint_fallback",
    "HAVE_SPGLIB",
    "write_config",
    "write_config_vacancy_only",
    "write_config_info",
    "parse_config_info",
    "run_general",
    "run_refine",
]

__version__ = "0.1.0"
