"""Common helpers and constants shared across the package."""

import hashlib
import sys


# Layer-group keywords accepted by --vac-layers. Anion layers are tagged
# surf_{top,bottom}, subsurf_{top,bottom}, inner_{N}; these keys select them.
VALID_LAYER_KEYS = {"surf", "subsurf", "inner"}


def die(msg):
    """Print an error and abort. Used for unrecoverable input or consistency errors."""
    print(f"\n[ERROR] {msg}")
    sys.exit(1)


def warn(msg):
    """Print a non-fatal warning (indented to stand out in the run log)."""
    print(f"  [WARNING] {msg}")


def stable_hash(obj):
    """Deterministic 16-hex-char sha1 of repr(obj).

    Used as the canonical-configuration key. Unlike Python's built-in hash(),
    sha1 is stable across interpreter sessions, so canonical hashes and
    --seed reproducibility hold between separate runs.
    """
    return hashlib.sha1(repr(obj).encode()).hexdigest()[:16]
