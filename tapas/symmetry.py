"""Space-group symmetry: build atom-index permutations from spglib
operations, canonical labels for exact deduplication, and a distance /
depth fingerprint fallback for the no-spglib case."""

import numpy as np
from ase.geometry import get_distances

try:
    import spglib
    HAVE_SPGLIB = True
except ImportError:
    HAVE_SPGLIB = False


def symmetry_permutations(atoms, symprec):
    """Build the index permutations induced by the pristine slab's symmetry
    operations.

    Each spglib symmetry operation (rotation R, translation t) is applied to
    the fractional coordinates and matched back to the original atoms
    (minimum image) to produce a permutation of atom indices. The set of
    permutations includes lateral supercell translations and, if the slab is
    mirror-symmetric, the top <-> bottom flip. These permutations are the
    group used to deduplicate configurations exactly.

    Returns a list of integer permutation arrays, or None if spglib is
    missing or finds nothing usable (caller then falls back to a distance
    fingerprint).
    """
    if not HAVE_SPGLIB:
        return None
    numbers = atoms.get_atomic_numbers()
    frac = atoms.get_scaled_positions(wrap=True)
    cell = (atoms.cell.array, frac, numbers)
    sym = spglib.get_symmetry(cell, symprec=symprec)
    if sym is None:
        return None
    perms = []
    L = atoms.cell.array
    for R, t in zip(sym["rotations"], sym["translations"]):
        new = (frac @ R.T + t) % 1.0
        # Match each transformed position to an original atom (minimum image).
        perm = np.full(len(atoms), -1, dtype=int)
        ok = True
        for i in range(len(atoms)):
            delta = frac - new[i]
            delta -= np.round(delta)
            dcart = np.linalg.norm(delta @ L, axis=1)
            j = int(dcart.argmin())
            # reject if no atom matches closely or species differ
            if dcart[j] > max(symprec * 10, 1e-3) or numbers[j] != numbers[i]:
                ok = False
                break
            perm[i] = j
        # keep only true bijections
        if ok and len(set(perm.tolist())) == len(atoms):
            perms.append(perm)
    return perms if perms else None


def canonical_label(vac, red, perms):
    """Canonical key of a configuration: the lex-min over all symmetry
    operations of (sorted vacancy indices, sorted reduced-cation indices).

    Two configurations related by any host symmetry operation map to the
    same canonical label, so hashing this label deduplicates them exactly.
    """
    best = None
    for p in perms:
        lab = (tuple(sorted(int(p[v]) for v in vac)),
               tuple(sorted(int(p[c]) for c in red)))
        if best is None or lab < best:
            best = lab
    return best


def fingerprint_fallback(atoms, vac_sites, red_sites, ndec=2):
    """Distance + depth fingerprint, used ONLY when spglib is unavailable
    or found no usable operations.

    Builds a rounded, sorted tuple of intra/inter-group distances
    (vac-vac, vac-red, red-red) plus the (min, max) depth of each site
    from the two surfaces. This is an approximate invariant; exact spglib
    deduplication is strongly preferred.
    """
    pos, cell, pbc = atoms.positions, atoms.cell, atoms.pbc
    z_top, z_bot = pos[:, 2].max(), pos[:, 2].min()
    groups = []
    for A, B in [(vac_sites, vac_sites), (vac_sites, red_sites), (red_sites, red_sites)]:
        if not A or not B:
            groups.append(())
            continue
        _, D = get_distances(pos[list(A)], pos[list(B)], cell=cell, pbc=pbc)
        vals = (
            sorted(round(float(x), ndec) for x in np.unique(np.triu(D, 1)) if x > 1e-6)
            if A is B
            else sorted(round(float(x), ndec) for x in D.flatten())
        )
        groups.append(tuple(vals))
    for sites in (vac_sites, red_sites):
        depths = sorted(
            (round(min(z_top - pos[s, 2], pos[s, 2] - z_bot), ndec),
             round(max(z_top - pos[s, 2], pos[s, 2] - z_bot), ndec)) for s in sites
        )
        groups.append(tuple(depths))
    return tuple(groups)
