"""Slab geometry: vacuum detection, centering, anion layer tagging,
radii auto-calibration, coordination checks and NN/NNN shell lookups."""

import numpy as np
from ase.geometry import get_distances

from .common import die, warn


def validate_and_center(atoms, host, proxy, anion, strict=True):
    """Validate the slab contains the expected species and re-center it so
    the vacuum gap lies at the cell boundary along z.

    Re-centering matters because layer detection and the top/bottom mirror
    rely on the slab being contiguous in z; a slab that wraps across the
    periodic boundary would be split into spurious layers.

    Parameters
    ----------
    atoms : ase.Atoms
    host, proxy, anion : str
        Chemical symbols for the host cation, proxy species and anion.
    strict : bool
        If True, abort when host or anion are missing. If False, tolerate
        either absence (used by the refinement path when the input POSCAR
        is pristine and we only need it to anchor the symmetry group and
        atomic positions).

    Returns
    -------
    atoms : ase.Atoms
        Re-centered copy of the input (the original is not modified).
    vacuum : float
        Detected vacuum thickness in Angstrom.
    """
    syms = set(atoms.get_chemical_symbols())
    if strict:
        if host not in syms or anion not in syms:
            die(f"The POSCAR must contain {host} and {anion}; it contains {sorted(syms)}")
    extra = syms - {host, anion, proxy}
    if extra:
        warn(f"unexpected species in the slab: {sorted(extra)} (ignored)")

    c = atoms.cell[2, 2]
    if (abs(atoms.cell[2, 0]) > 1e-6 or abs(atoms.cell[2, 1]) > 1e-6
            or abs(atoms.cell[0, 2]) > 1e-6 or abs(atoms.cell[1, 2]) > 1e-6):
        warn("the c axis is not orthogonal to the ab plane; assuming surface normal is ~z")

    # Largest cyclic gap in z = vacuum region.
    z = np.sort(atoms.positions[:, 2] % c)
    gaps = np.diff(np.append(z, z[0] + c))
    imax = int(gaps.argmax())
    vacuum = float(gaps[imax])
    if vacuum < 6.0:
        die(f"No clear vacuum gap detected along z (largest gap {vacuum:.2f} A < 6 A).")

    # Shift so the centre of the vacuum lands on the cell edge.
    gap_start = z[imax]
    shift = -(gap_start + vacuum / 2.0)
    atoms = atoms.copy()
    atoms.positions[:, 2] = (atoms.positions[:, 2] + shift) % c
    return atoms, vacuum


def classify_anion_layers(atoms, anion, layer_tol=None):
    """Group anion atoms into atomic layers along z.

    If `layer_tol` is None the threshold is chosen adaptively: anion atoms are
    sorted by z, consecutive separations are computed, and the threshold is
    placed at the largest jump in the sorted list of separations (i.e. the
    natural divide between "same plane" and "different plane" spacings). This
    avoids a brittle fixed tolerance.

    Layer tags encode position relative to BOTH surfaces:
      - outermost layers      -> surf_top / surf_bottom
      - second-outermost       -> subsurf_top / subsurf_bottom
      - everything in between  -> inner_{layer_index}

    Returns
    -------
    info : dict[int, tuple[int, str]]
        Maps each anion atom index to (layer_index, tag).
    layers : list[list[int]]
        List of layers, each a list of anion atom indices, ordered by
        increasing z.
    """
    o_idx = np.array([a.index for a in atoms if a.symbol == anion])
    zs = atoms.positions[o_idx, 2]
    order = np.argsort(zs)
    zs_sorted = zs[order]
    dz = np.diff(zs_sorted)

    if layer_tol is None:
        nz = np.sort(dz)
        if len(nz) == 0:
            die("Only one anion position found.")
        jumps = np.diff(nz)
        if nz[-1] < 0.15:
            die("All z separations between anion atoms are < 0.15 A.")
        k = int(jumps.argmax())
        layer_tol = float((nz[k] + nz[k + 1]) / 2.0)
        layer_tol = max(layer_tol, 0.15)  # floor to avoid splitting a buckled plane
    splits = np.where(dz > layer_tol)[0]

    layers, start = [], 0
    for s in splits:
        layers.append(o_idx[order[start:s + 1]].tolist())
        start = s + 1
    layers.append(o_idx[order[start:]].tolist())

    n = len(layers)
    info = {}
    for li, lay in enumerate(layers):
        if li == n - 1:
            tag = "surf_top"
        elif li == n - 2 and n >= 4:
            tag = "subsurf_top"
        elif li == 0:
            tag = "surf_bottom"
        elif li == 1 and n >= 4:
            tag = "subsurf_bottom"
        else:
            tag = f"inner_{li}"
        for idx in lay:
            info[idx] = (li, tag)

    return info, layers


def suggest_radii(atoms, host, anion):
    """Suggest r_nn and r_nnn from the real host-anion distance distribution.

    The first large gap after the first coordination shell defines r_nn;
    the next gap defines r_nnn. Returned values are only suggestions;
    explicit --r-nn/--r-nnn override them.
    """
    ce = [a.index for a in atoms if a.symbol == host]
    o = [a.index for a in atoms if a.symbol == anion]
    _, D = get_distances(atoms.positions[o], atoms.positions[ce],
                         cell=atoms.cell, pbc=atoms.pbc)
    d = np.sort(D[D < 6.5])
    gaps = np.diff(d)
    cand = np.where((d[:-1] > 2.0) & (gaps > 0.3))[0]
    r_nn = float(d[cand[0]] + gaps[cand[0]] / 2.0) if len(cand) else 2.7
    cand2 = np.where((d[:-1] > r_nn + 1.0) & (gaps > 0.2))[0]
    r_nnn = float(d[cand2[0]] + gaps[cand2[0]] / 2.0) if len(cand2) else r_nn + 2.2
    return round(r_nn, 3), round(r_nnn, 3)


def coordination_check(atoms, allowed_o, r_nn, host, allowed_coord, force=False):
    """Sanity-check that each candidate anion has a host coordination in the
    allowed set (e.g. {3, 4} for fluorite). Aborts on failure unless force.

    Returns
    -------
    coord : dict[int, int]
        anion index -> number of host cations within r_nn.
    """
    ce_idx = np.array([a.index for a in atoms if a.symbol == host])
    bad = []
    coord = {}
    for o in allowed_o:
        d = atoms.get_distances(o, ce_idx, mic=True)
        nnn = int((d < r_nn).sum())
        coord[o] = nnn
        if nnn not in allowed_coord:
            bad.append((o, nnn))
    if bad:
        msg = (f"{len(bad)} anion site(s) have an anomalous host coordination "
               f"(expected {sorted(allowed_coord)}).")
        if force:
            warn(msg + " Continuing anyway.")
        else:
            r1, r2 = suggest_radii(atoms, host, atoms[allowed_o[0]].symbol)
            die(msg + f"\nSuggested: --r-nn {r1} --r-nnn {r2}.")
    return coord


def shells_for_vacancy(atoms, o_index, r_nn, r_nnn, host):
    """Return (nn, nnn) host-atom indices for the anion site at `o_index`,
    using minimum-image distances. Host atoms beyond r_nnn are excluded by
    construction (never reduced-cation candidates)."""
    ce_idx = np.array([a.index for a in atoms if a.symbol == host])
    d = atoms.get_distances(o_index, ce_idx, mic=True)
    nn = ce_idx[d < r_nn].tolist()
    nnn = ce_idx[(d >= r_nn) & (d < r_nnn)].tolist()
    return nn, nnn
