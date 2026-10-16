"""Verified POSCAR writers and config_info.txt reader/writer.

POSCARs are written with ASE's sorted output and then re-read to:
    - verify composition matches the expectation,
    - map proxy indices in the written file back to the original positions
      by coordinate matching (ASE's sort would otherwise mis-label sites).

config_info.txt is both an output artefact of this script and the input
contract for the refinement mode. The writer here produces files that are
directly consumable by `parse_config_info`.
"""

import ast
import os
import re

from ase.io import read, write
from ase.geometry import get_distances

from .common import die


# ---------------------------------------------------------------------------
# POSCAR writers
# ---------------------------------------------------------------------------

def write_config_vacancy_only(atoms, combo, cfg_dir, anion):
    """Vacancy-only: write POSCAR with the selected anions removed, no proxy
    substitution. Used by --no-proxy in general mode."""
    new = atoms.copy()
    del new[[v for v in sorted(combo, reverse=True)]]
    os.makedirs(cfg_dir, exist_ok=True)
    poscar = os.path.join(cfg_dir, "POSCAR")
    write(poscar, new, format="vasp", sort=True, direct=True)
    chk = read(poscar)
    n_anion_exp = sum(1 for a in atoms if a.symbol == anion) - len(combo)
    if chk.get_chemical_symbols().count(anion) != n_anion_exp:
        die(f"Incorrect composition in {poscar}")


def write_config(atoms, combo, red, cfg_dir, host, proxy, anion):
    """Write POSCAR with vacancies removed and reduced cations substituted as
    proxy. Verifies composition and proxy indices by coordinate round-trip
    (ASE sorts on write and would otherwise mis-label the proxy sites).

    Returns
    -------
    proxy_written : list[int]
        Indices of the proxy atoms in the WRITTEN POSCAR (base 0).
    counts : dict
        Element counts in the written file (host, proxy, anion).
    (n_con_in, n_con_out) : tuple[int, int]
        Number of selective-dynamics constrained atoms before/after writing.
    """
    new = atoms.copy()
    for c in red:
        new[c].symbol = proxy
    red_pos = atoms.positions[list(red)].copy()
    del new[[v for v in sorted(combo, reverse=True)]]

    os.makedirs(cfg_dir, exist_ok=True)
    poscar = os.path.join(cfg_dir, "POSCAR")
    write(poscar, new, format="vasp", sort=True, direct=True)

    chk = read(poscar)
    n_host_exp = sum(1 for a in atoms if a.symbol == host) - len(red)
    n_anion_exp = sum(1 for a in atoms if a.symbol == anion) - len(combo)
    counts = {s: chk.get_chemical_symbols().count(s) for s in (host, proxy, anion)}
    if (counts[proxy] != len(red) or counts[host] != n_host_exp
            or counts[anion] != n_anion_exp):
        die(f"Incorrect composition in {poscar}: {counts}; "
            f"expected {host}={n_host_exp} {proxy}={len(red)} {anion}={n_anion_exp}")

    proxy_written = []
    for p in red_pos:
        _, D = get_distances(p.reshape(1, 3), chk.positions, cell=chk.cell, pbc=chk.pbc)
        j = int(D.argmin())
        if D[0, j] > 1e-4 or chk[j].symbol != proxy:
            die(f"Failed to map proxy cation indices in {poscar}")
        proxy_written.append(j)

    n_con_in = sum(len(getattr(c, "index", [])) for c in atoms.constraints)
    n_con_out = sum(len(getattr(c, "index", [])) for c in chk.constraints)
    return sorted(proxy_written), counts, (n_con_in, n_con_out)


# ---------------------------------------------------------------------------
# config_info.txt writer
# ---------------------------------------------------------------------------

def write_config_info(path, atoms, provenance, name, source_poscar, nvac,
                      mode_tag, pair_type, dedup_hash, counts, host, proxy, anion,
                      vac_combo, info_layers, face_fn, red_sites, red_labels,
                      shell_pools, proxy_written, constraints=None):
    """Write a config_info.txt in the format consumed by `parse_config_info`
    and the downstream pipeline. All atom indices are 0-based on the ORIGINAL
    centered POSCAR; proxy indices in the WRITTEN POSCAR are reported both in
    base 0 and base 1 (VASP) convention."""
    pos = atoms.positions
    lines = ["# === CONFIG INFO (control) ==="] + provenance + [""]
    lines += [
        f"config_id        : {name}",
        f"source_poscar    : {source_poscar}",
        f"n_vacancies      : {nvac}",
        f"mode             : {mode_tag}  |  pair_type: {pair_type}",
        f"dedup_label_hash : {dedup_hash}",
        f"POSCAR composition: {host}={counts[host]} {proxy}={counts[proxy]} {anion}={counts[anion]}",
    ]
    if constraints is not None and (constraints[0] or constraints[1]):
        lines.append(f"constraints      : in={constraints[0]} out={constraints[1]} fixed atoms")
    lines.append("")

    lines.append(f"# VACANCIES (indices in the ORIGINAL centered POSCAR, base 0)")
    for v in vac_combo:
        x, y, z = pos[v]
        if info_layers is not None and v in info_layers:
            li, tag = info_layers[v]
            face_s = face_fn(v) if face_fn else "?"
            nn_v = sorted(shell_pools[v][0]) if shell_pools and v in shell_pools else []
            lines.append(
                f"  vac {anion} idx={v:4d}  layer={li}  tag={tag}  face={face_s}  "
                f"coord_NN={len(nn_v)}  xyz=({x:8.4f},{y:8.4f},{z:8.4f})"
            )
        else:
            lines.append(f"  vac {anion} idx={v:4d}  xyz=({x:8.4f},{y:8.4f},{z:8.4f})")

    lines += ["", f"# {host}3+ (= {proxy} in the generated POSCAR)"]
    for c, lbl in zip(red_sites, red_labels):
        x, y, z = pos[c]
        dstr = ", ".join(
            f"d(vac{v})={float(atoms.get_distance(c, v, mic=True)):.3f}A" for v in vac_combo
        )
        lines.append(
            f"  {host}3+ orig_idx={c:4d}  shell={lbl:4s}  "
            f"xyz=({x:8.4f},{y:8.4f},{z:8.4f})  {dstr}"
        )

    if shell_pools:
        lines += ["", "# POOLS per vacancy (original indices)"]
        for v in vac_combo:
            nn_v, nnn_v = shell_pools[v]
            lines.append(f"  vac {v}: NN={sorted(nn_v)}  NNN={sorted(nnn_v)}")

    lines += [
        "",
        f"# {proxy} INDICES IN THE WRITTEN POSCAR (verified by coordinate matching)",
        f"  base 0: {proxy_written}",
        f"  base 1 (VASP): {[i + 1 for i in proxy_written]}",
        "",
        "# POST-DFT CONTROL",
        f"# 1. Resubstitute {proxy}->{host} keeping positions (la_to_ce / prep_vasp).",
        "# 2. After VASP: magnetization ~1 muB on the target sites, ~0 elsewhere.",
        "# 3. If any polaron migrated, reclassify NN/NNN using the distances above.",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


# ---------------------------------------------------------------------------
# config_info.txt parser (refinement mode input)
# ---------------------------------------------------------------------------

def parse_config_info(path, anion):
    """Parse a config_info.txt produced by this package.

    Extracts: source_poscar (optional), n_vacancies, vacancy indices
    (anion atom indices in the original POSCAR, 0-indexed), and per-vacancy
    NN/NNN pools. Vacancy indices are detected from lines of the form
    `vac {anion} idx=  NN`; pools from `vac N: NN=[...]  NNN=[...]`.
    """
    text = open(path).read()

    m = re.search(r"source_poscar\s*:\s*(\S+)", text)
    source_poscar = m.group(1) if m else None

    m = re.search(r"n_vacancies\s*:\s*(\d+)", text)
    if not m:
        die("n_vacancies not found in config_info.txt")
    n_vac = int(m.group(1))

    vac_indices = [int(m.group(1))
                   for m in re.finditer(rf"vac\s+{re.escape(anion)}\s+idx=\s*(\d+)", text)]
    if not vac_indices:
        # Backward compatibility: also try without the explicit species label
        vac_indices = [int(m.group(1)) for m in re.finditer(r"vac\s+\S+\s+idx=\s*(\d+)", text)]
    if len(vac_indices) != n_vac:
        die(f"Expected {n_vac} vacancy indices, found {len(vac_indices)}")

    pools = {}
    for m in re.finditer(r"vac (\d+):\s*NN=(\[.*?\])\s+NNN=(\[.*?\])", text):
        v = int(m.group(1))
        pools[v] = (set(ast.literal_eval(m.group(2))),
                    set(ast.literal_eval(m.group(3))))
    if not pools:
        die("No NN/NNN pools found in config_info.txt")

    missing = [v for v in vac_indices if v not in pools]
    if missing:
        die(f"Pools missing for vacancies: {missing}")

    return source_poscar, n_vac, vac_indices, pools
