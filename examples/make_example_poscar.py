#!/usr/bin/env python3
"""Generate a minimal CeO2 (110) slab POSCAR for the examples/ folder.

Fluorite CeO2 (a = 5.411 A), p(1x1), 3 CeO2 layers stacked along [110].
Each (110) layer is itself stoichiometric CeO2 (2 Ce + 4 O per layer),
so the slab is dipole-free and well-terminated without needing manual
atom removal. Total: Ce6 O12, >20 A vacuum along z.

This is a didactic slab for examples/ only; it is NOT the slab used in
any production study.
"""
import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import surface

# CeO2 fluorite conventional cubic cell
a = 5.411
Ce_frac = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
O_frac = np.array([[0.25, 0.25, 0.25], [0.75, 0.75, 0.25],
                   [0.75, 0.25, 0.75], [0.25, 0.75, 0.75],
                   [0.75, 0.75, 0.75], [0.25, 0.25, 0.75],
                   [0.25, 0.75, 0.25], [0.75, 0.25, 0.25]])
cell = a * np.eye(3)
positions_cart = np.vstack([Ce_frac, O_frac]) @ cell
symbols = ["Ce"] * 4 + ["O"] * 8
bulk = Atoms(symbols=symbols, positions=positions_cart, cell=cell, pbc=True)

slab = surface(bulk, (1, 1, 0), layers=3, vacuum=10.0)
slab.center(axis=2, vacuum=10.0)

syms = slab.get_chemical_symbols()
print(f"Composition: Ce={syms.count('Ce')} O={syms.count('O')}")

write("POSCAR_CeO2_110_demo", slab, format="vasp", direct=True, sort=True)
print(f"Wrote POSCAR_CeO2_110_demo: {slab.get_chemical_formula()}")
print(f"Cell c = {slab.cell[2,2]:.3f} A")
