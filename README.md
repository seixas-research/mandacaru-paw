# Mandacaru PAW-LCAO pseudopotentials

This repository holds the **pseudopotential datasets** used by
[Mandacaru](https://github.com/seixas-research/mandacaru), a framework for
simulating molecules with variational quantum algorithms (VQE, ADAPT-VQE) on
simulators and quantum hardware. It has one file per element for every
element from hydrogen to uranium (Z ≤ 92).

Mandacaru generated every dataset here from scratch, with its own
all-electron atomic solver. Nothing is converted from another code. The
datasets live in their own repository, not in the Mandacaru package, because
of their size.

## Why Mandacaru needs them

Mandacaru computes its Hamiltonian integrals on a real-space grid. A grid
fine enough for valence electrons (0.15–0.30 Å) cannot resolve the core
electrons, whose orbitals vary on a scale of hundredths of an ångström for
oxygen and heavier atoms. Removing the core and replacing the nuclear
potential by a smooth one solves this, and it also shrinks the problem: only
valence electrons are left to put on qubits.

**PAW-LCAO** is the projector augmented-wave (PAW) method of Blöchl in a form
built for a localized basis. The smooth pseudo-atomic orbitals of each
dataset are also the **basis functions** of the molecular calculation. PAW
keeps the information needed to reconstruct the true all-electron orbitals
near each nucleus, and adds the correction terms this reconstruction implies
to the molecular Hamiltonian.

## Contents

| Folder | Elements | What it is | Size |
|---|---|---|---|
| `lda-sr/` | 92 (H–U) | scalar-relativistic, LDA — **the default** | ~200 MB |
| `lda-fr/` | 91 (H–U, no Pa) | the same, plus spin-orbit coupling | ~250 MB |
| `pbe-sr/` | 92 (H–U) | scalar-relativistic, PBE | ~210 MB |

- **LDA** here is Slater exchange with Perdew–Zunger correlation, with the
  relativistic correction to exchange of MacDonald and Vosko.
- **Scalar-relativistic** means the reference atoms include mass-velocity
  and Darwin effects, but not spin-orbit coupling.
- **`lda-fr/`** adds spin-orbit coupling from the Dirac equation; see
  [FR.md](FR.md).
- **PBE** here is the functional of Perdew, Burke and Ernzerhof, with
  Perdew–Wang correlation and the same relativistic correction to exchange.
  [PBE.md](PBE.md) describes the set, its checks and its flagged elements.
  There is no PBE spin-orbit set.

## Quick start

Install Mandacaru, then clone this repository and tell Mandacaru where it is:

```bash
git clone https://github.com/seixas-research/mandacaru-paw
mandacaru --set-paw mandacaru-paw     # writes MANDACARU_PAW_PATH to ~/.zshrc or ~/.bashrc
# open a new terminal (or `source` that file), then check:
mandacaru --pseudo-status
```

PAW-LCAO is chosen as the **basis** of a calculation:

```python
from ase.build import molecule
from mandacaru import Mandacaru

atoms = molecule("H2O")
atoms.center(vacuum=4.0)          # the cell is the real-space box
atoms.calc = Mandacaru(method="adapt-vqe",
                       basis={"name": "PAW-LCAO", "size": "DZP"},
                       grid_spacing=0.25)    # Angstrom
energy = atoms.get_potential_energy()
```

The same from the command line:

```bash
mandacaru H2O --cell 8 --basis PAW-LCAO --basis-option size=DZP --h 0.25
```

If `MANDACARU_PAW_PATH` is not set, a PAW-LCAO calculation stops before it
starts, with a `LibraryPathError` naming the command above.

## Options

**Basis size.** `size` picks how many atomic orbitals each atom contributes.

| `size` | Orbitals per valence shell |
|---|---|
| `"SZ"` (default) | single zeta: one per occupied valence orbital |
| `"DZ"` | double zeta: two |
| `"DZP"` | double zeta plus polarization functions |
| `"TZP"` | triple zeta plus polarization functions |

`SZP`, `DZ2P`, `TZ`, `TZ2P`, `QZ`, `QZP` and `QZ2P` follow the same pattern.
A larger basis is more accurate and needs more qubits.

**Confinement.** By default each orbital is confined by an energy shift of
0.1 eV, which contracts the diffuse free-atom orbitals toward their size in a
molecule. It is a basis option:
`basis={"name": "PAW-LCAO", "energy_shift": None}` switches it off.

**Which set.** The calculator's `paw_dataset` option names a folder of this
repository. It defaults to `"lda-sr"`:

```python
atoms.calc = Mandacaru(method="adapt-vqe", pool="spin-orbit",
                       basis={"name": "PAW-LCAO", "size": "SZ"},
                       paw_dataset="lda-fr", grid_spacing=0.25)

atoms.calc = Mandacaru(method="dft", xc="pbe",
                       basis={"name": "PAW-LCAO", "size": "DZP"},
                       paw_dataset="pbe-sr", grid_spacing=0.20)
```

On the command line the same option is `--paw-dataset pbe-sr`. A Kohn-Sham
calculation whose functional differs from its datasets' warns, and names the
folder of the matching set.

Spin-orbit coupling breaks the conservation of the spin projection $S_z$, so
it needs a spin-orbit operator pool (or the generalized Hartree–Fock mean
field, `method="ghf"`). [FR.md](FR.md) explains what changes.

**Mixed bases.** A per-element basis may give each atom its own size, as long
as all atoms use PAW-LCAO:

```python
basis={"O": {"name": "PAW-LCAO", "size": "DZP"}, "H": {"name": "PAW-LCAO"}}
```

The full list of options is in Mandacaru's
[pseudopotential guide](https://mandacaru.readthedocs.io/en/latest/guide/pseudopotentials.html).

## How a dataset is built

The method is Blöchl's projector augmented-wave method, *Phys. Rev. B*
**50**, 17953 (1994). Mandacaru uses its frozen-core form. For each element:

1. **Reference atom.** The all-electron atom is solved self-consistently,
   scalar-relativistically, in the set's functional (LDA or PBE).
2. **Partial waves.** For each angular momentum `l` of the valence, the
   all-electron partial waves are taken at two energies: the bound level and
   one scattering energy above it.
3. **Smooth partial waves.** Inside a cutoff radius `r_c`, each partial wave
   is replaced by a smooth one, expanded in spherical Bessel functions.
4. **Local potential.** A smooth local potential replaces the singular
   nuclear attraction inside a radius that follows the largest cutoff.
5. **Projectors.** Projectors are made dual to the smooth partial waves.
6. **One-center terms.** These are the overlap correction `q`, the kinetic
   and potential differences, the coupling matrix `D`, the compensation
   charge, and a partial core density for the nonlinear core correction.

In a molecule, the smooth partial waves (and their multiple-zeta and
polarization companions) are the basis. The projectors and `D` add a
nonlocal term to the Hamiltonian, and `q` makes the basis overlap
`S + C q Cᵀ`.

**Two simplifications relative to the full PAW method:**

- **Linearized one-center terms.** The one-center energies are linearized
  around the reference atom, so `D` is a fixed matrix per element. This is
  the ultrasoft-pseudopotential form of PAW, and the error is second order
  in how far the atom's density differs from the reference atom's.
- **Frozen core.** The core electrons stay frozen.

**Where the functional enters.** It enters only through the dataset: the
reference atom, the unscreening of the local potential and the one-center
energies. The valence electrons of the molecule are treated with the exact
Coulomb interaction, by Hartree–Fock or a quantum algorithm; there is no
density functional in the molecular Hamiltonian.

**Frozen 4f for Tl–Rn.** For Tl through Rn the filled 4f shell is part of
the frozen core. An empty f channel takes its place, scattering at +0.25 Ha.

**Semicore valence.** Na, Mg, K, Ca, Rb, Sr, Sc, Ti, Y–Tc and Ta carry
their filled $(n-1)s\,(n-1)p$ shell in the valence (sodium: 9 electrons),
so their `p` channel has projectors; Li and Be have an empty `p` channel,
and K, Ca, Rb and Sr an empty `d` one. Cs, Ba, Fr and Ra are not semicore
yet and are flagged.

**Per-element settings.** Cutoffs, local-potential raises, Bessel counts and
second reference energies follow the generator's defaults. Some elements
carry repairs of their own; the generator's `DEFAULT_*_BY_DATASET` tables
list them. The repairs cover mainly:

- s channels that must represent a hydrogen 1s entering the sphere;
- the compact semicore d of Ga–Kr, I and Xe.

## Quality checks and flagged elements

Every dataset was checked when it was generated:

| Check | Requirement |
|---|---|
| Reference atom | the all-electron self-consistent field converged; otherwise the element is refused |
| Ghost states | no spurious bound state in any channel's spectrum, compared with and without projectors, against the all-electron atom |
| Scattering | the phase of the logarithmic derivative at the projector radius within 0.05 rad of the all-electron atom over ε ± 0.5 Ha, and 0.3 rad over ε ± 1 Ha |
| Channels without projectors | every angular momentum up to one above the highest channel, which sees the local potential alone, within 0.2 rad of the all-electron atom over the highest valence level ± 0.5 Ha |
| Intruding 1s | the `s` projectors reproduce a hydrogen 1s orbital entering the sphere: the miss, the norm of what the projected expansion leaves out divided by the norm of the function, below 1 |
| The file | written, read back and moved into place only if every channel's bound levels match the generated dataset to 1e-5 Ha (builds from 2026-10-11 on) |

When the default construction failed a check, the generator tried repairs:
zero norm deficit, a shorter local radius, a raised local potential,
rebalanced cutoffs and more Bessel functions. An element that no repair
cleans is still written, with its defect recorded in the file. Loading it
raises a `GhostStateWarning`, and its `repr` says `SCATTERING OFF` when the
phase is off.

The table below is generated from the files' own records by
`tools/readme_tables.py`; do not edit it by hand. A phase is given as
near/far window (rad); "no projectors" is the local potential alone.

<!-- generated: flagged -->
- `lda-sr/`: 22 of 92 flagged, none for a ghost state.
- `lda-fr/`: 21 of 91 flagged, none for a ghost state.
- `pbe-sr/`: 27 of 92 flagged, none for a ghost state.

| Element | `lda-sr/` | `lda-fr/` | `pbe-sr/` |
|---|---|---|---|
| Ag | p (no projectors) 0.42 rad | p (no projectors) 0.42 rad | p (no projectors) 0.30 rad |
| Cd | p (no projectors) 0.41 rad | p (no projectors) 0.41 rad | p (no projectors) 0.42 rad |
| Cs | s 0.07/0.07 rad; p (no projectors) 0.37 rad | s 0.07/0.07 rad; p (no projectors) 0.38 rad | p (no projectors) 0.67 rad |
| Ba | p (no projectors) 0.44 rad | p (no projectors) 0.44 rad | p (no projectors) 0.33 rad |
| La | f (no projectors) 0.51 rad | — | p (no projectors) 0.58 rad; f (no projectors) 1.01 rad |
| Ce | f 0.28/0.28 rad; p (no projectors) 1.18 rad; d (no projectors) 1.42 rad; g (no projectors) 0.21 rad | f 0.29/0.29 rad; p (no projectors) 1.19 rad; d (no projectors) 1.36 rad; g (no projectors) 0.21 rad | f 0.33/0.33 rad; p (no projectors) 1.27 rad; d (no projectors) 1.31 rad; g (no projectors) 0.22 rad |
| Pr | f 0.14/0.14 rad; p (no projectors) 0.90 rad; d (no projectors) 1.14 rad | f 0.05/0.05 rad | f 0.10/0.10 rad; p (no projectors) 0.72 rad; d (no projectors) 1.12 rad |
| Nd | — | — | f 0.05/0.05 rad; p (no projectors) 0.49 rad; d (no projectors) 1.30 rad |
| Pm | f 0.06/0.06 rad; p (no projectors) 0.55 rad; d (no projectors) 1.20 rad | f 0.06/0.06 rad; p (no projectors) 0.55 rad; d (no projectors) 1.28 rad | f 0.05/0.05 rad; p (no projectors) 0.29 rad; d (no projectors) 1.17 rad |
| Sm | f 0.11/0.11 rad; p (no projectors) 0.98 rad; d (no projectors) 1.37 rad | f 0.11/0.11 rad; p (no projectors) 0.99 rad; d (no projectors) 1.25 rad | f 0.13/0.13 rad; p (no projectors) 1.05 rad; d (no projectors) 1.29 rad |
| Eu | — | f 0.06/0.06 rad; p (no projectors) 0.29 rad; d (no projectors) 1.33 rad | f 0.06/0.06 rad; p (no projectors) 0.38 rad; d (no projectors) 1.28 rad |
| Gd | f 0.10/0.10 rad; p (no projectors) 0.72 rad; d (no projectors) 1.24 rad | f 0.09/0.09 rad; p (no projectors) 0.47 rad; d (no projectors) 1.28 rad | f 0.10/0.10 rad; p (no projectors) 0.55 rad; d (no projectors) 1.21 rad |
| Tb | f 0.17/0.17 rad; p (no projectors) 0.89 rad; d (no projectors) 1.36 rad | f 0.15/0.15 rad; p (no projectors) 0.89 rad; d (no projectors) 1.28 rad | f 0.08/0.08 rad; p (no projectors) 0.95 rad; d (no projectors) 1.29 rad |
| Dy | f 0.10/0.10 rad; p (no projectors) 0.24 rad; d (no projectors) 1.37 rad | f 0.06/0.06 rad | f 0.17/0.17 rad; p (no projectors) 0.92 rad; d (no projectors) 1.47 rad |
| Ho | f 0.15/0.15 rad; d (no projectors) 1.32 rad | f 0.15/0.15 rad; d (no projectors) 1.37 rad | f 0.17/0.17 rad; d (no projectors) 1.36 rad |
| Er | f 0.17/0.17 rad; p (no projectors) 0.62 rad; d (no projectors) 1.30 rad | f 0.17/0.17 rad; p (no projectors) 0.39 rad; d (no projectors) 1.33 rad | f 0.19/0.19 rad; p (no projectors) 0.46 rad; d (no projectors) 1.30 rad |
| Tm | f 0.27/0.27 rad; p (no projectors) 0.80 rad; d (no projectors) 1.34 rad | f 0.22/0.22 rad; p (no projectors) 0.80 rad; d (no projectors) 1.31 rad | p (no projectors) 0.85 rad; d (no projectors) 1.29 rad |
| Yb | f 0.24/0.24 rad; p (no projectors) 0.78 rad; d (no projectors) 1.43 rad | f 0.31/0.31 rad; p (no projectors) 0.78 rad; d (no projectors) 1.40 rad | f 0.23/0.23 rad; p (no projectors) 0.64 rad; d (no projectors) 1.44 rad |
| Lu | s 0.08/0.08 rad; f 0.13/0.13 rad; p (no projectors) 0.89 rad | s 0.10/0.10 rad; p (no projectors) 0.89 rad | 1s miss 2.25; s 0.06/0.06 rad; p (no projectors) 0.31 rad |
| Hf | — | — | p (no projectors) 0.65 rad |
| W | — | — | p (no projectors) 0.54 rad |
| Re | p (no projectors) 1.55 rad | p (no projectors) 0.47 rad | — |
| Ir | p (no projectors) 0.21 rad | p (no projectors) 0.20 rad | p (no projectors) 0.21 rad |
| Fr | — | — | p (no projectors) 0.39 rad |
| Ra | p (no projectors) 0.40 rad | p (no projectors) 0.44 rad | p (no projectors) 0.22 rad |
| Ac | p (no projectors) 0.75 rad; f (no projectors) 1.01 rad | — | p (no projectors) 0.79 rad; f (no projectors) 1.02 rad |
| Th | s 1.53/1.53 rad; d 0.58/0.58 rad; f 0.39/0.39 rad; p (no projectors) 0.93 rad; g (no projectors) 0.23 rad | f 0.08/0.08 rad; p (no projectors) 0.75 rad | — |
| Pa | — | not built | f 0.06/0.06 rad; p (no projectors) 0.61 rad; d (no projectors) 1.30 rad |
| U | — | — | p (no projectors) 0.36 rad; d (no projectors) 1.41 rad |
<!-- end generated: flagged -->

How to read the table:

- **1s miss:** a ratio of norms, not a small excess over a limit. A miss
  of 1 means the projected expansion misses as much as the function holds;
  1.76 is 176 %, an expansion farther from the function than zero is.
  Without the check, misses reached 37.
- **No projectors:** an angular momentum the dataset has no channel for
  (the `p` of Cs and Ba, the `p` and `d` under the lanthanide 4f) scatters
  off the local potential alone; a few tenths of a radian there puts
  lattice constants and bonds several percent off.
- **`f` phase:** the lanthanide 4f is very compact, which limits how well any
  smooth construction reproduces its scattering.
- **What to do:** for molecules containing these elements, check bond lengths
  against an energy scan rather than trusting a relaxed geometry alone.

A record with no defects means that the checks above passed when the file
was generated; it is not a certificate for every chemical environment.
A channel may also not bind more levels below −0.1 Ha than the
all-electron atom has (an excited ghost). The files built before that
check carried one just above the reference of a compact channel whose
second reference energy was still bound: the 4f of W, Re, Os, Ir, Pt, Au
and Hg (W −0.15 to Hg −2.5 Ha, every set), the 4p of Tc (−0.60 Ha) and
the LDA 5p of Ta (−0.46 Ha). They were rebuilt on 2026-10-11 with that
channel's second reference at +0.25 Ha; none of the 26 carries one now.

**Semicore levels.** A few compact semicore shells reproduce their
reference levels to 1–5 mHa rather than the 1 mHa most channels reach:

- the 3d of Se–Kr;
- the 4d of I and Xe;
- the 4f of W–Hg.

Most of that is the overlap, not the construction: the molecular overlap
carries the plain charge of each channel, while the scalar-relativistic
partial waves conserve a mass-weighted norm that differs from it by about
1e-4 to 1e-3. With the norm in the overlap the Au 4f level is exact
(−3.9 mHa with the charge), Pb and Bi 5d move from −0.3 to within 0.1 mHa,
and valence levels change by at most 0.08 mHa.

## Regenerating the datasets

With a Mandacaru development install and `MANDACARU_PAW_PATH` set:

```bash
mandacaru-build --pp PAW --relativistic --xc LDA --all --workers 7 --check --ghosts flag --install
mandacaru-build --pp PAW --dirac --xc LDA --all --workers 7 --check --ghosts flag --install
mandacaru-build --pp PAW --relativistic --xc PBE --all --workers 4 --ghosts flag --install
```

`--install` writes into `lda-sr/`, into `lda-fr/` with `--dirac`, or into
`pbe-sr/` with `--xc PBE`.
A scalar set takes 1.5–2 hours on 7 cores, and a Dirac set about 3 hours.
`--element Fe Cu` rebuilds single elements.

Build on an otherwise idle machine. A worker killed under memory pressure
leaves the process pool waiting forever.

## File format

Each `<Symbol>.parquet` is a self-describing Mandacaru record:

- format `mandacaru-pseudopotential`, version 2, family `paw-lcao`;
<!-- generated: grid -->
- a uniform radial grid starting one step from the origin and ending at 30 Bohr, finer for heavier atoms: 1500 points 0.0200 Bohr apart (hydrogen) to 34500 points 0.00087 Bohr apart (uranium) -- every fourth point of the grid the reference atom was solved on -- with every quantity in atomic units.
<!-- end generated: grid -->

A record holds:

- the all-electron and smooth partial waves;
- the projectors and the one-center matrices;
- the local potential;
- the core and partial core densities;
- the compensation charge and the frozen one-center energy;
- the generation record: functional, relativity, relativistic exchange and
  core correction;
- any recorded defects;
- in `lda-fr/`, the spin-orbit term.

Read one with:

```python
from mandacaru.pseudopotentials.io import load_pseudopotential
pp = load_pseudopotential("lda-sr/O.parquet")
```

## Further reading

- [FR.md](FR.md): how spin-orbit coupling is built and used, and its
  limits.
- [PBE.md](PBE.md): the PBE set, how it is built and checked, its flagged
  elements and its comparison with all-electron PBE.
- Mandacaru's
  [pseudopotential guide](https://mandacaru.readthedocs.io/en/latest/guide/pseudopotentials.html):
  the method in detail, validation, and the unitary variant UPAW-LCAO.

## References

- P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994) — the projector
  augmented-wave method.
- D. D. Koelling and B. N. Harmon, *J. Phys. C* **10**, 3107 (1977) — the
  scalar-relativistic equation.
- J. P. Perdew and A. Zunger, *Phys. Rev. B* **23**, 5048 (1981) — LDA
  correlation.
- J. P. Perdew, K. Burke and M. Ernzerhof, *Phys. Rev. Lett.* **77**, 3865
  (1996) — the PBE functional.
- A. H. MacDonald and S. H. Vosko, *J. Phys. C* **12**, 2977 (1979) — the
  relativistic correction to exchange.
- S. G. Louie, S. Froyen and M. L. Cohen, *Phys. Rev. B* **26**, 1738 (1982)
  — the nonlinear core correction.

## License

MIT, see [LICENSE](LICENSE).
