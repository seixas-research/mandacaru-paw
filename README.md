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
| `lda-dirac/` | 91 (H–U, no Pa) | the same, plus spin-orbit coupling | ~250 MB |

- **LDA** here is Slater exchange with Perdew–Zunger correlation, with the
  relativistic correction to exchange of MacDonald and Vosko.
- **Scalar-relativistic** means the reference atoms include mass-velocity
  and Darwin effects, but not spin-orbit coupling.
- **`lda-dirac/`** adds spin-orbit coupling from the Dirac equation; see
  [DIRAC.md](DIRAC.md).
- **PBE sets** are not shipped. [PBE.md](PBE.md) explains how they are built
  and what is still wrong with them.

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
                       h=0.25)    # grid spacing, Angstrom
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

**Which set.** The calculator's `directory` option names a folder of this
repository. It defaults to `"lda-sr"`:

```python
atoms.calc = Mandacaru(method="adapt-vqe", pool="spin-orbit",
                       basis={"name": "PAW-LCAO", "size": "SZ"},
                       directory="lda-dirac", h=0.25)
```

Spin-orbit coupling breaks the conservation of the spin projection $S_z$, so
it needs a spin-orbit operator pool (or the generalized Hartree–Fock mean
field, `method="ghf"`). [DIRAC.md](DIRAC.md) explains what changes.

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
   scalar-relativistically, in the LDA.
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
| Intruding 1s | the `s` projectors reproduce a hydrogen 1s orbital entering the sphere, with a miss below 1 (in units of its norm); an `s` channel that fails it gives wrong molecules even when every atomic check passes |

When the default construction failed a check, the generator tried repairs:
zero norm deficit, a raised local potential, rebalanced cutoffs and more
Bessel functions. An element that no repair cleans is still written, with
its defect recorded in the file. Loading it raises a `GhostStateWarning`,
and its `repr` says `SCATTERING OFF` when the phase is off.

**No dataset holds a ghost state.** Thirteen elements are flagged for other
checks:

| Element | `lda-sr/` 1s miss | `lda-sr/` `f` phase | `lda-dirac/` 1s miss | `lda-dirac/` `f` phase |
|---|---|---|---|---|
| Tc | 1.01 | — | 1.01 | — |
| Ce | 1.11 | 0.068 rad | 1.17 | 0.069 rad |
| Pr | passes | — | 1.04 | — |
| Pm | 1.05 | 0.055 rad | 1.11 | — |
| Gd | 1.29 | 0.071 rad | 1.28 | 0.068 rad |
| Tb | 1.42 | 0.085 rad | 1.44 | 0.081 rad |
| Dy | 1.76 | 0.081 rad | 1.76 | 0.084 rad |
| Ho | 1.60 | 0.069 rad | 1.60 | 0.071 rad |
| Er | 1.28 | 0.071 rad | 1.45 | 0.061 rad |
| Tm | 1.66 | 0.079 rad | 1.64 | 0.077 rad |
| Yb | 1.64 | 0.054 rad | 1.63 | 0.055 rad |
| Lu | 1.18 | — | 1.21 | — |
| Th | 1.19 | — | 1.18 | 0.050 rad |

How to read the table:

- **1s miss:** the misses are small, 1.0–1.8, just above the limit. Without
  the check they reached 37.
- **`f` phase:** the lanthanide 4f is very compact, which limits how well any
  smooth construction reproduces its scattering.
- **What to do:** for molecules containing these elements, check bond lengths
  against an energy scan rather than trusting a relaxed geometry alone.

**Semicore levels.** A few compact semicore shells reproduce their
reference levels to 1–5 mHa rather than the 1 mHa most channels reach:

- the 3d of Se–Kr;
- the 4d of I and Xe;
- the 4f of W–Hg.

## Regenerating the datasets

With a Mandacaru development install and `MANDACARU_PAW_PATH` set:

```bash
mandacaru-build --pp PAW --relativistic --xc LDA --all --workers 7 --check --ghosts flag --install
mandacaru-build --pp PAW --dirac --xc LDA --all --workers 7 --check --ghosts flag --install
```

`--install` writes into `lda-sr/`, or into `lda-dirac/` with `--dirac`.
A scalar set takes 1.5–2 hours on 7 cores, and a Dirac set about 3 hours.
`--element Fe Cu` rebuilds single elements.

Build on an otherwise idle machine. A worker killed under memory pressure
leaves the process pool waiting forever.

## File format

Each `<Symbol>.parquet` is a self-describing Mandacaru record:

- format `mandacaru-pseudopotential`, version 2, family `paw-lcao`;
- a radial grid of 3000 points from 0.01 to 30 Bohr, with every quantity in
  atomic units.

A record holds:

- the all-electron and smooth partial waves;
- the projectors and the one-center matrices;
- the local potential;
- the core and partial core densities;
- the compensation charge and the frozen one-center energy;
- the generation record: functional, relativity, relativistic exchange and
  core correction;
- any recorded defects;
- in `lda-dirac/`, the spin-orbit term.

Read one with:

```python
from mandacaru.pseudopotentials.io import load_pseudopotential
pp = load_pseudopotential("lda-sr/O.parquet")
```

## Further reading

- [DIRAC.md](DIRAC.md): how spin-orbit coupling is built and used, and its
  limits.
- [PBE.md](PBE.md): the PBE datasets, and why they are not shipped.
- Mandacaru's
  [pseudopotential guide](https://mandacaru.readthedocs.io/en/latest/guide/pseudopotentials.html):
  the method in detail, validation, and the other families (ONCVPSP,
  UPAW-LCAO).

## References

- P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994) — the projector
  augmented-wave method.
- D. D. Koelling and B. N. Harmon, *J. Phys. C* **10**, 3107 (1977) — the
  scalar-relativistic equation.
- J. P. Perdew and A. Zunger, *Phys. Rev. B* **23**, 5048 (1981) — LDA
  correlation.
- A. H. MacDonald and S. H. Vosko, *J. Phys. C* **12**, 2977 (1979) — the
  relativistic correction to exchange.
- S. G. Louie, S. Froyen and M. L. Cohen, *Phys. Rev. B* **26**, 1738 (1982)
  — the nonlinear core correction.

## License

MIT, see [LICENSE](LICENSE).
