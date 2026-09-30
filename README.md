# Mandacaru PAW-LCAO pseudopotentials

Projector augmented-wave datasets in the PAW-LCAO form for
[Mandacaru](https://github.com/seixas-research/mandacaru): one file per
element for every element with **Z ≤ 92** (H through U), in the LDA, in two
sets: scalar-relativistic and Dirac. Mandacaru generates them from scratch with
its own all-electron radial solver and its `mandacaru.pseudopotentials.paw`
module; nothing here comes from another code.

The datasets live here, not in the Mandacaru package, because of their size
(about 200 MB per set).

## Contents

| Folder | Datasets | Reference atom |
|---|---|---|
| `lda-sr/` | 92 (H–U) | LDA, scalar-relativistic, nonlinear core correction |
| `lda-dirac/` | 91 (H–U, no Pa) | LDA, Dirac: the scalar set plus a j-resolved spin-orbit term |

Both sets use the same construction, cutoffs and checks; the Dirac set adds
the spin-orbit term (see *Spin-orbit coupling* below). PBE sets are not
shipped; [PBE.md](PBE.md) records how to build them and what was still wrong
with them.

## Using the datasets

```bash
git clone https://github.com/seixas-research/mandacaru-paw
mandacaru --set-paw /path/to/mandacaru-paw   # writes MANDACARU_PAW_PATH to ~/.zshrc or ~/.bashrc
# open a new terminal, then
mandacaru --pseudo-status
```

The family is selected as a basis, and reads `$MANDACARU_PAW_PATH/lda-sr/` by
default:

```python
from ase.build import molecule
from mandacaru import Mandacaru

atoms = molecule("H2O")
atoms.center(vacuum=4.0)          # the cell is the real-space box
atoms.calc = Mandacaru(method="adapt-vqe",
                       basis={"name": "PAW-LCAO", "size": "DZP"},
                       h=0.25)
atoms.get_potential_energy()
```

The calculator's `directory` option names the set: `"lda-sr"` (the default)
or `"lda-dirac"` (see *Spin-orbit coupling* below).

Without `MANDACARU_PAW_PATH` a PAW-LCAO calculation stops before it starts,
with a `LibraryPathError` that names the command above. The basis option
`directory="..."` points a single run at any folder of datasets.

## Construction

Following P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994), in a frozen-core,
one-center form with a localized (LCAO) basis:

1. a scalar-relativistic all-electron LDA atom, with the relativistic correction to its exchange (MacDonald and
   Vosko, *J. Phys. C* **12**, 2977 (1979)); the valence partial waves at two
   energies per `l` (the bound level and one scattering energy above it).
   The 4f14 of Tl–Rn is in the frozen core, with an empty f channel
   scattering at +0.25 Ha in its place;
2. smooth partial waves as spherical-Bessel expansions inside each channel's
   cutoff `r_c`;
3. a smooth local potential inside `r_cl`, which follows the **largest**
   cutoff; each projector reaches out to `r_cl`, where it is
   `(v_AE − v_loc) φ`, so no channel sits in the bare all-electron well;
4. projectors dual to the smooth waves, the one-center matrices (overlap
   correction `q`, kinetic and potential differences, coupling `D`), a
   compensation charge, and a partial core density for the nonlinear core
   correction.

Cutoffs, local-potential shifts, Bessel counts and second reference energies
follow the generator's defaults, which include per-dataset repairs for most of
the d block (an s channel that has to represent a hydrogen 1s entering the
sphere, and the compact semicore d of Ga–Kr, I and Xe). The generator's
`DEFAULT_*_BY_DATASET` tables list them.

## Checks, and the flagged elements

Every dataset was checked when it was generated:

- **reference atom:** the all-electron SCF must have converged, or the element
  is refused rather than built;
- **ghost states:** the spectrum of every channel, with and without
  projectors, against the all-electron reference;
- **scattering:** the phase `arctan L(E)` of the logarithmic derivative,
  compared with the all-electron atom at the projector radius, within
  0.05 rad over ε ± 0.5 Ha and 0.3 rad over ε ± 1 Ha;
- **intruding 1s:** how badly the `s` projectors miss a hydrogen 1s orbital
  entering the sphere, below 1 (in units of its norm). An `s` channel that
  fails it cannot represent a neighbor's orbital, and a molecule built on it
  goes wrong although every atomic check passes.

When the default construction failed a check, the generator tried a zero
norm deficit, raised local potentials and balanced cutoffs. **No element in
either set holds a ghost state, and none is missing a valence channel.** An
element that no repair cleans is still written, with its defect recorded in
the file; loading it raises a `GhostStateWarning`, and its `repr` says
`SCATTERING OFF`.

The flagged elements, as recorded in the files (built 2026-09-28/30): the
intruding-1s miss, and the largest `f`-channel phase error near the
reference energy where it exceeds 0.05 rad.

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

The misses are all between 1.0 and 1.8; before the check was part of
generation, 29 datasets per set missed by 1 or more, up to 37 (Yb). For the lanthanides, the compact 4f is the same limit
it always was. Treat these elements' bond lengths with care, and check
them against an energy scan.

A few reference levels are reproduced less tightly than the 1 mHa most
channels reach, in the same elements in every set: the compact semicore 3d/4d
of Se–Kr and I–Xe and the 4f of W–Hg (about 1–5 mHa).

## Spin-orbit coupling: the Dirac sets

`lda-dirac/` is the scalar set with, for each `l ≥ 1`,
two more unitary partial-wave branches built from the Dirac atom, one per
`j = l ∓ 1/2`. They are stored as their (2j+1) average and their L·S
difference on the union of the two branches' projectors, so each `j` is
exact; the scalar channels, the overlap and the compensation charges are those
of the scalar construction. Against the Dirac atom: per-`j` levels within
0.08 mHa at the median, most bound channels within 1 mHa (the rest are the
compact semicore shells above, both `j` together); splittings within 0.07 %
at the median and 2.5 % at worst (the 4f of
the light lanthanides). Protactinium is missing: one of its 5f `j` branches is
not bound in the Dirac atom.

```python
atoms.calc = Mandacaru(method="adapt-vqe", pool="spin-orbit",
                       basis={"name": "PAW-LCAO", "size": "SZ"},
                       directory="lda-dirac", h=0.25)
```

The spin-orbit term breaks `S_z`: only ADAPT-VQE with `pool="spin-orbit"`
(and the `"ghf"` mean field) takes it, and the orbitals are generalized
Hartree–Fock spinors. See Mandacaru's pseudopotential and active-space guides.

## Generating the datasets

With a Mandacaru development install and `MANDACARU_PAW_PATH` set:

```bash
mandacaru-build --pp PAW --relativistic --xc LDA --all --workers 7 --check --ghosts flag --install
mandacaru-build --pp PAW --dirac --xc LDA --all --workers 7 --check --ghosts flag --install
```

`--install` writes into `$MANDACARU_PAW_PATH/lda-sr/` and, with `--dirac`,
into `$MANDACARU_PAW_PATH/lda-dirac/`. The radial kernels run in C; a
scalar set takes about an hour and a half on 7 cores, a Dirac set about three
hours.

## File format

Each `<Symbol>.parquet` is a self-describing Mandacaru pseudopotential record
(format `mandacaru-pseudopotential`, version 2, family `paw-lcao`) on a
3000-point radial grid (0.01 Bohr out to 30 Bohr). It holds the all-electron
and smooth partial waves, the projectors, the one-center matrices, the local
potential, the core and partial core densities, the compensation charge and
the frozen one-center energy, plus the generation record (functional,
relativity, relativistic exchange, core correction), any recorded defects and,
in the Dirac sets, the `spin_orbit` term. All quantities are in
atomic units. Read one with
`mandacaru.pseudopotentials.io.load_pseudopotential(path)`.

## License

MIT, see [LICENSE](LICENSE).
