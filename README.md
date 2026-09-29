# Mandacaru PAW-LCAO pseudopotentials

Projector augmented-wave datasets in the PAW-LCAO form for
[Mandacaru](https://github.com/seixas-research/mandacaru): one file per
element for every element with **Z ≤ 92** (H through U), for two
exchange-correlation functionals. Mandacaru generates them from scratch with
its own all-electron radial solver and its `mandacaru.pseudopotentials.paw`
module; nothing here comes from another code.

The datasets live here, not in the Mandacaru package, because of their size
(about 200 MB per functional).

## Contents

| Folder | Datasets | Reference atom |
|---|---|---|
| `lda/` | 92 (H–U) | LDA, scalar-relativistic, nonlinear core correction |
| `pbe/` | 92 (H–U) | PBE, scalar-relativistic, nonlinear core correction |
| `lda-dirac/` | 91 (H–U, no Pa) | LDA, Dirac: the scalar set plus a j-resolved spin-orbit term |
| `pbe-dirac/` | 91 (H–U, no Pa) | PBE, Dirac: the scalar set plus a j-resolved spin-orbit term |

All sets use the same construction, cutoffs and checks; only the functional
of the reference atom (and of the unscreening) differs, and the Dirac sets
add the spin-orbit term (see *Spin-orbit coupling* below).

## Using the datasets

```bash
git clone https://github.com/seixas-research/mandacaru-paw
mandacaru --set-paw /path/to/mandacaru-paw   # writes MANDACARU_PAW_PATH to ~/.zshrc or ~/.bashrc
# open a new terminal, then
mandacaru --pseudo-status
```

The family is selected as a basis, and reads `$MANDACARU_PAW_PATH/lda/` by
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

To use the PBE datasets, name the `pbe/` folder with the calculator's
`directory` option (in Mandacaru releases after v26.9.52; `"lda"` is the default):

```python
atoms.calc = Mandacaru(method="adapt-vqe",
                       basis={"name": "PAW-LCAO", "size": "DZP"},
                       directory="pbe",   # $MANDACARU_PAW_PATH/pbe/
                       h=0.25)
```

Without `MANDACARU_PAW_PATH` a PAW-LCAO calculation stops before it starts,
with a `LibraryPathError` that names the command above. The basis option
`directory="..."` points a single run at any folder of datasets.

## Construction

Following P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994), in a frozen-core,
one-center form with a localized (LCAO) basis:

1. a scalar-relativistic all-electron atom in the set's functional (LDA or
   PBE), with the relativistic correction to its exchange (MacDonald and
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

For PBE, the gradient correction differentiates the density twice. Its radial
derivatives are therefore taken on grid points spaced `max(h, 0.01 r)`: every
point near the nucleus, logarithmic spacing further out. On the fine uniform
grid of a heavy atom, differentiating at every point amplified the orbitals'
round-off until the reference atom no longer converged and the local
potential, matched through its fourth derivative at `r_cl`, broke.

## Checks, and the flagged elements

Every dataset was checked when it was generated:

- **reference atom:** the all-electron SCF must have converged, or the element
  is refused rather than built;
- **ghost states:** the spectrum of every channel, with and without
  projectors, against the all-electron reference;
- **scattering:** the phase `arctan L(E)` of the logarithmic derivative,
  compared with the all-electron atom at the projector radius, within
  0.05 rad over ε ± 0.5 Ha and 0.3 rad over ε ± 1 Ha.

When the default construction failed a check, the generator tried a zero
norm deficit, raised local potentials and balanced cutoffs. **No element in
either set holds a ghost state, and none is missing a valence channel.** An
element that no repair cleans is still written, with its defect recorded in
the file; loading it raises a `GhostStateWarning`, and its `repr` says
`SCATTERING OFF`.

| Elements | `lda/` | `pbe/` |
|---|---|---|
| Ce, Pr, Pm, Sm, Eu, Gd, Tb, Dy | `f`-channel phase off by 0.06–0.12 rad | `f`-channel phase off by 0.05–0.10 rad |
| Nd | clean | `f`-channel phase off by 0.16 rad |
| Ho, Er, Yb | clean | `f`-channel phase off by 0.05–0.07 rad |

The Dirac sets flag the same compact lanthanide 4f channels (LDA: Ce, Pr, Pm,
Sm, Eu, Gd, Tb, Dy; PBE: Ce, Sm, Eu, Gd, Tb, Dy, Ho, Er), 0.06–0.13 rad.

A few reference levels are reproduced less tightly than the 1 mHa most
channels reach, in the same elements in every set: the compact semicore 3d/4d
of Se–Kr and I–Xe and the 4f of W–Hg (about 1–5 mHa).

## Spin-orbit coupling: the Dirac sets

`lda-dirac/` and `pbe-dirac/` are the scalar sets with, for each `l ≥ 1`,
two more unitary partial-wave branches built from the Dirac atom, one per
`j = l ∓ 1/2`. They are stored as their (2j+1) average and their L·S
difference on the union of the two branches' projectors, so each `j` is
exact; the scalar channels, the overlap and the compensation charges are those
of the scalar construction. Against the Dirac atom (LDA / PBE): per-`j`
levels within 0.08 mHa at the median, 94 of 106 bound channels within 1 mHa
(the rest are the compact semicore shells above, both `j` together);
splittings within 0.07 / 0.08 % at the median and 2.5 % at worst (the 4f of
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
mandacaru-build --pp PAW --relativistic --xc PBE --all --workers 7 --check --ghosts flag --install
mandacaru-build --pp PAW --dirac --xc LDA --all --workers 7 --check --ghosts flag -o <dir>
mandacaru-build --pp PAW --dirac --xc PBE --all --workers 7 --check --ghosts flag -o <dir>
```

`--install` writes into `$MANDACARU_PAW_PATH/<xc>/`; the Dirac sets are
copied into `lda-dirac/` and `pbe-dirac/`. The radial kernels run in C; a
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
