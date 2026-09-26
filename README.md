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

Both sets use the same construction, cutoffs and checks; only the functional
of the reference atom (and of the unscreening) differs.

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

To use the PBE datasets, point the basis at the `pbe/` folder:

```python
import os

pbe = os.path.join(os.environ["MANDACARU_PAW_PATH"], "pbe")
atoms.calc = Mandacaru(method="adapt-vqe",
                       basis={"name": "PAW-LCAO", "size": "DZP",
                              "directory": pbe},
                       h=0.25)
```

Without `MANDACARU_PAW_PATH` a PAW-LCAO calculation stops before it starts,
with a `LibraryPathError` that names the command above. The basis option
`directory="..."` points a single run at any folder of datasets.

## Construction

Following P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994), in a frozen-core,
one-center form with a localized (LCAO) basis:

1. a scalar-relativistic all-electron atom in the set's functional (LDA or
   PBE); the valence partial waves at two energies per `l` (the bound level
   and one scattering energy above it);
2. smooth partial waves as spherical-Bessel expansions inside each channel's
   cutoff `r_c`;
3. a smooth local potential inside `r_cl`, which follows the **largest**
   cutoff; each projector reaches out to `r_cl`, where it is
   `(v_AE − v_loc) φ`, so no channel sits in the bare all-electron well;
4. projectors dual to the smooth waves, the one-center matrices (overlap
   correction `q`, kinetic and potential differences, coupling `D`), a
   compensation charge, and a partial core density for the nonlinear core
   correction.

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

| Elements | LDA | PBE |
|---|---|---|
| Ce, Pr, Sm, Eu, Gd, Tb, Dy | `f`-channel phase off by 0.06–0.11 rad | `f`-channel phase off by 0.06–0.11 rad |
| Nd, Pm | `f`-channel phase off by 0.06–0.07 rad | clean |
| Ho, Er | clean | `f`-channel phase off by 0.06–0.07 rad |
| Tl, Pb, Bi, Po, At, Rn | semicore 4f phase off by 0.49–0.65 rad | semicore 4f phase off by 0.49–0.65 rad |

For the compact 4f channel this phase is measured at a radius where the 4f
wave has all but vanished, so the last row may overstate the error.

A few reference levels are reproduced less tightly than the 1 mHa most
channels reach, in the same elements in both sets: the semicore 3d/4d of
Zn–Kr and I–Xe and the 4f of Pm, Gd and W–Hg (about 1–6 mHa), and the flagged
Tl–Rn semicore 4f (about 12–40 mHa).

## Generating the datasets

With a Mandacaru development install and `MANDACARU_PAW_PATH` set:

```bash
mandacaru-build --pp PAW --relativistic --xc LDA --all --workers 7 --check --ghosts flag --install
mandacaru-build --pp PAW --relativistic --xc PBE --all --workers 7 --check --ghosts flag --install
```

`--install` writes into `$MANDACARU_PAW_PATH/<xc>/`. The radial kernels run
in C; each set takes roughly an hour and a half on 7 cores.

## File format

Each `<Symbol>.parquet` is a self-describing Mandacaru pseudopotential record
(format `mandacaru-pseudopotential`, version 2, family `paw-lcao`) on a
3000-point radial grid (0.01 Bohr out to 30 Bohr). It holds the all-electron
and smooth partial waves, the projectors, the one-center matrices, the local
potential, the core and partial core densities, the compensation charge and
the frozen one-center energy, plus the generation record (functional,
relativity, core correction) and any recorded defects. All quantities are in
atomic units. Read one with
`mandacaru.pseudopotentials.io.load_pseudopotential(path)`.

## License

MIT, see [LICENSE](LICENSE).
