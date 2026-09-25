# Mandacaru PAW-LCAO pseudopotentials

Projector augmented-wave datasets in the PAW-LCAO form for
[Mandacaru](https://github.com/seixas-research/mandacaru): one file per
element for every element with **Z ≤ 92** (H through U). Mandacaru generates
them from scratch with its own all-electron radial solver and its
`mandacaru.pseudopotentials.paw` module; nothing here comes from another code.

The datasets live here, not in the Mandacaru package, because of their size
(about 200 MB).

## Contents

| Folder | Datasets | Reference atom |
|---|---|---|
| `lda/` | 92 (H–U) | LDA, scalar-relativistic, nonlinear core correction |

Mandacaru reads `<checkout>/<xc>/`, so a PBE set would go in `pbe/` beside
`lda/`.

## Using the datasets

```bash
git clone https://github.com/seixas-research/mandacaru-paw
mandacaru --set-paw /path/to/mandacaru-paw   # writes MANDACARU_PAW_PATH to ~/.zshrc or ~/.bashrc
# open a new terminal, then
mandacaru --pseudo-status
```

The family is selected as a basis:

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

Without `MANDACARU_PAW_PATH` a PAW-LCAO calculation stops before it starts,
with a `LibraryPathError` that names the command above. The basis option
`directory="..."` points a single run at another folder of datasets.

## Construction

Following P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994), in a frozen-core,
one-center form with a localized (LCAO) basis:

1. a scalar-relativistic LDA atom; the valence partial waves at two energies
   per `l` (the bound level and one scattering energy above it);
2. smooth partial waves as spherical-Bessel expansions inside each channel's
   cutoff `r_c`;
3. a smooth local potential inside `r_cl`, which follows the **largest**
   cutoff; each projector reaches out to `r_cl`, where it is
   `(v_AE − v_loc) φ`, so no channel sits in the bare all-electron well;
4. projectors dual to the smooth waves, the one-center matrices (overlap
   correction `q`, kinetic and potential differences, coupling `D`), a
   compensation charge, and a partial core density for the nonlinear core
   correction.

## Checks, and the flagged elements

Every dataset was checked when it was generated:

- **ghost states:** the spectrum of every channel, with and without
  projectors, against the all-electron reference;
- **scattering:** the phase `arctan L(E)` of the logarithmic derivative,
  compared with the all-electron atom at the projector radius, within
  0.05 rad over ε ± 0.5 Ha and 0.3 rad over ε ± 1 Ha.

When the default construction failed a check, the generator tried a zero
norm deficit, raised local potentials and balanced cutoffs. **No element
holds a ghost state.** An element that no repair cleans is still written,
with its defect recorded in the file; loading it raises a
`GhostStateWarning`, and its `repr` says `SCATTERING OFF`.

| Elements | Defect |
|---|---|
| Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy | `f`-channel phase off by 0.06–0.11 rad |
| Tl, Pb, Bi, Po, At, Rn | phase of the semicore 4f channel off by 0.5–0.65 rad |

For the compact 4f channel this phase is measured at a radius where the 4f
wave has all but vanished, so the second row may overstate the error.

## Generating the datasets

With a Mandacaru development install and `MANDACARU_PAW_PATH` set:

```bash
mandacaru-build --pp PAW --relativistic --xc LDA --all --workers 7 --check --ghosts flag --install
```

`--install` writes into `$MANDACARU_PAW_PATH/lda/`. The radial kernels run
in C; the full set takes about an hour and a half on 7 cores.

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
