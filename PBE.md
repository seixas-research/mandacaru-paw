# PBE PAW-LCAO datasets

`pbe-sr/` holds a scalar-relativistic PAW-LCAO dataset generated with the
PBE functional for every element from hydrogen to uranium (Z ≤ 92). This note
describes how the set is built, how it was checked, what is still flagged,
and how it compares with all-electron PBE. There is no PBE Dirac
(spin-orbit) set.

## Selecting it

The LDA set `lda-sr/` stays the default. A calculation reads the PBE set
with the calculator's `directory` option:

```python
atoms.calc = Mandacaru(method="dft", xc="pbe",
                       basis={"name": "PAW-LCAO", "size": "DZP"},
                       directory="pbe-sr")
```

On the command line it is `--directory pbe-sr`, and from Python
`get_paw("O", xc="pbe")` loads `pbe-sr/O.parquet`. A Kohn-Sham calculation
whose functional is not the datasets' warns that the two differ, and the
warning names the folder of the matching set.

## What the functional changes

In Mandacaru the functional belongs to the **dataset**. The functional
enters the dataset in three places:

- the all-electron **reference atom** the partial waves are cut from;
- the **unscreening**, which removes the valence Hartree and
  exchange-correlation potential from the atom's potential to leave the
  ionic local potential;
- the **one-center energies**, and the frozen core–valence interaction,
  linearized around the reference atom.

A correlated calculation on these datasets (Hartree–Fock, VQE, ADAPT-VQE)
treats the valence electrons with the exact Coulomb interaction; only the
core and the augmentation spheres carry PBE. A Kohn-Sham calculation with
`xc="pbe"` is PBE throughout.

## How a PBE dataset is built

The construction is the LDA one (see the [README](README.md)), with the
functional swapped. Everything below is in Mandacaru
(`src/mandacaru/pseudopotentials/paw.py`, `partial_waves.py`,
`src/mandacaru/basis/xc.py`).

**The functional.** PBE is the functional of Perdew, Burke and Ernzerhof.
Its correlation is built on the Perdew–Wang (1992) uniform gas. A
relativistic reference atom scales PBE's exchange by the MacDonald–Vosko
factor, as for LDA; the potential gets the exact derivative.

**Radial derivatives.** A gradient functional differentiates the density
twice. The radial derivatives are taken on points spaced `max(h, 0.01 r)`
(`GGA_DERIVATIVE_STEP`) and splined back, so the PBE potential of the
reference atom is a cubic spline through those points.

**The local potential is matched over a wide window.** The local potential
continues the atom's potential inside `r_cl` with a polynomial matched
through the fourth derivative at `r_cl`. For LDA those derivatives are read
off a fit over ±25 grid points. For PBE that window spans two or three
spline points, and their kinks read as a fourth derivative: aluminum's came
out −51 instead of about −0.3, its local potential sank to −13 Ha at the
origin (LDA −1.5), and at shorter radii it bound ghost states up to 112 Ha
deep. The PBE sets built before this was fixed (October 2026) needed a
repair for almost every element for that reason alone. PBE derivatives are
now read over ±0.15 `r_cl` (`GGA_LOCAL_FIT_FRACTION`), and a PBE dataset is
built like its LDA counterpart: the same light elements come out with the
same settings, and aluminum, oxygen, sodium and magnesium no longer need a
repair.

**Semicore shells.** Na, Mg, K, Ca, Rb, Sr and the early d block Sc, Ti,
Y–Tc and Ta carry their `(n−1)s (n−1)p` shell in valence, in both
functionals. Without it these elements have no `p` projectors, and the
local potential alone scattered `p` 0.2–1.5 rad off for every repair: PBE
ScH came out 11 % long, bcc Nb 4 %, Mo 12 %, hcp Ti 6 %. With it each
builds clean and comes within 0.7 % of all-electron PBE. Sc and Ti also take
larger `d` spheres. Li and Be carry an empty `p` channel (PBE BeO was 11 %
long without it, 1 % with it), K–Sr an empty `d`.

Semicore Lu, Hf, W and Re passed every atomic check but failed in the solid
(bcc W and hcp Re 4–5 % short, no converged SCF for hcp Hf and Lu); semicore
Ag, Cd and Ir lost their filtered-`d` projection. These keep their valence
and their flags (below).

**Per-dataset settings.** A few elements carry their own settings, the
`("<element>", "pbe")` entries of the `DEFAULT_*_BY_DATASET` tables in
`paw.py`: s-channel repairs of Ac, Ag, Au, Bi, Cd, Cr, Hg, Ir, La, Mn, Os,
Pd, Pt, Re, Rh, Ru, Sn and V; d-channel repairs of Ga, Ge, As, Se, Br, Kr,
I and Xe; lithium at 9 Bessel functions. Ga, Ge, Sn and Bi now take the LDA
entries; the others were re-scanned where it mattered (Li, V, As–Kr) and
kept.

## Building the set

```bash
mandacaru-build --pp PAW --relativistic --xc pbe --all --workers 4 --ghosts flag --install
```

`--install` writes into `pbe-sr/`. The October 2026 build took 12.2 hours on
4 workers on a machine shared with two LDA builds (173,000 worker-seconds;
actinium, thorium, protactinium and uranium 2–2.4 hours each, hydrogen 4 s). Build on an
otherwise idle machine if possible: a worker killed under memory pressure
leaves the process pool waiting forever.

## Checks

Every dataset was checked when it was generated, and the set was audited
again from the files:

| Check | Result |
|---|---|
| Reference atom | converged for every element |
| Ghost states | none, in any channel or in the local potential alone |
| Intruding hydrogen 1s (`s` projectors) | below 1 for every element but Lu (2.25, flagged) |
| Scattering of the projected channels | within 0.05 rad near the reference and 0.3 rad over ±1 Ha, but the lanthanide `4f` and Pa (below) |
| Scattering of the channels without projectors | within 0.2 rad, but the flagged elements below |
| Filtered `d` projection (h = 0.20 Å) | within 0.13 of exact for every `d` channel (worst Kr 0.12, Se 0.11) |
| Loading | every file loads through `get_paw(symbol, xc="pbe")` |

**Flagged elements.** 27 datasets are written with their defects recorded,
and warn when loaded. None holds a ghost:

| Element | Channels | What is flagged |
|---|---|---|
| Cs, Ba, Fr, Ra | `s` | the empty `p` (0.22–0.67 rad) |
| La, Ac | `s`, `d` | the empty `p` (0.58, 0.79 rad) and `f` (1.01, 1.02 rad) |
| Ce–Yb | `s`, `f` | the empty `p` (0.29–1.27 rad) and `d` (1.12–1.47 rad); the `4f` itself 0.05–0.33 rad |
| Pa, U | `s`, `f` | the empty `p` (0.61, 0.36 rad) and `d` (1.30, 1.41 rad) |
| Lu | `s`, `d`, `f` | the empty `p` (0.31 rad); the `s` projectors miss an intruding hydrogen 1s by 2.25 |
| Hf, W, Ir | `s`, `d`, `f` | the empty `p` (0.65, 0.54, 0.21 rad); PBE bcc W comes out 12 % long |
| Ag, Cd | `s`, `d` | the empty `p` (0.30, 0.42 rad) |

The lanthanides and actinides need projectors in channels below their
highest one (an empty `d` under the `4f`), which the generator does not
build yet; the LDA set has the same flags. Check a bond to any flagged
element against an energy scan.

The LDA set staged before these fixes flags 36 elements. It lacks the
semicore early d block and beryllium's `p` channel, which a rebuild with the
current generator gives it; thorium, flagged there, is clean in PBE.

## Against all-electron PBE

References: all-electron PBE with pyscf, def2-QZVPP; crystals from all-electron
(LAPW) PBE in the literature. PAW-LCAO with `pbe-sr/`: DZP basis, h = 0.20 Å
for molecules, 0.25 Å (and 0.18 Å) for crystals.

Molecules (bond lengths in Å; DZP, h = 0.20 Å):

| | All-electron PBE | `pbe-sr/` | error |
|---|---|---|---|
| H₂ | 0.7501 | 0.7840 | +4.5 % |
| LiH | 1.6053 | 1.6341 | +1.8 % |
| Li₂ | 2.7286 | 2.7571 | +1.0 % |
| N₂ | 1.1021 | 1.1558 | +4.9 % |
| CO | 1.1353 | 1.1832 | +4.2 % |
| HCl | 1.2883 | 1.3226 | +2.7 % |
| AlH | 1.6771 | 1.7135 | +2.2 % |
| LiF | 1.5786 | 1.6093 | +1.9 % |
| NaH | 1.8935 | 1.9085 | +0.8 % |
| Na₂ | 3.0891 | 3.0823 | −0.2 % |
| NaCl | 2.3740 | 2.3873 | +0.6 % |
| ScH | 1.7582 | 1.7456 | −0.7 % |
| BeO | 1.3408 | 1.3539 | +1.0 % |

The same generator's LDA datasets with LDA against all-electron LDA are as
far off (+1.9 % on average against +2.2 % here). The PBE − LDA change of
the bond, which is what the datasets' functional decides, follows the
all-electron change for N₂, CO and LiF to 0.004 Å; PBE run on the LDA
datasets misses it there by 0.008–0.020 Å. Sodium compounds come out
0.02 Å longer than all-electron PBE.

Crystals (lattice constants in Å; DZP, 6³ k-points, h = 0.25 Å unless
noted):

| | All-electron PBE | `pbe-sr/` | error |
|---|---|---|---|
| diamond Si | 5.469 | 5.563 | +1.7 % |
| fcc Al | 4.040 | 4.138 | +2.4 % |
| bcc Na | 4.200 | 4.185 | −0.4 % |
| rocksalt MgO | 4.258 | 4.265 (h = 0.18) | +0.2 % |
| rocksalt NaCl | 5.697 | 5.680 (h = 0.18) | −0.3 % |
| bcc Nb | 3.31 | 3.289 | −0.6 % |
| bcc Mo | 3.16 | 3.163 | +0.1 % |
| bcc Ta | 3.32 | 3.320 | 0.0 % |
| hcp Ti | 2.94 | 2.919 | −0.7 % |
| hcp Re | ~2.77 | 2.795 | +0.9 % |
| fcc Ag | 4.15 | 4.067 | −2.0 % |
| bcc W (flagged) | 3.17 | 3.548 | +12 % |

The DZP basis dominates the absolute bond errors (H₂ is 4–5 % long in every
column, and QZ2P removes most of it). The crystals of semicore Na and Mg
need a finer grid than 0.25 Å.

## References

- J. P. Perdew, K. Burke and M. Ernzerhof, *Phys. Rev. Lett.* **77**, 3865
  (1996) — the PBE functional.
- J. P. Perdew and Y. Wang, *Phys. Rev. B* **45**, 13244 (1992) — the
  uniform-gas correlation PBE is built on.
- A. H. MacDonald and S. H. Vosko, *J. Phys. C* **12**, 2977 (1979) — the
  relativistic correction to exchange.
- P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994) — the projector
  augmented-wave method.
- P. Haas, F. Tran and P. Blaha, *Phys. Rev. B* **79**, 085104 (2009) —
  all-electron LDA and PBE lattice constants.
