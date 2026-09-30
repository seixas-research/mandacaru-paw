# PBE PAW-LCAO datasets

This repository ships only LDA datasets (`lda-sr/` and `lda-dirac/`).
Mandacaru's generator can also build PAW-LCAO datasets with the PBE
functional. This note describes how those datasets are constructed, what
works, and what is still wrong with them — the reasons they are not shipped.

## What the functional changes

In Mandacaru the functional belongs to the **dataset**, not to the molecular
calculation. A PAW-LCAO calculation treats the valence electrons with an
exact Coulomb interaction (Hartree–Fock, VQE, ADAPT-VQE, exact
diagonalization); there is no exchange-correlation functional in the
molecular Hamiltonian. The functional enters in three places:

- the all-electron **reference atom** the partial waves are cut from;
- the **unscreening**, which removes the valence Hartree and
  exchange-correlation potential from the atom's potential to leave the
  ionic local potential;
- the **one-center energies**, and the frozen core–valence interaction,
  linearized around the reference atom.

A PBE dataset therefore describes the same atom with a different core and a
different core–valence exchange-correlation. It does not make the molecular
calculation "PBE". Differences between LDA and PBE datasets in a molecule are
differences in the pseudopotential. They are not the functional-level
differences a DFT user would expect.

## How a PBE dataset is built

The construction is the LDA one (see the [README](README.md)), with the
functional swapped. Everything below is encoded in Mandacaru
(`src/mandacaru/pseudopotentials/paw.py`, `src/mandacaru/basis/xc.py`) and
applies automatically to `--xc PBE`.

**The functional.** PBE is the functional of Perdew, Burke and Ernzerhof.
Its correlation is built on the Perdew–Wang (1992) uniform gas, which is
what PBE is defined against, not on the Perdew–Zunger parameterization the
LDA sets use.

**Radial derivatives.** A gradient functional differentiates the density
twice. The radial derivatives are taken on points spaced `max(h, 0.01 r)`
(`GGA_DERIVATIVE_STEP`) rather than at every grid point. On a heavy atom's
fine inner grid, point-by-point differencing amplifies round-off until the
reference atom does not converge. The local potential, which is matched
through its fourth derivative, then breaks.

**Relativistic exchange.** A relativistic reference atom scales PBE's
exchange by the MacDonald–Vosko factor $\Phi_E(k_F/c)$, as for LDA. The
potential gets the exact derivative, which for PBE includes an extra term.
The correlation is unchanged.

**Per-dataset settings.** The generic construction does not pass every check
for every element, so some elements carry their own settings. These are the
`("<element>", "pbe")` entries of the `DEFAULT_*_BY_DATASET` tables in
`paw.py`. Each entry is the full parameter set a scan found, not a single
adjusted knob: changing one channel's cutoff and leaving the rest to the
generic rule has produced new ghost and phase errors in other channels.

| Problem | Elements | Setting |
|---|---|---|
| a 2s pair with a node too far out | Li | 9 Bessel functions (at 8, LiH binds at 2.4 Å instead of 1.6 Å) |
| a d-channel expansion that is too short | Fe, Co, Ni, Sn | 9–10 Bessel functions |
| a norm deficit that passed every atomic check but collapsed the molecule | Fe, Co, Ni, Cu, Zn and every repaired dataset | norm deficit 0 (unitary partial waves) |
| an s channel that cannot represent a hydrogen 1s entering the sphere | Ac, Ag, Au, Cd, Cr, Hg, Ir, La, Mn, Pd, Pt, Rh, Ru, Y, Zr | 10 Bessel functions and a raised local potential (0–50 Ha) |
| the same, when that is not enough | Os, Re, Ta, V | a shorter s sphere (2.98–3.22 Bohr) |
| a semicore d that the grid's basis filter distorts | As, Se, Br, Kr, I, Xe | smaller d spheres (1.0–1.6 Bohr) with the second reference energy at 2–3 Ha |

## Building a PBE set

```bash
mandacaru-build --pp PAW --relativistic --xc PBE --all --workers 7 --check --ghosts flag -o <dir>
mandacaru-build --pp PAW --dirac --xc PBE --all --workers 7 --check --ghosts flag -o <dir>
```

A scalar set takes about 2 hours on 7 cores, a Dirac set about 3.

Build on an otherwise idle machine. If a worker is killed under memory
pressure, the process pool waits forever.

**Checks.** Every dataset gets the same checks as the LDA sets, and a
dataset that fails one is written with its defect recorded:

- no ghost state;
- the scattering phase within 0.05 rad near the reference energy and 0.3 rad
  over ε ± 1 Ha;
- the intruding-1s miss below 1.

## What works

- **Reference atoms converge for every element**, with the derivative
  spacing above.
- **No ghost states.** The last complete builds of both PBE sets had none.
- **The repaired hydrides bind.**
  - Once copper is built at norm deficit 0, PBE CuH has its minimum near
    1.6 Å, like LDA; at the default deficit it collapsed with no minimum.
  - With the 9-function lithium, LiH binds near 1.6 Å.
  - FeH stays within 16 eV over 1.4–2.9 Å.
- **Spin-orbit datasets build.** The Dirac PBE set matched the Dirac atom as
  closely as the LDA one:
  - per-j levels within 0.08 mHa at the median;
  - splittings within 0.08 % at the median and 2.5 % at worst.

## What does not work

**Incomplete s channels.** The last complete PBE sets were built before the
intruding-1s check was part of generation. In the scalar set, 23 datasets
missed an intruding hydrogen 1s by more than its norm: Na, Mg, Ca, Sc, Ga,
Nb, Tc, Ce, Sm, Eu, Dy, Ho, Er, Tm, Yb, Lu, Hf (12.2), W (3.6), Pt, Pb
(1.9), Th, Pa and U. A molecule built on such an s channel goes wrong even
though every atomic check passes. The generator now repairs most such cases
during the build, as it did for the LDA sets, but no PBE set has been
rebuilt with the check.

**No repair found.** The scans found no clean construction for:

- Hf, Lu, Nb, Sc, Tc, Th and W;
- the semicore d of Ga and Ge.

A rebuild writes these with their defects recorded.

**Lanthanide 4f phases.** In the scalar set, Ce, Pr, Nd (0.16 rad), Pm, Sm,
Eu, Gd, Tb, Dy, Ho, Er and Yb miss the phase window by 0.05–0.10 rad. The 4f
is compact, and the same limit appears in LDA.

**Protactinium (Dirac).** One 5f j-branch is unbound in the Dirac atom
(ε ≈ +0.001 Ha), so the element cannot be built. The same happens in LDA.

**Not loadable by name.** Mandacaru's table of library folders
(`environment.LIBRARY_FOLDERS`) has no PBE entry, and there are no PBE
PAW-LCAO tests. A PBE set built into a checkout would have to be added
there before calculations can select it.

## To ship a PBE set

1. Build both sets on an idle machine, with the commands above.
2. Check the result: no ghosts, the flagged elements listed, the
   intruding-1s miss of every dataset, and the filtered-d projection for the
   d block, as done for `lda-sr/`.
3. For the Dirac set, audit the per-j levels and splittings against the
   Dirac atom.
4. Add the folders to Mandacaru's library table (`LIBRARY_FOLDERS`, for
   example `pbe-sr` and `pbe-dirac`).
5. Add PBE PAW-LCAO tests, and document the folders in the README.

## References

- J. P. Perdew, K. Burke and M. Ernzerhof, *Phys. Rev. Lett.* **77**, 3865
  (1996) — the PBE functional.
- J. P. Perdew and Y. Wang, *Phys. Rev. B* **45**, 13244 (1992) — the
  uniform-gas correlation PBE is built on.
- A. H. MacDonald and S. H. Vosko, *J. Phys. C* **12**, 2977 (1979) — the
  relativistic correction to exchange.
- P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994) — the projector
  augmented-wave method.
