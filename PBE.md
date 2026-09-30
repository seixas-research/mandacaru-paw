# PBE PAW-LCAO datasets: shelved (2026-09-29)

This repository ships two PAW-LCAO sets, both LDA: `lda-sr/`
(scalar-relativistic) and `lda-dirac/` (with the j-resolved spin-orbit
term). The PBE sets (`pbe/` and `pbe-dirac/`) were removed on 2026-09-29 to
keep one functional well characterized rather than two half-finished. This
note keeps what was learned about building them, so the work can be picked up
again without rediscovering it. Nothing built for PBE is kept.

The generator still builds PBE datasets. Every rule below is encoded in
Mandacaru (`src/mandacaru/pseudopotentials/paw.py`, `basis/xc.py`) and applies
automatically to `--xc PBE`; nothing here has to be re-entered by hand.

## How to build them

```bash
mandacaru-build --pp PAW --relativistic --xc PBE --all --workers 7 --check --ghosts flag -o <dir>
mandacaru-build --pp PAW --dirac --xc PBE --all --workers 7 --check --ghosts flag -o <dir>
```

A scalar set took about 2 h on 7 cores (2026-09-29); a Dirac set about 3 h.
Build on an otherwise idle machine: one Dirac build hung when a worker was
killed under memory pressure from other jobs, and the process pool then waited
forever.

## The rules the generator applies to PBE

**Functional.** PBE differentiates the density twice. Its radial derivatives
are taken on points spaced `max(h, 0.01 r)` (`basis/xc.py`,
`GGA_DERIVATIVE_STEP`); differentiating at every point of a heavy atom's fine
grid amplified round-off until the reference atom did not converge and the
local potential, matched through its fourth derivative, broke. A relativistic
reference atom scales PBE's exchange by the relativistic factor of
MacDonald and Vosko, `Phi_E(k_F/c)`; the correlation is unchanged.

**Per-dataset settings** (the `("<element>", "pbe")` entries of the tables in
`paw.py`), each the combination a scan measured:

- *Bessel functions*: Li 9 (at 8 its 2s pair has a node at 1.22 Bohr and LiH
  binds at 2.4 Angstrom instead of 1.6), Fe, Co, Ni 9, Sn 10.
- *Norm deficit 0* (unitary) for Fe, Co, Ni, Cu, Zn and every repaired
  dataset below. Cu-PBE at deficit 0.1 passed every atomic check and still
  collapsed CuH.
- *s channel* (the projectors had to represent a hydrogen 1s entering the
  sphere): 10 Bessel functions and a local shift of 0-50 Ha for Ac, Ag, Au,
  Cd, Cr, Hg, Ir, La, Mn, Pd, Pt, Rh, Ru, Y, Zr, and a shorter s sphere for
  Os (3.041 Bohr), Re (2.981), Ta (3.223) and V (3.046).
- *Semicore d* (the d projectors amplified the basis filter at
  h = 0.20 Angstrom): d spheres of As 1.612, Se 1.531, Br 1.009, Kr 1.070,
  I 1.113, Xe 1.070 Bohr with the second reference energy at 2-3 Ha.
- The full cutoffs of each repaired dataset are in
  `DEFAULT_CUTOFFS_BY_DATASET`: a repair is its whole parameter set, and
  leaving the other channels to the generic rule broke Os-PBE's f channel.

**Generation checks** (every dataset, both functionals): no ghost state;
scattering phase within 0.05 rad near the reference and 0.3 rad over
eps +- 1 Ha; and, since 2026-09-29, the s channel's miss of an intruding
hydrogen 1s below 1 (`paw.intruding_s_miss`), which the repair search
enforces like a phase error.

## What was still wrong when it was shelved

Measured on the last PBE builds (2026-09-28/29), before the s-channel check
was in the generator:

- **s channel incomplete** (intruding-1s miss >= 1) in 23 scalar datasets:
  Na, Mg, Ca, Sc, Ga, Nb, Tc, Ce, Sm, Eu, Dy, Ho, Er, Tm, Yb, Lu, Hf (12.2),
  W (3.6), Pt, Pb (1.9), Th, Pa, U. The check now repairs most such cases
  (Hf-LDA went from 6.8 to 0.43), but PBE has not been rebuilt with it.
- **f-channel phase** flagged for lanthanides: Ce, Pr, Nd (0.16 rad), Pm, Sm,
  Eu, Gd, Tb, Dy, Ho, Er, Yb (0.05-0.10 rad); the compact 4f is the same
  limit in LDA.
- **No repair found** in the scans: Hf, Lu, Nb, Sc, Tc, Th, W; the semicore d
  of Ga and Ge.
- **Dirac:** Pa cannot be built (one 5f j-branch is unbound in the Dirac atom,
  eps = +0.0010 Ha; Mandacaru TODO 1.15). Tl-PBE-Dirac and Bi-PBE-Dirac had
  s-channel misses of 3.1 and 2.2 before the check; the check repairs Tl
  (0.24).

## To bring PBE back

1. Build both sets with the commands above, on an idle machine.
2. Verify: no ghost, flags listed, the intruding-1s miss over every dataset,
   and the filtered-d projection for the d block (as done for `lda-sr/`).
3. For the Dirac set, audit per-j levels and splittings against the Dirac
   atom (Mandacaru's phase-6 audit).
4. Add PBE folders to Mandacaru's library-folder table
   (`environment.LIBRARY_FOLDERS`, e.g. `pbe-sr` and `pbe-dirac`), bring
   back PBE PAW tests, and document the folders in the README.
