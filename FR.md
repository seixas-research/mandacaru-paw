# Fully relativistic PAW-LCAO: the `lda-fr/` and `pbe-fr/` sets

The datasets in `lda-fr/` and `pbe-fr/` carry **spin-orbit coupling** on top
of the scalar-relativistic construction of `lda-sr/` and `pbe-sr/`. Sections
1 to 6 describe the construction, which is the same for both functionals;
section 7 describes `pbe-fr/`. This note describes how the
spin-orbit term is built, how a Mandacaru calculation uses it, what the
generator checks, and what can go wrong. Code references are to the Mandacaru
repository (`src/mandacaru/...`).

## 1. The reference atom

Every dataset starts from an all-electron atom. For `lda-fr/` it is solved
with the radial **Dirac** equation (`basis/relativity.py`,
`basis/atomic_solver.py`, `solve_atom(relativity="dirac")`). The small
component is eliminated exactly, not perturbatively. With

$$
M(r) = 1 + \frac{\varepsilon - V(r)}{2c^2},
$$

the large component $P = rg$ obeys one second-order equation, and its only
$j$-dependent term, $\kappa M' P/(M r)$, is the spin-orbit interaction.
Each $(n, l, \kappa)$ is solved separately:

- $\kappa = l$ is $j = l - 1/2$;
- $\kappa = -(l+1)$ is $j = l + 1/2$.

The occupations are divided between the two $j$ in proportion to $2j+1$.
Replacing $\kappa$ by $-1$ for every $l$ gives the scalar-relativistic
(Koelling–Harmon) equation, which is what `lda-sr/` is built from.

The solver reproduces the hydrogenic Dirac spectrum in closed form. It also
reproduces measured atomic splittings: Ar 3p gives 0.1786 eV against a
measured 0.178, and Kr 4p gives 0.648 against 0.666.

**Exchange.** A relativistic reference atom uses the relativistic correction
to LDA exchange of MacDonald and Vosko. With $\beta = k_F/c$, the exchange
energy per electron is scaled by

$$
\Phi_E = 1 - \tfrac{3}{2}F^2,\qquad
F = \frac{\sqrt{1+\beta^2}}{\beta} - \frac{\mathrm{asinh}(\beta)}{\beta^2},
$$

and the potential by the exact derivative of that energy
(`basis/xc.py`, `relativistic_exchange_factors`). The correction acts in:

- the reference atom;
- the unscreening of the local potential;
- the one-center energies.

A dataset records that it used the correction (`relativistic_exchange`), so
it is rescreened the same way it was unscreened. The correction is large in
the core (Au +41.5 Ha, Bi +48.3 Ha in total energy) and small for valence
levels (6p ±0.1 mHa, 5d about −2 mHa). It changes spin-orbit splittings by
at most about 1 %.

## 2. The spin-orbit term in a dataset

`pseudopotentials/paw.py`, `j_resolved_spin_orbit`.

**One branch per $j$.** For every channel with $l \ge 1$, the generator
builds two complete PAW branches, one for each $j = l \mp 1/2$. Each branch
has:

- its reference waves from the Dirac atom, with that branch's $\kappa$;
- its own smooth partial waves, matched at the channel's cutoff;
- its own projectors and one-center Hamiltonian matrices, in the same
  screened local potential as the scalar channel.

**Two terms reproduce both branches.** Any $j$-dependent separable operator
can be written as a $j$-average plus a coupling term, because
$\mathbf{L}\cdot\mathbf{S}$ takes the value $l/2$ on $j = l+1/2$ and
$-(l+1)/2$ on $j = l-1/2$:

$$
V_j = V^{\mathrm{avg}} + V^{\mathrm{SO}}\,\mathbf{L}\cdot\mathbf{S},\qquad
V^{\mathrm{avg}} = \frac{(l+1)\,V_{l+1/2} + l\,V_{l-1/2}}{2l+1},\qquad
V^{\mathrm{SO}} = \frac{2}{2l+1}\left(V_{l+1/2} - V_{l-1/2}\right).
$$

Both terms are stored on the **union** of the two branches' projectors, the
$j = l-1/2$ projectors first. On that union they are block diagonal:

- $V^{\mathrm{avg}}$ has the blocks $\frac{l}{2l+1}V_{l-1/2}$ and
  $\frac{l+1}{2l+1}V_{l+1/2}$;
- $V^{\mathrm{SO}}$ has the blocks $-\frac{2}{2l+1}V_{l-1/2}$ and
  $+\frac{2}{2l+1}V_{l+1/2}$.

Each $j$ is therefore recovered exactly, not to first order. The same
two-term storage is how ONCVPSP (Hamann) stores its $j$ channels
(`oncv.py`, `_combine_j_channels`).

**Stored fields.** Each $l$ entry of the file's `spin_orbit` field holds:

- the union projectors;
- the unscreened and screened average and coupling matrices;
- the same split of the branches' overlap corrections;
- the reference energies of each branch;
- each branch's level error against the Dirac atom.

**The branches are unitary** (norm deficit 0). Their overlap corrections are
about 1e-4, so the overlap operator $S = 1 + \sum |\tilde p\rangle q \langle\tilde p|$
stays spin-free. A $j$-dependent $q$ would give $S$ an
$\mathbf{L}\cdot\mathbf{S}$ structure. Every consumer of $S$ would then have
to become spin-aware: the Löwdin orthogonalization, the basis, the forces.

**Why not a first-order term.** The simpler alternative evaluates
$\int \xi(r)\,\varphi_i\varphi_j$ on the scalar partial waves, with
$\xi = (1/2c^2M^2r)\,dV/dr$. That fails where it matters: it is 7 % short
for 5p splittings and 19–20 % short for 6p. A single $j$-averaged partial
wave cannot represent the contraction of p1/2 relative to p3/2 (6p1/2 has
27–29 % more of its norm inside the sphere than 6p3/2). With the $j$-resolved
branches the 6p splitting of Tl, Pb and Bi is within 0.01 %. The loader
refuses a file carrying a first-order term.

## 3. How a calculation uses it

Select the folder, and a spin-orbit pool or the generalized mean field:

```python
atoms.calc = Mandacaru(method="adapt-vqe", pool="spin-orbit",
                       basis={"name": "PAW-LCAO", "size": "SZ"},
                       paw_dataset="lda-fr", grid_spacing=0.25)
```

**The Hamiltonian.** The $\mathbf{L}\cdot\mathbf{S}$ coupling acts through
its own union projectors (`paw_spin_orbit_projectors`) and is added to the
scalar Hamiltonian (`core/spin_orbit.py`, `MolecularIntegrals`). Its matrix
is complex and Hermitian, with genuine α–β blocks, and it is transformed to
the Löwdin-orthonormal basis like every other one-body term. The scalar
channel stands in for the $j$-average $V^{\mathrm{avg}}$, so:

- the overlap, the compensation charges, the basis functions and the scalar
  forces are those of `lda-sr/`;
- the price is the gap between the scalar-relativistic level and the Dirac
  atom's $(2j+1)$-weighted average (section 5).

**Consequences for the quantum problem.** Spin-orbit coupling breaks $S_z$;
only the particle number $N$ is conserved. Therefore:

- **Orbitals:** the molecular orbitals are **generalized Hartree–Fock
  spinors** (`method="ghf"`, `hartree_fock.GHF`): one determinant of complex
  two-component spinors. Every spin-orbit run starts from the GHF
  determinant.
- **Kramers pairs:** a closed shell comes out paired without being forced to
  (within 4e-15 Ha for Pb); pairing is measured, not imposed.
- **Frozen core and active space:** both are chosen over **Kramers pairs**,
  ranked by spinor energy. MP2 and threshold selection are refused.
- **Pools and reductions:** only `pool="spin-orbit"` (with spin-flip and
  complex generators) and the `"ghf"` mean field accept the term. The
  $(n_\alpha, n_\beta)$ sector, the parity reduction and the $S_z$ tapering
  symmetries do not exist here.
- **Forces:** they include the derivative of the spin-orbit term. The atom's
  union projectors move (Hellmann–Feynman) and the basis functions move
  (Pulay). On HI the spin-orbit part is right to about 1 % of itself.
- **Densities and charges:** they are computed from the full spin-orbital
  reduced density matrices.

## 4. What the generator checks

Every dataset in `lda-fr/` passed these checks, or carries its defect in
the file (a `GhostStateWarning` on loading):

| Check | What is compared | Requirement |
|---|---|---|
| Reference atom | the all-electron Dirac self-consistent field | converged, or the element is refused |
| Ghost states | each channel's spectrum with and without projectors, against the all-electron atom | no spurious bound state |
| Scattering | the logarithmic-derivative phase against the all-electron atom at the projector radius | 0.05 rad over ε ± 0.5 Ha, 0.3 rad over ε ± 1 Ha |
| Intruding 1s | the `s` projectors' miss of a hydrogen 1s entering the sphere | below 1 |
| Per-$j$ levels | each branch's lowest level against the Dirac atom | reported (`level_errors`) |
| Splittings | the splitting of every $l \ge 1$ channel against the Dirac atom | reported |

At the last audit:

- **Per-$j$ levels:** the error is 0.08 mHa at the median, and 94 of 106
  bound spin-orbit channels are within 1 mHa.
- **Splittings:** the error is 0.07 % at the median and 2.5 % at worst.

## 5. What can go wrong

**Missing states: an unbound $j$ branch.** A branch is built from a bound
state of the Dirac atom. If one $j$ of a valence shell is not bound, that
channel has no reference, and the generator refuses the element rather than
invent one. Protactinium is missing from `lda-fr/` for this reason: one
of its 5f branches is unbound in the Dirac atom (ε ≈ +0.001 to +0.005 Ha).
Building it needs a reference at a chosen energy for an unbound branch, or a
different reference configuration. `pbe-fr/` has Pa, built on the second
of these (section 7).

**Ghost states.** The two branches share the scalar channel's local
potential. A local potential deep enough to bind a spurious state produces
the ghost in both $j$. The repairs are the scalar ones:

- raise the local potential at the origin;
- rebalance the cutoffs;
- add Bessel functions;
- use norm deficit 0.

Each repair is the exact parameter set a scan found, kept per dataset in
`paw.py`'s `DEFAULT_*_BY_DATASET` tables. No dataset in the set holds a
ghost. A repair must be the whole parameter set: changing one channel's
cutoff and leaving the rest to the generic rule has produced new ghosts and
phase errors in other channels.

**Scattering errors in compact shells.** The lanthanide 4f is compact
($r_c$ well under 1 Bohr) and its phase misses the near window by 0.05–0.09
rad. This is the same limit as in `lda-sr/`, not a spin-orbit error. The
flagged elements are listed in the [README](README.md).

**Level errors of compact semicore shells.** In these shells both $j$ levels
miss the Dirac atom by about the same amount:

- 3d of Se–Kr: 1.1–2.4 mHa;
- 4d of I and Xe: 1.0–1.3 mHa;
- 4f of W–Hg: 1.2–4.9 mHa.

Because both $j$ move together, this is the scalar channel's error, and the
splittings stay within 2.5 %.

**The intruding 1s.** An `s` channel whose projectors cannot represent a
neighbor's 1s gives wrong molecules even when every atomic check passes. The
check is part of generation. The elements still above the limit (misses
1.0–1.8, mostly lanthanides) are listed in the README.

**The scalar channel as the $j$-average.** The molecular Hamiltonian uses
the scalar-relativistic channel in place of $V^{\mathrm{avg}}$. Both $j$
levels therefore shift by the difference between the scalar level and the
Dirac $(2j+1)$ average:

- about 0.001 mHa for O 2p;
- 1.5–1.9 mHa for 5p;
- 3.7–6.8 mHa for the 6p of Tl–Bi.

The splitting is not exactly unaffected: on the atom's own radial operator
(scalar channel plus the coupling term, as the molecule assembles it) the
6p splitting of Pb and Bi comes out 2.0–2.4 % short and the 6p1/2 level up
to 2.6 mHa high (I 5p: 0.4 %; the d shells within 0.15 %). The stored
average through the union projectors, with a spin-free metric, reproduces
both branches to 0.02 mHa and the splitting to 0.02 %; it is not
implemented, because the averaged channel would also need its own partial
waves (compensation multipoles, one-center terms), which the files do not
store.

**The grid.**

- *Compact semicore d shells:* a real-space grid mixes their $m$
  components. At $h = 0.25$ Å, $[h_{SO}, J_z]$ is broken in the I 4d and
  Bi 5d blocks, while the valence p blocks conserve $J_z$ to 1e-8.
  Spin-orbit on such a shell needs a finer grid; check $J_z$ before
  trusting it.
- *Near the nucleus:* a p1/2 state goes as $r^\gamma$ with
  $\gamma = \sqrt{\kappa^2 - (Z\alpha)^2} < 1$ (0.80 for bismuth). The
  uniform radial grid of the generators resolves this only approximately.
  The per-$j$ audit compares against a Dirac atom solved on the same grid,
  so it does not measure that error.

**The overlap is spin-free only approximately.** The unitary branches leave
overlap corrections of about 1e-4. Their $j$-difference is stored but not
used, so the metric is the scalar one.

**Cost.** Without $S_z$ there is no $(n_\alpha, n_\beta)$ sector, no parity
reduction and fewer tapering symmetries. The spin-orbit pool is larger than
a spin-conserving one, and the Hamiltonian is complex. The same molecule
therefore needs more qubits and deeper circuits than with `lda-sr/`. Budget
this before any hardware run.

**Not supported.** These are refused rather than approximated:

- spin-orbit forces with a truncated virtual space (deleted Kramers pairs);
- the orbital-response residual of spinor forces, which is not measured.

## 6. What is not implemented

- **Fully relativistic PAW.** A $j$-resolved augmentation sphere, as in Dal
  Corso's fully relativistic PAW, would make the overlap, the compensation
  charges and the one-center densities $j$-dependent. It would also need
  two-component spinor basis functions, a spin-coupled Löwdin
  orthogonalization and spin structure in every force term. The
  $j$-resolved Hamiltonian term above reproduces the valence per-$j$ levels
  to 0.08 mHa (median), so this has not been needed for valence shells.
- **Non-collinear magnetism.** The frozen one-center densities are
  spin-free.
- **Spin-orbit Kohn-Sham.** `method="dft"` has no spin-orbit term in its
  Kohn-Sham operator and refuses a spin-orbit dataset (`NotImplementedError`),
  so `pbe-fr/` and `lda-fr/` serve `method="ghf"`, the spin-orbit pools and
  the other correlated methods; a Kohn-Sham calculation reads `pbe-sr/` or
  `lda-sr/`.

## 7. The PBE set: `pbe-fr/`

`pbe-fr/` is the Dirac construction of sections 1 to 6 with the PBE
functional of [PBE.md](PBE.md): the reference atom is the radial Dirac atom
in PBE with the relativistic exchange factor of section 1, and the
unscreening and the one-center energies are PBE. It holds all 92 elements,
H to U, 282 MB. The generator's per-dataset tables are keyed by element and
functional, so the repairs of `pbe-sr/` apply to it unchanged.

**Protactinium.** `pbe-fr/` has Pa, which `lda-fr/` lacks. The scalar
reference configuration of Pa, $5f^3 7s^2$, leaves one 5f branch of the Dirac
atom unbound. The check that the relaxed configuration binds
every valence orbital (`bound_valence`, `basis/atomic_solver.py`) now looks
at each $j$ branch of a Dirac atom, not only at the $(2j+1)$ average, and
the relaxed configuration of the Dirac atom becomes $5f^2 6d^1 7s^2$, the observed ground
configuration. All its branches are bound (5f −0.089 and −0.064 Ha, 6d −0.094
and −0.077 Ha, 7s −0.151 Ha) and the dataset has `s`, `d` and `f` channels
(the scalar `pbe-sr/` file holds `s` and `f`). Pa's reference configuration
therefore differs between `pbe-sr/` and `pbe-fr/`. The same change would let
`lda-fr/` build Pa; that set has not been rebuilt.

**Checks.** Every file was checked when it was generated (section 4) and
again from the file by an independent audit (92 files, 128 bound
$j$-resolved channels): the reference atom solved
again, the scalar checks of the generator, and for every $l \ge 1$ channel
the spectrum and the scattering phase of each $j$ branch against the Dirac
atom (the audit script is not part of the generator; it reads the stored
union projectors and the $j$ average as section 2 describes).

| Check | Result |
|---|---|
| Reference atom | converged for every element |
| Ghost states | none in any scalar channel, none in either $j$ branch, in 92 files |
| Missing states | every occupied reference level of both branches is reproduced: 111 of 128 bound channels within 1 mHa (median 0.10 mHa), the rest 1.0–4.8 mHa (compact 3d of Se–Kr, 4d of I and Xe, 4f of Yb and Ta–Hg, as in `lda-fr/`), except the empty bound 3d of Ca and the empty 4d of Rb and Sr, which the dataset places 67–84 mHa too high (also in `lda-fr/` and `pbe-sr/` for Sr; Ca and Rb not compared) |
| Scalar scattering | within 0.05 rad over ε ± 0.5 Ha and 0.3 rad over ε ± 1 Ha, and the channels without projectors within 0.2 rad, for 68 files; the other 24 are flagged |
| Per-$j$ scattering | the median over 128 channels is 0.002 rad; the largest unflagged: Kr 3d 0.054, Br 3d 0.051, Re 4f 0.048, Se 3d 0.041; the flagged 4f of Pr–Yb 0.11–0.60 rad (below) |
| Splittings | 123 channels, median error 0.07 %; light-element $p$ (B, C, Ne) 1.3–1.6 % of a splitting of a few meV; the 4f of the lanthanides 1.8–3.4 % (Gd 3.4, Tm and Yb 3.3, Er 3.2, Pm 3.0), Ta 4f 2.4 %, Cr 3d 2.2 % |
| Loading | all 92 files load, record PBE, `relativity="dirac"` and the relativistic exchange |

The `lda-fr/` files for Ho and Yb, audited the same way, give per-$j$ phases
of 0.37 and 0.53 rad (pbe-fr: 0.41 and 0.60) and Kr 3d 0.052 (0.054):
the $j$ phase of the lanthanide 4f is a limit of the construction, not of
the functional. `lda-fr/` Yb's splitting error is 3.9 %, not the 2.5 % of
section 4.

**Splittings against the all-electron atom** (eV; the all-electron value is
the Dirac atom's, in the set's own functional):

| shell | PBE all-electron | `pbe-fr/` | LDA all-electron | `lda-fr/` |
|---|---|---|---|---|
| Ar 3p | 0.1752 | 0.1751 | 0.1789 | 0.1788 |
| Kr 4p | 0.6335 | 0.6335 | 0.6497 | 0.6498 |
| Xe 5p | 1.2243 | 1.2244 | 1.2608 | 1.2610 |
| Au 5d | 1.6114 | 1.6119 | 1.6241 | 1.6247 |
| Pb 6p | 1.3579 | 1.3581 | 1.4028 | 1.4030 |
| Bi 6p | 1.8360 | 1.8362 | 1.8795 | 1.8797 |
| U 5f | 0.8226 | 0.8205 | 0.7274 | 0.7257 |

PBE splittings are 1–4 % smaller than LDA ones for the p shells and the file
reproduces its own atom to 0.07 % or better, as in `lda-fr/` (U 5f: −0.25 %).
These are splittings of the file's own radial operator, as in section 2; the
molecular operator uses the scalar channel for the $j$ average (section 5).

**Flagged elements.** 24 of 92 files record a defect and warn when loaded;
none is flagged for a ghost (`pbe-sr/`: 27, `lda-fr/`: 21). The list and the
defects are in the [README](README.md). All come from the same structural
limit as in the other sets: no projectors under the 4f (empty $p$ and $d$ of
Ce–Yb, 0.3–1.4 rad), none for the $p$ of Ag, Cd, Cs, Ba, Fr, Hf, W and Ir
(0.21–0.67 rad), Lu's intruding 1s (miss 2.26), and the 4f phase of the
lanthanides (scalar 0.05–0.23 rad). Compared with `pbe-sr/`: Pa and U are clean
here (they have a $d$ channel, and Pa a different reference, above), Ra is
clean (its unprojected $p$ is 0.16 rad, inside the 0.2 limit), Th is clean
after the replacement below, and La's unprojected $f$ is 0.43 rad rather
than 1.01.

**Choosing among constructions.** The flagged elements were compared with the
construction of the same element in `pbe-sr/` (its channel cutoffs, local
radius, local raise and norm deficit, rebuilt as a Dirac atom with the
generator's own checks, each evaluated against a freshly solved atom and
ranked as the generator ranks its repairs). Two files changed. Thorium: the
default Dirac build was flagged (unprojected $p$ 0.55 rad, $f$ 0.056 rad);
the `pbe-sr/` construction with ten Bessel functions in the $s$ channel is
clean (unprojected $p$ 0.11 rad, intruding-1s miss 0.99), and is the
installed file. Lutetium: the intruding-1s miss falls from 2.52 to 2.26.
For Ac, Gd, Nd and Pm the `pbe-sr/` construction is no better, for Cs, Ba,
Fr, Hf, W and Pr the two agree to 3 %, and those files were kept as
generated. The Th and Lu files are therefore not reproduced by a default
`mandacaru-build`; a default rebuild gives the flagged thorium.

**Storage.** `mandacaru-build` writes a file only when its bound levels read
back within 1e-5 Ha of the generated dataset (the checks table of
[README.md](README.md)). At the library stride (every fourth point) the
check refused nine of the 92 files, whose compact channels lose accuracy on
the thinned grid. They are stored at a finer stride, the coarsest that
passes (Eu, Ga, Yb, B every third point; Ho, Mn, Pb, V every second; Ne
every point), and are larger than the `lda-fr/` files (Pb 15.5 MB against 7.8, Ho 7.0
against 3.5, Ne 1.9 against 0.5). `lda-fr/` and `pbe-sr/` were built before
the check existed (the README dates it) and were not tested against it
here.

**Validation.** Spin-orbit Kohn-Sham does not exist, so `pbe-fr/` was
validated through `method="ghf"` against its scalar sibling `pbe-sr/`
(`method="rhf"`), with `lda-fr/` against `lda-sr/` as the control, at
SZ or DZP with $h$ = 0.25 Å. The shift between a Dirac and a scalar dataset
is the same in both functionals: H₂O +0.075 eV and N₂ +0.134 eV in PBE
(+0.075 and +0.135 eV in LDA), Xe +1.845 eV (+1.808), HI +3.914 eV (+4.013),
Hg +16.74 eV (+16.71). Bond lengths agree with `pbe-sr/`: N₂ 1.1252 Å
against 1.1249 and the O–H bond of water 0.9788 against 0.9786 Å (`lda-fr/`
against `lda-sr/`: 1.1233 against 1.1230 and 0.9754 against 0.9753 Å). The Pb atom (6p²) gains 2.238 eV
from the 6p1/2² determinant over the closed 6p² one, `lda-fr/` 2.296 eV,
a ratio of 0.975 against 0.968 for their all-electron 6p splittings
(1.358 and 1.403 eV); the Kramers pairing is 1e-14 Ha in both.

## Building the set

```bash
mandacaru-build --pp PAW --dirac --xc LDA --all --workers 7 --check --ghosts flag --install
mandacaru-build --pp PAW --dirac --xc PBE --all --workers 6 --check --ghosts flag --install
```

An LDA set takes about three hours on 7 cores: twice the scalar work per
$l \ge 1$, one branch per $j$. `pbe-fr/` took 77 worker-hours (15.5 hours on
six cores, one element per process with one thread each, 1.4 GB at most per
worker): thorium 4.0 hours, iridium 3.8, tungsten 3.8, lutetium 3.6, hafnium
3.5, protactinium 3.5, uranium 3.0, hydrogen 9 seconds. Build on an otherwise
idle machine; a worker killed under memory pressure leaves the process pool
waiting forever.

## References

- P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994) — the projector
  augmented-wave method.
- D. D. Koelling and B. N. Harmon, *J. Phys. C* **10**, 3107 (1977) — the
  scalar-relativistic radial equation.
- A. H. MacDonald and S. H. Vosko, *J. Phys. C* **12**, 2977 (1979) — the
  relativistic correction to LDA exchange.
- L. Kleinman, *Phys. Rev. B* **21**, 2630 (1980) — relativistic
  norm-conserving pseudopotentials, the $j$-average plus
  $\mathbf{L}\cdot\mathbf{S}$ form.
- G. B. Bachelet and M. Schlüter, *Phys. Rev. B* **25**, 2103 (1982) —
  $j$-averaged and spin-orbit pseudopotentials.
- D. R. Hamann, *Phys. Rev. B* **88**, 085117 (2013) — optimized
  norm-conserving Vanderbilt pseudopotentials, whose $j$ channels are stored
  the same way.
- A. Dal Corso and A. Mosca Conte, *Phys. Rev. B* **71**, 115106 (2005) —
  spin-orbit coupling with ultrasoft pseudopotentials.
- A. Dal Corso, *Phys. Rev. B* **82**, 075116 (2010) — fully relativistic
  PAW.
- C. A. Jiménez-Hoyos, T. M. Henderson and G. E. Scuseria,
  *J. Chem. Theory Comput.* **7**, 2667 (2011) — generalized Hartree–Fock.
