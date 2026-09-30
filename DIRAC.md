# Fully relativistic PAW-LCAO: the `lda-dirac/` set

The datasets in `lda-dirac/` carry **spin-orbit coupling** on top of the
scalar-relativistic construction of `lda-sr/`. This note describes how the
spin-orbit term is built, how a Mandacaru calculation uses it, what the
generator checks, and what can go wrong. Code references are to the Mandacaru
repository (`src/mandacaru/...`).

## 1. The reference atom

Every dataset starts from an all-electron atom. For `lda-dirac/` it is solved
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
F = \frac{\sqrt{1+\beta^2}}{\beta} - \frac{\operatorname{asinh}\beta}{\beta^2},
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
                       directory="lda-dirac", h=0.25)
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

Every dataset in `lda-dirac/` passed these checks, or carries its defect in
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
invent one. Protactinium is missing from `lda-dirac/` for this reason: one
of its 5f branches is unbound in the Dirac atom (ε ≈ +0.001 to +0.005 Ha).
Building it needs a reference at a chosen energy for an unbound branch, or a
different reference configuration.

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

The splitting is unaffected. Using the stored average through the union
projectors instead would remove the shift; that is not implemented.

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
- **Other functionals.** Only LDA ships. The generator builds PBE Dirac
  datasets (`--xc PBE`), but they are shelved ([PBE.md](PBE.md)).

## Building the set

```bash
mandacaru-build --pp PAW --dirac --xc LDA --all --workers 7 --check --ghosts flag --install
```

A full set takes about three hours on 7 cores: twice the scalar work per
$l \ge 1$, one branch per $j$. Build on an otherwise idle machine; a worker
killed under memory pressure leaves the process pool waiting forever.

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
