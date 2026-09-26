# Dirac-relativistic PAW-LCAO: analysis and plan

*Draft, 2026-09-26.* What it would take for the PAW-LCAO datasets in this
repository, and the Mandacaru calculations that use them, to carry
**spin-orbit coupling** (SOC), not only the scalar-relativistic
(Koelling-Harmon) correction the `lda/` and `pbe/` sets carry today. Nothing
here is implemented yet beyond what section 1 lists. Code references are to
the Mandacaru repository (`src/mandacaru/...`).

---

## 1. Where things stand

### 1.1 What already exists

**The radial Dirac atom** (`basis/relativity.py`, `basis/atomic_solver.py`).
The small component is eliminated exactly, not perturbatively. With
`M = 1 + (eps - V)/2c^2`, the large component `P = rg` obeys one second-order
equation in which the only j-dependent term is `kappa M'P/(Mr)`, and that term
*is* the spin-orbit interaction. `kappa = -1`, the (2j+1) average for every
`l`, is the scalar (Koelling-Harmon) equation. `solve_atom(relativity="dirac")`
solves every `(n, l, kappa)` separately, with jj-proportional occupations
(`split_configuration`), and stores `orbitals_j` / `eigenvalues_j` beside the
(2j+1) averages. It is validated against the closed-form hydrogenic Dirac
spectrum and reproduces measured splittings (Ar 3p 0.1786 eV against 0.178;
Kr 4p 0.648 against 0.666).

**ONCVPSP with SOC** (`pseudopotentials/oncv.py`, `_combine_j_channels`).
There, `relativity="dirac"` builds one full channel per j and stores their
union: 4 projectors per `l`, a (2j+1)-averaged coupling as the ordinary
channel, and `spin_orbit[l]` as the difference. Because `L.S` takes one
value per j, the pair reproduces each j exactly.

**PAW-LCAO with a spin-orbit term** (`pseudopotentials/paw.py`,
`spin_orbit_blocks`). `generate_paw(relativity="dirac")` is **scalar
partial waves plus a first-order one-center term**:

    D_SO_ij = int_0^rc [ xi(r) phi_i phi_j - xi~(r) phi~_i phi~_j ] r^2 dr,
    xi = (1 / 2c^2 M^2 r) dV/dr

One partial-wave set per `l` is kept (the j average), so `q`, `Delta T`, the
compensation charges and the overlap `S` are unchanged. The design reason,
from `paw.py`, is still valid: a j-dependent `q` would give the *metric*
`S = 1 + sum |p~> q <p~|` an `L.S` structure, and every consumer of `S`
(Loewdin orthogonalization above all) would have to know about spin. Only
the reference atom is solved with Dirac. Oxygen 2p: `D_SO` gives 0.03647 eV
against the Dirac atom's 0.03674 eV (99 %); the hydrogenic 2p fine structure
is reproduced at ratio 1.000000 for Z = 1, 2, 5.

**The build command.** `mandacaru-build --pp PAW --dirac` exists (PAW, UPAW,
ONCV; refused for NCPP, which has one projector per `l`).

**The Hamiltonian side** (`core/spin_orbit.py`, `core/hamiltonian.py`).
`ls_matrix(l)` and `spin_orbit_one_body` assemble a complex Hermitian
`(2M, 2M)` spin-orbital matrix with genuine alpha-beta blocks from the
projector blocks, `MolecularIntegrals(spin_orbit_coupling=...)` accepts them
for projector families (PAW-LCAO, UPAW-LCAO, ONCVPSP; all-electron and
Gaussian bases have nothing for the term to act on), and `L.S` is checked by
its spectrum, Hermiticity, `[L.S, J_z] = 0` and `[L.S, S_z] != 0`.

### 1.2 What does not work

**There is no working end-to-end path.** Both standard Hamiltonian builders
(`algorithms/_hamiltonian_from_atoms.py` and `build_valence_hamiltonian` in
`pseudopotentials/families.py`) call `molecular_hamiltonian(mo_basis=True,
num_particles=(n_alpha, n_beta))`, and `_refuse_spin_orbit_without_sz`
refuses exactly that: with SOC, `S_z` is not conserved and there is no
`(n_alpha, n_beta)` sector. A calculation that loaded a Dirac dataset would
stop with `NotImplementedError` at Hamiltonian assembly. The guard is
deliberate: the downstream machinery assumes a fixed `(n_alpha, n_beta)`.

**No Dirac dataset ships.** All 92 files in `lda/` and in `pbe/` load as
`relativity="scalar"`, `has_spin_orbit=False`.

**The PAW term is first order.** It is evaluated on j-averaged waves with
`M` at `energy=0`. That is accurate where SOC is a small perturbation of the
radial function (oxygen: 1 %), and increasingly wrong where it is not: the
6p of Tl-Bi and Po-Rn, the 5d of Au and Hg, the 5p of I and Xe, and every 5f.
There the p1/2 orbital contracts relative to p3/2 (a second-order effect, and
large), and a single j-averaged partial wave cannot represent both radial
shapes. No test measures the error for a heavy element yet.

**The uniform radial grid** cannot converge a point-nucleus `s` or `p1/2`
state at second order: they go as `r^gamma` with
`gamma = sqrt(kappa^2 - (Z alpha)^2) < 1`, so gamma = 0.80 for bismuth's p1/2.
The log grid (`grid="log"`) resolves it but is not used by the generators,
which need waves that satisfy the uniform-grid equation to fourth order.
`D_SO` integrates `xi ~ 1/r^3` weighted by the waves near the nucleus, which is
exactly where the uniform grid is weakest.

**Tests.** `test/core/test_spin_orbit.py` uses synthetic projectors; nothing
drives a real Dirac dataset through `MolecularIntegrals`, then
`molecular_hamiltonian`, then an energy. `test_oncv.py` has no dedicated test
of the j-resolved ONCV recombination.

---

## 2. Two ways to make PAW-LCAO Dirac-relativistic

### Option A: keep the term, make it better (recommended first)

Keep one scalar partial-wave set per `l`, so the overlap `S` stays
spin-independent, and improve only the Hamiltonian term:

1. **j-resolved one-center Hamiltonian.** Solve the reference partial waves
   for both j (the Dirac atom already provides `phi_{l,j}`), and build the
   one-center *Hamiltonian* correction as the exact per-j difference,
   projected onto the scalar projectors:

       D_l,j = <phi_l,j | H_AE | phi_l,j> - <phi~ | H~ | phi~>  (per j)
       D_avg = [(2l+2) D_l,l+1/2 + 2l D_l,l-1/2] / (4l+2)
       D_SO  = 2/(2l+1) (D_l,l+1/2 - D_l,l-1/2)

   This is ONCV's exact two-term decomposition applied to `D` only. It
   captures the p1/2 contraction through the partial waves, not through a
   first-order `xi` integral.
2. **The energy argument of `M`.** Evaluate `xi` (or the per-j `D`) at the
   reference energies instead of at `energy=0`.
3. **Measure what the fixed `S` costs.** With the overlap j-averaged, the
   per-j norm deficit is not reproduced. Quantify it (per-j eigenvalue and
   scattering phase against the Dirac atom) before deciding whether option B
   is needed.

Advantages: the dataset layout, the overlap, the Loewdin orthogonalization,
the basis functions and forces all stay as they are; only `spin_orbit` in
the file changes meaning (and gains a version). Expected accuracy: exact per
j at the reference energies, degrading away from them like any separable
form.

### Option B: a j-resolved augmentation sphere (fully relativistic PAW)

Partial waves, projectors, `q`, `Delta T`, compensation charges and the
one-center densities per `(l, j)`, as in fully relativistic PAW
(A. Dal Corso, *Phys. Rev. B* **82**, 075116 (2010)) and fully relativistic
ultrasoft potentials (A. Dal Corso and A. Mosca Conte, *Phys. Rev. B* **71**,
115106 (2005)).

What it costs, and why it is a different construction rather than an option:

- `S` gains an `L.S` structure. The basis becomes two-component spinors, the
  generalized eigenproblem and the Loewdin orthogonalization become complex
  and spin-coupled, and every consumer of `S` must change.
- The one-center densities and compensation charges become spin densities
  (for a non-collinear magnetization, a 2x2 density matrix).
- The PAW-LCAO basis (numerical atomic orbitals from the pseudo atom) would
  need j-resolved or spinor orbitals.
- Force and stress expressions all gain spin structure.

Recommendation: do option A, measure (phase 3 below), and start option B only
if the measured per-j errors for 6p and 5f elements are unacceptable.

---

## 3. The Hamiltonian side: the real blocker

A correct dataset is useless until a calculation can use it. Everything below
is in Mandacaru, not in this repository, and is needed for either option.

1. **A register without an `(n_alpha, n_beta)` sector.** The builders must
   pass `num_particles` as a total `N` when SOC is on. Particle number is still
   conserved, so the Jordan-Wigner register and the variational loop are
   untouched. What assumes `S_z` and must be generalized or refused:
   - the particle-number sector reduction (`core/sector.py`): use the
     `N`-electron sector without the `S_z` split;
   - the parity two-qubit reduction (`core/mapping.py`), which reads
     `(n_alpha, n_beta)` directly: refuse, or reduce to one parity qubit;
   - Z2 tapering: the `S_z` parity symmetries disappear, and time reversal is
     antiunitary, so it is not a Z2 qubit symmetry. Rely on the symmetry
     finder (it reads symmetries off the Hamiltonian) and test that it finds
     fewer;
   - the excitation pools (`circuits/pools.py`) and UCCSD: add spin-flip
     excitations (alpha to beta singles and mixed doubles), and complex
     generators or both real and imaginary parts, since the Hamiltonian is
     complex Hermitian.
2. **A mean-field reference with SOC.** RHF and UHF refuse SOC
   (`algorithms/mean_field.py`). Needed: generalized (two-component) Hartree-
   Fock with a complex spinor Fock matrix, giving a spinor MO basis and a
   reference determinant (`TODO.md` 4.5). Kramers-restricted GHF is the
   natural choice for closed shells.
3. **Frozen core and active spaces over spinors** (`core/hamiltonian.py`,
   `algorithms/active_space.py`, `algorithms/mp2.py`): select Kramers pairs of
   spinors, not spatial orbitals.
4. **Observables.** Spin-summed spatial RDMs (`algorithms/pseudo_forces.py`,
   `algorithms/volumetric.py`, charges, cubes) must accept a spin-orbital RDM
   with alpha-beta blocks, and the spin-orbit term needs its own force
   contribution (the derivative of `D_SO` projector blocks with the atoms, as
   for `D`).
5. **Selecting the relativity at run time.** Datasets are chosen by folder
   (`library_directory(family, xc)` reads `<checkout>/<xc>/`). Add a relativity
   level to the layout and a basis option to select it (see section 4).
6. **Guards stay until each piece lands.** `_refuse_spin_orbit_without_sz`
   is the right behavior today; lift it piece by piece, with a test per piece.

---

## 4. This repository

- **Layout.** Keep `lda/` and `pbe/` as the scalar sets. Add `lda-dirac/`
  and `pbe-dirac/` (or `<xc>/dirac/`; decide once, in
  `pseudopotentials/environment.py`, together with the run-time option).
  A Dirac dataset carries the scalar part too, so it can serve a run that does
  not want SOC, but a separate folder keeps provenance explicit.
- **Format.** `spin_orbit` is already serialized as `{l: D_SO}`. Option A
  changes its content (exact per-j difference), so bump the payload version
  and record `spin_orbit_method` (`"first-order"` or `"j-resolved"`) so an old
  file is never read as a new one.
- **Build.** `mandacaru-build --pp PAW --dirac --xc LDA --all --workers 7
  --check --ghosts flag --install`, and the same for PBE. Expect up to twice
  the scalar build time (both j per `l`).
- **Checks per dataset**, beyond the current ghost and phase checks: each j's
  reference level and scattering phase against the Dirac atom (not only the
  j average), and the splitting of every channel with `l >= 1`.
- **README.** Document the folders, the method, and the flagged elements, as
  for the scalar sets.

---

## 5. Phased plan

| Phase | Work | Acceptance |
|---|---|---|
| 0. Baseline | Generate Dirac PAW (current first-order term) and Dirac ONCV for O, Ar, Kr, I, Xe, Au, Tl, Pb, Bi, U. Tabulate `D_SO` splittings against the Dirac atom's `spin_orbit_splitting(n, l)` and against ONCV's exact per-j result. Add the missing ONCV recombination test. | A measured table of the first-order error by element and channel. |
| 1. Option A | Per-j reference partial waves; exact per-j `D` difference; `M` at the reference energies; versioned payload with `spin_orbit_method`. | Per-j reference levels within 1 mHa of the Dirac atom and phases within the existing 0.05 / 0.3 rad windows, for the phase-0 elements. |
| 2. Hamiltonian, no reductions | Builders pass total `N` under SOC; sector without `S_z`; spin-flip and complex pool generators; parity reduction and tapering refused under SOC. | A heavy-atom and a molecular test (Tl or Bi atom; HI or TlH) reach the exact diagonalization of the same SOC Hamiltonian; `J_z` conserved. |
| 3. Measure option B's need | Compare option A per-j errors for 6p, 5d and 5f elements with a j-resolved reference. | A decision: A is sufficient, or B is scheduled. |
| 4. Mean field | Kramers-restricted GHF with complex spinors; spinor MO basis; `as_quantum_problem()` for SOC. | GHF reproduces RHF with SOC off; with SOC, the GHF energy is variationally below the scalar RHF. |
| 5. Reductions and observables | Frozen core and active spaces over Kramers pairs; spin-orbital RDMs in forces, densities and charges; the `D_SO` force term. | Finite-difference force check with SOC at `tol=1e-8`; active-space energy converges to the full-space one. |
| 6. Libraries | Build `lda-dirac/` and `pbe-dirac/` (92 each); audit against the Dirac atom; README. | No ghosts, per-j phases within windows or flagged, as for the scalar sets. |
| 7. Validation | Molecular splittings and bond lengths: HI, I2, TlH, Bi2, Au2, PbO against published relativistic references. | Documented agreement, and a guide section in `docs/source/guide/pseudopotentials.md`. |

Phases 0-1 touch only the generators and this repository. Phase 2 is where
SOC first becomes usable; it and phases 4-5 are Mandacaru work. Phase 3 is a
decision point, not code.

---

## 6. Risks and open questions

- **Accuracy of a fixed j-averaged overlap.** Option A reproduces the
  Hamiltonian per j but not the per-j norm. Phase 3 measures whether that
  matters; for 6p and 5f it may.
- **The uniform grid near the nucleus.** `D_SO` and any per-j `D` are
  dominated by `r -> 0`, where p1/2 goes as `r^gamma`. Moving the generators
  to the log grid is a larger change (`TODO.md`, core-vs-NIST offset); measure
  the grid convergence of `D_SO` first.
- **Qubit cost.** SOC removes the `S_z` sector and the parity reduction, and
  spin-flip excitations enlarge the pools: the same molecule needs more qubits
  and deeper circuits than without SOC. Budget this before any hardware run.
- **Complex Hamiltonians on hardware.** Measurement is unchanged (Pauli
  expectation values), but the ansatz must represent complex amplitudes.
- **Magnetism.** Non-collinear spin densities are out of scope unless option
  B is chosen; option A keeps the frozen one-center densities spin-free.

## References

- P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994) — the PAW method.
- D. D. Koelling and B. N. Harmon, *J. Phys. C* **10**, 3107 (1977) — the
  scalar-relativistic equation.
- G. B. Bachelet and M. Schlüter, *Phys. Rev. B* **25**, 2103 (1982) —
  j-averaged and spin-orbit pseudopotentials.
- A. Dal Corso and A. Mosca Conte, *Phys. Rev. B* **71**, 115106 (2005) —
  spin-orbit coupling with ultrasoft pseudopotentials.
- A. Dal Corso, *Phys. Rev. B* **82**, 075116 (2010) — fully relativistic PAW.
