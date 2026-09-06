# Computation reference

The 38 scripts in this directory reproduce calculations in the three
manuscripts. Most use exact integer, rational and symbolic arithmetic;
the explicitly identified numerical period/wall calculations are exceptions.
They supplement the geometric and analytic proofs, not replace them.

## Source map and execution

- [Paper 1](../papers/01-global-constraints/): general constraints.
- [Paper 2](../papers/02-quantum-transition/): the compact quantum transition.
- [Paper 3](../papers/03-magnetic-sheaves/): magnetic-sheaf technical note.
- `check_x9_weighted_abel.sage` independently checks the weighted transverse
  finite-difference identities in main Appendix C: 672 exact polynomial
  comparisons through spectator degree six. The all-degrees uniform Abel
  bound is an analytic proof in the paper, not a conclusion of this cutoff.
  Run `sage -c "load('checks/check_x9_weighted_abel.sage')"`.

Run from the repository root:

```bash
make verify       # ten Python checks and all three LaTeX builds
make verify-full  # all 38 scripts and all three LaTeX builds
make papers       # all three LaTeX builds, without computations
```

The Python checks require SymPy 1.14.0 and mpmath 1.3.0, specified in
`requirements.txt`. SageMath 10.9 supplies these dependencies and the
additional algebra, NumPy and SciPy used by the Sage calculations. Individual
Python checks run as `sage -python checks/<name>.py`, or with a Python
environment containing the pinned requirements. Individual Sage checks run
as `sage -c "load('checks/<name>.sage')"`; this avoids generated `.sage.py`
files. No network access or external data is needed once dependencies exist.

The full runner uses total nef degree eight for the standalone mirror/GV
check and spectator degree four for the regular magnetic-period check.
The other scripts use their documented defaults. Several scripts independently
recompute shared geometry; the full suite can take several minutes.

## Paper 1: finite algebra

The first six Python scripts are `check_nodal_charge.py`,
`check_e2_gauging.py`, `check_e3_gauging.py`,
`check_e3_magnetic_quivers.py`, `check_xcirc_e3_profiles.py` and
`check_x9_spectrum.py`. They cover, respectively, the nodal charge kernel,
the $E_2$ invariant scheme and $X_9$ matrix, the $E_3$ branchwise invariant
scheme, the two $E_3$ magnetic-quiver Hilbert series, all four $X^{\circ}$
profiles, and the $X_9$ divisor/Hodge/multiplet calculation. The related
geometric papers and their computations are available separately in
[calabi-yau-smoothability](https://github.com/bwuebben/calabi-yau-smoothability).
The scripts here do not depend on that checkout.

## Geometry, fibrations and anomalies

- `check_x9_sixd_anomaly.py` — exact 6d anomaly bookkeeping of the genus-one
  fibration of X̂₉ (T=2, su(2) on H−E₃−E₄ with 10 fundamentals, H=219=123+20+76);
  solves the Abelian conditions within its stated half-integral charge bound
  (absolute charge at most two), finding one spectrum up to normalization
  (doublets 1, singlets 2, b = 22H − 6E₃ − 6E₄), and verifies
  that b·H/4 = 11/2 equals the Shioda-corrected height pairing on X̂₉.
  Paper 2 derives the center quotient from the geometric Shioda
  class; it does not claim uniqueness for an unrestricted anomaly problem.
- `x9_toric_base.sage` — the projection of the fan onto the last two coordinates
  is the fan of S₇ = Bl₂ℙ²; the flop of C₁ makes it a toric morphism (flat
  elliptic fibration over the 6d base); base dictionary for D₃…D₈.
- `check_x9_prepotential.sage` — resolved fan (29 unimodular cones), κ_IJK,
  c₂·D, linear equivalences, local-scale plane and its nef quadrant, norm
  formula, annihilators/charge rows of the three rigid limits, explicit limits
  along D₀ and D₁, the cubic (83,147,81,14), the formal quotient XY=β⁴
  in physical current operators.
- `x9_fibrations.sage` — K3 fibration (fibre D₅, lattice U⊕⟨−4⟩), the surfaces
  D₅|_X and D₆|_X, genus-one fibration over ℙ² with sections/multisections.
- `x9_weierstrass.sage` — Morrison–Park form and Weierstrass model of the
  genus-one fibration; discriminant data ((0,0,2) on L₀, (4,6,12) at q₃,q₄,
  3×III + 1×IV, 76 isolated I₂).

- `check_x9_smooth_side.py` — the smooth side X_t of the X₉ transition: neutral
  classes = H²(X_t) mod E, the full cubic form (83,49,27,14|0,0,1,0,−3,9) and
  c₂ = (146,80,−6) on (D₀,D₁,ν) from κ plus the local correction on V₇ = Bl_pℙ³
  (ν³: 8→9, c₂·ν: −4→−6), E-shift invariance, the ℙ² invariants, RR integrality.

## Integral charges, current exchange and sheaf counts

- `check_x9_rank_one_index.py` — centered Hilbert-scheme Poincaré products
  through q^12, independent colored-partition Euler counts, refined blowup
  ratio, eta vacuum shift, GRR charge formula and half-canonical flux
  lattices; verifies D8/D5/D6 Dirac pairings with the two nodes and their
  shared root. Imports the smooth-side and M5 anomaly checks.
  Standard-library only. The geometric identification of the compact
  moduli component is proved in the manuscript, not by this finite check.
  Neither the full compact-condensate index nor protection of the refined
  spin character is asserted.
- `check_x9_integral_screening.sage` — explicit unimodular divisor–curve
  pairing proving the full free resolved H² lattice; E2 surface restriction
  and index-two local charge lattice; center-compatible compact cocharacters
  and nilpotent shear; integral local gluing; primitive smooth basis and
  cubic, checked against the earlier rational basis and all RR residues
  modulo 12; local P² line charge −3 and compact D1² charge +1. This
  establishes screening of the local E0 Z3 line charge, not BPS multiplicities
  or the absence of torsion. Run with
  `sage -c "load('checks/check_x9_integral_screening.sage')"`.
- `check_x9_dressed_instantons.sage` — primitive compact current charges;
  physical-current elimination XY=β⁴; completeness of invariant generators;
  independent integral four-hyper quotient with Smith invariants (1,1,1);
  exact rational Hilbert-series identity and refined Molien coefficients
  through degree 16. Removing J± gives (X,Y,β²), whose reduction is a point.
  This checks a formal reduced ring, not a compact index. Run with
  `sage -c "load('checks/check_x9_dressed_instantons.sage')"`.
- `check_x9_heterotic_transition.sage` — independent flat intersection ring
  and flop correction; the E8 roots r_i=s_i−o_i−F_i and equality of their
  compact images; complete Shioda class and height; both zero-section/KK
  frames; `7/24 <= Vol(X)/(J.F)^3 <= 83/162`; affine charge relation;
  residual instanton bookkeeping (11,12;1). No BPS multiplicity is asserted.
- `check_x9_heterotic_k3_candidate.sage` — integral subdivision of the
  162-point Newton polytope by m4+2m1+m2=0; reflexive 30-point middle slice;
  fixed rational quartic with discriminant s² times a squarefree degree-22
  factor, coprime to f,g, and smooth infinity. Establishes a candidate K3,
  not a global semistable degeneration or a heterotic bundle.
- `check_x9_current_exchange.sage` — reconstructs the full six-vector Hodge
  metric and inverse from the resolved fan; contracts the two nodal charge
  columns; proves positivity on Λ ≥ 0 for the D0, D1 and D0+D1 rays; checks
  exact asymptotic matrices, the root-exchange coefficients 2/3, 1, 2/5,
  graviphoton subtraction, non-orthogonal integral basis covariance, and
  the difference between limiting kernels and kernels of the limit. This
  is a smooth-chamber two-derivative calculation, not a full SCFT potential.
- `check_x9_m5_string_anomaly.py` — independent surface restriction lattices
  of F1 and P2, compared with the smooth-side cubic; transverse SU(2) and
  gravitational anomalies; verifies that the jump is minus the anomaly
  of one unit-charge left-moving chiral lattice boson. Standard-library
  only. No full BPS index or interacting tensionless-string central charges.

Run from the project root (the `load` form avoids generated `.sage.py` files):

```bash
sage -c "load('checks/check_x9_heterotic_transition.sage')"
sage -c "load('checks/check_x9_heterotic_k3_candidate.sage')"
sage -c "load('checks/check_x9_current_exchange.sage')"
python3 checks/check_x9_m5_string_anomaly.py
python3 checks/check_x9_rank_one_index.py
```

## Quantum periods and magnetic stability

- `check_x9_d6_period.sage` — constructs the integral K-theory charge
  basis, including curvature/flux shifts, and checks its symplectic
  pairing and the surviving quadratic matrix. Verifies the smooth-side
  torsion exclusion and index-ten divisor-intersection curve sublattice;
  the full D2 basis is topological, not a list of BPS states. Checks the
  one-dimensional annihilator of the seven fixed periods, the Euler
  characteristics (-232,-240), and the four exact cubic Gamma tails
  giving -8*zeta(3). The proof fixes the last period by the polarization
  and reality of the marked continuation, and obtains the instanton
  prepotential shift -4*zeta(3)/(2*pi*i)^3. Run:
  `sage -c "load('checks/check_x9_d6_period.sage')"`.
- `check_x9_neutral_period_descent.sage` — verifies the primitive
  rank-three vanishing lattice, its rank-eight symplectic quotient and
  the index-four distinction from the combined monodromy image. Checks
  independence of the eight logarithmic polynomials (quadratic minor 10),
  matching cubic contractions and zero constant corrections, and replays
  the independent neutral/GV check. The all-orders proof additionally
  uses conifold Hodge-theoretic descent, the residue-preserving mutation
  and the full Frobenius basis. It identifies seven periods and all three
  instanton gradients. The D6 marking and prepotential constant require
  the separate calculation above. Run:
  `sage -c "load('checks/check_x9_neutral_period_descent.sage')"`.
- `check_x9_neutral_quantum_periods.sage` — constructs both smooth-side
  ambient subdivisions and their Mori cones; independently matches all
  cubic/c2 data and the plane-line charge. Compares every electric and
  neutral magnetic coefficient through spectator degree four (35 triples)
  with the fully summed old local degrees. Computes smooth-side genus-zero
  GV invariants through degree eight: 164 nonzero degree vectors, 31
  nonzero invariants, with three-component integrability and integrality.
  An independent local P2 calculation checks the diagonal sequence
  (3,-6,27,-192) and its mirror logarithm. This is a finite-order neutral
  comparison; the all-orders period argument is checked separately above.
  Neither check is a full integral period-frame proof. The smooth Mori-cone
  bound supplies the separate all-orders circle-suppression argument.
  Run: `sage -c "load('checks/check_x9_neutral_quantum_periods.sage')"`.
- `check_x9_circle_limit.sage` — reconstructs the intersection data;
  checks the del Pezzo restrictions and the extra primitive beta2 charge,
  rescaled cubic volume, graviphoton-orthogonal bulk-gauging matrix and its
  positive eigenvalues, and vanishing local Coulomb norm. Checks the
  bounded E structure-sheaf period's polynomial terms, surviving-plane
  volume, correlated mirror weights and KK mode-count scaling. The
  analytic period bounds and circle dictionary are supplied in the text.
  No uniform 5d Higgs action or magnetic BPS survival is asserted.
  Run: `sage -c "load('checks/check_x9_circle_limit.sage')"`.
- `check_x9_mirror_mutation.sage` — independently reconstructs the
  smooth-side mirror polytope, its polar and Hodge counts; verifies the
  full Laurent mutation, inverse and actual small-resolution chart,
  log-volume preservation, root coefficient gauge, and unimodular torus
  and parameter maps. An exact sample checks all 50 positive-dimensional
  faces, including the full torus. Replays 56 period coefficients under
  the printed all-orders index bijection. The comparison is generic and
  up to flops, not a selected crepant chamber or a 5d-limit calculation.
  Run: `sage -c "load('checks/check_x9_mirror_mutation.sage')"`.
- `check_x9_moving_probes.sage` — reconstructs the resolved intersection
  data via the prepotential check; counts the two toric sections of D5/D6;
  verifies the base-curve pairings, the explicit ample polarization and
  twelve proper-subunion stability numerators; checks the GRR seed
  charges, Euler/Behrend arithmetic and support-charge restrictions.
  It does not certify the reduced-pencil hypothesis needed for the
  conditional seed components, or a physical stability condition. Run:
  `sage -c "load('checks/check_x9_moving_probes.sage')"`.
- `check_x9_glued_probe.sage` — verifies the complete-intersection
  normal degrees of C2, its double-curve incidence, a two-term local
  presentation of I_(C2,T) and its dual, and the GRR charges from that
  resolution. Checks all six strict stability inequalities for each
  dual sheaf (and the original ideal), plus the P1 cohomology inputs
  to the rigidity proof. The two isolated +1 contributions do not use
  the reduced-pencil hypothesis. The check does not classify the full
  moduli or certify a physical stability chamber. Run:
  `sage -c "load('checks/check_x9_glued_probe.sage')"`.
- `check_x9_probe_transport.sage` — verifies coefficientwise positive
  stability bounds for J=aD0+bD1+epsilon*J*, the first B=-sD5 slope
  wall for each sheaf and its limiting rational function. Checks the
  spectral-flow transformation of the core and relative charge, and
  both magnetic extensions and their Dirac pairings. It does not
  establish exact quantum stability or a full DT/particle index. Run:
  `sage -c "load('checks/check_x9_probe_transport.sage')"`.

- `check_x9_probe_phase_wall.sage` — verifies the constituent volumes,
  charges and four strict quotient stability comparisons. Derives the
  finite-volume polynomial phase walls including the D0 terms, checks
  their slope limits, transversality and same-phase condition. The text
  separately proves spherical rigidity and asymptotic persistence of
  alignment under small quantum corrections. This is not a numerical
  quantum-period computation or a full wall-crossing index. Run:
  `sage -c "load('checks/check_x9_probe_phase_wall.sage')"`.
- `check_x9_mirror_periods.sage` — reconstructs the integral nef basis,
  the six-parameter GLSM system and the 24-cone coarser fan, checking
  that removal of the facet-interior ray leaves the same nef chamber.
  Computes rational Frobenius coefficients including negative Gamma
  arguments, the mirror map and genus-zero GV invariants. Independently
  verifies 168 Gamma recurrences, the quintic degree-one/two/three
  benchmark, integrality and agreement of all six magnetic equations,
  the electric indices and the printed quadratic magnetic periods.
  Default degree six: 923 nonzero lattice classes, 41 nonzero invariants.
  Degree eight: 3002 classes, 80 nonzero invariants. Run:
  `MIRROR_DEGREE=8 sage -c "load('checks/check_x9_mirror_periods.sage')"`.
  This is a large-volume series, not analytic continuation to the face.
- `check_x9_quantum_wall.sage` — computes the degree-eight mirror data,
  evaluates the full Li2 covers of every retained GV term and reproduces
  both numerical alignment tables at 50-digit working precision.
  Checks root residuals, transversality, the same-phase condition and
  the observed degree-six/eight agreement at J*. Comparisons of
  cutoffs are numerical diagnostics, not rigorous tail bounds or tests
  of constituent stability at the contraction. Run:
  `sage -c "load('checks/check_x9_quantum_wall.sage')"`.

- `check_x9_mirror_continuation.sage` — verifies the del Pezzo face
  lattice, quartic discriminant and Horn parametrization, the two
  integral symplectic nodal monodromies, their ranks and probe charge
  shifts, the D0-tower endpoint condition, and the Li2 monodromy sign
  by direct numerical continuation. It does not identify all compact
  discriminant components or stable vanishing branes. Run:
  `sage -c "load('checks/check_x9_mirror_continuation.sage')"`.
- `check_x9_local_elliptic_period.sage` — derives the local elliptic
  Gauss-Manin connection by exact differential reduction, checks the
  Gamma and hypergeometric coefficients through w^24, the positive-real
  discriminant locus and the I2 fiber, and continues the A-period to
  w=0.9999. Independent 50-digit hypergeometric quadrature evaluates
  the endpoint integral at w=1. Numerical agreement is not an interval
  bound; exact vanishing is proved separately by Chu--Vandermonde. This calculation
  is restricted to z3=z4=0, u=1, v=w, not the finite-volume compact wall.
  Run: `sage -c "load('checks/check_x9_local_elliptic_period.sage')"`.
- `check_x9_local_vanishing_charges.sage` — checks the factorized I2
  fiber, independent nodal smoothing parameters, root-mass splitting,
  the match to the conventional circle E2 curve, the terminating
  Vandermonde polynomial and Weyl symmetry through v^20, the integral
  compact charge lifts, saturation, and probe monodromies. The text's
  all-orders identity and nonsingular continuation argument prove the
  exact endpoint; this finite check does not replace them or establish
  finite-bulk stability. Run:
  `sage -c "load('checks/check_x9_local_vanishing_charges.sage')"`.
- `check_x9_finite_bulk_nodes.sage` — independent exact compact-mirror
  calculation, with both bulk parameters retained. Verifies the coefficient
  gauge, smooth toric charts, four-node factorization, Hessians, rank-three
  smoothing covectors, unit relative residue coefficients, and the four
  divisor intersections. At an exact rational sample the saturated
  critical ideal has length four, and all other toric strata are checked.
  No global integral period transport or magnetic stability is inferred.
  Run: `sage -c "load('checks/check_x9_finite_bulk_nodes.sage')"`.
- `check_x9_integral_vanishing_marking.sage` — checks the fixed-spectator
  first-derivative bounds and binomial convolution proof, resumming every
  local degree for all 56 spectator coefficients through degree five.
  Verifies the exact restricted fundamental period, the three zero electric
  corrections, unit limiting residue vectors, the full C1,C2,e1,e2
  zero-D0/D4/D6 charge embedding, saturation, and the commuting monodromies
  in the 14-dimensional symplectic lattice. Unimodular original and polar
  vertex spans exclude homological torsion. The printed convergence and
  Picard--Lefschetz proof establishes the all-orders result; finite tests
  alone do not. The marking is fixed on a sufficiently small, nonzero
  bulk branch. Regular magnetic periods and stability are not computed.
  Run: `sage -c "load('checks/check_x9_integral_vanishing_marking.sage')"`.
- `check_x9_regular_magnetic_periods.sage` — resums every local degree
  at fixed spectator degree for the six contracted magnetic periods.
  Checks bounded second-derivative sectors, the four rational-tail
  finite-difference formulas, exact sums in Q+Q*zeta(2), and independent
  radial Abel limits retaining all transverse weights. Also checks the
  regular-period constants and linear terms, constituent charges, and
  the exact negative phase sign on the real small-bulk branch. The
  printed proof establishes convergence; the finite tests alone do not.
  This does not establish constituent BPS survival or a full magnetic index.
  Run: `MAGNETIC_BULK_DEGREE=4 sage -c "load('checks/check_x9_regular_magnetic_periods.sage')"`
  for all 35 spectator coefficients through degree four; the default is two.

- `check_x9_magnetic_higgs.sage` — checks the exact light-emission phase
  identity on the final positive-imaginary electric ray, the five probe
  pairing/flux vectors, the saturated A3 flux image and its primitive
  magnetic kernel. The explicit E-shift identifies that kernel with
  the old neutral lattice modulo E. All ten cubic and three second-Chern
  pairings agree with the independent smooth-side surgery formula.
  This does not establish global BPS survival or compute string tensions.
  Run: `sage -c "load('checks/check_x9_magnetic_higgs.sage')"`.

Verified with Sage 10.9. If the optional local venv is absent, the Python
checks can also be run with `sage -python checks/<name>.py`; Sage supplies
the required SymPy dependency.
