# Compact Higgs transitions

Research manuscripts and reproducibility code by **Bernd Johannes Wuebben**.

The common question is how local Higgs branches change when their singular
regions belong to the same compact Calabi–Yau threefold. In M-theory, vector
fields of the compactification can gauge flavor symmetries shared by several
local sectors. Higgs directions that exist in separate noncompact models
can then be obstructed, or become available only through a joint transition.

## Read the papers

| | Manuscript | Read | Sources and guide |
|---|---|---|---|
| **1** | Global Constraints on Higgsing Localized Five-Dimensional Sectors in Compact M-Theory Vacua | [PDF · 38 pages](papers/01-global-constraints/paper.pdf) | [Paper 1](papers/01-global-constraints/) |
| **2** | Bulk Gauging and Quantum Couplings across a Compact Higgs Transition | [PDF · 54 pages](papers/02-quantum-transition/paper.pdf) | [Paper 2](papers/02-quantum-transition/) |
| **3** | Magnetic Sheaves and Stability Walls in a Compact Calabi–Yau Transition | [PDF · 24 pages](papers/03-magnetic-sheaves/paper.pdf) | [Technical note](papers/03-magnetic-sheaves/) |

Start with **paper 2** for the explicit compact transition and its quantum
vector couplings. **Paper 1** develops the general deformation/gauging
framework. **Paper 3** studies magnetic sheaves and partial stability
constraints; it does not determine the physical BPS spectrum at the endpoint.
Each paper's directory has its PDF, source, and a short guide to its results
and limitations.

## Paper 1 · Global constraints on local Higgs flows

**Question.** Which combinations of local Higgs directions survive a compact
embedding, and how are the physical gauging conditions related to the
geometry of smoothing singularities?

The paper treats conifold points together with rank-one interacting E₂
and E₃ sectors. Wrapped-M2 charges and compact divisor–curve pairings
determine a common Abelian charge and flavor-gauging matrix. Its kernel
correlates deformation coordinates belonging to otherwise separate sectors.

The main results are:

- An identification, on each chosen E₃ branch, of the invariant-coordinate
  image of the rigid triplet-moment-map zero locus with the kernel of the
  compact gauging matrix.
- For the stated admissible toric hypersurfaces and Hodge-theoretic
  hypotheses, an identification of that same kernel with the image of
  global first-order deformations. Combining this with the exact mixed
  smoothing theorem gives a necessary and sufficient existence criterion:
  the kernel must contain a vector outside every local discriminant.
  This includes both E₃ branches and arbitrary analytic smoothing arcs.
- Smoothness of the selected global deformation base when the relation
  space projects nontrivially to every selected E₃ summand, by the
  companion integrability theorem. Every tangent on that base integrates
  to a convergent deformation. Nodes and E₂ points need no additional
  projection hypothesis for this smoothness assertion.
- Explicit examples of three outcomes: an obstructed isolated E₂ sector;
  the X₉ model, in which two conifold sectors and one E₂ sector can move
  only together; and a compact model with thirty locally smoothable sectors
  but no trajectory that activates all thirty. A separate thirty-node,
  two-E₃ example has nonzero local deformations on mixed profiles but is
  forced onto an E₃ discriminant line: integrability does not imply that
  the resulting fibers are smooth.

The result concerns invariant deformation coordinates, not the entire
hypermultiplet space or its metric. The geometric existence and integrability
results are imported from
[CY paper 5](https://github.com/bwuebben/calabi-yau-smoothability/tree/main/paper5)
with their hypotheses; the contribution here identifies their relation map
with compact physical gauging. On the E₃ plane, all three discriminant
coefficients must be nonzero, not merely the plane vector. The selected-base
theorem does not settle all zero-projection E₃ profiles, unrestricted
nilpotent deformation schemes, or the finite-Planck scalar geometry.

## Paper 2 · An explicit transition and its quantum couplings

**Question.** What does the shared gauging do in an actual compact model,
does it survive decoupling gravity, and what are the vector-multiplet
couplings after the transition?

The model is X₉. If C₁ and C₂ are its two exceptional nodal curves and R
is the compact image of the E₂ flavor root, the integral relation
**R = C₁ + C₂** forces the three sectors to Higgs together. The constructed
smoothing changes the Hodge numbers from **(6,122) to (3,123)**.
Its existence also follows from the exact mixed smoothing criterion;
the explicit construction identifies the smooth phase needed for the
charge and quantum-period calculations.

The paper develops this example in three directions:

- **Integral charges and classical couplings.** It determines the primitive
  surviving charge lattice, its cubic and curvature couplings, and the
  local replacement of a contractible Hirzebruch surface F₁ by a rigid
  projective plane P². A wrapped-M5 anomaly comparison supplies an
  independent check of the change.
- **Gravity-decoupling limits.** It classifies the relevant nef directions
  with divergent total volume and bounded local Kähler scales. The shared
  E₂ flavor gauging becomes non-dynamical in every such limit. The conclusion
  includes vector mixing through the full inverse kinetic matrix and is
  restricted to these local-scale assumptions.
- **Quantum vector geometry.** At finite auxiliary-circle radius, four
  conifold states span a primitive rank-three electric lattice. A Laurent
  mutation and analytic continuation of all eight regular periods identify
  the surviving genus-zero prepotential with that of the smooth phase,
  including its integral normalization, quadratic ambiguity and constant.
  This is an all-orders period comparison on the marked branch, not only
  agreement of a finite set of expansion coefficients.

K3 and elliptic fibrations give a complementary interpretation: a dilaton
gauging, an anomaly-consistent six-dimensional U(2) phase, equal compact
images of two E-string roots, and a candidate heterotic K3. An exact volume
bound excludes the six-dimensional F-theory limit on the transition face.

The global heterotic bundle construction and its massless U(1), as well as
the finite-Planck hypermultiplet metric, are not supplied. Nor does the
four-state finite-radius description become a complete five-dimensional
theory: the Kaluza–Klein gap closes in the decompactification limit.

## Paper 3 · Magnetic sheaves and stability walls

**Question.** What can magnetic objects reveal about the same transition,
and which conclusions require stability information beyond charge and
period calculations?

This supporting technical note studies surface-supported sheaves in the
large-volume description and then uses paper 2's continued periods to
constrain their possible physical behavior near the quantum transition.

Its concrete results include:

- **Specified sheaf-moduli contributions.** Fixed-support Hilbert schemes
  on the rigid F₁ and P² surfaces give numerical Donaldson–Thomas counts.
  Their vacuum-normalized, fixed-flux oscillator series are η⁻⁴ and η⁻³,
  where η is the Dedekind eta function. The ratio is an oscillator-count
  comparison, not an identification of equal-charge states across the
  transition.
- **Two isolated stable sheaves.** On a reducible moving divisor, two
  explicitly constructed sheaves each contribute **+1** to the numerical
  invariant. Their Gieseker stability persists along ample polarizations
  approaching the contraction, without the additional hypothesis needed
  for the full-pencil counts.
- **Persistence within the smooth resolved phase.** Their first and second
  self-Ext groups vanish. The relative stable-sheaf moduli is locally étale
  over a smooth polarized deformation base, so the specified isolated
  contributions remain +1 nearby. This is not transport across the
  singular, topology-changing transition.
- **Partial wall and confinement results.** An identified magnetic split
  has one extension mode and a transverse phase alignment near large volume.
  Exact continued periods fix a constituent phase sign near the marked
  transition and exclude emission of the known light electric charges on
  a specified final segment. The probe charges carry confined magnetic
  flux after Higgsing.

These statements do not establish that the constituents survive every wall
along the continuation, classify every moduli component, or determine the
physical BPS spectrum at the endpoint. Counts for complete moving pencils
require a separate reduced-member hypothesis, isolated in Appendix A.

Together, the papers move from an exact geometric existence criterion with
a compact-gauging interpretation, to the
classical and quantum vector geometry of one realized transition, to a
separate investigation of its magnetic objects. Paper 3 uses the geometry
and periods of paper 2; paper 2's period proofs do not depend on the note.

## Organization

```text
papers/
  01-global-constraints/   Paper 1: PDF, main.tex, sections and figures
  02-quantum-transition/   Paper 2: PDF, main.tex and sections
  03-magnetic-sheaves/     Technical note: PDF, main.tex and sections
shared/                   Common notation, bibliographies and journal styles
checks/                   39 computational scripts, indexed by topic
scripts/                  Build and check runners
```

## Build and verify

Run from the repository root. Builds require TeX Live and `latexmk`;
the full computational suite also requires SageMath (tested with 10.9).

```bash
make papers       # build all three manuscripts
make paper2       # build just paper 2 (also: paper1, paper3)
make checks       # run the eleven Python checks
make check-full   # run all 39 Python/Sage scripts
make verify-full  # run all checks and build all papers
```

Generated PDFs and LaTeX auxiliaries go to the ignored `build/` directory.
The versioned `papers/*/paper.pdf` files are never overwritten by a build.
See the [computation guide](checks/README.md) for dependencies, individual
commands and the distinction between exact checks and numerical diagnostics.

## Citation and related work

[BibTeX entries](CITATION.bib) are provided for all three manuscripts;
[CITATION.cff](CITATION.cff) supplies machine-readable citation metadata.
Cite the paper relevant to the result used. The geometric companion papers
and their computations are in
[calabi-yau-smoothability](https://github.com/bwuebben/calabi-yau-smoothability).

Code is MIT licensed; the manuscripts and repository prose are CC BY 4.0.
The journal support files retain their upstream LPPL license. See [LICENSE](LICENSE).
