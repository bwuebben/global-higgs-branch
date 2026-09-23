# Paper 4 · The Higgs branch at finite Planck mass

*The Higgs Branch of a Compact Calabi–Yau Transition at Finite Planck Mass*
— Bernd Johannes Wuebben.

**[Read the paper · 45 pages](paper.pdf)** · [LaTeX entry point](main.tex)
· [Back to the collection](../../README.md)

With gravity decoupled, the Higgs branch opened by light wrapped-M2
hypermultiplets is locally a hyperkähler quotient. In a compact
compactification the hypermultiplet moduli space is quaternionic Kähler, and
its metric is not known. The paper determines the local structure of the Higgs
branch anyway. The underlying fact is that a quaternionic Kähler moment map has no
additive constant, so supergravity has no Fayet–Iliopoulos parameter that could
resolve or deform the singular point.

The main results are:

- **A local model for quaternionic Kähler quotients.** Near a fixed point at
  which the moment map vanishes, and under a freeness condition on the normal
  representation, the quotient is locally the product of the fixed set with the
  flat hyperkähler quotient of that representation. The metric tangent cone is
  the flat product, with relative corrections of order r² at distance r. The
  result holds for either sign of the scalar curvature and without
  completeness, and also for hyperkähler and Kähler quotients. When the flat
  quotient is four-dimensional, the leading correction is the curvature
  restricted to a quaternionic line.
- **Faces with one relation.** Suppose exactly k disjoint (−1,−1)-curves have
  zero area on a face of the nef cone with finite gauge couplings, and their
  classes satisfy one primitive relation with nonzero coefficients aᵢ. Then,
  under a stated regularity hypothesis, the Higgs branch at the root is locally
  a smooth factor times C²/Z_m, m = Σ|aᵢ|, as without gravity.
- **The X₉ transition.** A face of the nef cone of the resolved X₉ of
  [paper 2](../02-quantum-transition/) has exactly four such curves, giving
  C²/Z₄. The leading correction there lies in a three-dimensional space. The
  same local structure holds for the metric completion of the smooth-threefold
  vacua where the del Pezzo surface collapses.
- **Type IIA in four dimensions.** Under the corresponding hypotheses, the
  A_{k−1} point near k homologous vanishing three-cycles, expected by Greene,
  Morrison and Vafa, persists at finite Planck mass.
- **A census.** Calabi–Yau hypersurfaces from favourable Kreuzer–Skarke
  polytopes with h¹¹ ≤ 5, together with a sample at h¹¹ = 6, give 5605 faces of
  toric Kähler cones whose light spectrum has one relation. All of them have
  unit coefficients, so m = k throughout. Whether faces with m ≠ k exist
  remains open.

The results about quaternionic Kähler quotients are proved unconditionally.
The physical theorems assume the regularity hypothesis, whose main content is
a smooth Wilsonian hypermultiplet manifold containing the light states. The
four-dimensional corollary also assumes an identification of the Higgs vacua
with type IIA compactifications on the smoothings. The paper makes no
complex-analytic claim, and it does not identify the completion at the del
Pezzo point with the Higgs branch of the interacting E₂ theory.

## Files and build

- `main.tex` — manuscript entry point.
- `sections/` — body and appendix.
- `paper.pdf` — versioned reading copy.
- References specific to this paper are in
  [shared/bibliography/finite_planck.bib](../../shared/bibliography/finite_planck.bib).

From the repository root:

```bash
make paper4
```

The result is `build/04-finite-planck-vacua/main.pdf`. The build runner
locates the [shared notation, bibliographies and styles](../../shared/README.md).
The census and curvature computations are not part of this repository's check
suite; the scripts are available from the author.
