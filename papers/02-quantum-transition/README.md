# Paper 2 · The compact quantum transition

*Bulk Gauging and Quantum Couplings across a Compact Higgs Transition*
— Bernd Johannes Wuebben.

**[Read the paper · 56 pages](paper.pdf)** · [LaTeX entry point](main.tex)
· [Back to the collection](../../README.md)

The X₉ model has a mixed transition from Hodge numbers (6,122) to (3,123),
with an explicit primitive surviving charge lattice and cubic couplings.
The shared flavor gauging cannot remain dynamical when gravity decouples
with bounded local Kähler scales. On a marked finite-radius branch, four
conifold states span a primitive rank-three electric lattice. Laurent
mutation and continuation of all eight regular periods identify the
surviving genus-zero prepotential with the smooth-phase prepotential,
including its constant and integral normalization.

The four-hypermultiplet one-loop thresholds give the same integral
monodromy as the geometric vanishing cycles. Their neutral projections
vanish; determining the remaining regular couplings still requires the
full period comparison. The latter fixes both the scalar metric and
the electromagnetic coupling matrix, including graviphoton mixing.

The exact mixed smoothing criterion and selected-base integrability theorem
of [CY paper 5](https://github.com/bwuebben/calabi-yau-smoothability/tree/main/paper5)
now supply the general existence argument. The explicit X₉ construction
is retained because it identifies the smooth phase used in the quantum
comparison. This geometric integrability does not supply a
finite-Planck hypermultiplet metric.

The fibration and heterotic K3 data do not supply the global heterotic
bundle construction or massless U(1). At fixed five-dimensional Planck scale, the circle
limit retains finite bulk gauging but closes the Kaluza–Klein gap, so the
four-state description has no uniform energy range. This is distinct
from decoupling gravity. The period results refer to the stated marking.

## Files and build

- `main.tex` — manuscript entry point.
- `sections/` — body, lattice data and analytic period proofs.
- `paper.pdf` — versioned reading copy.

From the repository root:

```bash
make paper2
```

The result is `build/02-quantum-transition/main.pdf`. The build runner
locates the [shared notation, bibliographies and styles](../../shared/README.md).
See the [geometry](../../checks/README.md#compact-geometry-and-gauging) and
[period](../../checks/README.md#quantum-periods) check groups. Magnetic
sheaves are treated separately in [paper 3](../03-magnetic-sheaves/).
