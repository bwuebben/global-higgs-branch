# Paper 2 · The compact quantum transition

*Bulk Gauging and Quantum Couplings across a Compact Higgs Transition*
— Bernd Johannes Wuebben.

**[Read the paper · 54 pages](paper.pdf)** · [LaTeX entry point](main.tex)
· [Back to the collection](../../README.md)

The X₉ model has a mixed transition from Hodge numbers (6,122) to (3,123),
with an explicit primitive surviving charge lattice and cubic couplings.
The shared flavor gauging cannot remain dynamical when gravity decouples
with bounded local Kähler scales. On a marked finite-radius branch, four
conifold states span a primitive rank-three electric lattice. Laurent
mutation and continuation of all eight regular periods identify the
surviving genus-zero prepotential with the smooth-phase prepotential,
including its constant and integral normalization.

The fibration and heterotic K3 data do not supply the global heterotic
bundle construction or massless U(1). The finite-Planck hypermultiplet
metric remains open. The four-state description is not uniform in the
five-dimensional circle limit, and the period results refer to the stated
continuation marking.

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
