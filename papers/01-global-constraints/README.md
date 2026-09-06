# Paper 1 · Global constraints

*Global Constraints on Higgsing Localized Five-Dimensional Sectors in
Compact M-Theory Vacua* — Bernd Johannes Wuebben.

**[Read the paper · 37 pages](paper.pdf)** · [LaTeX entry point](main.tex)
· [Back to the collection](../../README.md)

Compact divisor pairings determine the Abelian flavor gaugings shared by
conifold, E₂ and E₃ sectors. Their branchwise charge kernel describes both
the invariant-coordinate image of the rigid moment-map zero locus and,
for the stated toric class, the global first-order deformation image.
Examples include the joint transition in X₉ and a thirty-sector obstruction.

The separate all-orders necessary condition applies at nodes, E₂ points
and line-smoothed E₃ points. General integrability and the all-orders
two-dimensional E₃ extension are not asserted.

## Files and build

- `main.tex` — manuscript entry point.
- `sections/` — body and appendices.
- `figures/` — native TikZ figures.
- `paper.pdf` — versioned reading copy.

From the repository root:

```bash
make paper1
```

The result is `build/01-global-constraints/main.pdf`. The build runner
locates the [shared notation, bibliography and styles](../../shared/README.md).
The relevant calculations are indexed in the
[computation guide](../../checks/README.md#general-constraints).
