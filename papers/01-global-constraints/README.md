# Paper 1 · Global constraints

*Global Constraints on Higgsing Localized Five-Dimensional Sectors in
Compact M-Theory Vacua* — Bernd Johannes Wuebben.

**[Read the paper · 38 pages](paper.pdf)** · [LaTeX entry point](main.tex)
· [Back to the collection](../../README.md)

Compact divisor pairings determine the Abelian flavor gaugings shared by
conifold, E₂ and E₃ sectors. Their branchwise charge kernel describes both
the invariant-coordinate image of the rigid moment-map zero locus and,
for the stated toric class, the global first-order deformation image.
Examples include the joint transition in X₉ and a thirty-sector obstruction.

Combined with the companion mixed smoothing theorem, the kernel gives an
exact existence criterion after all local discriminants are excluded.
On the E₃ plane this means H₁H₂(H₁−H₂) ≠ 0, not merely a nonzero vector.
The selected deformation base is smooth when its relation space projects
nontrivially to each selected E₃ summand, and every tangent on it integrates.
A separate thirty-node/two-cone example shows that nonzero local deformations
can be forced to retain a singularity even on smooth deformation bases.

These geometric results do not identify the compact hypermultiplet metric,
settle all zero-projection E₃ profiles, or determine unrestricted nilpotent
deformation schemes. The theorem and its analytic proof are credited to
[CY paper 5](https://github.com/bwuebben/calabi-yau-smoothability/tree/main/paper5).

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
