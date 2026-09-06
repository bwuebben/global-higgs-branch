# Paper 3 · Magnetic sheaves

*Magnetic Sheaves and Stability Walls in a Compact Calabi–Yau Transition*
— Bernd Johannes Wuebben. Supporting technical note.

**[Read the note · 23 pages](paper.pdf)** · [LaTeX entry point](main.tex)
· [Back to the collection](../../README.md)

This note computes specified rank-one sheaf-moduli contributions and
constructs two isolated stable sheaves on a reducible moving divisor.
It studies polarization transport, constituent phase alignments,
final-segment electric-emission constraints and confined magnetic flux
after Higgsing. It uses the analytically continued periods of
[paper 2](../02-quantum-transition/).

These are component counts and partial stability results, not a
determination of the physical BPS spectrum at the quantum endpoint.
Complete-pencil counts require the explicit reduced-member hypothesis
stated in Appendix A.

## Files and build

- `main.tex` — note entry point.
- `sections/` — sheaf constructions, stability arguments and appendices.
- `paper.pdf` — versioned reading copy.

From the repository root:

```bash
make paper3
```

The result is `build/03-magnetic-sheaves/main.pdf`. The build runner locates
the [shared notation, bibliographies and styles](../../shared/README.md).
The supporting scripts are indexed under
[magnetic sheaves and stability](../../checks/README.md#magnetic-sheaves-and-stability).
