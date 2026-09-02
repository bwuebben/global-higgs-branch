# Global constraints on Higgsing localized five-dimensional sectors

Source, figures, and exact finite checks for

> **Global Constraints on Higgsing Localized Five-Dimensional Sectors in
> Compact M-Theory Vacua**
> Bernd Johannes Wuebben, 2026.

**[Read the current manuscript](global-higgs-branch.pdf)**

The paper determines which local Higgs-branch flows survive when conifold
points and rank-one E₂ and E₃ sectors occur together in a compact
Calabi--Yau threefold. Wrapped M2-branes supply electric charges; compact
divisor pairings determine the flavor cocharacters gauged by dynamical
five-dimensional vectors. On each reduced E₃ branch, the resulting compact
matrix has the same kernel as the projected rigid moment-map locus and, for
the stated toric class, the image of global first-order deformations.

The examples include a cooperative transition on X₉, where two conifold
coordinates and one E₂ coordinate can move only along one common direction,
and a compact model with thirty locally smoothable sectors but no trajectory
that activates all thirty.

## Verify and build

Requirements are Python 3 and a current TeX Live installation with
`latexmk`. From the repository root, one command creates a local Python
environment, runs all six exact checks, and compiles the manuscript:

```bash
./verify.sh
```

The checks use integer, rational, polynomial, or symbolic arithmetic. They
test the finite algebra printed in the manuscript: the nodal and E₂/E₃
invariant-coordinate calculations, magnetic-quiver Hilbert series, all four
thirty-sector branch profiles, and the X₉ divisor, Hodge, and multiplet
calculation. They supplement rather than replace the geometric proofs.

The related geometric papers and their source computations are available in
[calabi-yau-smoothability](https://github.com/bwuebben/calabi-yau-smoothability).

## Repository layout

| Path | Contents |
|---|---|
| `global-higgs-branch.pdf` | Current 37-page US-Letter manuscript. |
| `main.tex`, `macros.tex`, `references.bib` | Manuscript entry point, notation, and bibliography. |
| `sections/` | Paper body and appendices. |
| `figures/` | Native TikZ figures. |
| `checks/` | Six exact finite checks for the calculations printed in the paper. |
| `jheppub.sty`, `JHEP.bst` | JHEP support files under their upstream license. |

The manuscript is prepared for submission to JHEP. An arXiv identifier has
not yet been assigned.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Until an
arXiv identifier is assigned, cite the manuscript by author, title, and year.

## License

The code in `checks/` and `verify.sh` is MIT licensed. The paper, figures,
bibliography, and repository prose are CC BY 4.0. The JHEP support files retain
the LaTeX Project Public License stated in their headers. See [`LICENSE`](LICENSE).
