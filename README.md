# Compact Higgs transitions

Research manuscripts and reproducibility code by **Bernd Johannes Wuebben**.

## Read the papers

| | Manuscript | Read | Sources and guide |
|---|---|---|---|
| **1** | Global Constraints on Higgsing Localized Five-Dimensional Sectors in Compact M-Theory Vacua | [PDF · 37 pages](papers/01-global-constraints/paper.pdf) | [Paper 1](papers/01-global-constraints/) |
| **2** | Bulk Gauging and Quantum Couplings across a Compact Higgs Transition | [PDF · 54 pages](papers/02-quantum-transition/paper.pdf) | [Paper 2](papers/02-quantum-transition/) |
| **3** | Magnetic Sheaves and Stability Walls in a Compact Calabi–Yau Transition | [PDF · 23 pages](papers/03-magnetic-sheaves/paper.pdf) | [Technical note](papers/03-magnetic-sheaves/) |

Start with **paper 2** for the explicit compact transition and its quantum
vector couplings. **Paper 1** develops the general deformation/gauging
framework. **Paper 3** studies magnetic sheaves and partial stability
constraints; it does not determine the physical BPS spectrum at the endpoint.
Each paper's directory has its PDF, source, and a short guide to its results
and limitations.

## Organization

```text
papers/
  01-global-constraints/   Paper 1: PDF, main.tex, sections and figures
  02-quantum-transition/   Paper 2: PDF, main.tex and sections
  03-magnetic-sheaves/     Technical note: PDF, main.tex and sections
shared/                   Common notation, bibliographies and journal styles
checks/                   38 computational scripts, indexed by topic
scripts/                  Build and check runners
```

## Build and verify

Run from the repository root. Builds require TeX Live and `latexmk`;
the full computational suite also requires SageMath (tested with 10.9).

```bash
make papers       # build all three manuscripts
make paper2       # build just paper 2 (also: paper1, paper3)
make checks       # run the ten Python checks
make check-full   # run all 38 Python/Sage scripts
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
