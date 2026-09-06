# Compact Higgs transitions: constraints, quantum couplings and magnetic sheaves

Three manuscripts by Bernd Johannes Wuebben, with LaTeX sources, figures
and reproducibility code. The common setting is M-theory on compact
Calabi–Yau threefolds with several localized five-dimensional sectors.

## The papers

| Paper | Manuscript | PDF | Source |
|---|---|---|---|
| 1 | **Global Constraints on Higgsing Localized Five-Dimensional Sectors in Compact M-Theory Vacua** | [37 pages](global-higgs-branch.pdf) | [main.tex](main.tex) |
| 2 | **Bulk Gauging and Quantum Couplings across a Compact Higgs Transition** | [54 pages](bulk-gauging-quantum-couplings.pdf) | [main2.tex](main2.tex) |
| 3 | **Magnetic Sheaves and Stability Walls in a Compact Calabi–Yau Transition** — supporting technical note | [23 pages](magnetic-sheaves-stability-walls.pdf) | [magnetic_probes.tex](magnetic_probes.tex) |

**Paper 1: general constraints.** Compact divisor pairings determine the
Abelian flavor gaugings shared by conifold, E₂ and E₃ sectors. Their
branchwise charge kernel describes the invariant-coordinate image of the
rigid moment-map zero locus and, for the stated toric class, the global
first-order deformation image. A separate all-orders necessary condition
applies at nodes, E₂ points and line-smoothed E₃ points. Examples include
the joint transition in X₉ and a thirty-sector obstruction. General
integrability is not asserted.

**Paper 2: one transition and its quantum vector geometry.** The X₉ model
has a mixed transition from Hodge numbers (6,122) to (3,123), with an
explicit primitive surviving charge lattice and cubic couplings. The
shared flavor gauging cannot stay dynamical when gravity decouples while
the local Kähler scales remain bounded. On a marked finite-radius branch,
four conifold states span a primitive rank-three electric lattice.
Laurent mutation and continuation of all eight regular periods identify
the surviving genus-zero prepotential with that of the smooth phase,
including its constant and integral normalization. The paper also gives
fibration, six-dimensional anomaly and heterotic K3 data.

**Paper 3: magnetic sheaves and partial stability results.** The technical
note computes specified rank-one sheaf-moduli contributions and constructs
two isolated stable sheaves on a reducible moving divisor. It studies
polarization transport, constituent phase alignments, final-segment
electric-emission constraints and confined magnetic flux after Higgsing.
Complete-pencil counts are conditional on an explicit reduced-member
hypothesis. These results do not determine the physical BPS spectrum at
the quantum endpoint.

For the explicit compact physics, start with paper 2. Paper 1 supplies the
general deformation/gauging framework; paper 3 develops the magnetic
questions separately. Each manuscript builds independently. The geometric
companion papers cited throughout are available in
[calabi-yau-smoothability](https://github.com/bwuebben/calabi-yau-smoothability).

## Scope

These are research manuscripts, not a claim of journal acceptance. The
finite-Planck hypermultiplet metric, global heterotic bundle construction
and massless U(1), and physical magnetic stability remain open. The
finite-radius four-state description is not a uniform five-dimensional
effective theory: its Kaluza–Klein gap closes in the circle limit. The
all-orders period results refer to the specified continuation marking.

## Reproduce the calculations and build the papers

Requirements: Python 3, a TeX Live installation with `latexmk`, and
SageMath 10.9 for the full computational suite. The Python dependencies
are pinned in [requirements.txt](requirements.txt). If Sage is available,
the runner uses its Python; otherwise it creates `.venv` and installs the
pinned Python requirements. The full suite requires Sage and never silently
skips its calculations.

```bash
./verify.sh                # ten Python checks and all three builds
./verify.sh --full         # all 38 scripts and all three builds
./verify.sh --build-only   # build all three manuscripts without checks
```

`SAGE=/path/to/sage` selects a Sage executable. The full suite includes
degree-eight mirror/GV data and spectator-degree-four regular magnetic
periods. It can take several minutes. See [checks/README.md](checks/README.md)
for individual commands, mathematical scope and numerical qualifications.
The calculations are self-contained; no private checkout or external data
is required after dependencies are installed. Finite checks supplement,
but do not replace, the proofs.

To build a single manuscript:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error main2.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error magnetic_probes.tex
```

These produce `main.pdf`, `main2.pdf` and `magnetic_probes.pdf` locally.
The descriptively named PDFs linked above are the versioned release copies.
Building does not overwrite those release copies.

## Repository layout

| Path | Contents |
|---|---|
| `main.tex`, `sections/` | Paper 1 and its appendices. |
| `main2.tex`, `sections2_restructured/` | Paper 2 and its analytic and lattice appendices. |
| `magnetic_probes.tex`, `sections_probes/` | Supporting technical note. |
| `macros.tex`, `references.bib`, `references_new.bib` | Shared notation and bibliographies. |
| `figures/` | Native TikZ figures. |
| `checks/` | Ten Python and 28 Sage scripts, with a computation guide. |
| `verify.sh`, `requirements.txt` | Build/check runner and Python dependencies. |
| `jheppub.sty`, `JHEP.bst` | Journal support files under their upstream license. |

## Citation and license

Cite the individual manuscript relevant to the result used.
[CITATION.bib](CITATION.bib) provides all three BibTeX entries;
[CITATION.cff](CITATION.cff) provides repository metadata and manuscript
references, with paper 2 as the preferred citation for the compact-transition
calculations.

The code in `checks/` and `verify.sh` is MIT licensed. The manuscripts,
figures, bibliographies and repository prose are CC BY 4.0. The JHEP support
files retain the LaTeX Project Public License stated in their headers.
See [LICENSE](LICENSE).
