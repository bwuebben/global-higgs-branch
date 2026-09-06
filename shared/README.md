# Shared manuscript resources

The three papers use a single copy of these resources:

- [macros.tex](macros.tex) — notation.
- [bibliography/foundations.bib](bibliography/foundations.bib) — references
  used by paper 1 and the later manuscripts.
- [bibliography/transition.bib](bibliography/transition.bib) — additional
  references for the quantum-transition paper and magnetic-sheaf note.
- [tex/jheppub.sty](tex/jheppub.sty) and [tex/JHEP.bst](tex/JHEP.bst) —
  upstream JHEP support files, preserved under their original license.

Use `make papers` or `make paper1`, `make paper2`, `make paper3` from the
repository root. The build runner supplies the TeX and bibliography search
paths and keeps generated files under `build/`.
