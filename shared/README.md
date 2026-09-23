# Shared manuscript resources

The four papers use a single copy of these resources:

- [macros.tex](macros.tex) — notation.
- [bibliography/foundations.bib](bibliography/foundations.bib) — references
  used by paper 1 and the later manuscripts.
- [bibliography/transition.bib](bibliography/transition.bib) — additional
  references for the quantum-transition paper and magnetic-sheaf note.
- [bibliography/finite_planck.bib](bibliography/finite_planck.bib) — further
  references for paper 4.
- [tex/jheppub.sty](tex/jheppub.sty) and [tex/JHEP.bst](tex/JHEP.bst) —
  upstream JHEP support files, preserved under their original license.

Use `make papers` or `make paper1` … `make paper4` from the
repository root. The build runner supplies the TeX and bibliography search
paths and keeps generated files under `build/`.
