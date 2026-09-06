#!/usr/bin/env bash
set -euo pipefail

task_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
task_selection="${1:-all}"
if [[ $# -gt 1 ]]; then
  printf 'Usage: %s [1|2|3|all]\n' "$0" >&2
  exit 2
fi
case "$task_selection" in
  all) task_papers=(01-global-constraints 02-quantum-transition 03-magnetic-sheaves) ;;
  1) task_papers=(01-global-constraints) ;;
  2) task_papers=(02-quantum-transition) ;;
  3) task_papers=(03-magnetic-sheaves) ;;
  *) printf 'Usage: %s [1|2|3|all]\n' "$0" >&2; exit 2 ;;
esac

if ! command -v latexmk >/dev/null 2>&1; then
  printf 'latexmk is required; install TeX Live before building.\n' >&2
  exit 1
fi

# A trailing separator retains the standard TeX distribution search paths.
export TEXINPUTS="$task_root/shared/tex:${TEXINPUTS:-}"
export BIBINPUTS="$task_root/shared/bibliography:${BIBINPUTS:-}"
export BSTINPUTS="$task_root/shared/tex:${BSTINPUTS:-}"

for task_paper in "${task_papers[@]}"; do
  task_output="$task_root/build/$task_paper"
  mkdir -p -- "$task_output"
  printf 'Building %s\n' "$task_paper"
  (
    cd -- "$task_root/papers/$task_paper"
    latexmk -pdf -silent -interaction=nonstopmode -halt-on-error \
      -outdir="$task_output" main.tex
  )
  printf 'PDF: build/%s/main.pdf\n' "$task_paper"
done
