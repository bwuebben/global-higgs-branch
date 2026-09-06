#!/usr/bin/env bash
set -euo pipefail

cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
task_mode="${1:-default}"
task_sage="${SAGE:-sage}"

if [[ $# -gt 1 ]]; then
  printf 'Usage: %s [--full|--help]\n' "$0" >&2
  exit 2
fi
case "$task_mode" in
  default|--full) ;;
  --help|-h)
    printf '%s\n' \
      'Usage: ./scripts/check.sh [--full]' \
      'Default: twelve Python checks. --full: all 40 scripts; requires SageMath.' \
      'Set SAGE to select Sage; PYTHON selects venv creation when Sage is absent.'
    exit 0 ;;
  *) printf 'Unknown option: %s\n' "$task_mode" >&2; exit 2 ;;
esac

if [[ "$task_mode" == --full ]] && ! command -v "$task_sage" >/dev/null 2>&1; then
  printf 'The full suite requires SageMath (tested with 10.9). Set SAGE if needed.\n' >&2
  exit 1
fi

if command -v "$task_sage" >/dev/null 2>&1; then
  task_python=("$task_sage" -python)
else
  if [[ ! -x .venv/bin/python ]]; then
    "${PYTHON:-python3}" -m venv .venv
  fi
  .venv/bin/python -m pip install --quiet -r requirements.txt
  task_python=(.venv/bin/python)
fi

task_count=0
for task_check in checks/check_*.py; do
  printf 'Running %s\n' "$task_check"
  "${task_python[@]}" "$task_check"
  task_count=$((task_count + 1))
done

if [[ "$task_mode" == --full ]]; then
  for task_check in checks/*.sage; do
    printf 'Running %s\n' "$task_check"
    case "$task_check" in
      checks/check_x9_mirror_periods.sage)
        MIRROR_DEGREE=8 "$task_sage" -c "load('$task_check')" ;;
      checks/check_x9_regular_magnetic_periods.sage)
        MAGNETIC_BULK_DEGREE=4 "$task_sage" -c "load('$task_check')" ;;
      *) "$task_sage" -c "load('$task_check')" ;;
    esac
    task_count=$((task_count + 1))
  done
fi
printf '%s computational scripts completed successfully.\n' "$task_count"
