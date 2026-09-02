#!/usr/bin/env bash
set -euo pipefail

task_python="${PYTHON:-python3}"

if [[ ! -x .venv/bin/python ]]; then
  "$task_python" -m venv .venv
fi

.venv/bin/python -m pip install --quiet -r requirements.txt

.venv/bin/python checks/check_nodal_charge.py
.venv/bin/python checks/check_e2_gauging.py
.venv/bin/python checks/check_e3_gauging.py
.venv/bin/python checks/check_e3_magnetic_quivers.py
.venv/bin/python checks/check_xcirc_e3_profiles.py
.venv/bin/python checks/check_x9_spectrum.py

latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
