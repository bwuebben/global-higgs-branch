# Exact checks

The six scripts in this directory verify the finite algebra printed in the
paper. They use only exact integer, rational, and symbolic arithmetic.

Run them from the project root with the commands printed in the manuscript's
reproducibility appendix:

```bash
./venv/bin/python checks/check_nodal_charge.py
./venv/bin/python checks/check_e2_gauging.py
./venv/bin/python checks/check_e3_gauging.py
./venv/bin/python checks/check_e3_magnetic_quivers.py
./venv/bin/python checks/check_xcirc_e3_profiles.py
./venv/bin/python checks/check_x9_spectrum.py
```

The checks cover, respectively, the nodal charge kernel, the $E_2$ invariant
scheme and $X_9$ matrix, the $E_3$ branchwise invariant scheme, the two $E_3$
magnetic-quiver Hilbert series, all four $X^{\circ}$ profiles, and the $X_9$
divisor/Hodge/multiplet calculation. Source geometric verification files remain in
`../../cy_smoothing/paper4/certificates/` and
`../../cy_smoothing/paper5/certificates/`.
