# Computational checks

The 41 scripts below support the three papers. Most are exact integer,
rational or symbolic computations. The numerical period and wall diagnostics
are identified separately. None replaces the geometric or analytic proofs.

## Run the checks

From the repository root:

```bash
make checks       # thirteen Python scripts
make check-full   # all 41 Python/Sage scripts
```

The full suite requires SageMath (tested with 10.9). Set `SAGE=/path/to/sage`
if needed. For Python checks, the runner uses Sage's Python when available;
otherwise it creates `.venv` with the pinned [requirements](../requirements.txt).
There are no external datasets or private-checkout dependencies.

Individual checks can be run with:

```bash
sage -python checks/check_x9_spectrum.py
sage -c "load('checks/check_x9_integral_screening.sage')"
```

The full runner uses degree eight for the standalone mirror/GV calculation
and spectator degree four for the regular magnetic-period calculation.
See [REFERENCE.md](REFERENCE.md) for all cutoffs, commands and qualifications.
The scripts share computations, so they are kept together rather than
duplicated across paper directories.

## General constraints

Primarily [paper 1](../papers/01-global-constraints/).

| Script | Calculation |
|---|---|
| [check_nodal_charge.py](check_nodal_charge.py) | Nodal charge kernels. |
| [check_e2_gauging.py](check_e2_gauging.py) | E₂ invariant ring and X₉ gauging. |
| [check_e3_gauging.py](check_e3_gauging.py) | Branchwise E₃ invariant rings. |
| [check_e1_e4_dictionary.py](check_e1_e4_dictionary.py) | 45 exact checks of integral ruling/pencil markings, moment-map pairings and balanced lifts. |
| [check_smoothing_criterion.py](check_smoothing_criterion.py) | E₃ discriminant factors, the exact row-space test and four counterexample restrictions; standard library only. |
| [check_e3_magnetic_quivers.py](check_e3_magnetic_quivers.py) | E₃ magnetic-quiver Hilbert series. |
| [check_xcirc_e3_profiles.py](check_xcirc_e3_profiles.py) | Four thirty-sector branch profiles. |
| [check_x9_spectrum.py](check_x9_spectrum.py) | Divisor, Hodge and multiplet counts. |

## Compact geometry and gauging

Primarily [paper 2](../papers/02-quantum-transition/).

| Script | Calculation |
|---|---|
| [check_x9_prepotential.sage](check_x9_prepotential.sage) | Resolved fan, intersection ring and nef limits. |
| [check_x9_smooth_side.py](check_x9_smooth_side.py) | Smooth-side cubic and second-Chern pairings. |
| [check_x9_integral_screening.sage](check_x9_integral_screening.sage) | Primitive lattices and compact screening. |
| [check_x9_dressed_instantons.sage](check_x9_dressed_instantons.sage) | Formal invariant ring and dressed currents. |
| [check_x9_current_exchange.sage](check_x9_current_exchange.sage) | Vector mixing and current exchange. |
| [check_x9_m5_string_anomaly.py](check_x9_m5_string_anomaly.py) | Wrapped-M5 anomaly comparison. |
| [x9_fibrations.sage](x9_fibrations.sage) | K3 and genus-one fibrations. |
| [x9_toric_base.sage](x9_toric_base.sage) | Flop and flat elliptic base. |
| [x9_weierstrass.sage](x9_weierstrass.sage) | Weierstrass discriminant and fiber counts. |
| [check_x9_sixd_anomaly.py](check_x9_sixd_anomaly.py) | Six-dimensional anomaly arithmetic. |
| [check_x9_heterotic_transition.sage](check_x9_heterotic_transition.sage) | E-string roots, Shioda class and volume bound. |
| [check_x9_heterotic_k3_candidate.sage](check_x9_heterotic_k3_candidate.sage) | Candidate heterotic K3; no bundle construction. |

## Quantum periods

Primarily [paper 2](../papers/02-quantum-transition/), with period inputs to
the magnetic-sheaf note.

| Script | Calculation |
|---|---|
| [check_x9_mirror_periods.sage](check_x9_mirror_periods.sage) | Frobenius series, mirror map and genus-zero GV data. |
| [check_x9_mirror_continuation.sage](check_x9_mirror_continuation.sage) | Nodal monodromies and local continuation. |
| [check_x9_local_elliptic_period.sage](check_x9_local_elliptic_period.sage) | Exact connection and numerical endpoint quadrature. |
| [check_x9_local_vanishing_charges.sage](check_x9_local_vanishing_charges.sage) | Local exceptional-curve charge marking. |
| [check_x9_finite_bulk_nodes.sage](check_x9_finite_bulk_nodes.sage) | Compact four-node degeneration. |
| [check_x9_integral_vanishing_marking.sage](check_x9_integral_vanishing_marking.sage) | Primitive compact vanishing charges. |
| [check_x9_thresholds.py](check_x9_thresholds.py) | Hypermultiplet thresholds, unit monodromies, neutral projections and special-geometry identities. |
| [check_x9_mirror_mutation.sage](check_x9_mirror_mutation.sage) | Smooth-side Laurent mutation. |
| [check_x9_regular_magnetic_periods.sage](check_x9_regular_magnetic_periods.sage) | Resummed regular magnetic periods. |
| [check_x9_neutral_quantum_periods.sage](check_x9_neutral_quantum_periods.sage) | Neutral-period coefficients and smooth-side GV data. |
| [check_x9_neutral_period_descent.sage](check_x9_neutral_period_descent.sage) | Symplectic quotient and period-descent data. |
| [check_x9_d6_period.sage](check_x9_d6_period.sage) | D6 marking and prepotential constant. |
| [check_x9_circle_limit.sage](check_x9_circle_limit.sage) | Auxiliary-circle scaling. |
| [check_x9_weighted_abel.sage](check_x9_weighted_abel.sage) | Weighted transverse finite-difference identities. |

## Magnetic sheaves and stability

Primarily [paper 3](../papers/03-magnetic-sheaves/).

| Script | Calculation |
|---|---|
| [check_x9_rank_one_index.py](check_x9_rank_one_index.py) | Hilbert-scheme component counts and flux shifts. |
| [check_x9_moving_probes.sage](check_x9_moving_probes.sage) | Support restrictions and conditional pencil counts. |
| [check_x9_glued_probe.sage](check_x9_glued_probe.sage) | Two isolated sheaves and stability inequalities. |
| [check_x9_probe_transport.sage](check_x9_probe_transport.sage) | Polarization transport and slope walls. |
| [check_x9_probe_phase_wall.sage](check_x9_probe_phase_wall.sage) | Constituent phases near large volume. |
| [check_x9_quantum_wall.sage](check_x9_quantum_wall.sage) | Numerical alignment tables at finite degree. |
| [check_x9_magnetic_higgs.sage](check_x9_magnetic_higgs.sage) | Final-segment emission constraints and confined flux. |

Component counts are not a full physical BPS index. Numerical cutoff
agreement is not a certified instanton-tail bound. These qualifications
are retained in the individual scripts, detailed reference and manuscripts.
