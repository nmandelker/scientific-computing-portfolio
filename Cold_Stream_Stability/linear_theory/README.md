# Linear theory — Mandelker et al. 2016

A representative sample of the MATLAB analysis behind the linear theory of
Mandelker et al. 2016, presented as-is. The full analysis is considerably
larger — the paper works out the marginal-stability and resonance behaviour of
the slab in detail across several appendices, in many more scripts. These are a
curated illustration of how the numerically-solved dispersion relations became
the paper's figures, not the complete set.

| File | Role |
|---|---|
| `nir_test_adiabatic.nb` | Mathematica: numerically solves the linear KHI dispersion relation for the planar-slab geometry — the central case of the analysis. The single-interface *sheet* is simpler (a quartic polynomial, solved directly in MATLAB in `compressible_sheet_growth_times.m`); the *cylinder* is harder, so the paper solves the slab first and carries the results over to the cylinder by analogy (the tanh → Bessel-function correspondence, §2.4). The cylinder's dispersion relation was solved in a separate notebook, not included here. Exports the slab solution grids as CSV (see `sample_output_ImP_00.csv`). |
| `compressible_sheet_growth_times.m` | Solves and analyzes the sheet (single-interface) dispersion relation directly. (M16, Fig. 1). |
| `compressible_slab_growth_times.m` | Growth-time analysis of the slab solutions from the notebook. |
| `KHI_phase_diagrams_analytic.m` | Stability phase diagrams over the (M_b, δ) parameter space (M16, Fig. 11). |
| `azimuthal_modes.m` | Structure of the azimuthal modes (m = 0, 1, 2, ...) of the cylinder — the cleanest cylinder example, and a paper figure (M16, Fig. 6). |
| `P1_cyl.m`, `cyl_mode_structure.m` | Pressure-perturbation mode structure and growth rates for the cylindrical stream — the cylinder analogue of the slab analysis (paper §2.4). |

## `verification/`

The RAMSES patch used to verify the linear theory in simulations, plus the
measurement scripts:

| File | Role |
|---|---|
| `condinit.f90` | Initial conditions: sheet/slab setup with a perturbation that is either a single perturbed variable (density, pressure, velocity, interface shape) or a *full analytic eigenmode*, selected by namelist parameters. |
| `hydro_parameters.f90`, `read_hydro_params.f90`, `init_flow_fine.f90` | Parameter definitions, namelist reading, and flow initialization for the patch. |
| `kh_eigenmode.nml` | Configuration for eigenmode-seeded verification runs, including the complex eigenfrequency `pert_omega` computed by the Mathematica notebook. |
| `kh_production.nml` | Configuration for the runs with a single-variable perturbation. |
| `sample_submit.sh` | Example cluster submission script. |
| `Eigen_growth.m` | Measures perturbation growth rates in the simulation outputs for comparison with theory (M16, Fig. 8). |

The measured growth rates were verified to converge with resolution (stable
across `levelmax` 12–14); those convergence tests are not shown in the paper,
whose rigor rested on the exact analytic dispersion relation that the
simulations were checked against.

The patch files extend the public RAMSES code (Teyssier 2002), which is not
included here.
