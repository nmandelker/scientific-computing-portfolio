# Linear theory — Mandelker et al. 2016

The analytic engine of the cold-streams program, presented as-is.

| File | Role |
|---|---|
| `nir_test_adiabatic.nb` | Mathematica: numerically solves the linear KHI dispersion relation for the planar-slab geometry — the hardest of the three cases — exporting solution grids as CSV (see `sample_output_ImP_00.csv`). (The single-interface sheet reduces to a quartic polynomial solved directly in MATLAB; the cylinder was solved in a separate notebook, not included here.) |
| `compressible_sheet_growth_times.m` | Solves and analyzes the sheet (single-interface) dispersion relation directly — the quartic case (M16, Fig. 1). |
| `compressible_slab_growth_times.m` | Growth-time analysis of the slab solutions from the notebook. |
| `KHI_phase_diagrams_analytic.m` | Stability phase diagrams over the (M_b, δ) parameter space (M16, Fig. 11). |
| `marginally_stable_cylinder.m` | Marginal-stability boundaries for the cylindrical stream. |
| `azimuthal_modes.m` | Structure of the azimuthal (m = 0, 1, 2, ...) modes of the cylinder (M16, Fig. 6). |
| `P1_cyl.m`, `cyl_mode_structure.m` | Eigenmode pressure structure and growth rates in the cylinder. |

## `ramses_verification/`

The RAMSES patch used to verify the linear theory in simulations, plus the
measurement scripts:

| File | Role |
|---|---|
| `condinit.f90` | Initial conditions: sheet/slab setup with a perturbation that is either a single perturbed variable (density, pressure, velocity, interface shape) or a *full analytic eigenmode*, selected by namelist parameters. |
| `hydro_parameters.f90`, `read_hydro_params.f90`, `init_flow_fine.f90` | Parameter definitions, namelist reading, and flow initialization for the patch. |
| `kh_eigenmode.nml` | Configuration for eigenmode-seeded verification runs, including the complex eigenfrequency `pert_omega` computed by the Mathematica notebook. |
| `kh_production.nml` | Configuration for the paper's production runs with a single-variable perturbation. |
| `sample_submit.sh` | Example cluster submission script. |
| `Eigen_growth.m` | Measures perturbation growth rates in the simulation outputs for comparison with theory (M16, Fig. 8). |

The measured growth rates were verified to converge with resolution (stable
across `levelmax` 12–14); those convergence tests are not shown in the paper,
whose rigor rested on the exact analytic dispersion relation that the
simulations were checked against.

The patch files extend the public RAMSES code (Teyssier 2002), which is not
included here.
