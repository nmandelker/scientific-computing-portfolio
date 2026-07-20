# Linear theory — Mandelker et al. 2016

The analytic engine of the cold-streams program, presented as-is.

| File | Role |
|---|---|
| `nir_test_adiabatic.nb` | Mathematica: numerically solves the linear KHI dispersion relations (sheet, slab, cylinder geometries), exporting solution grids as CSV (see `sample_output_ImP_00.csv`). |
| `compressible_sheet_growth_times.m`, `compressible_slab_growth_times.m` | Growth-time analyses for the single-interface and slab cases. |
| `KHI_phase_diagrams.m`, `KHI_phase_diagrams_analytic.m`, `launch_KHI_phase_diagrams.m` | Stability phase diagrams over the (M_b, δ) parameter space. |
| `marginally_stable_cylinder.m` | Marginal-stability boundaries for the cylindrical stream. |
| `azimuthal_modes.m` | Structure of the azimuthal (m = 0, 1, 2, ...) modes of the cylinder (M16, Fig. 6). |
| `P1_cyl.m`, `cyl_mode_structure.m` | Pressure-perturbation mode structure in the cylinder. |

## `ramses_verification/`

The RAMSES patch used to verify the linear theory in simulations, plus the
measurement scripts:

| File | Role |
|---|---|
| `condinit.f90` | Initial conditions: sheet/slab setup with a perturbation that is either a single perturbed variable (density, pressure, velocity, interface shape) or a *full analytic eigenmode*, selected by namelist parameters. |
| `hydro_parameters.f90`, `read_hydro_params.f90`, `init_flow_fine.f90` | Parameter definitions, namelist reading, and flow initialization for the patch. |
| `kh_eigenmode.nml` | Configuration for eigenmode-seeded verification runs — including the complex eigenfrequency `pert_omega` computed by the Mathematica notebook. |
| `kh_production.nml` | Configuration for the paper's production runs (single-variable perturbation, higher resolution). |
| `sample_submit.sh` | Example cluster submission script. |
| `Eigen_growth.m` | Measures perturbation growth rates in the simulation outputs for comparison with theory (M16, Fig. 8). |
| `convergence.m`, `plot_convergence_resolution.m` | Resolution-convergence study of the measured growth rates. |
| `read_ramses_2d.m` | Reader for the raw 2D RAMSES outputs (these runs were small enough to analyze without a conversion step). |

The patch files extend the public RAMSES code (Teyssier 2002), which is not
included here.
