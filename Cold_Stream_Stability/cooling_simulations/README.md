# Stage 4 (M20a): cooling simulations — the full vertical slice

This directory is the complete pipeline behind
[Mandelker et al. 2020a](https://arxiv.org/abs/1910.05344): from running the
idealized cold-stream simulations to the figures in the paper. It is the one
place in the cold-streams project shown end to end — designing and running
the numerical experiments, converting the raw outputs, measuring stream
properties, and producing the analytic estimates and plots. Presented as-is,
research code.

The simulations use the public AMR code
[RAMSES](https://arxiv.org/abs/astro-ph/0111367); RAMSES itself is not
included — only my patches to it and my analysis code.

## The pipeline, in order

### 1. `ramses_patch/` — set up and run the simulations

A RAMSES "patch" replaces or extends specific source files of the public
code. This one configures an idealized cold stream (a dense, cold cylinder
in a hot wind) with radiative cooling:

| File | Role |
|---|---|
| `condinit.f90` | Initial conditions: the cold cylindrical stream, hot background, shear velocity, and the perturbation (single-variable or eigenmode). |
| `cooling_module.f90`, `cooling_fine.f90` | The radiative cooling and UV-heating routines. The stream and background are given different metallicities (`met_s`, `met_b` — passive scalars set in `condinit.f90`), with UV-background heating (Haardt-Madau) and self-shielding as runtime options. |
| `hydro_parameters.f90`, `read_hydro_params.f90`, `read_params.f90` | Parameter definitions and namelist reading for the added stream/cooling parameters. |
| `adaptive_loop.f90`, `init_time.f90` | Modifications to the main time loop and initialization for the idealized-box setup and custom outputs. |
| `kh_example.nml` | A representative namelist. Many runs span a grid of Mach number, density contrast, stream density, and metallicity, so no single namelist captures the whole suite — this is one example; individual parameters (e.g. the cooling temperature floor) vary per run. |

### 2. `conversion/` — raw outputs → compact analysis format

RAMSES writes the full AMR hierarchy. For analysis I convert each snapshot to
a compact format keeping only the AMR *leaf* cells, laid out similarly (though
not identically) to the post-processed VELA cosmological outputs used in the
giant-clumps projects — much smaller and simpler to analyze.

| File | Role |
|---|---|
| `make_ART_format.f90` | Reads a raw RAMSES snapshot and writes the leaf-only binary (density, velocity, pressure, passive scalar / "colour", cell size). |
| `loop_output_ART.f90` | Driver that runs the conversion over all snapshots of a simulation. |
| `submit_loop_output_ART.sh` | Example HPC batch submission script. |

### 3. `analysis/` — measure stream properties

Fortran programs that read the leaf-format snapshots and measure the physical
quantities behind the paper's figures. These began as several separate
scripts (one per quantity) that I later consolidated into a single executable;
this is that consolidated version.

| File | Role |
|---|---|
| `stream_analysis_combined.f90` | The main analysis: cold-gas mass, stream volume, centre-of-mass velocity (deceleration), turbulent velocity dispersion in the mixing layer, kinetic/thermal energy budgets, radial profiles, and density PDFs. |
| `stream_analysis_collate.f90` | Aggregates the per-snapshot outputs into time series for plotting. |
| `compute_cooling_rates.f90` | Computes the net radiative cooling / UV-heating rates in the simulations, for the energy-budget figures. |

A numerical note: the turbulent velocity dispersion is accumulated as the
mass-weighted `(v − ⟨v⟩)²` in double precision, rather than the
`⟨v²⟩ − ⟨v⟩²` form — the latter suffers catastrophic cancellation and gave
unreliable dispersions in early single-precision versions.

### 4. `matlab/` — analytic estimates and final plots

| File | Role |
|---|---|
| `Rs_crit_panels.m` | Computes and plots the critical stream radius where `t_cool,mix = t_shear` — the paper's Fig. 1 (the survival criterion). |
| `Rs_over_cs_tcool.m`, `tcool_mix_over_tsc.m`, `find_tcool_mix_over_tshear.m` | The governing dimensionless numbers: mixing-layer cooling time vs. shear and sound-crossing times, as functions of stream properties. |
| `timescale_comparison.m` | Compares the relevant timescales (cooling, shear, sound-crossing, disruption) across the parameter space. |
| `simulation_cooling_check.m` | Validates the cooling implemented in the simulations against the analytic cooling function. |
| `find_cooling_equilibrium.m`, `cooling_equilibrium.m`, `dlnL_dlnT.m` | Helpers: thermal-equilibrium temperature and local slope of the cooling curve, used by the estimates above. |

Raw simulation outputs are not included (too large); the figures show what the
pipeline produces.
