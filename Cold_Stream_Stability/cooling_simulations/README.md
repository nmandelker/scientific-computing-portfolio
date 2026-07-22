# Stage 4 (M20a): cooling simulations — the full vertical slice

This directory is the complete pipeline behind
[Mandelker et al. 2020a](https://arxiv.org/abs/1910.05344): from running the
idealized cold-stream simulations to the figures in the paper. It is the one
place in the cold-streams project shown end to end — designing and running
the numerical experiments, converting the raw outputs, and measuring stream
properties. Presented as-is, research code.

The simulations use the public AMR code
[RAMSES](https://arxiv.org/abs/astro-ph/0111367); RAMSES itself is not
included, only my patches to it and my analysis code.

## The pipeline, in order

### 1. `ramses_patch/` — set up and run the simulations

A RAMSES "patch" replaces or extends specific source files of the public
code. This one configures an idealized cold stream (a dense, cold cylinder
pressure confined by and flowing through a static, hot background) with 
radiative cooling and heating:

| File | Role |
|---|---|
| `condinit.f90` | Initial conditions: the cold cylindrical stream, hot background, shear velocity, and the perturbation (single-variable or eigenmode). |
| `cooling_module.f90`, `cooling_fine.f90` | RAMSES's stock radiative-cooling and UV-heating routines with targeted modifications for the idealized (non-cosmological) setup: control over the redshift of the UV background, and a **temperature ceiling `Tmax_cool` above which cooling is switched off** (keeping the hot background from cooling). The stream and background are assigned different metallicities (`met_s`, `met_b`) which impacts their cooling rates, and different passive scalars which allow us to trace phase mixing, in `condinit.f90`. |
| `hydro_parameters.f90`, `read_hydro_params.f90`, `read_params.f90` | Parameter definitions and namelist reading for the added stream/cooling parameters. |
| `adaptive_loop.f90`, `init_time.f90` | Modifications to the main time loop and initialization for the idealized-box setup and custom outputs. |
| `kh_example.nml` | A representative namelist. Many runs span a grid of Mach number, density contrast, stream density, and metallicity, so no single namelist captures the whole suite. This is one example; individual parameters vary per run. |

### 2. `conversion/` — raw outputs → compact analysis format

RAMSES writes the full AMR hierarchy and outputs one file per core used during runtime for every snapshot. 
For analysis, I convert each snapshot to a single file per snapshot, keeping only the AMR *leaf* cells in 
a compact format. These are much smaller and simpler to analyze.

| File | Role |
|---|---|
| `make_ART_format.f90` | Reads a raw RAMSES snapshot and writes the leaf-only binary (density, velocity, pressure, passive scalar / "colour", cell size). |
| `loop_output_ART.f90` | Driver that runs the conversion over all snapshots of a simulation. |
| `time.f90` | Reads each snapshot's `info_*.txt` header and writes a catalog of the exact output times (in code units), so the analysis stage can attach a timestamp to each compacted snapshot without re-parsing the raw outputs. |
| `submit_loop_output_ART.sh` | Example HPC batch submission script. |

### 3. `analysis/` — measure stream properties

Fortran programs that read the leaf-format snapshots and measure the physical
quantities behind the paper's figures. These began as several separate
scripts (one per quantity) that I later consolidated into a single executable;
this is that consolidated version.

| File | Role |
|---|---|
| `stream_analysis_combined.f90` | The main analysis: cold-gas mass, stream volume, centre-of-mass velocity (deceleration), turbulent velocity dispersion in the mixing layer, kinetic/thermal energy budgets, radial profiles, density PDFs, and other diagnostics of phase mixing and stream evolution. |
| `stream_analysis_collate.f90` | Aggregates the per-snapshot outputs into time series for plotting. |
| `compute_cooling_rates.f90` | Computes the net radiative cooling / UV-heating rates in the simulations, for the energy-budget figures. |
| `Sightlines.f90` | A *forward model*, not a property measurement: ray-traces synthetic quasar absorption sightlines through the simulated stream and hot halo, predicting the ion column densities and line-of-sight velocities a telescope would observe. Written for the code-comparison study of Hafen et al. 2024 (MNRAS 528, 39), which used simulations as a known ground truth to quantify the biases of common observational methods (not to reproduce specific real systems) — see the "From model to telescope" section of the [project README](../). |

The analytic estimates behind the cooling criterion (the governing timescale
ratios and the survival radius) were computed in a small MATLAB layer, not
included here; the more substantial analytic modelling is in the
[`cosmological_model/`](../cosmological_model/) directory of this project.

Raw simulation outputs are not included (too large); the figures show what the
pipeline produces.
