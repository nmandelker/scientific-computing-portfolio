# MATLAB post-processing — turbulence decomposition

Post-processing of the binary grid outputs produced by the Fortran
decomposition pipeline in `../src/`, presented as-is.

| Script | Role |
|---|---|
| `sol_over_comp.m` | Time series of the compressive-to-solenoidal power ratio across a simulation's history — the machinery behind Mandelker et al. 2025, Fig. 2. |
| `turbulence_density_maps.m` | Maps of gas surface density alongside the local converging-flow energy fraction — the machinery behind Mandelker et al. 2025, Fig. 3. |

The proto-clump statistics (Mandelker et al. 2025, Fig. 4) combine these
outputs with the proto-clump positions traced by the clump finder project's
`thin_timesteps/proto_clump_positions*.f90`.
