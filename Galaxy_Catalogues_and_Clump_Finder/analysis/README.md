# Selected MATLAB analysis scripts

The statistical analysis and figure production for the clump papers was done
in MATLAB, on top of the catalogs produced by the Fortran pipeline. This is a
small representative selection from a much larger collection (~60 scripts) —
chosen to show the analysis layer end to end, presented as-is.

| Script | Role |
|---|---|
| `load_data.m` + the eight scripts it calls (`load_gen2.m`, `load_gen3.m`, `gen2_vs_gen3_galaxy_discs_spheres.m`, `post_referee_properties.m`, `add_properties.m`, `common_sample_gen2_3.m`, `nis2_common_and_friends.m`, `load_Behroozi.m`) | The complete data-ingestion chain: loads the in-situ / ex-situ / bulge clump catalogs (raw and disc-normalized) and galaxy/disc catalogs for all 35 galaxies of both simulation generations, harmonizes them into a common comparison sample, attaches derived properties, and loads an external observational calibration (the Behroozi et al. stellar-to-halo-mass relation). `load_gen3.m` also serves as a working reference for parsing the files in `../sample_output/`. The underlying catalog files are not included here (samples are in `../sample_output/` and `../galaxy_catalog_sample/`). |
| `face_on_clumpy_images.m` | Renders the disc images with detected clumps overlaid, encoding clump class, mass, and size in the marker symbols — the script behind Fig. 1 of Mandelker et al. 2017, i.e. the header image of the main README. |
| `launch_gradients_paper.m`, `clump_gradients_paper.m` | Driver + analysis function for radial gradients: how clump properties vary with position in the host disc, with configurable sample cuts (mass, lifetime, height, redshift) — the machinery behind the gradient figures of Mandelker et al. 2017 (Fig. 15 and their Appendix B counterparts). |
| `launch_classification_v4_parameters.m`, `classification_v4_parameters.m` | Driver + analysis function behind the clump mass–size plane analysis of Mandelker et al. 2017 (Fig. 6). The smoothing scale F_W and detection threshold δ_min of the catalog under analysis are input parameters, so the same machinery could be pointed at catalogs built with different finder settings. |
| `launch_comparative_clump_mass_function.m`, `comparative_clump_mass_function_multiplot.m` | Driver + analysis function for the clump mass and SFR functions compared across samples — Fig. 8 of Mandelker et al. 2017. |

A note on style: these scripts are working research code from the same era as
the Fortran pipeline — long functions, positional arguments, and column
indices documented in comments rather than named constants. The column
conventions match the table in `../sample_output/README.md`.
