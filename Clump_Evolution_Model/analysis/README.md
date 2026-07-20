# MATLAB analysis scripts — clump evolution model

The complete analysis layer of Dekel, Mandelker et al. 2022 (D22), presented
as-is.

| Script | Role |
|---|---|
| `clump_evolution_exact.m` | **The model itself**: the three coupled ODEs (gas mass, stellar mass, migration radius) in eight lines. Everything below calls this. |
| `fig_clump_bathtub.m` | Integrates the model to produce the model-solution curves — D22 Fig. 1. |
| `launch_prop_vs_t_over_tdyn.m`, `prop_vs_t_over_tdyn_rectangle.m` | Driver + the central analysis function: applies sample cuts to the tracked clump histories, aligns and stacks them, smooths, and integrates the model inside the plotting routine to overlay model curves on the simulation stacks — D22 Figs. 4–7. |
| `read_Guo2.m` | The observational test: ingests the Guo et al. HST/CANDELS clump catalog, applies selection cuts, and compares observed trends against the model — D22 Fig. 10. |
| `launch_specific_rates.m`, `specific_rates.m` | Direct measurement of the model's rate parameters (accretion, star formation, outflow, stripping) from the simulations. |
| `launch_clump_property_histograms.m`, `clump_property_histograms.m` | Clump property distributions. |
| `load_data.m`, `load_gen3.m`, `add_properties.m` | Catalog ingestion (the format is documented in the clump finder project's `sample_output/`). |
| `clump_evolution_exact2.m` | Unpublished exploration: maps final clump properties, interpolated to the moment of arrival at the disc center, across the feedback-strength parameter space. |

Same era and style as the rest of the MATLAB in this repository: long
functions, positional arguments, column indices documented in comments.
