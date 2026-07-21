# Cosmological forward model — Mandelker et al. 2020b

The analytic forward model that turns the idealized-simulation results into a
prediction of an observable signal (Lyman-α emission) as a function of halo
mass and redshift. MATLAB, presented as-is. Unlike the other stages, this one
ran no new simulations — it is pure analytic modelling on top of the
scaling relations measured from the cooling simulations.

## The model, from inputs to observable

| File | Role |
|---|---|
| `params_vs_mass.m` | Sets the stream's properties on entering the halo (density contrast δ, radius, velocity) as a function of halo mass and redshift, from cosmological accretion. |
| `halo_decel_cooling_NFW.m` | The core model: the coupled equations for how a stream decelerates and grows (by cooling-driven entrainment) as it falls through an NFW dark-matter halo potential. |
| `solve_halo_decel_cooling_panels_NFW.m` | Driver: integrates the model over a grid of halo mass and redshift and produces the multi-panel figures (velocity, mass, luminosity). |
| `Ltot.m`, `Ltot_panels_NFW.m` | Convert the dissipated energy into the total emitted luminosity (mostly Lyman-α) and plot it vs. halo mass and redshift. |
| `Rscrit_vs_mass.m` | The survival criterion: the ratio of stream radius to the critical radius for cooling to beat disruption — i.e. which streams reach the galaxy. |
| `Mstream_Daddi.m` | Comparison of the modelled stream mass against the observational constraints of Daddi et al. |
| `cooling_equilibrium.m`, `find_cooling_equilibrium2.m`, `xHI.m` | Helpers: thermal-equilibrium temperature, cooling-curve evaluation, and neutral-hydrogen fraction. |

## Why this piece matters beyond astrophysics

This is a compact example of the operation at the heart of a lot of applied
scientific and engineering work: take a physical model of a state you cannot
directly see, forward-model it into a quantity an instrument *can* measure, and
hand that prediction to others to test against data. The model here was applied
by independent observational teams (Daddi et al. 2022b; Emonts et al. 2024,
*Science*) — the prediction survived contact with real measurements.
