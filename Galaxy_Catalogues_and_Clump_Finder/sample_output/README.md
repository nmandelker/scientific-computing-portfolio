# Sample clump catalogs

Example output of the clump finder & tracker for one simulation (VELA V07,
thin-snapshot-spacing run). Each row is one clump detection in one snapshot;
tracked clumps keep the same ID across snapshots, so an object's history is
the set of rows sharing its ID.

| File | Contents |
|---|---|
| `in_situ.out` | Clumps classified as formed within the disc (in-situ) |
| `ex_situ.out` | Clumps classified as accreted (ex-situ; remnants of infalling satellites) |
| `bulge.out` | Detections associated with the central bulge, kept separate from the disc clump populations |
| `normalized_*.out` | Same rows, with clump properties divided by the host disc's corresponding property at that snapshot (mass by disc mass, velocity dispersion by disc dispersion, etc.), for comparisons across galaxies and cosmic times |

## Column format

Written by `src/thin_timesteps/clump_finder_combined.f90` (Fortran format
`(1x,i5.5,1x,f6.3,23(1x,es10.2),1x,i3,1x,i1,20(1x,es10.2))`):

| # | Quantity | Units / notes |
|---|---|---|
| 1 | Clump ID | persistent across snapshots for tracked clumps |
| 2 | Redshift | dimensionless (cosmic time stamp; larger = earlier) |
| 3 | Clump radius | kpc |
| 4 | Gas mass | M☉ |
| 5 | Stellar mass | M☉ |
| 6 | Total baryonic mass | M☉ |
| 7 | Gas fraction | of baryonic mass |
| 8 | Dark-matter mass fraction | of total mass |
| 9 | Gas surface-mass density | M☉/pc² |
| 10 | Stellar surface-mass density | M☉/pc² |
| 11 | Baryonic surface-mass density | M☉/pc² |
| 12 | Mean stellar age | Myr |
| 13 | Gas metallicity | log [O/H]+12 |
| 14 | Stellar metallicity | log [O/H]+12 |
| 15 | Star-formation rate (SFR) | M☉/yr |
| 16 | SFR surface density | M☉/yr/kpc² |
| 17 | Specific SFR | 1/Gyr |
| 18 | Gas depletion time | Gyr |
| 19 | Radial position | units of disc radius R_d |
| 20 | Vertical position | units of disc half-height H_d (signed) |
| 21 | Radial position | kpc |
| 22 | Vertical position | kpc (signed) |
| 23 | Mean overdensity δρ/ρ | detection significance of the clump |
| 24 | η (shape parameter) | min/max eigenvalue of inertia tensor |
| 25 | Clump/background dark-matter density ratio | key in-situ/ex-situ discriminant |
| 26 | Ex-situ flag | integer |
| 27 | Merger flag | 1 if the clump underwent a clump-clump merger |
| 28 | Internal Clump Free-fall time | Myr |
| 29 | Disc dynamical time at clump position | Myr |
| 30 | Global disc dynamical time | Myr |
| 31–33 | Gas mass inflow rate | three measurement radii |
| 34–43 | Gas mass outflow rate | ten velocity/radius bins |
| 44 | Stellar mass inflow | M☉ |
| 45 | Stellar mass outflow | M☉ |
| 46 | Stellar mass formed in place | M☉ |
| 47 | Virial parameter α_vir | ~1 for self-gravitating bound objects |

In the `normalized_*` files, columns 4–18 are expressed relative to the host
disc's value at the same snapshot instead of in physical units.
