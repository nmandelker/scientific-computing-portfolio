# Scientific Computing Portfolio — Nir Mandelker

I'm a computational astrophysicist (PhD, Hebrew University; postdoctoral
fellowships at Yale and KITP/UC Santa Barbara; most recently faculty at the
Hebrew University of Jerusalem) with 14 years of experience building analytical 
and numerical models, simulation pipelines, and analysis tools for terabyte-scale 
scientific data. My work has led to over 60 peer-reviewed publications with over 5000 
citations. This repository is a curated showcase of the code behind that career, 
as I begin my transition from academia to industry research/data-science roles. 

**What this is:** research code, presented honestly. These tools were written
to answer scientific questions under deadline, by a physicist, for himself
and his collaborators, not as software products. What they demonstrate is
algorithm design from scratch, performance-critical programming against real
memory and compute limits, rigorous validation habits, and the ability to
explain hard technical material. Each project README is written for a
technical reader with no astrophysics background.

**Two kinds of work are showcased here.** In the giant-clumps projects I
built the measurement instruments, including detection, tracking, and analysis
pipelines applied to ~100 TB of simulations run by collaborators. In the
cold-streams project I designed and ran the numerical experiments myself,
from analytic theory to supercomputer simulations to analysis. Different
skills; both represented. In the cosmic web refinement project, I combined 
both skills, designing and running, together with my collaborators, cosmological 
simulations with novel refinement methods and building analysis pipelines 
from scratch to deal with ~10s of TB of data. 

## Projects

### The giant-clumps suite
Three connected projects studying dense, transient, star-forming structures
in galaxies as they were ~10 billion years ago:

| Project | What it is | Status |
|---|---|---|
| [Clump finder & tracker](Galaxy_Catalogues_and_Clump_Finder/) | End-to-end pipeline in Fortran 90: unsupervised object detection, multi-object tracking, and classification in 3D adaptive-resolution simulation data. Its catalogs later served as training data for neural-network studies of real Hubble Space Telescope images. | ✅ |
| [Clump evolution model](Clump_Evolution_Model/) | An analytic "bathtub" model of how clumps grow, migrate, and die, solved with 
8 lines of differential equations and validated against high-cadence tracking from the pipeline above and compared to HST observations (Dekel, Mandelker et al. 2022). | ✅ |
| [Turbulence decomposition](Turbulence_Decomposition/) | Helmholtz decomposition of turbulent velocity fields (Fortran + MKL FFTs), validated on synthetic turbulence with known ground truth, showing that clumps form where compressive turbulence is anomalously strong, resolving a puzzle where classical stability theory said clumps shouldn't form (Mandelker et al. 2025). | ✅ |

### Cold streams (design + execution of numerical experiments)
Analytic stability theory for cold gas streams feeding early galaxies,
tested with idealized high-resolution simulations I designed and ran with
the RAMSES code. Adding one physical ingredient at a time (gravity,
radiative cooling, magnetic fields), I identified at each step the
dimensionless number that determines whether the new physics changes the
answer. Included here is also a forward model of an observable signal 
(Lyman-alpha emission) later applied by independent observational teams, 
including in *Science*. — 🚧 in preparation

### AREPO / moving-mesh projects (Python)
Analysis pipelines for a different class of simulations, plus
adaptive-refinement techniques achieving ~100× effective resolution gain for
~5× cost in targeted regions, developed in collaboration with the code's
author and since adopted by multiple simulation groups. — 🚧 in preparation

## Languages and tools

Fortran 90 (performance-critical pipelines), MATLAB and Python/Matplotlib
(analysis and visualization), Mathematica (symbolic/analytic work), MPI-era
HPC cluster computing, Intel MKL, RAMSES / ART / AREPO simulation codes.

## More

- Academic site: [nirmandelker.com](https://nirmandelker.com) — research
  overview, publications, movies
- Publications: [Google Scholar](https://scholar.google.com/citations?user=iOHVGEoAAAAJ&hl=en)
