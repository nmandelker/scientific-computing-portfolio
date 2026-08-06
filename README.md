# Scientific Computing Portfolio — Nir Mandelker

I'm a computational astrophysicist (PhD, Hebrew University; postdoctoral
fellowships at Yale, the Heidelberg Institute for Theoretical Studies, 
and KITP/UC Santa Barbara; Assistant Professor at the Hebrew University 
of Jerusalem, Visiting Scholar at UW Seattle) with 14 years of experience 
building analytical and numerical models, simulation pipelines, and analysis 
tools for terabyte-scale scientific data. My work has led to over 60 peer-reviewed 
publications with over 6000 citations and an h-index of 36. This repository 
is a curated showcase of the code behind that career. 

**What this is:** In this repo, I present my code from my various research 
projects, as I wrote it at the time. These tools were written to answer scientific 
questions under deadline, by a physicist, for himself and his collaborators, 
not as software products. What they demonstrate is algorithm design from scratch, 
performance-critical programming against real memory and compute limits, rigorous 
validation habits, and the ability to explain hard technical material. Each project 
README is written for a technical reader with no astrophysics background.

**Two kinds of work are showcased here.** In the giant-clumps projects I
built the measurement instruments, including detection, tracking, and analysis
pipelines applied to ~100 TB of simulations run by a collaborator. In the
cold-streams project I designed and ran the numerical experiments myself,
from analytic theory to supercomputer simulations to analysis. This entire 
pipeline is represented in the attached code and project outcomes. In the 
cosmic web refinement project, I combined both skills. Together with my collaborators, 
I designed and ran cosmological simulations with novel refinement methods 
and built analysis pipelines from scratch to deal with ~10s of TB of data. 

## Projects

### The giant-clumps suite
Three connected projects studying dense, transient, star-forming structures
in disk galaxies as they were ~10 billion years ago, roughly 2-3 billion years 
after the Big Bang:

| Project | What it is |
|---|---|
| [Clump finder & tracker](Galaxy_Catalogues_and_Clump_Finder/) | End-to-end pipeline in Fortran 90: unsupervised object detection, multi-object tracking, and classification in 3D adaptive-resolution simulation data. Its catalogs later served as training data for neural-network studies of real Hubble Space Telescope images. |
| [Clump evolution model](Clump_Evolution_Model/) | An analytic model of how clumps grow, migrate, and die, solved with 8 lines of differential equations and validated against high-cadence tracking from the pipeline above and compared to HST observations (Dekel, Mandelker et al. 2022). |
| [Turbulence decomposition](Turbulence_Decomposition/) | Helmholtz decomposition of turbulent velocity fields in Fortran 90, validated on synthetic turbulence with known ground truth, showing that clumps form where compressive turbulence is anomalously strong, resolving a puzzle where classical stability theory said clumps shouldn't form (Mandelker et al. 2025). |

### [Cold streams](Cold_Stream_Stability/)
A multi-year, multi-step, international collaboration I initiated and spearheaded. 
I developed an analytic theory for the evolution and stability properties of cold 
gas streams feeding early massive galaxies, and tested this with idealized high-resolution 
simulations I designed and ran with the RAMSES code. Adding one physical ingredient 
at a time (hydrodynamics, gravity, radiative cooling, magnetic fields), I identified 
at each step the dimensionless number that determines whether the new physics changes the
answer. Included here is also a forward model of an observable signal (Lyman-alpha emission) 
later applied by independent observational teams, including in a *Science* paper I co-authored.

### [Cosmic web refinement](Cosmic_Web_Refinement/)
I developed novel refinement methods for the moving-mesh code AREPO, in collaboration 
with the code's author, that target resolution in the diffuse gas around galaxies that 
every previous simulation left under-resolved. We gained up to ~100× more resolution 
elements for ~5× the compute. The new resolution revealed physics that lower-resolution
simulations miss entirely. Included here is the Python analysis behind that
result, and additional studies made possible by the enhanced resolution. This includes 
Lagrangian tracking of gas through a mesh with no persistent cells, and a 
study where an interpolation scheme was rejected because it failed validation tests 
in a way that became worse as resolution increased. This scheme was then replaced and 
the alternative scheme's error was bounded by an independent check. 

### [LLM confidence probing](LLM_Confidence_Probing/)

This is not from my academic career, and not a contribution to the literature. This 
is the result of three days spent reading up on transformers and working through the 
PyTorch tutorials, then attempting something real with what I had learned. Four experiments 
on a small open-weights model, asking where inside the network a rephrased question
changes the answer. I included this because it shows how my research methodology transferred 
to a domain I had started learning that week. I applied the same habits as the projects above, 
isolated the confounding variables, distrusted the encouraging first result, and reported the weaker 
one. I wrote the code with the help of an AI coding assistant, but all of the analysis and 
interpretation are entirely my own. This is the only project here that anyone can run from start to finish.

## Languages and tools

Fortran 90 (performance-critical pipelines), MATLAB and Python/Matplotlib
(analysis and visualization), Mathematica (symbolic/analytic work), PyTorch and
Hugging Face transformers, MPI-era HPC cluster computing, Intel MKL,
RAMSES / ART / AREPO simulation codes.

## More

- Academic site: [nirmandelker.com](https://nirmandelker.com) — research
  overview, publications, movies
- Publications: [Google Scholar](https://scholar.google.com/citations?user=iOHVGEoAAAAJ&hl=en)
