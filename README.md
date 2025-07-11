# IGE-Stochastic
Repository for my 6 months internship at the IGE on the simulation of uncertainties in an oceanic circulation model. Study focused on the South-West Indian ocean and used the model **CROCO** (*Coastal and Regional Ocean COmmunity model*).

We looked into the impact of intrinsic variability, as well as the wind stress and vertical turbulent mixing perturbation. For that we used stochastic processes that acted as the uncertainties on the terms we wanted to study. I personally develop the turbulent mixing scheme perturbation by modifying the Fortran source code of the model and i made some diagnostics.

## Structure of the repo
- [bibliography](bibliography/) contains documents and main papers used for bibliography.
- [book](book/) Contains some weekly reports.
- [documentation](documentation/) contains diverse code notes.
- [report](report/) contains the report.
- [scripts](scripts/) contains all the scripts.
- [slides](slides/) contains the slides.

## ADASTRA login
To connect to ADASTRA 
```bash
ssh adastra
```

Then I enter the password and activate my python environment.
```bash
source ./python_environment/bin/activate
``` 

## Links
- [Adastra documentation](https://dci.dci-gitlab.cines.fr/webextranet/)
- [Adastra access documentation](https://dci.dci-gitlab.cines.fr/webextranet/user_support/index.html#adastra-accessing-account-opening)
- [Croco documentation](https://croco-ocean.gitlabpages.inria.fr/croco_doc/)
- 