## SUMMER (python)
Scalable, Universal Mathematical Modelling in R.

This repository is a Python-coded version of the jtrauer/summer repository, coded in R.

# Philosophical approach
The purpose of this repository and that of SUMMER is to provide a codebase to facilitate the construction of
mathematical models simulating infectious disease transmission.
Many aspects of models of infectious disease transmission are built using higher-level packages.
These include numeric integrators (e.g. ODE solvers), packages to manipulate data structures, Bayesian calibration tools
and graphing packages to visualise outputs.
However, it is our group's experience that pre-built packages are less frequently used for defining the model structures
themselves (although this is not to say that this is never done).
That is, it remains common modelling practice to define compartmental models of infectious disease transmission by hand-
coding a series of ordinary differential equations representing the rates of change of each modelled population
compartment.

