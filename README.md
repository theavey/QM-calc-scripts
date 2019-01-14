This repository contains some scripts and a package related to running and 
analyzing Gaussian calculations (and possibly a little for Q-Chem). 
There is also some older code for looking at the convergence of metadynamics
calculations run with PLUMED.

Two dependencies of gautools that might be harder to find and will likely not
install automatically are:
  * [ParaTemp](https://github.com/theavey/ParaTemp)
  * [thtools](https://github.com/theavey/thtools)

There was previously a script to help setup parallel tempering calculations,
but that was removed because a newer version is available in 
[ParaTemp](https://github.com/theavey/ParaTemp).
