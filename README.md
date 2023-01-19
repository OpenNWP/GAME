
# GAME

The Geophysical Fluids Modeling Framework (GAME) is a non-hydrostatic hexagonal C-grid dynamical core with the possibility to advect a variable number of constituents. For radiation, it is coupled to RTE+RRTMGP.

## Overview

It is known that the forecast skill of a NWP model depends more on physics and data assimilation than on the dynamical core. However, all dynamical cores have inconsistencies even in the most fundamental dynamical quantities (mass, energy forms and entropy). That is why the aim of this project is to develop a next generation dynamical core with the following properties:

* stability
* the numerical dispersion relation shall contain no unphysical branches
* the numerical dispersion relation shall contain all relevant physical branches, including a geostrophic mode
* strong scalability on massively parallel computer architectures
* mass conservation to machine precision
* energy conversions shall be based on a spatial discretization of Poisson brackets
* energy conservation in the spatial discretization
* consistent friction and entropy production
* absence of unphysical numerical stabilizers like divergence damping, fixers and filters
* no stability problems with terrain following coordinates
* ellipsoidal grid geometry
* a capable and flexible framework for coupling to physics and to other components of an Earth system model
* good self-consistency also in the presence of multiple constituents and radiation

No discretization is able to achieve all these goals simultaneously. However, a hexagonal C-grid is probably one of the best choices when it comes to finding a reasonable balance of all these aspects.

## Documents

### Scientific derivations

The derivations of both the continuous equations as well as of the discretization techniques can be found in this textbook on theoretical meteorology (in German): [Kompendium Theoretische Meteorologie](https://raw.githubusercontent.com/MHBalsmeier/kompendium/main/kompendium.pdf). The most fundamental numerical techniques have firstly been published in the literature cited below.

### Handbook

The handbook of the model can be found in the subdirectory handbook. It is to be understood as a user manual and contains information on how to generate grid files and explains how to configure, compile and run the model.

## Installing the model

It is recommended to run the model on Linux. These installation instructions are tested for Ubuntu, for other Linux distributions they might have to be modified.

### Dependencies

Everything is easy and quick to install.

	sudo apt-get install gfortran make cmake wget python3-pip libnetcdff-dev

* Clone the RTE+RRTMGP repository: `git clone https://github.com/earth-system-radiation/rte-rrtmgp`

#### For using the plotting routines

The following packages are additionally required if you want to make use of the plotting routines:

* Python visualization library scitools-iris (installation manual: https://scitools-iris.readthedocs.io/en/latest/installing.html#installing-from-source-without-conda-on-debian-based-linux-distros-developers)
* FFMPEG (`sudo apt-get install ffmpeg`)

### Download

	git clone https://github.com/OpenNWP/GAME
	cd GAME
	./create_directories.sh

### Compiling the code

	./compile.sh -f

Before being able to run the model, however, you need to create a grid file. The handbook gives further information on that. A first try could be

	cd grid_generator
	./compile.sh
	./run_script.sh

## Fundamental literature

* Staniforth, A. and Thuburn, J. (2012), Horizontal grids for global weather and climate prediction models: a review. Q.J.R. Meteorol. Soc., 138: 1-26. doi:10.1002/qj.958
* Thuburn, John. (2008). Numerical wave propagation on the hexagonal C-grid. Journal of Computational Physics. 227. 5836-5858. 10.1016/j.jcp.2008.02.010. 
* Gassmann, Almut & Herzog, Hans-Joachim. (2008). Towards a consistent numerical compressible non‐hydrostatic model using generalized Hamiltonian tools. Quarterly Journal of the Royal Meteorological Society. 134. 1597 - 1613. 10.1002/qj.297.
* Thuburn, John et al. “Numerical representation of geostrophic modes on arbitrarily structured C-grids.” J. Comput. Phys. 228 (2009): 8321-8335.
* Ringler, Todd & Thuburn, John & Klemp, J. & Skamarock, W.C.. (2010). A unified approach to energy conservation and potential vorticity dynamics on arbitrarily structured C-grids. J. Comput. Physics. 229. 3065-3090. 10.1016/j.jcp.2009.12.007.
* Gassmann, A. (2013), A global hexagonal C‐grid non‐hydrostatic dynamical core (ICON‐IAP) designed for energetic consistency. Q.J.R. Meteorol. Soc., 139: 152-175. doi:10.1002/qj.1960
* Gassmann, A. Discretization of generalized Coriolis and friction terms on the deformed hexagonal C‐grid. Q J RMeteorol Soc. 2018; 144: 2038– 2053. https://doi.org/10.1002/qj.3294



