# AnalysisTools

A set of utilities to analyse trajectories from particle-based molecular
simulations. The utilities could be roughly divided into four types:

* utilities to calculate system-wide properties; e.g., pair
  correlation functions or particle densities along a simulation box axis
* utilities to calculate per-molecule or per-aggregate properties
  (where aggregate stands for any supramolecular structure); e.g., shape
  descriptors for individual molecules or whole micelles
* utilities to manipulate a configuration; e.g., create initial
  configuration for a simulation from scratch or by adding molecules to
  an existing one
* helper utilities to analyse text files; e.g., calculate averages
  and standard deviations of a data series

Installation
============

Requirements: GSL (GNU Scientific Library), cmake, C and Fortran compilers

Download from [tags](https://github.com/KaGaSi/AnalysisTools/tags) or
clone this repository.

To compile all utilities in a `build` directory,

`cd /path/to/build; cmake /path/to/AnalysisTools/; make`
