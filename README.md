# AnalysisTools

A set of utilities to analyse trajectories from particle-based molecular
simulations. The 34 utilities (including one from the [DL_MESO simulation
package](https://www.scd.stfc.ac.uk/Pages/DL_MESO.aspx)) could be roughly
divided into four types:

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

Either download the latest zipped/tarred version
[here](https://github.com/KaGaSi/AnalysisTools/tags), or clone this repository:

`git clone https://github.com/KaGaSi/AnalysisTools`

AnalysisTools requires `C` and `Fortran` compilers and `cmake`.  The
compilation should be done in a separate directory, such as `<root>/build`
(where `<root>` is the AnalysisTools root directory). The following
generates a `Makefile` within `build`:

`mkdir <root>/build; cd <root>/build; cmake ../`

Then run `make` to compile all utilities or `make <utility name>` to
compile a single utility.

The binaries are located in the `<root>/build/bin` directory.

A reference manual is included in the repository (manual.pdf).
