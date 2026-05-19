# AddToSystem Usage Examples

The examples are ordered roughly from simple to complex:

- **01-ManualExamples**
  Minimal, single-invocation examples showing the required inputs and default
  behaviour as presented in the AnalysisTools manual (Fig. 3.1 TODO:?).
- **02-IncrementalBuild**
  Illustrates how systems can be grown incrementally by chaining coordinate
  files and repeatedly using `--add` combined with other options.
- **03-AssymetricBilayer**
  A more realistic, multistep workflow that builds an asymmetric bilayer using
  spatial constraints, mixed molecule types, and composition gradients.
- **04-Cylinder**
  Another multistep workflow that builds a cylindrical aggregate spanning the
  simulation box.
- **05-LibraryBuild**
  Builds a system from scratch using `-lib`/`-mol`/`-ntot`, then generates a
  complete FIELD file with DPD interactions using `Info -lib`.  Also shows
  `-sys` for writing and re-using a system info file when extending an existing
  system with `-lib`.  Uses the shared molecule library in `Examples/library/`.
- **06-ExchangeBeads**
  Shows exchange mode (`-xb`), distance constraints relative to bonded beads
  (`--bonded`), placement anchored to the molecule's tail bead (`--tail`), and
  reproducible placement (`-s`).
- **07-RealCoordinates**
  Shows `--real` with negative box coordinates (LAMMPS-style centred box) and
  `-ebt` to reserve extra bead type slots in a LAMMPS data file output.

## Scope and intent

These examples show how different pieces fit together when building up a system
with repeated calls to the utility.

They are deliberately not concerned with:
- defending particular parameter choices or molecule counts,
- validating physical realism,
- optimising performance,
- recommending "best" simulation setups.

Those decisions are left to the user and their use cases. Note that the examples
02 through 04 give valid configurations for classic dissipative particle
dynamics simulations as they assume particle density of 3.
