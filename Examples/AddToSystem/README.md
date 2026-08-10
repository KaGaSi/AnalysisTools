# AddToSystem Usage Examples

The examples are ordered roughly from simple to complex:

- **01-ManualExamples**
  Minimal, single-invocation examples showing the required inputs and default
  behaviour as presented in the AnalysisTools manual (Fig. 3.1 TODO:?).
- **02-LibraryBuild**
  Builds a system from scratch using `-lib`/`-mol`/`-ntot`, then extends it with
  more molecules.  Also shows `-sys` for writing and re-using a system info file
  when working with `-lib`, and the percentage form of `-mol`.  Uses the shared
  molecule library in `Examples/library/`.
- **03-ExchangeBeads**
  Shows exchange mode (`-xb`), distance constraints relative to bonded beads
  (`--bonded`), placement anchored to the molecule's tail bead (`--tail`), and
  reproducible placement (`-s`).
- **04-RealCoordinates**
  Contrasts the two coordinate spaces of the `-cx`/`-cy`/`-cz` constraints,
  fractions of the box and `--real`, which are interchangeable for an orthogonal
  box (here a LAMMPS-style centred one) but not for a tilted one.
- **05-IncrementalBuild**
  Illustrates how systems can be grown incrementally by chaining coordinate
  files and repeatedly using `--add` combined with other options.
- **06-AssymetricBilayer**
  A more realistic, multistep workflow that builds an asymmetric bilayer using
  spatial constraints, mixed molecule types, and composition gradients.
- **07-Cylinder**
  Another multistep workflow that builds a cylindrical aggregate spanning the
  simulation box.

## Scope and intent

These examples show how different pieces fit together when building up a system
with repeated calls to the utility.

They are deliberately not concerned with:
- defending particular parameter choices or molecule counts,
- validating physical realism,
- optimising performance,
- recommending "best" simulation setups.

Those decisions are left to the user and their use cases. Note that the examples
05 through 07 give valid configurations for classic dissipative particle
dynamics simulations as they assume particle density of 3.
