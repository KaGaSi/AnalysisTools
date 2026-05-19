# AddToSystem — Usage Examples

This directory contains examples showing how `AddToSystem` can be used in
practice: from minimal, single-command invocations to more realistic,
multistep system construction workflows.

These examples are not meant as tutorials. They assume you have already looked
through the manual and are comfortable with the basic command-line interface,
including bash scripting (somewhat advanced bash scripting for the more complex
examples 03 and 04).

## What you’ll find here

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

## Assumptions

A few practical assumptions are made throughout:
- Examples should be run from within their own directories because intermediate
  files (e.g., coordinate and FIELD files) are created and left in place.
- FIELD files are often treated as simple templates with placeholders (e.g.
  particle counts) replaced using standard shell tools.

## How to approach the examples

If you’re new to the utility, start with the minimal example and then look at
the incremental build. The bilayer and cylinder examples are best read as a
reference for a realistic workflow, rather than something to type through line
by line.
