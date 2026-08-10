# Workflow Scripts

Where the per-utility folders each demonstrate a single utility, the scripts
here combine several AnalysisTools utilities, and often outside programs, into a
complete analysis or visualisation workflow.

- **01-VisualizeAggregates**
  Renders a VMD snapshot per timestep from an `Aggregates` output file, giving
  each aggregate its own colour, and crops the results with ImageMagick.

## Scope and intent

Because these workflows reach outside AnalysisTools, they need more than a built
`AnalysisTools` to run - VMD, ImageMagick, gnuplot and the like. Each script
states its own requirements in its `README.txt` and exits with a clear message
when something is missing, rather than failing halfway through.

Shared shell helpers live in `Examples/func.sh`; scripts source it relative to
their own location, so they can be run from any directory.

`snap.tcl` in this folder is the shared VMD helper: it renders the current scene
and converts the result with ImageMagick. Both `01-VisualizeAggregates` and
`AddToSystem/05-IncrementalBuild` source it, so the VMD-based examples need
nothing outside the repository beyond `vmd` and `magick` themselves.
