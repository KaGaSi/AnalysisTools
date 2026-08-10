This script renders one VMD snapshot per timestep from an aggregate file,
giving every aggregate its own colour, and crops the results with ImageMagick.

Unlike the per-utility examples, scripts in this folder combine AnalysisTools
with outside programs, so they need more than a built AnalysisTools to run.

Inputs

  VisAgg.sh [<in.agg> [<in.vtf>]]

  <in.agg>  aggregate file produced by the Aggregates utility
            (default: red.agg)
  <in.vtf>  matching structure and coordinate file
            (default: join.vtf)

  Both default to names from the original run; pass your own on the command
  line.

  Outputs: snapshots are numbered per frame and written to the current
  directory (0001.tga, then cropped to 0001.jpg), so run the script from
  wherever you want them. The temporary tcl script handed to vmd is named
  after the aggregate file's stem and removed after each frame.

  The .agg parsing follows the AnalysisTools v4.0 format: a 'Step:' line, a
  count, then two lines per aggregate (core and border molecules). The final
  'Last Step: <num>' line gives the number of frames to render.

  The two kinds of molecule get a representation each, sharing the aggregate's
  colour: core molecules opaque and border ones transparent. Aggregates of a
  single molecule are drawn as glass instead. Note that with the Aggregates
  utility's default of -nbr 1, every molecule is a core one, so the border
  representation only appears for files made with a higher -nbr.

External requirements

  vmd     must be in $PATH
  magick  ImageMagick, must be in $PATH; snap.tcl uses it to convert and crop
          each rendered frame
  vmdrc   optional VMD startup file, taken from $VMDRC or ~/.config/vmdrc. It
          is only passed to vmd when readable, so the script runs without one.

  The 'snap' proc that renders and writes each frame comes from ../snap.tcl,
  which ships with AnalysisTools, so nothing outside the repository is needed
  beyond vmd and magick. Set VMD_SCRIPT_DIR to use a different copy:

    VMD_SCRIPT_DIR=~/my/vmd/scripts ./VisAgg.sh run.agg run.vtf

Things worth adjusting

  Near the top of the script: canvas resolution (xres/yres), whether snapshots
  are trimmed, the physical snapshot size and dpi, and the range of frames.
  The molecule selections use 'name A C D' to accentuate a hydrophobic core;
  change those in the add_rep function to match your own bead names, along with
  the drawing style shared by all representations. The materials picked for the
  core and the border molecules are set in the loop over aggregates.

  The final magick call crops a fixed region (1600x1600+2160+2125) that suits
  the original canvas size. If you change xres/yres, that crop needs redoing.
