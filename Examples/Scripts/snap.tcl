################################################################################
# Render the current VMD scene to an image file.
#
# Sourced by the example scripts that drive VMD (Scripts/01-VisualizeAggregates
# and AddToSystem/05-IncrementalBuild). Requires ImageMagick's 'magick' in the
# $PATH; the .tga that VMD renders is converted to .jpg and then removed.
#
#   snap <name> <trim> [<size_cm>] [<dpi>] [<dir>]
#     name     output file name, ending in .tga
#     trim     1 to trim surrounding whitespace (leaving a 2px border), else 0
#     size_cm  physical size of the image in centimetres; empty to skip resizing
#     dpi      resolution used together with size_cm (default: 300)
#     dir      which dimension size_cm refers to: 'x' (default) or 'y'
#
#   Returns 0 on success, -1 on failure.
################################################################################
proc snap {name trim_flag {size_cm ""} {dpi 300} {dir "x"}} {
  # render the snapshot using TachyonInternal
  render TachyonInternal $name
  # prepare the ImageMagick command as a list
  set cmd [list magick $name]
  # trim whitespace if requested, leaving 2px border
  if {$trim_flag == 1} {
    lappend cmd -trim -bordercolor white -border 2
  }
  # if size_cm is provided (not empty), resize image accordingly
  if {$size_cm ne ""} {
    set pixels [expr {$size_cm / 2.54 * $dpi}]
    if {$dir eq "x"} {
      # width specified, height auto-adjusted
      lappend cmd -units PixelsPerInch -resize ${pixels}x -density $dpi
    } elseif {$dir eq "y"} {
      # height specified, width auto-adjusted
      lappend cmd -units PixelsPerInch -resize x${pixels} -density $dpi
    } else {
      puts stderr "Invalid direction '$dir'. Use 'x' or 'y'."
      return -1
    }
  }
  # convert output filename from .tga to .jpg
  set outputname [string map {".tga" ".jpg"} $name]
  lappend cmd $outputname
  # execute the command and catch errors
  if {[catch {eval exec $cmd} err]} {
    puts stderr "Error executing ImageMagick command: $err"
    return -1
  }
  # remove the original .tga file
  if {[catch {exec rm $name} err]} {
    puts stderr "Warning: could not remove $name: $err"
  }
  return 0
}
