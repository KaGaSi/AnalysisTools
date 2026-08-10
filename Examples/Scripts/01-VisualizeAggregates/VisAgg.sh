#!/usr/bin/env zsh

################################################################################
# Render VMD snaphsots of aggregates (one colour per aggregate) for each
# timestep from .agg and .vtf files. The script's parsing of the .agg file is
# based on the AnalysisTools v4.0 format
#
# Usage: ./VisAgg.sh [<in.agg> [<in.vtf>]]
#
# See README.txt for the external programs and files this needs.
################################################################################

# exit on any command failure, error on unset variables
set -euo pipefail

# make the script runnable from anywhere
ROOT="${0:A:h}"
# shared helper functions (zpad)
source "${ROOT}/../../func.sh"

################################################################################
# file names and other variables
################################################################################
# input files; override on the command line
agg="${1:-red.agg}" # aggregate file (generated via Aggregates)
vtf="${2:-join.vtf}" # structure and coordinate file
# expand to full paths (necessary if a path starts with ~)
agg=${agg:A}
vtf=${vtf:A}
# derived names: strip the extension and reuse the stem
name=${agg:r}
tcl=${name}.tcl # tcl script for vmd
################################################################################
# external dependencies
#
# snap.tcl provides the 'snap' proc that renders and saves each frame; the copy
# in Examples/Scripts/ is used by default. Set VMD_SCRIPT_DIR to point at a
# different one. The VMD startup file is optional - it is only passed on if it
# is readable.
################################################################################
vmd_dir="${VMD_SCRIPT_DIR:-${ROOT}/..}"
# absolute, because vmd resolves it from wherever the script was run
vmd_dir=${vmd_dir:A}
vmdrc="${VMDRC:-${HOME}/.config/vmdrc}"
# basic error checks
[[ -r ${agg} ]] || { print "Error: File not readable: ${agg}" >&2; exit 1 }
[[ -r ${vtf} ]] || { print "Error: File not readable: ${vtf}" >&2; exit 1 }
[[ -r ${vmd_dir}/snap.tcl ]] || {
  print "Error: Missing snap.tcl in ${vmd_dir}; set VMD_SCRIPT_DIR" >&2
  exit 1
}
command -v vmd >/dev/null || { print "Error: vmd not in \$PATH" >&2; exit 1 }
command -v magick >/dev/null || {
  print "Error: magick (ImageMagick) not in \$PATH" >&2
  exit 1
}
# resolution (size of vmd canvas in pixels)
xres="5000"
yres="5000"
# should resultant snapshot be trimmed?
trim="0"
# size of the snapshot in centimetres
snap_size_cm="50"
# final snapshot resolution
dpi="300"
# first snapshot
first=1
# last snapshot - extract <num> from the last line: 'Last Step: <num>'
last=$(awk 'END{print $3}' "${agg}")
[[ ${last} == <-> ]] || { print "Invalid last timestep: ${last}" >&2; exit 1 }
# append one vmd representation to the tcl script //{{{
# $1 ... selection (e.g., 'resid 1 2 3')
# $2 ... vmd ColorID
# $3 ... vmd material (Opaque, Transparent, Glass1, ...)
add_rep() {
  tcl_script+=("set rep [expr \$rep + 1]")
  tcl_script+=("mol addrep \${mol}")
  tcl_script+=("mol modselect   \${rep} \${mol} ${1} and name A C D")
  tcl_script+=("mol modstyle    \${rep} \${mol} cpk 0.7 0.5")
  tcl_script+=("mol modcolor    \${rep} \${mol} ColorID ${2}")
  tcl_script+=("mol modmaterial \${rep} \${mol} ${3}")
} #}}}
for ((i = first; i <= last; i++)); do
  # numbered snapshot names
  pic=$(zpad "${i}" 4).tga
  echo "${pic}"
  ##############################################################################
  # start of the tcl script for vmd
  ##############################################################################
  step=$(( ${i} - 1 ))
  tcl_script=()
  tcl_script+=("set dir ${vmd_dir}")
  # 'pbc box' comes from this plugin; vmd does not load it by default when
  # running headless (-dispdev text) without a startup file
  tcl_script+=("package require pbctools")
  tcl_script+=("pbc box -width 1")
  tcl_script+=("axes location Off")
  tcl_script+=("display depthcue on")
  # tcl_script+=("display projection orthographic")
  tcl_script+=("display projection perspective")
  tcl_script+=("scale by 0.7")
  # tcl_script+=("rotate x by 10")
  # tcl_script+=("rotate y by 30")
  # tcl_script+=("translate by 0 0 0.5")
  tcl_script+=("animate goto ${step}")
  tcl_script+=("set mol 0")
  tcl_script+=("set rep 0")
  tcl_script+=("mol modselect \${rep} \${mol} none")
  ##############################################################################
  # find number of aggregates for each timestep in the agg file and pick the
  # timestep to visualize
  ##############################################################################
  # the '|| ...' catches a failing grep here; without it, 'set -e' would abort
  # the script on the assignment itself, never reaching the message
  steps=$(grep -n "Step:" "${agg}") || {
    print "No steps found in ${agg}" >&2
    exit 1
  }
  # ${(@f)...} split the grep output into lines
  lines1=("${(@f)steps}")
  # ${(@)...%%:*} trims off the text after the line number (123:Step: -> 123)
  lines1=("${(@)lines1%%:*}")
  # indices in arrays - zsh start from 1
  i1=${i}
  i2=$(( i + 1 ))
  if (( i2 > ${#lines1[@]} )); then
    print "Missing second 'Step:' for i=${i}" >&2
    continue
  fi
  # line numbers of aggreagetes in the step
  # first line; +2 to ignore the 'Step:' and number of aggregates lines
  l1=$(( ${lines1[${i1}]} + 2 ))
  # last line; -1 to stop before the next 'Step:' line
  l2=$(( ${lines1[${i2}]} - 1 ))
  ##############################################################################
  # loop over all aggregate lines in step ${i}, creating vmd representation (or
  # more) for each aggregate
  ##############################################################################
  n=-1 # ColorID
  for ((j = l1; j <= l2; j+=2)); do
    # each aggregate takes two lines: core molecules first, border ones second
    pair=("${(@f)$(sed -n "${j},$(( j + 1 ))p" "${agg}")}")
    # drop the '<count> :' prefix of each line; ${=...} splits the ids into
    # words, getting rid of the leading space in the process
    core=(${=pair[1]#*:})
    border=(${=pair[2]#*:})
    # the aggregate's size counts both kinds of molecule
    size=$(( ${#core[@]} + ${#border[@]} ))
    n=$(( (n + 1) % 33 )) # vmd ColorID: 0 to 32
    # a lone molecule is barely an aggregate, so make it hardly visible
    if [[ "${size}" -gt "1" ]]; then
      mat_core="Opaque"
      mat_border="Transparent"
    else
      mat_core="Glass1"
      mat_border="Glass1"
    fi
    ############################################################################
    # core and border molecules get a representation each, sharing the
    # aggregate's colour, so that the two can be styled differently (change the
    # materials above to draw them alike)
    ############################################################################
    if [[ ${#core[@]} -gt 0 ]]; then
      add_rep "resid ${core}" "${n}" "${mat_core}"
    fi
    if [[ ${#border[@]} -gt 0 ]]; then
      add_rep "resid ${border}" "${n}" "${mat_border}"
    fi
  done
  tcl_script+=("source \$dir/snap.tcl")
  tcl_script+=("set st ${step}")
  tcl_script+=("set rc [snap ${pic} ${trim} ${snap_size_cm} ${dpi}]")
  tcl_script+=("if {\$rc != 0} { puts stderr \"Snapshot failed\" }")
  tcl_script+=("exit")
  printf "%s\n" "${tcl_script[@]}" > "${tcl}"
  ##############################################################################
  # run vmd with the new script and delete the script
  ##############################################################################
  vmd_args=(-size "${xres}" "${yres}" -dispdev text)
  # the startup file is optional; only pass it on if there is one to read
  if [[ -r ${vmdrc} ]]; then
    vmd_args+=(-startup "${vmdrc}")
  fi
  vmd "${vtf}" "${vmd_args[@]}" -e "${tcl}"
  # remove the used vmd tcl file
  rm "${tcl}"
  magick ${pic/tga/jpg} -crop 1600x1600+2160+2125 ${pic/tga/jpg}
done
