#!/bin/env zsh

################################################################################
# Render VMD snaphsots of aggregates (one colour per aggregate) for each
# timestep from .agg and .vtf files. The script's parsing of the .agg file is
# based on the AnalysisTools v4.0 format
################################################################################

# exit on any command failure, error on unset variables
set -euo pipefail

source ~/Scripting/shell/func.sh

################################################################################
# file names and other variables
################################################################################
name=red
# expand to full path (neceessary if $name starts with ~)
name=${name:A}
tcl=${name}.tcl # tcl script for vmd
agg=${name}.agg # aggregate file (generated via Aggregates)
vtf=join.vtf # structure and coordinate file
pic=${name}.tga # un-numbered snapshot name
# basic error checks
[[ -r ${agg} ]] || { print "Error: File not readable: ${agg}" >&2; exit 1 }
[[ -r ${vtf} ]] || { print "Error: File not readable: ${vtf}" >&2; exit 1 }
[[ -r ~/Scripting/vmd/snap.tcl ]] || { print "Missing snap.tcl" >&2; exit 1 }
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
for ((i = first; i <= last; i++)); do
  # numbered snapshot names
  pic=$(zpad "${i}" 4).tga
  echo "${pic}"
  ##############################################################################
  # start of the tcl script for vmd
  ##############################################################################
  step=$(( ${i} - 1 ))
  tcl_script=()
  tcl_script+=("set dir ~/Scripting/vmd")
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
  # ${(@f)...} split the grep output into lines
  lines1=("${(@f)$(grep -n "Step:" "${agg}")}")
  # check grep actually found something
  (( ${#lines1[@]} > 0 )) || { print "No steps found in ${agg}" >&2; exit 1 }
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
    # line=$(head -n ${j} ${agg} | tail -n 1) # read the appropriate line
    line=$(sed -n "${j}p" "${agg}") # read the appropriate line
    # remove <size> : from the line
    size=${line%% :*}
    line="resid ${(j: :)${line#*:}}"
    n=$(( (n + 1) % 33 )) # vmd ColorID: 0 to 32
    if [[ "${size}" -gt "1" ]]; then
      ##########################################################################
      # print all aggregates
      ##########################################################################
      # # whole aggregate via small balls
      # tcl_script+=("set rep [expr \$rep + 1]")
      # tcl_script+=("mol addrep \${mol}")
      # tcl_script+=("mol modselect \${rep} \${mol} ${line}")
      # tcl_script+=("mol modstyle  \${rep} \${mol} cpk 0.1 0.1")
      # tcl_script+=("mol modcolor  \${rep} \${mol} ColorID ${n}")
      # accentuate hydrophobic core via big balls
      tcl_script+=("set rep [expr \$rep + 1]")
      tcl_script+=("mol addrep \${mol}")
      tcl_script+=("mol modselect \${rep} \${mol} ${line} and name A C D")
      tcl_script+=("mol modstyle  \${rep} \${mol} cpk 0.7 0.5")
      tcl_script+=("mol modcolor  \${rep} \${mol} ColorID ${n}")
    else
      tcl_script+=("set rep [expr \$rep + 1]")
      tcl_script+=("mol addrep \${mol}")
      tcl_script+=("mol modselect \${rep} \${mol} ${line} and name A C D")
      tcl_script+=("mol modstyle  \${rep} \${mol} cpk 0.7 0.5")
      tcl_script+=("mol modcolor  \${rep} \${mol} ColorID ${n}")
      tcl_script+=("mol modmaterial \${rep} 0 Glass1")
    fi
  done
  tcl_script+=("source \$dir/snap.tcl")
  tcl_script+=("set st ${step}")
  tcl_script+=("set rc [snap ${pic} ${trim} ${snap_size_cm} ${dpi}]")
  tcl_script+=("if {\$rc != 0} { puts stderr 'Snapshot failed' }")
  tcl_script+=("exit")
  printf "%s\n" "${tcl_script[@]}" > "${tcl}"
  ##############################################################################
  # run vmd with the new script and delete the script
  ##############################################################################
  vmd \
    "${vtf}" \
    -size "${xres}" "${yres}" \
    -dispdev text \
    -startup ~/.config/vmdrc \
    -e "${tcl}"
  # remove the used vmd tcl file
  rm "${tcl}"
  magick ${pic/tga/jpg} -crop 1600x1600+2160+2125 ${pic/tga/jpg}
done
