################################################################################
# Command to start an interactive vmd session:
#   vmd <file>.vtf -e vmd.tcl
#
# Colourscheme:
#   A4B6 molecules ... pink and cyan balls
#   W solvent beads ... transparent grey balls
################################################################################

package require pbctools

color Display Background white
axes location Off
display projection orthographic
display depthcue off
scale by 1.5

display resetview
rotate x by 20
rotate y by 30
pbc box

set mol 0
set rep 0
mol modselect   ${rep} ${mol} resname A4B6
mol modstyle    ${rep} ${mol} CPK 0.8 0.5
mol modcolor    ${rep} ${mol} Name
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name W
mol modstyle    ${rep} ${mol} CPK 0.3
mol modcolor    ${rep} ${mol} ColorID 6
mol modmaterial ${rep} ${mol} Transparent
