################################################################################
# Command to start an interactive vmd session:
#   vmd <file>.vtf -e vmd.tcl
#
# Colourscheme:
#   A bead ... small yellow balls (the original system)
#   W bead ... big magenta balls (the added beads)
################################################################################

package require pbctools
# general settings
color Display Background white
axes location LowerLeft
axes location Off
display projection orthographic
display depthcue off
display resetview
translate by 0.5 0 0
pbc box
# visualize original beads
set mol 0
set rep 0
mol modselect   ${rep} ${mol} name A
mol modstyle    ${rep} ${mol} CPK 0.5 0.0 12.0 12.0
mol modcolor    ${rep} ${mol} ColorID 4
mol modmaterial ${rep} ${mol} Opaque
# visualize new beads
set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name W
mol modstyle    ${rep} ${mol} CPK 1.0 0.0 12.0 12.0
mol modcolor    ${rep} ${mol} ColorID 13
mol modmaterial ${rep} ${mol} Opaque
