################################################################################
# Command to start an interactive vmd session:
#   vmd <file>.vtf -e vmd.tcl
#
# Colourscheme:
#   S1 and S2     ... small blue and red balls
#   AB molecule   ... bigger green balls and sticks
#   C bead        ... big magenta balls
#   D1E9 molecule ... big magenta and yellow (the D bead; molecule's head) balls
#                     and stics
################################################################################

package require pbctools
# general settings
color Display Background white
axes location LowerLeft
display projection orthographic
display depthcue off
display resetview
translate by 0.5 0 0
pbc box
# visualize original beads
set mol 0
set rep 0
mol modselect   ${rep} ${mol} name S1
mol modstyle    ${rep} ${mol} CPK 0.1 0.0 12.0 12.0
mol modcolor    ${rep} ${mol} ColorID 0
mol modmaterial ${rep} ${mol} Opaque
set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name S2
mol modstyle    ${rep} ${mol} CPK 0.1 0.0 12.0 12.0
mol modcolor    ${rep} ${mol} ColorID 1
mol modmaterial ${rep} ${mol} Opaque
set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} resname AB
mol modstyle    ${rep} ${mol} CPK 0.4 0.4 12.0 12.0
mol modcolor    ${rep} ${mol} ColorID 7
mol modmaterial ${rep} ${mol} Opaque
# visualize new beads
set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name C
mol modstyle    ${rep} ${mol} CPK 0.5 0.0 12.0 12.0
mol modcolor    ${rep} ${mol} ColorID 13
mol modmaterial ${rep} ${mol} Opaque
set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} resname D1E9
mol modstyle    ${rep} ${mol} CPK 0.5 0.4 12.0 12.0
mol modcolor    ${rep} ${mol} ColorID 13
mol modmaterial ${rep} ${mol} Opaque
set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name D
mol modstyle    ${rep} ${mol} CPK 0.5 0.4 12.0 12.0
mol modcolor    ${rep} ${mol} ColorID 4
mol modmaterial ${rep} ${mol} Opaque
