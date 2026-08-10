################################################################################
# Command to start an interactive vmd session:
#   vmd <file>.vtf -e vmd.tcl
#
# Colourscheme:
#   A5B1 molecules ... pink and cyan balls
#   E5D1 molecules ... magenta and green balls
#   P and N ions ... blue and red balls, respectively
#   W beads ... transparent grey balls
################################################################################

package require pbctools

color Display Background white
axes location Off
display projection orthographic
display depthcue off

display resetview
rotate y by 110
rotate x by 20
scale by 1.4
pbc box

set mol 0
set rep 0
mol modselect   ${rep} ${mol} resname A5B1
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} Name
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} resname E5D1
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} Name
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name W
mol modstyle    ${rep} ${mol} CPK 0.3 0.5
mol modcolor    ${rep} ${mol} ColorID 6
mol modmaterial ${rep} ${mol} Opaque
mol modmaterial ${rep} ${mol} Transparent

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name P
mol modstyle    ${rep} ${mol} CPK 0.3 0.5
mol modcolor    ${rep} ${mol} ColorID 0
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name N
mol modstyle    ${rep} ${mol} CPK 0.3 0.5
mol modcolor    ${rep} ${mol} ColorID 1
mol modmaterial ${rep} ${mol} Opaque
