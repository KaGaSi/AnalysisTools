################################################################################
# Command to generate pictures in manual (generates snap.tga image):
#   vmd <file>.vtf -e vmd.tcl -size 500 500 -dispdev text
#
# For interactive vmd session, comment out the last two commands and run:
#   vmd <file>.vtf -e vmd.tcl
#
# Note that this vmd script shows only a cross-section of the beads added to the
# simulation box; to inspect the whole system, remove the z-axis constraints for
# 'w' beads.
################################################################################

package require pbctools
# general settings
color Display Background white
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
mol modstyle    ${rep} ${mol} CPK 1.000000 0.000000 12.000000 12.000000
mol modcolor    ${rep} ${mol} ColorID 0
mol modmaterial ${rep} ${mol} Opaque
set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name B
mol modstyle    ${rep} ${mol} CPK 1.000000 0.000000 12.000000 12.000000
mol modcolor    ${rep} ${mol} ColorID 1
mol modmaterial ${rep} ${mol} Opaque
# visualize new beads
set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name w and z > 6 and z < 7
mol modstyle    ${rep} ${mol} CPK 1.000000 0.000000 12.000000 12.000000
mol modcolor    ${rep} ${mol} ColorID 7
mol modmaterial ${rep} ${mol} Opaque
# generate snapshot
render TachyonInternal "snap.tga"
# terminate vmd
exit
