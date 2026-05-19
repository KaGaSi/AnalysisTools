# VMD for LINUXAMD64, version 1.9.3 (November 30, 2016)

color Display Background white
axes location Off
display projection orthographic
display depthcue off

display resetview
pbc box
scale by 3

set mol 0
set rep 0
mol modselect   ${rep} ${mol} name A
mol modstyle    ${rep} ${mol} CPK 0.5
mol modcolor    ${rep} ${mol} ColorID 4
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name B
mol modstyle    ${rep} ${mol} CPK 0.5
mol modcolor    ${rep} ${mol} ColorID 6
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} resname CD
mol modstyle    ${rep} ${mol} CPK 0.5 0.8
mol modcolor    ${rep} ${mol} Name
mol modmaterial ${rep} ${mol} Opaque

set dir ~/Scripting/vmd
# cli arguments to pass to snap.tcl:
# 1: snapshot name
# 2: trim(1)/notrim(0)
# 3: xsize in cm (optional)
# 4: dpi (optional)
source $dir/snap.tcl
set rc [snap snap.tga 1 10 300]
if {$rc != 0} {
  puts stderr "Snapshot failed for frame $i"
}
exit
