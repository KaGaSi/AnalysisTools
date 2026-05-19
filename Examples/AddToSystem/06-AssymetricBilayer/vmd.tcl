# VMD for LINUXAMD64, version 1.9.3 (November 30, 2016)

color Display Background white
axes location Off
display projection orthographic
display depthcue off

display resetview
rotate y by 110
rotate x by 20
translate by 0 0 -5
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
mol modcolor    ${rep} ${mol} ColorID 8
mol modmaterial ${rep} ${mol} Opaque

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
