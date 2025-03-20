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
mol modselect   ${rep} ${mol} resname FA_C16 FA_C18 BTAC
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} ColorID 5
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name pos
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} ColorID 1
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name co
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} ColorID 10
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name CI
mol modstyle    ${rep} ${mol} CPK 0.5 0.5
mol modcolor    ${rep} ${mol} ColorID 0
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name w
mol modstyle    ${rep} ${mol} CPK 0.1 0.5
mol modcolor    ${rep} ${mol} ColorID 8
mol modmaterial ${rep} ${mol} Opaque
