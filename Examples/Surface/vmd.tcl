# VMD for LINUXAMD64, version 1.9.3 (November 30, 2016)

color Display Background white
axes location Off
display depthcue off

pbc box

set mol 0
set rep 0
mol modselect   ${rep} ${mol} resname FA_C16 CTAC
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} Name
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name c
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} ColorID 12
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name co
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} ColorID 9
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name pos
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} ColorID 27
mol modmaterial ${rep} ${mol} Opaque

set rep [expr $rep + 1]
mol addrep ${mol}
mol modselect   ${rep} ${mol} name CI
mol modstyle    ${rep} ${mol} CPK 1.0 0.5
mol modcolor    ${rep} ${mol} ColorID 8
mol modmaterial ${rep} ${mol} Opaque

animate goto start

display resetview
rotate x by 90
rotate y by 35
