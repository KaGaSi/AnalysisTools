Examples from the figure in the manual section for AddToSystem:

Commands to generate the new systems (System1.lammpstrj is the same as in Examples/JoinSystems):
Figure (a): AddToSystem System1.lammpstrj FIELD fig_a.vtf --add -ld 3 -bt A
Figure (b): AddToSystem System1.lammpstrj FIELD fig_b.vtf --add -hd 4 -bt A
Figure (c): AddToSystem System1.lammpstrj FIELD fig_c.vtf --add -ld 3 -hd 4 -bt A
Figure (d): AddToSystem System1.lammpstrj FIELD fig_d.vtf --add -ld 3 -hd 4 -bt A -cx 0.5 1
Figure (e): AddToSystem System1.lammpstrj FIELD fig_e.vtf --add -ld 3 -hd 4 -bt A -cx 0.5 1 -b 30 20 25
Figure (f): AddToSystem System1.lammpstrj FIELD fig_f.vtf --add -ld 3 -hd 4 -bt A -cx 0.5 1 -b 30 20 25 -off -0.2 0.2 0

The vmd.tcl file is a vmd script that was used to generate the figures inside
the manual.
