Example of getting distributions of angles, both specified in a structure file
and specified by user; the traj.nfo contains system information generated via
the Info utility. Gnuplot scripts were used to generate the attached plots.

Command:
AngleMolecules traj.vcf 1 angles.txt -i traj.data -n angles_extra.txt 1 2 4 1 2 5 --all

traj.vcf ... input coordinate file (containing only the molecules)
traj.data ... input structure file
angles.txt ... output file with distributions and averages of angles specified in
               the structure file
angle_extra.txt ... output file with distributions and average of angles
                    specified in the AngleMolecules command
