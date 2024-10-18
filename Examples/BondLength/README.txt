Example of getting distributions of bond lengths and distances between the first
and last beads for molecules in the system in the ../AngleMolecules directory.
Gnuplot scripts were used to generate the attached plots.

Command:
BondLength ../AngleMolecules/traj.vcf 0.01 bonds.txt -i ../AngleMolecules/traj.data -d end_to_end.txt

traj.vcf ... input coordinate file (containing only the molecules)
traj.data ... input structure file
bonds.txt ... output file with distributions and averages of bonds specified in
              the structure file
end_to_end.txt ... output file with distributions and average of distances
                   between the first and last bead in each molecule
