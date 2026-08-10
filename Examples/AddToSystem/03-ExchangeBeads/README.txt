This example demonstrates four options that are not shown in the earlier
examples: -xb, --bonded, --tail, and --head

The initial system is a 20x20x20 with 3000 beads: 100 two-bead AB amphiphiles
plus two-component solvent (1500 S1 and 1300 S2 beads)

Exchange mode:
  Without --add, AddToSystem switches beads from the input system rather than
  appending, leaving the total bead count unchanged. By default, it targets the
  most numerous bead type (here S1). The -xb option names the type explicitly:

  default behaviour       ... S1 bead type chosen automatically
  relevant option: -xb S2 ... S2 bead type chosen instead

Distance check from bonded beads:
  --bonded replaces -bt for the -ld/-hd distance check: all bonded beads in
  the system are used as the reference set, regardless of their type.

  relevant options: --bonded -ld 0.1 -hd 0.2

Tail/Head-anchored placement:
  -cx/-cy/-cz/-hd/-ld constrain the molecule's geometric centre (cog) by
  default; --head uses the first bead instead (D in this case); conversely,
  --tail uses the last one (the last E bead in this case):

  relevant options: -cx 0 0.01              ... cog at the box edge
                    --head -cx 0 0.01       ... head bead at the box edge
                    --tail -cx 0 0.01       ... tail bead the box edge
                    --bonded -hd 0.1        ... cog near any bonded bead
                    --head --bonded -hd 0.1 ... head bead near any bonded bead
                    --tail --bonded -hd 0.1 ... tail bead near any bonded bead

Run all steps:
  ./run.sh

The new systems' composition should be investigated using the Info utility, and
the positions of added species through vmd <file>.vtf -e vmd.tcl command.
