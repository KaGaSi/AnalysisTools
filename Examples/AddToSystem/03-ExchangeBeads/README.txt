This example demonstrates four options that are not shown in the earlier
examples: -xb, --bonded, --tail, and -s.

The initial system is a 10x10x10 box at DPD density 3 (3000 beads): 100
two-bead AB amphiphiles plus two-component solvent (1300 S1 and 1500 S2 beads)

Exchange mode:
  Without --add, AddToSystem switches beads from the input system rather than
  appending, leaving the total bead count unchanged. By default, it targets the
  most numerous bead type (here S1). The -xb option names the type explicitly:

    AddToSystem 1.vtf E.FIELD 2.vtf        # S1 chosen automatically
    AddToSystem 1.vtf E.FIELD 3.vtf -xb S2 # S2 chosen instead of the default S1

Distance check from bonded beads:
  --bonded replaces -bt for the -ld/-hd distance check: all bonded beads in
  the system are used as the reference set, regardless of their type.

    AddToSystem 1.vtf E.FIELD 4.vtf --add --bonded -ld 1.5 -hd 3

Tail/Head-anchored placement:
  -cx/-cy/-cz constrain the molecule's geometric centre by default; --head
  uses the first bead (D in this case); --tail uses the last one (the last E
  bead in this case):

    AddToSystem 1.vtf D1E9.FIELD 5.vtf --add --tail -cx 0.5 1
    AddToSystem 1.vtf D1E9.FIELD 6.vtf --add --head -cx 0.5 1
