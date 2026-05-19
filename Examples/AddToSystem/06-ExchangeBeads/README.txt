This example demonstrates four options that are not shown in the earlier
examples: -xb, --bonded, --tail, and -s.

The initial system is a 10x10x10 box at DPD density 3 (3000 beads): 100
two-bead AB amphiphiles (A=head, B=tail) plus 2800 W solvent beads.

Exchange mode:
  Without --add, AddToSystem switches beads from the input system rather than
  appending. By default it targets the most numerous bead type (here S1).
  -xb names the type explicitly:

    AddToSystem 1.vtf E.FIELD 2a.vtf           # W chosen automatically
    AddToSystem 1.vtf E.FIELD 2b.vtf -xb W     # same, explicit

  Both replace 100 W beads with 100 E beads; total bead count stays at 3000.

Distance check from bonded beads (--bonded):
  --bonded replaces -bt for the -ld/-hd distance check: all bonded beads in
  the system are used as the reference set, regardless of their type.

    AddToSystem 1.vtf E.FIELD 3.vtf --add --bonded -ld 1.5

  Places 100 E beads at least 1.5 units from any bead in an AB molecule.
  Useful when bead type names are not known in advance or change between runs.

Tail-anchored placement (--tail):
  -cx/-cy/-cz constrain the molecule's geometric centre by default.  --head
  uses the first bead; --tail uses the last bead.

    AddToSystem 1.vtf AB.FIELD 4.vtf --add --tail -cx 0.5 1

  Places AB molecules with their B (tail) bead in x in [5, 10].  Swapping
  --tail for --head would constrain the A (head) bead to the same region
  instead.

Reproducible placement (-s):
  The random seed defaults to the system clock.  -s fixes it:

    AddToSystem 1.vtf E.FIELD 5.vtf -s 42

  Two runs with the same -s and the same input produce identical output.
