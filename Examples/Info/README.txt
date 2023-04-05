Files in.FIELD and in.vtf contain j
.nfo files show outputs of Info command with various options

1) simple case: Info <file>
  * files: vtf.nfo
           FIELD.nfo
  * this command just lists the contents of the given system
  * use vtf.nfo to compare with the following results

2) using --detailed option: Info in.vtf --detailed
  * file: vtf_detailed.nfo
  * this command identifies bead types (and, consequently, molecule types)
    according to not just their names, but also according to their mass, charge,
    and/or radius
    * compare with vtf.nfo

3) using -i option: Info in.vtf -i[i] in.FIELD
  * files: opt_i.nfo
           opt_i-excl.nfo
  * this command adds information from in.FIELD to the information from in.vtf,
    adding charge and mass to beads with unspecified values and enriching the
    molecules with bonds, angles, etc.
    * using just -i changes only molecules that share bead types and bead order
    * using -i! first exchanges beads from in.vtf molecule types with those from
      in.FIELD molecule types (for types with the same number of beads)
      before adding the extra information
    * in both cases, charge is added to the B beads

4) using -c option: Info in.vtf -c in.vtf
  * file: vtf_coor.nfo
  * this commands takes the first timestep from in.vtf, removing from the
    structure information all beads not present in that timestep
