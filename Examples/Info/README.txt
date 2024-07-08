Files in.FIELD and in.vtf contain a small sandbox system, and .nfo files show
outputs of Info command with various options.

1) simple case: Info <file>
  * files: vtf.nfo
           FIELD.nfo
  * this command lists the contents of the given system in the most simplest way
  * use vtf.nfo to compare with the following output files

2) using --detailed option: Info in.vtf --detailed
  * file: vtf_detailed.nfo
  * this command identifies bead types (and, consequently, molecule types)
    according to not just bead names, but also according to beads' mass, charge,
    and radius

3) using -i option: Info in.vtf -i in.FIELD
  * files: opt_i.nfo
  * this command supplements information from in.vtf by information from
    in.FIELD, adding charge and mass to beads with unspecified values and
    enriching the molecules with bonds, angles, etc. (for the molecules share a
    name and a bead order)

4) using -c option: Info in.vtf -c in.vtf
  * file: vtf_coor.nfo
  * this commands takes the first timestep from in.vtf, removing from the
    structure information all beads not present in that timestep
