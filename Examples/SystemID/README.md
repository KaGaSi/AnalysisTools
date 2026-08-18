# Identifying a system

Unlike the other example folders, this one is not about a single utility. It
covers the options every utility shares for working out what system it has been
given: `-i` and `-ft` for which file is read and how, and `-lib` and `-sys` for
what the bead and molecule types in it are called and what they weigh.

- **01-FilesAndOptions**
  What each file type carries and what it leaves out, and which option supplies
  the rest - shown on one small system written out as vtf, xyz, and lammps data
  (both as these tools write it and as another program would). Ends with a
  mass-weighted calculation that returns `-nan` until the masses are supplied.

The example uses `Selected` as a neutral read-and-write-back utility and `Info`
as the inspector. Note that `Info` and `AddToSystem` have their own, different
`-lib`: see their sections in the manual.
