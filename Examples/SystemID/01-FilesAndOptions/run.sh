#!/usr/bin/env bash
# make the script runnable from anywhere
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

################################################################################
# How a system gets identified: what each file type carries, and which of the
# -i/-ft/-lib/-sys options supplies the rest.
#
# Selected is used throughout as a neutral read-it-and-write-it-back utility, so
# that nothing but the identification is going on, and Info (with no options of
# its own) as the microscope showing what came out. Note that Info's and
# AddToSystem's own -lib do more than the common one shown here.
#
# Everything is generated from one source system, so this directory starts out
# with nothing but the script and the README.
################################################################################
lib="${ROOT}/../../library/"

################################################################################
# 0) the source system, and copies of it that have lost something
################################################################################
# 4 CTAC (9 beads each) + the 4 Cl- counterions they declare + 160 water beads
# -s fixes the placement, so that this example gives the same numbers each run
AddToSystem - system.vtf -lib "${lib}" -mol CTAC 4 -ntot 200 water \
  -b 10 10 10 -sys system.nfo -s 1

Selected system.vtf system.xyz   # names and coordinates, nothing else
Selected system.vtf system.data  # a lammps data file as these tools write it

# a) what another program's data file looks like: same data, no name comments
awk '$1 !~ /^[0-9.-]/ && NF > 0 {section = $1}
     (section == "Masses" || section == "Atoms") && $1 ~ /^[0-9]+$/ {
       sub(/[ \t]*#.*$/, "")
     } {print}' system.data > foreign.data
# b) the same file with the Masses section gone as well
awk '/^Masses/ {skip = 1; next} skip && /^[A-Z]/ {skip = 0} !skip' \
  foreign.data > nomass.data
# c) a vtf that has lost its per-bead-type masses
sed -E 's/ mass +[0-9.]+//' system.vtf > nomass.vtf
# d) a file whose extension says nothing about its format
cp system.xyz coords.out

################################################################################
# 1) what each file says on its own - look at what is missing
################################################################################
Info system.vtf   > 01-vtf.nfo      # names, molecules, mass/charge/radius
Info system.xyz   > 01-xyz.nfo      # names only: no molecules, no bead data
Info system.data  > 01-data.nfo     # everything, via the '#' name comments
Info foreign.data > 01-foreign.nfo  # bead types b0..b4, molecule type 'm'

################################################################################
# 2) the options that fill in the rest
################################################################################
# the coordinate file has no topology - take it from another file
Selected system.xyz 11-xyz-with-i.vtf -i system.vtf
Info 11-xyz-with-i.vtf > 11-xyz-with-i.nfo

# the extension says nothing - name the format
Selected coords.out 12-forced-type.vtf -ft xyz
Info 12-forced-type.vtf > 12-forced-type.nfo

# names are already right, the data is missing - the library fills it in and
# renames nothing (an xyz has no charges to match on anyway)
Selected system.xyz 13-xyz-lib.vtf -lib "${lib}"
Info 13-xyz-lib.vtf > 13-xyz-lib.nfo

# a foreign data file with -lib alone: only the free beads can be matched (H2O
# by charge, Cl- by charge and count), so it warns and the molecule's own bead
# types are left as b1..b3
Selected foreign.data 14-lib-only.vtf -lib "${lib}"
Info 14-lib-only.vtf > 14-lib-only.nfo

# -sys alone names the molecule types but touches no bead type: the two options
# do different jobs
Selected foreign.data 15-sys-only.vtf -sys system.nfo
Info 15-sys-only.vtf > 15-sys-only.nfo

# together: molecule types get their names, so the library can match them and
# every bead type is identified
Selected foreign.data 16-lib-sys.vtf -lib "${lib}" -sys system.nfo
Info 16-lib-sys.vtf > 16-lib-sys.nfo

# and what it refuses to do: with no Masses section the file cannot tell its
# neutral beads apart, so they arrive as one bead type that stands for three
# library beads at once. That type is left alone, with a warning, instead of
# taking whichever name came last
Selected nomass.data 17-ambiguous.vtf -lib "${lib}" -sys system.nfo
Info 17-ambiguous.vtf > 17-ambiguous.nfo

################################################################################
# 3) why it is worth the trouble: a mass-weighted result needs the masses
################################################################################
Aggregates system.vtf system.agg -d 1.0 -c 1
# nomass.vtf keeps its names, so -lib matches them and fills the masses back in
GyrationAggregates nomass.vtf system.agg rg-nolib.txt --joined
GyrationAggregates nomass.vtf system.agg rg-lib.txt --joined -lib "${lib}"
# rg-nolib.txt has -nan in every mass-weighted column, rg-lib.txt has numbers
