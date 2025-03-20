#include "ReadWriteConfig.h"

void WriteConfig(const SYSTEM System, const char *file) { //{{{
  FILE *out = OpenFile(file, "w");
  // TODO: check triclinic box in dl_meso
  // print CONFIG file initial stuff
  fprintf(out, "NAME\n 0 1\n"); // not sure what 0 1 is...
  fprintf(out, "%.3f 0.000 0.000\n", System.Box.Length[0]);
  fprintf(out, "0.000 %.3f 0.000\n", System.Box.Length[1]);
  fprintf(out, "0.000 0.000 %.3f\n", System.Box.Length[2]);

  // bead coordinates
  // unbonded beads must be first (dl_meso requirement)
  for (int i = 0; i < System.Count.BeadCoor; i++) {
    int id = System.BeadCoor[i];
    int btype = System.Bead[id].Type;
    fprintf(out, "%s %d\n", System.BeadType[btype].Name, id + 1);
    fprintf(out, "%lf %lf %lf\n", System.Bead[id].Position[0],
                                  System.Bead[id].Position[1],
                                  System.Bead[id].Position[2]);
  }
  fclose(out);
} //}}}
