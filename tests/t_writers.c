/*
 * Writer round-trip tests.
 *
 * Every utility ends in a writer, but the write path had no coverage at all.
 * The shape of each test is the same: read a committed fixture, write it back
 * out with the real entry points (WriteStructure / WriteTimestep /
 * WriteAggregates), read the result again, and assert that the two SYSTEM
 * structs agree on everything the target format is able to carry.
 *
 * That makes each test a genuine round-trip: a writer that drops a bond, a
 * reader that mis-maps an index, and any disagreement between the two about
 * the on-disk layout all show up as a mismatch. Formats are lossy to different
 * degrees (xyz keeps only names and positions, FIELD carries no box, vtf
 * carries no angles), so each comparison is told what to check.
 *
 * Under -DSANITIZE=ON the writers are also exercised for memory errors, which
 * is worth having on its own - they build temporary index arrays that no other
 * test touches.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <fcntl.h>
#include <sys/wait.h>

#ifndef FIXTURES_DIR
#define FIXTURES_DIR "."
#endif
#define FIX(name) FIXTURES_DIR "/" name

// scratch directory for everything this executable writes
static char g_dir[] = "/tmp/writer_test_XXXXXX";

// build a path inside the scratch directory //{{{
/*
 * The returned buffer is reused, so a caller must not hold on to two paths at
 * once; every use here passes it straight into a reader or writer.
 */
static const char *tpath(const char *name) {
  static char buf[256];
  snprintf(buf, sizeof buf, "%s/%s", g_dir, name);
  return buf;
} //}}}

// delete everything this executable wrote, then the directory itself //{{{
static void cleanup_scratch(void) {
  DIR *d = opendir(g_dir);
  if (!d) {
    return;
  }
  struct dirent *e;
  while ((e = readdir(d)) != nullptr) {
    if (strcmp(e->d_name, ".") == 0 || strcmp(e->d_name, "..") == 0) {
      continue;
    }
    char path[512];
    snprintf(path, sizeof path, "%s/%s", g_dir, e->d_name);
    unlink(path);
  }
  closedir(d);
  rmdir(g_dir);
} //}}}

// argv the writers stamp into the output file's byline
static char *g_argv[] = {(char *)"t_writers", nullptr};
static const int g_argc = 1;

// read a structure file of the given type //{{{
static SYSTEM read_struct(const char *path, int type) {
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, path, LINE);
  in.stru.type = type;
  in.coor = in.stru;
  return ReadStructure(in, false);
} //}}}
// read a structure file and its first coordinate frame //{{{
static SYSTEM read_struct_frame(const char *path, int type) {
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, path, LINE);
  in.stru.type = type;
  in.coor = in.stru;
  SYSTEM S = ReadStructure(in, false);
  FILE *fr = OpenFile(path, "r");
  int lc = 0;
  CHECK(ReadTimestep(in, fr, &S, &lc));
  fclose(fr);
  return S;
} //}}}
// write System out as 'type' and read it straight back in //{{{
static SYSTEM write_and_reread(const SYSTEM S, const char *name, int type,
                               bool with_coor) {
  FILE_TYPE out = InitFile;
  s_strcpy(out.name, tpath(name), LINE);
  out.type = type;
  WriteStructure(out, S, -1, true, g_argc, g_argv);
  if (with_coor) {
    WriteTimestepAll(out, S, 1, g_argc, g_argv);
  }
  if (with_coor) {
    return read_struct_frame(out.name, type);
  } else {
    return read_struct(out.name, type);
  }
} //}}}

// ---- comparison helpers ---------------------------------------------------

// what a given format is able to carry through a round trip
enum {
  CMP_PROPS = 1 << 0, // bead charge/mass/radius
  CMP_BONDS = 1 << 1, // bonds
  CMP_EXTRA = 1 << 2, // angles, dihedrals, impropers
  CMP_BOX   = 1 << 3, // box size and angles
};

// bead types must survive by name, with their properties if the format keeps //{{{
static void cmp_bead_types(const char *label, const SYSTEM a, const SYSTEM b,
                           unsigned flags) {
  if (a.Count.BeadType != b.Count.BeadType) {
    g_failures++;
    g_checks++;
    fprintf(stderr, "  FAIL %s: BeadType count %d -> %d\n", label,
            a.Count.BeadType, b.Count.BeadType);
    return;
  }
  for (int i = 0; i < a.Count.BeadType; i++) {
    BEADTYPE *bt = &a.BeadType[i];
    // types may be reordered by the round trip, so match on the name
    int j = FindBeadType(bt->Name, b);
    if (j == -1) {
      g_failures++;
      g_checks++;
      fprintf(stderr, "  FAIL %s: bead type '%s' lost\n", label, bt->Name);
      continue;
    }
    BEADTYPE *bt2 = &b.BeadType[j];
    CHECK(bt->Number == bt2->Number);
    if (flags & CMP_PROPS) {
      CHECK_CLOSE(bt->Charge, bt2->Charge, 1e-5);
      CHECK_CLOSE(bt->Mass, bt2->Mass, 1e-5);
    }
  }
} //}}}
// molecule types must survive by name, with their topology //{{{
static void cmp_mol_types(const char *label, const SYSTEM a, const SYSTEM b,
                          unsigned flags) {
  if (a.Count.MoleculeType != b.Count.MoleculeType) {
    g_failures++;
    g_checks++;
    fprintf(stderr, "  FAIL %s: MoleculeType count %d -> %d\n", label,
            a.Count.MoleculeType, b.Count.MoleculeType);
    return;
  }
  for (int i = 0; i < a.Count.MoleculeType; i++) {
    MOLECULETYPE *mt = &a.MoleculeType[i];
    int j = FindMoleculeName(mt->Name, b);
    if (j == -1) {
      g_failures++;
      g_checks++;
      fprintf(stderr, "  FAIL %s: molecule type '%s' lost\n", label, mt->Name);
      continue;
    }
    MOLECULETYPE *mt2 = &b.MoleculeType[j];
    CHECK(mt->Number == mt2->Number);
    CHECK(mt->nBeads == mt2->nBeads);
    // the bead sequence inside the molecule must be preserved
    if (mt->nBeads == mt2->nBeads) {
      for (int k = 0; k < mt->nBeads; k++) {
        const char *n1 = a.BeadType[mt->Bead[k]].Name;
        const char *n2 = b.BeadType[mt2->Bead[k]].Name;
        CHECK(strcmp(n1, n2) == 0);
      }
    }
    if (flags & CMP_BONDS) {
      CHECK(mt->nBonds == mt2->nBonds);
      if (mt->nBonds == mt2->nBonds) {
        // SortAll() puts both tables in the same canonical order
        for (int k = 0; k < mt->nBonds; k++) {
          CHECK(mt->Bond[k][0] == mt2->Bond[k][0]);
          CHECK(mt->Bond[k][1] == mt2->Bond[k][1]);
        }
      }
    }
    if (flags & CMP_EXTRA) {
      CHECK(mt->nAngles == mt2->nAngles);
      CHECK(mt->nDihedrals == mt2->nDihedrals);
      CHECK(mt->nImpropers == mt2->nImpropers);
      if (mt->nAngles == mt2->nAngles) {
        for (int k = 0; k < mt->nAngles; k++) {
          for (int dd = 0; dd < 3; dd++) {
            CHECK(mt->Angle[k][dd] == mt2->Angle[k][dd]);
          }
        }
      }
    }
  }
} //}}}
// the whole structure: counts, box, bead types, molecule types //{{{
static void cmp_structure(const char *label, const SYSTEM a, const SYSTEM b,
                          unsigned flags) {
  CHECK(a.Count.Bead == b.Count.Bead);
  CHECK(a.Count.Molecule == b.Count.Molecule);
  CHECK(a.Count.Unbonded == b.Count.Unbonded);
  CHECK(a.Count.Bonded == b.Count.Bonded);
  if (flags & CMP_BONDS) {
    CHECK(a.Count.Bond == b.Count.Bond);
  }
  if (flags & CMP_EXTRA) {
    CHECK(a.Count.Angle == b.Count.Angle);
    CHECK(a.Count.Dihedral == b.Count.Dihedral);
    CHECK(a.Count.Improper == b.Count.Improper);
  }
  if (flags & CMP_BOX) {
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(a.Box.Length.v[dd], b.Box.Length.v[dd], 1e-4);
    }
    CHECK_CLOSE(a.Box.alpha, b.Box.alpha, 1e-3);
    CHECK_CLOSE(a.Box.beta, b.Box.beta, 1e-3);
    CHECK_CLOSE(a.Box.gamma, b.Box.gamma, 1e-3);
  }
  cmp_bead_types(label, a, b, flags);
  cmp_mol_types(label, a, b, flags);
} //}}}
// coordinates, bead for bead //{{{
/*
 * Writers print a fixed number of decimals, so the tolerance is the write
 * precision rather than a floating-point epsilon.
 */
static void cmp_coordinates(const char *label, const SYSTEM a, const SYSTEM b) {
  if (a.Count.BeadCoor != b.Count.BeadCoor) {
    g_failures++;
    g_checks++;
    fprintf(stderr, "  FAIL %s: BeadCoor %d -> %d\n", label,
            a.Count.BeadCoor, b.Count.BeadCoor);
    return;
  }
  for (int i = 0; i < a.Count.Bead; i++) {
    for (int dd = 0; dd < 3; dd++) {
      CHECK_CLOSE(a.Bead[i].Position.v[dd], b.Bead[i].Position.v[dd], 1e-3);
    }
  }
} //}}}

// ---- structure round trips ------------------------------------------------

// vtf -> vtf: bonds and bead properties survive, angles are not a vtf concept //{{{
static void test_vtf_structure_roundtrip(void) {
  SYSTEM a = read_struct(FIX("struct.vtf"), VTF_FILE);
  SYSTEM b = write_and_reread(a, "out.vtf", VTF_FILE, false);
  cmp_structure("vtf->vtf", a, b, CMP_PROPS | CMP_BONDS | CMP_BOX);
  FreeSystem(&a);
  FreeSystem(&b);
} //}}}
// data -> data: the one format that carries the full topology //{{{
static void test_data_structure_roundtrip(void) {
  SYSTEM a = read_struct(FIX("system.data"), LDATA_FILE);
  SYSTEM b = write_and_reread(a, "out.data", LDATA_FILE, false);
  cmp_structure("data->data", a, b,
                CMP_PROPS | CMP_BONDS | CMP_EXTRA | CMP_BOX);
  FreeSystem(&a);
  FreeSystem(&b);
} //}}}
// FIELD -> FIELD: full topology, but no box //{{{
static void test_field_structure_roundtrip(void) {
  SYSTEM a = read_struct(FIX("FIELD"), FIELD_FILE);
  SYSTEM b = write_and_reread(a, "FIELD_out", FIELD_FILE, false);
  cmp_structure("field->field", a, b, CMP_PROPS | CMP_BONDS | CMP_EXTRA);
  FreeSystem(&a);
  FreeSystem(&b);
} //}}}

// ---- cross-format conversions ---------------------------------------------

// data -> vtf -> data: topology must survive a trip through the leaner format //{{{
/*
 * vtf carries no angles/dihedrals/impropers, so only the bond topology can
 * come back; that is exactly the part this checks. Catches a writer that
 * renumbers molecule beads on the way out.
 */
static void test_data_via_vtf(void) {
  SYSTEM a = read_struct(FIX("system.data"), LDATA_FILE);
  SYSTEM b = write_and_reread(a, "cross.vtf", VTF_FILE, false);
  cmp_structure("data->vtf", a, b, CMP_PROPS | CMP_BONDS | CMP_BOX);
  FreeSystem(&a);
  FreeSystem(&b);
} //}}}
// vtf -> data -> read: the bond table must land in the data file intact //{{{
/*
 * A frame is loaded first because WriteLmpData writes its Atoms section from
 * BeadCoor, so it needs coordinates; without them it refuses outright, which
 * test_data_write_needs_coordinates below covers.
 */
static void test_vtf_via_data(void) {
  SYSTEM a = read_struct_frame(FIX("struct.vtf"), VTF_FILE);
  SYSTEM b = write_and_reread(a, "cross.data", LDATA_FILE, false);
  cmp_structure("vtf->data", a, b, CMP_PROPS | CMP_BONDS | CMP_BOX);
  FreeSystem(&a);
  FreeSystem(&b);
} //}}}

// writing a data file with no coordinates must be refused, not botched //{{{
/*
 * WriteLmpData builds its Atoms section from BeadCoor, so a structure-only
 * system used to yield "Atoms # full" followed by nothing - a file no reader
 * accepts, written without complaint. It now errors out the way the LTRJ_FILE
 * branch of WriteStructure always has. The writer exits, so run it in a child
 * and classify how that child died.
 */
static void test_data_write_needs_coordinates(void) {
  const char *path = tpath("nocoor.data");
  fflush(nullptr);
  pid_t pid = fork();
  if (pid == 0) { // child
    TestChildNoLeakReports(); // the writer exits without unwinding
    int devnull = open("/dev/null", O_WRONLY);
    if (devnull >= 0) {
      dup2(devnull, STDERR_FILENO);
      dup2(devnull, STDOUT_FILENO);
    }
    SYSTEM S = read_struct(FIX("struct.vtf"), VTF_FILE); // no frame loaded
    FILE_TYPE out = InitFile;
    s_strcpy(out.name, path, LINE);
    out.type = LDATA_FILE;
    WriteStructure(out, S, -1, true, g_argc, g_argv);
    _exit(0); // reached only if the writer accepted a system with no coords
  }
  int status = 0;
  waitpid(pid, &status, 0);
  CHECK(WIFEXITED(status));            // a clean error exit, not a crash
  if (WIFEXITED(status)) {
    CHECK(WEXITSTATUS(status) != 0);   // and it must actually report failure
  }
  // nothing may be left behind for a later read to trip over
  CHECK(access(path, F_OK) != 0);
} //}}}

// ---- coordinate round trips -----------------------------------------------

// vtf: structure and coordinates in a single file //{{{
static void test_vtf_coor_roundtrip(void) {
  SYSTEM a = read_struct_frame(FIX("struct.vtf"), VTF_FILE);
  SYSTEM b = write_and_reread(a, "coor.vtf", VTF_FILE, true);
  cmp_structure("vtf coor", a, b, CMP_PROPS | CMP_BONDS | CMP_BOX);
  cmp_coordinates("vtf coor", a, b);
  FreeSystem(&a);
  FreeSystem(&b);
} //}}}
// xyz: names and positions only, and the bead order must be preserved //{{{
static void test_xyz_coor_roundtrip(void) {
  SYSTEM a = read_struct_frame(FIX("traj.xyz"), XYZ_FILE);
  FILE_TYPE out = InitFile;
  s_strcpy(out.name, tpath("coor.xyz"), LINE);
  out.type = XYZ_FILE;
  // xyz has no structure block; the writer appends, so start from a new file
  FILE *fresh = OpenFile(out.name, "w");
  fclose(fresh);
  WriteTimestepAll(out, a, 1, g_argc, g_argv);
  SYSTEM b = read_struct_frame(out.name, XYZ_FILE);
  CHECK(a.Count.Bead == b.Count.Bead);
  cmp_bead_types("xyz", a, b, 0);
  cmp_coordinates("xyz", a, b);
  // xyz is an ordered format: bead i must still be bead i, of the same type
  if (a.Count.Bead == b.Count.Bead) {
    for (int i = 0; i < a.Count.Bead; i++) {
      const char *n1 = a.BeadType[a.Bead[i].Type].Name;
      const char *n2 = b.BeadType[b.Bead[i].Type].Name;
      CHECK(strcmp(n1, n2) == 0);
    }
  }
  FreeSystem(&a);
  FreeSystem(&b);
} //}}}
// lammpstrj: positions plus the box //{{{
static void test_ltrj_coor_roundtrip(void) {
  SYSTEM a = read_struct_frame(FIX("traj.lammpstrj"), LTRJ_FILE);
  FILE_TYPE out = InitFile;
  s_strcpy(out.name, tpath("coor.lammpstrj"), LINE);
  out.type = LTRJ_FILE;
  FILE *fresh = OpenFile(out.name, "w");
  fclose(fresh);
  WriteTimestepAll(out, a, 1, g_argc, g_argv);
  SYSTEM b = read_struct_frame(out.name, LTRJ_FILE);
  CHECK(a.Count.Bead == b.Count.Bead);
  for (int dd = 0; dd < 3; dd++) {
    CHECK_CLOSE(a.Box.Length.v[dd], b.Box.Length.v[dd], 1e-4);
  }
  cmp_coordinates("ltrj", a, b);
  FreeSystem(&a);
  FreeSystem(&b);
} //}}}
// two frames written in sequence must both read back, in order //{{{
/*
 * WriteTimestep appends, so a multi-frame file exercises a path a single
 * write does not: frame 2 must be found where frame 1 left off.
 */
static void test_multi_frame_vtf(void) {
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, FIX("struct.vtf"), LINE);
  in.stru.type = VTF_FILE;
  in.coor = in.stru;
  SYSTEM a = ReadStructure(in, false);
  FILE *fr = OpenFile(in.coor.name, "r");
  int lc = 0;

  FILE_TYPE out = InitFile;
  s_strcpy(out.name, tpath("multi.vtf"), LINE);
  out.type = VTF_FILE;
  WriteStructure(out, a, -1, true, g_argc, g_argv);
  // remember frame 1's first position so we can tell the frames apart
  CHECK(ReadTimestep(in, fr, &a, &lc));
  vec3d first_frame = a.Bead[0].Position;
  WriteTimestepAll(out, a, 1, g_argc, g_argv);
  CHECK(ReadTimestep(in, fr, &a, &lc));
  vec3d second_frame = a.Bead[0].Position;
  WriteTimestepAll(out, a, 2, g_argc, g_argv);
  fclose(fr);

  // the two frames must actually differ, or this test proves nothing
  CHECK(fabs(first_frame.x - second_frame.x) > 1e-6);

  SYS_FILES back = InitSysFiles;
  s_strcpy(back.stru.name, out.name, LINE);
  back.stru.type = VTF_FILE;
  back.coor = back.stru;
  SYSTEM b = ReadStructure(back, false);
  FILE *fr2 = OpenFile(out.name, "r");
  int lc2 = 0;
  CHECK(ReadTimestep(back, fr2, &b, &lc2));
  CHECK_CLOSE(b.Bead[0].Position.x, first_frame.x, 1e-3);
  CHECK(ReadTimestep(back, fr2, &b, &lc2));
  CHECK_CLOSE(b.Bead[0].Position.x, second_frame.x, 1e-3);
  CHECK(!ReadTimestep(back, fr2, &b, &lc2)); // exactly two frames
  fclose(fr2);
  FreeSystem(&a);
  FreeSystem(&b);
} //}}}

// ---- aggregate file round trip --------------------------------------------

// build a two-aggregate partition of the fixture's two molecules //{{{
static void make_aggregates(SYSTEM *S, AGGREGATE **Agg) {
  InitAggregate(*S, Agg);
  S->Count.Aggregate = 2;
  for (int i = 0; i < 2; i++) {
    (*Agg)[i].nMolecules = 1;
    (*Agg)[i].nCore = 1;
    (*Agg)[i].nBorder = 0;
    (*Agg)[i].Core[0] = i;
    S->Molecule[i].Aggregate = i;
  }
} //}}}
// WriteAggregates -> ReadAggregates must reproduce the partition //{{{
/*
 * The agg file stores molecules by resid (Molecule[].Index), not by the
 * compact molecule index, so the reader has to map back. This is the
 * regression guard for that mapping.
 */
static void agg_roundtrip(const char *label, bool renumber) {
  SYSTEM S = read_struct(FIX("struct.vtf"), VTF_FILE);
  CHECK(S.Count.Molecule == 2);
  if (renumber) {
    // non-contiguous resids: the compact index and the resid diverge
    S.Molecule[0].Index = 7;
    S.Molecule[1].Index = 41;
    if (S.Count.HighestResid < 41) {
      S.Count.HighestResid = 41;
    }
  }
  AGGREGATE *Agg = nullptr;
  make_aggregates(&S, &Agg);

  const char *path = tpath("agg.txt");
  // the reader expects the two-line header the utilities write
  PrintByline(path, g_argc, g_argv);
  bool use[2] = {true, true};
  WriteAggregates(1, path, S, Agg, use);

  // read it back into a fresh aggregate array
  AGGREGATE *back = nullptr;
  InitAggregate(S, &back);
  FILE *fr = OpenFile(path, "r");
  int lc = 0;
  for (int i = 0; i < 2; i++) { // skip the byline header
    SkipLine(fr);
    lc++;
  }
  int ret = ReadAggregates(fr, path, &S, back, &lc);
  fclose(fr);
  CHECK(ret == 1);
  CHECK(S.Count.Aggregate == 2);
  for (int i = 0; i < 2; i++) {
    CHECK(back[i].nMolecules == 1);
    CHECK(back[i].nCore == 1);
    CHECK(back[i].nBorder == 0);
    // the molecule must come back as the compact index, not the raw resid
    if (back[i].nCore == 1) {
      CHECK(back[i].Core[0] == i);
    }
  }
  if (ret != 1) {
    fprintf(stderr, "  (note: %s ReadAggregates returned %d)\n", label, ret);
  }
  FreeAggregate(S.Count, Agg);
  FreeAggregate(S.Count, back);
  FreeSystem(&S);
  unlink(path);
} //}}}
static void test_agg_roundtrip_contiguous(void) {
  agg_roundtrip("contiguous", false);
}
static void test_agg_roundtrip_sparse_resids(void) {
  agg_roundtrip("sparse", true);
}

// ---- documented writer side effect ----------------------------------------

// VtfWriteStruct renumbers the caller's resids //{{{
/*
 * SimplifyResid() shifts Molecule[].Index so the lowest becomes 1, and it does
 * so on the caller's system - SYSTEM is passed by value but Molecule is a
 * pointer into the same allocation. Any code that writes a vtf and then uses
 * resids (an aggregate file, say) sees the shifted numbering, so pin it here
 * rather than leave it as a surprise.
 */
static void test_vtf_write_shifts_resids(void) {
  SYSTEM S = read_struct(FIX("struct.vtf"), VTF_FILE);
  CHECK(S.Count.Molecule == 2);
  S.Molecule[0].Index = 7;
  S.Molecule[1].Index = 41;
  FILE_TYPE out = InitFile;
  s_strcpy(out.name, tpath("resid.vtf"), LINE);
  out.type = VTF_FILE;
  WriteStructure(out, S, -1, true, g_argc, g_argv);
  // lowest resid becomes 1, the gap to the next one is preserved
  CHECK(S.Molecule[0].Index == 1);
  CHECK(S.Molecule[1].Index == 35);
  FreeSystem(&S);
} //}}}

int main(void) {
  if (!mkdtemp(g_dir)) {
    fprintf(stderr, "cannot create scratch directory\n");
    return 1;
  }
  RUN(test_vtf_structure_roundtrip);
  RUN(test_data_structure_roundtrip);
  RUN(test_field_structure_roundtrip);
  RUN(test_data_via_vtf);
  RUN(test_vtf_via_data);
  RUN(test_data_write_needs_coordinates);
  RUN(test_vtf_coor_roundtrip);
  RUN(test_xyz_coor_roundtrip);
  RUN(test_ltrj_coor_roundtrip);
  RUN(test_multi_frame_vtf);
  RUN(test_agg_roundtrip_contiguous);
  RUN(test_agg_roundtrip_sparse_resids);
  RUN(test_vtf_write_shifts_resids);
  cleanup_scratch();
  return test_main_end();
}
