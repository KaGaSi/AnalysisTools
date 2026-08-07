/*
 * Reader tests (VTF, LAMMPS data, LAMMPS trajectory, XYZ).
 *
 * Two halves:
 *   1. Valid parse - the committed fixtures in tests/fixtures/ (all the same
 *      13-bead system) are read with the real entry points (ReadStructure /
 *      ReadTimestep) and their counts, box, and a known coordinate asserted.
 *   2. Malformed input - a fixture is mutated in memory, written to a temp
 *      file, and fed to the reader in a forked child. The readers call exit()
 *      on bad input, so we classify the child's termination:
 *        - clean exit(non-zero)  -> graceful rejection (good)
 *        - killed by a signal    -> crash / ASan abort (bad)
 *      Every malformed case must not crash; the VTF bond-index-overflow case
 *      (regression for commit 88031eb) must additionally be rejected.
 *
 * Under -DSANITIZE=ON each mutation doubles as a fuzz case: an overflow that
 * slips past validation becomes a SIGABRT and fails the test.
 */
#include "../src/AnalysisTools.h"
#include "test_util.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>

#ifndef FIXTURES_DIR
#define FIXTURES_DIR "."
#endif
#define FIX(name) FIXTURES_DIR "/" name

// data and FIELD are single-config structure files; the others hold frames
static bool has_frames(int ftype) {
  return ftype != LDATA_FILE && ftype != FIELD_FILE;
}

// read a whole file into a malloc'd NUL-terminated buffer //{{{
static char *slurp(const char *path) {
  FILE *f = fopen(path, "rb");
  if (!f) {
    return nullptr;
  }
  fseek(f, 0, SEEK_END);
  long n = ftell(f);
  fseek(f, 0, SEEK_SET);
  char *buf = malloc(n + 1);
  if (buf) {
    size_t got = fread(buf, 1, n, f);
    buf[got] = '\0';
  }
  fclose(f);
  return buf;
} //}}}

// return a new string with the first occurrence of needle replaced //{{{
static char *replace_first(const char *hay, const char *needle,
                           const char *repl) {
  const char *at = strstr(hay, needle);
  CHECK(at != nullptr); // the mutation target must exist in the fixture
  if (!at) {
    return nullptr;
  }
  size_t nlen = strlen(needle), rlen = strlen(repl), pre = at - hay;
  char *out = malloc(strlen(hay) - nlen + rlen + 1);
  memcpy(out, hay, pre);
  memcpy(out + pre, repl, rlen);
  strcpy(out + pre + rlen, at + nlen);
  return out;
}
// truncate the fixture at the first occurrence of marker
static char *truncate_at(const char *path, const char *marker) {
  char *buf = slurp(path);
  CHECK(buf != nullptr);
  char *at = nullptr;
  if (buf) {
    at = strstr(buf, marker);
  }
  CHECK(at != nullptr);
  if (at) {
    *at = '\0';
  }
  return buf;
} //}}}

// write content to a fresh temp file; returns the path in a static buffer //{{{
static const char *write_temp(const char *content) {
  static char path[64];
  strcpy(path, "/tmp/reader_test_XXXXXX");
  int fd = mkstemp(path);
  CHECK(fd >= 0);
  size_t len = strlen(content);
  CHECK(write(fd, content, len) == (ssize_t)len);
  close(fd);
  return path;
} //}}}

// run the reader on a file in a child process, classify how it terminates //{{{
struct outcome { bool exited; int code; bool signaled; int sig; };
static struct outcome read_in_child(const char *path, int ftype) {
  fflush(nullptr);
  pid_t pid = fork();
  if (pid == 0) { // child
    int devnull = open("/dev/null", O_WRONLY);
    if (devnull >= 0) {
      dup2(devnull, STDERR_FILENO); // readers are noisy on the error path
      dup2(devnull, STDOUT_FILENO);
    }
    SYS_FILES in = InitSysFiles;
    s_strcpy(in.stru.name, path, LINE);
    in.stru.type = ftype;
    in.coor = in.stru; // self-contained: structure and coords in one file
    SYSTEM s = ReadStructure(in, false);
    if (has_frames(ftype)) {
      FILE *fr = fopen(path, "r");
      if (fr) {
        int lc = 0;
        while (ReadTimestep(in, fr, &s, &lc)) {
          ;
        }
        fclose(fr);
      }
    }
    _exit(0); // reached only if the reader accepted the input
  }
  int status = 0;
  waitpid(pid, &status, 0);
  struct outcome o = {0};
  if (WIFEXITED(status)) {
    o.exited = true;
    o.code = WEXITSTATUS(status);
  }
  if (WIFSIGNALED(status)) {
    o.signaled = true;
    o.sig = WTERMSIG(status);
  }
  return o;
}
// a mutated fixture must not crash the reader (may accept or reject)
static void expect_no_crash(const char *label, int ftype, char *content) {
  if (!content) {
    return;
  }
  const char *path = write_temp(content);
  struct outcome o = read_in_child(path, ftype);
  if (o.signaled) {
    g_failures++;
    g_checks++;
    fprintf(stderr, "  FAIL %s: reader crashed (signal %d) on '%s'\n",
            __FILE__, o.sig, label);
  } else {
    CHECK(o.exited);
  }
  unlink(path);
  free(content);
}
// a mutated fixture must be rejected gracefully (clean non-zero exit)
static void expect_rejected(const char *label, int ftype, char *content) {
  if (!content) {
    return;
  }
  const char *path = write_temp(content);
  struct outcome o = read_in_child(path, ftype);
  CHECK(!o.signaled);             // must not crash
  CHECK(o.exited && o.code != 0); // must report an error and exit
  if (o.signaled) {
    fprintf(stderr, "  (note: '%s' crashed with signal %d)\n", label, o.sig);
  }
  unlink(path);
  free(content);
} //}}}

// ---- valid-parse assertions ----------------------------------------------

// VTF: structure + two indexed frames //{{{
static void test_vtf_valid(void) {
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, FIX("struct.vtf"), LINE);
  in.stru.type = VTF_FILE;
  in.coor = in.stru;

  SYSTEM S = ReadStructure(in, false);
  CHECK(S.Count.Bead == 13);
  CHECK(S.Count.BeadType == 3);     // A, B, C(default)
  CHECK(S.Count.Molecule == 2);
  CHECK(S.Count.MoleculeType == 1);
  CHECK(S.Count.Bond == 8);
  CHECK_CLOSE(S.Box.Length.x, 3.0, 1e-9);
  CHECK_CLOSE(S.Box.Length.y, 4.0, 1e-9);
  CHECK_CLOSE(S.Box.Length.z, 5.0, 1e-9);
  CHECK_CLOSE(S.Box.Volume, 60.0, 1e-6);

  FILE *fr = fopen(in.coor.name, "r");
  CHECK(fr != nullptr);
  int lc = 0;
  CHECK(ReadTimestep(in, fr, &S, &lc));
  CHECK(S.Count.BeadCoor == 13);
  // frame 1, indexed line "3  1.8996 1.3909 2.1082" (vtf ids are 0-based)
  CHECK_CLOSE(S.Bead[3].Position.x, 1.8996, 1e-4);
  CHECK_CLOSE(S.Bead[3].Position.y, 1.3909, 1e-4);
  CHECK_CLOSE(S.Bead[3].Position.z, 2.1082, 1e-4);
  CHECK(ReadTimestep(in, fr, &S, &lc));  // frame 2
  CHECK(!ReadTimestep(in, fr, &S, &lc)); // EOF
  fclose(fr);
  FreeSystem(&S);
} //}}}

// LAMMPS data: full topology in one config //{{{
static void test_data_valid(void) {
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, FIX("system.data"), LINE);
  in.stru.type = LDATA_FILE;
  in.coor = in.stru;

  SYSTEM S = ReadStructure(in, false);
  CHECK(S.Count.Bead == 13);
  CHECK(S.Count.BeadType == 3);
  CHECK(S.Count.Molecule == 2);
  CHECK(S.Count.MoleculeType == 1);
  CHECK(S.Count.Bond == 8);
  CHECK(S.Count.Angle == 6);
  CHECK(S.Count.Dihedral == 2);
  CHECK(S.Count.Improper == 2);
  CHECK_CLOSE(S.Box.Length.x, 3.0, 1e-9);
  CHECK_CLOSE(S.Box.Length.y, 4.0, 1e-9);
  CHECK_CLOSE(S.Box.Length.z, 5.0, 1e-9);
  CHECK_CLOSE(S.Box.Volume, 60.0, 1e-6);
  FreeSystem(&S);
} //}}}

// DL_MESO FIELD: topology only, no box, no coordinates //{{{
static void test_field_valid(void) {
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, FIX("FIELD"), LINE);
  in.stru.type = FIELD_FILE;
  in.coor = in.stru;

  SYSTEM S = ReadStructure(in, false);
  CHECK(S.Count.Bead == 13);        // 5 unbonded (A=1, C=4) + 2 x 4 in mols
  CHECK(S.Count.BeadType == 3);     // A, B, C
  CHECK(S.Count.Molecule == 2);
  CHECK(S.Count.MoleculeType == 1);
  CHECK(S.Count.Bond == 8);
  CHECK(S.Count.Angle == 6);
  CHECK(S.Count.Dihedral == 2);
  CHECK(S.Count.Improper == 2);
  CHECK(S.Count.Unbonded == 5);
  CHECK(S.Count.Bonded == 8);
  CHECK(S.Box.Volume == -1);        // FIELD carries no box
  FreeSystem(&S);
} //}}}

// LAMMPS trajectory used standalone (structure inferred, no topology) //{{{
static void test_ltrj_valid(void) {
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, FIX("traj.lammpstrj"), LINE);
  in.stru.type = LTRJ_FILE;
  in.coor = in.stru;

  SYSTEM S = ReadStructure(in, false);
  CHECK(S.Count.Bead == 13);
  CHECK(S.Count.BeadType == 3); // from the element column: C, A, B
  CHECK(S.Count.Molecule == 0); // no bonds without a structure file
  CHECK(S.Count.Bond == 0);
  CHECK_CLOSE(S.Box.Length.x, 3.0, 1e-9);
  CHECK_CLOSE(S.Box.Length.y, 4.0, 1e-9);
  CHECK_CLOSE(S.Box.Length.z, 5.0, 1e-9);

  FILE *fr = fopen(in.coor.name, "r");
  CHECK(fr != nullptr);
  int lc = 0;
  CHECK(ReadTimestep(in, fr, &S, &lc));
  CHECK(S.Count.BeadCoor == 13);
  // frame 1, "4  C  1.8996 1.3909 2.1082" (ltrj ids are 1-based -> Bead[3])
  CHECK_CLOSE(S.Bead[3].Position.x, 1.8996, 1e-4);
  CHECK_CLOSE(S.Bead[3].Position.y, 1.3909, 1e-4);
  CHECK_CLOSE(S.Bead[3].Position.z, 2.1082, 1e-4);
  CHECK(ReadTimestep(in, fr, &S, &lc));  // frame 2
  CHECK(!ReadTimestep(in, fr, &S, &lc)); // EOF
  fclose(fr);
  FreeSystem(&S);
} //}}}

// XYZ used standalone (structure inferred, no box) //{{{
static void test_xyz_valid(void) {
  SYS_FILES in = InitSysFiles;
  s_strcpy(in.stru.name, FIX("traj.xyz"), LINE);
  in.stru.type = XYZ_FILE;
  in.coor = in.stru;

  SYSTEM S = ReadStructure(in, false);
  CHECK(S.Count.Bead == 13);
  CHECK(S.Count.BeadType == 3);
  CHECK(S.Count.Molecule == 0);
  CHECK(S.Box.Volume == -1); // xyz carries no box information on read

  FILE *fr = fopen(in.coor.name, "r");
  CHECK(fr != nullptr);
  int lc = 0;
  CHECK(ReadTimestep(in, fr, &S, &lc));
  CHECK(S.Count.BeadCoor == 13);
  // frame 1 first line "C  0.0256 0.3966 1.6244" (xyz is ordered -> Bead[0])
  CHECK_CLOSE(S.Bead[0].Position.x, 0.0256, 1e-4);
  CHECK_CLOSE(S.Bead[0].Position.y, 0.3966, 1e-4);
  CHECK_CLOSE(S.Bead[0].Position.z, 1.6244, 1e-4);
  CHECK(ReadTimestep(in, fr, &S, &lc));  // frame 2
  CHECK(!ReadTimestep(in, fr, &S, &lc)); // EOF
  fclose(fr);
  FreeSystem(&S);
} //}}}

// every valid fixture must also be accepted by the forked reader (harness) //{{{
static void test_children_accept_valid(void) {
  const struct { const char *path; int type; } ok[] = {
    { FIX("struct.vtf"),     VTF_FILE },
    { FIX("system.data"),    LDATA_FILE },
    { FIX("FIELD"),          FIELD_FILE },
    { FIX("traj.lammpstrj"), LTRJ_FILE },
    { FIX("traj.xyz"),       XYZ_FILE },
  };
  for (size_t i = 0; i < sizeof ok / sizeof *ok; i++) {
    char *content = slurp(ok[i].path);
    CHECK(content != nullptr);
    if (!content) {
      continue;
    }
    const char *path = write_temp(content);
    struct outcome o = read_in_child(path, ok[i].type);
    if (!(o.exited && o.code == 0)) {
      g_failures++;
      fprintf(stderr, "  FAIL: valid fixture rejected: %s (exit %d, sig %d)\n",
              ok[i].path, o.code, o.sig);
    }
    g_checks++;
    unlink(path);
    free(content);
  }
} //}}}

// ---- malformed input ------------------------------------------------------

static void test_vtf_malformed(void) {
  char *v = slurp(FIX("struct.vtf"));
  CHECK(v != nullptr);
  if (!v) {
    return;
  }
  // bond to atom id 99 (> highest atom 12): regression for commit 88031eb
  expect_rejected("vtf bond-index-overflow", VTF_FILE,
                  replace_first(v, "# resid 2\n",
                                "bond      5:     99\n# resid 2\n"));
  expect_no_crash("vtf negative-atom-id", VTF_FILE,
                  replace_first(v, "atom       1 name", "atom      -1 name"));
  expect_no_crash("vtf non-numeric-bond", VTF_FILE,
                  replace_first(v, "bond      5:      6", "bond      x:      6"));
  expect_no_crash("vtf short-pbc", VTF_FILE,
                  replace_first(v, "pbc 3.000 4.000 5.000", "pbc 3.000 4.000"));
  expect_no_crash("vtf truncated", VTF_FILE,
                  truncate_at(FIX("struct.vtf"), "atom       7"));
  free(v);
}

static void test_data_malformed(void) {
  char *v = slurp(FIX("system.data"));
  CHECK(v != nullptr);
  if (!v) {
    return;
  }
  expect_no_crash("data atom-count-inflated", LDATA_FILE,
                  replace_first(v, "13 atoms", "40 atoms"));
  expect_no_crash("data bond-count-inflated", LDATA_FILE,
                  replace_first(v, "8 bonds", "40 bonds"));
  expect_no_crash("data non-numeric-mass", LDATA_FILE,
                  replace_first(v, "1 1.000000 # C", "1 xxxxx # C"));
  expect_no_crash("data truncated", LDATA_FILE,
                  truncate_at(FIX("system.data"), "Bonds"));
  // atom type outside 1..(atom types) would index name_mass out of bounds
  // (too high: 0-based type == atom_types; zero: 0-based type == -1)
  expect_rejected("data atom-type-too-high", LDATA_FILE,
                  replace_first(v, "     0     1        0.000000",
                                "     0     4        0.000000"));
  expect_rejected("data atom-type-zero", LDATA_FILE,
                  replace_first(v, "     0     1        0.000000",
                                "     0     0        0.000000"));
  free(v);
}

static void test_field_malformed(void) {
  char *v = slurp(FIX("FIELD"));
  CHECK(v != nullptr);
  if (!v) {
    return;
  }
  expect_no_crash("field species-count-inflated", FIELD_FILE,
                  replace_first(v, "species 3", "species 6"));
  expect_no_crash("field beads-count-inflated", FIELD_FILE,
                  replace_first(v, "beads 4", "beads 9"));
  expect_no_crash("field non-numeric-mass", FIELD_FILE,
                  replace_first(v, "A   0.8", "A   xxx"));
  expect_no_crash("field truncated", FIELD_FILE,
                  truncate_at(FIX("FIELD"), "bonds 4"));
  free(v);
}

static void test_ltrj_malformed(void) {
  char *v = slurp(FIX("traj.lammpstrj"));
  CHECK(v != nullptr);
  if (!v) {
    return;
  }
  // atom id 99 in frame 2 (> 13): regression for the commit 88031eb ltrj fix
  // (id is used as an array index; without the guard this overflows)
  expect_no_crash("ltrj atom-id-overflow", LTRJ_FILE,
                  replace_first(v, "       2        A   1.5253",
                                "      99        A   1.5253"));
  expect_no_crash("ltrj atom-count-mismatch", LTRJ_FILE,
                  replace_first(v, "OF ATOMS\n13", "OF ATOMS\n99"));
  expect_no_crash("ltrj non-numeric-coord", LTRJ_FILE,
                  replace_first(v, "1.8996", "xxxxx"));
  expect_no_crash("ltrj truncated", LTRJ_FILE,
                  truncate_at(FIX("traj.lammpstrj"), "      12        B"));
  // 'ITEM: ATOMS' without the mandatory 'id' column must be rejected, not
  // dereference split[-1] in LtrjReadCoorLine()
  expect_rejected("ltrj missing-id-column", LTRJ_FILE,
                  replace_first(v, "ITEM: ATOMS id element x y z",
                                "ITEM: ATOMS element x y z"));
  free(v);
}

static void test_xyz_malformed(void) {
  char *v = slurp(FIX("traj.xyz"));
  CHECK(v != nullptr);
  if (!v) {
    return;
  }
  expect_no_crash("xyz count-mismatch", XYZ_FILE,
                  replace_first(v, "13\npbc", "99\npbc"));
  expect_no_crash("xyz non-numeric-coord", XYZ_FILE,
                  replace_first(v, "0.0256", "xxxxx"));
  expect_no_crash("xyz truncated", XYZ_FILE,
                  truncate_at(FIX("traj.xyz"), "1.2162"));
  // a line with >=32 whitespace tokens must not overflow the split[] array
  // (the xyz comment line is split but otherwise ignored)
  expect_no_crash("xyz long-token-line", XYZ_FILE,
                  replace_first(v, "pbc 3.000 4.000 5.000",
                                "x x x x x x x x x x x x x x x x x x x x "
                                "x x x x x x x x x x x x x x x x x x x x"));
  free(v);
}

int main(void) {
  RUN(test_vtf_valid);
  RUN(test_data_valid);
  RUN(test_field_valid);
  RUN(test_ltrj_valid);
  RUN(test_xyz_valid);
  RUN(test_children_accept_valid);
  RUN(test_vtf_malformed);
  RUN(test_data_malformed);
  RUN(test_field_malformed);
  RUN(test_ltrj_malformed);
  RUN(test_xyz_malformed);
  return test_main_end();
}
