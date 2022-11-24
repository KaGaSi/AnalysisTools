#include "../AnalysisTools.h"

#define e 1.60217646e-19         // electron charge
#define N_A 6.0221367e23         // Avogadro
#define k_B 1.380658e-23         // Boltzman
#define eVA3_to_Pa (e * 1e30)    // eV/Angstrom^3 to Pa
#define Jmol_to_eV 1.03642755e-5 // J/mol to eV/particle

void Help(char cmd[50], bool error) { //{{{
  FILE *ptr;
  if (error) {
    ptr = stderr;
  } else {
    ptr = stdout;
    fprintf(ptr, "\
SuttenChenParameters calculates parameters for simulations using Sutten-Chen \
potentials (https://doi.org/10.1080/09500839008206493).\n\n");
  }

  fprintf(ptr, "Usage:\n");
  fprintf(ptr, "   %s <input> <element> [options]\n\n", cmd);

  fprintf(ptr, "      <input>       input coordinate file \
(vcf or vtf format)\n");
  fprintf(ptr, "   [options]\n");
  fprintf(ptr, "      -m <int>      m parameter for Sutton-Chen potential\n");
  fprintf(ptr, "      -n <int>      n parameter for Sutton-Chen potential\n");
  fprintf(ptr, "      -mu <float>   atomic mass (kg/mol)\n");
  fprintf(ptr, "      -B <float>    bulk modulus (Pa)\n");
  fprintf(ptr, "      -E <float>    cohesion energy (eV)\n");
  fprintf(ptr, "      -rho <float>  density (kg/m3)\n");
  fprintf(ptr, "      -a <float>    lattice constant (Å)\n");
  fprintf(ptr, "      --fcc/bcc/hcp crystal lattice\n");
  CommonHelp(error);
} //}}}

typedef struct element {
  char Symbol[3];
  int State, Radius, Crystal;
  double Mass, IonizationEnergy, MeltingPoint, BoilingPoint, Density,
      BulkModulus, CohesionEnergy;
} ELEMENT;

// fill element array with info //{{{
void FillElements(ELEMENT element[]) {
  // Crystal: fcc = 1; bcc = 2; hcp = 3
  // State: gas = 0; solid = 1; liquid = 2
  // Mass in kg/mol
  // Radius in pm
  // IonizationEnergy in eV
  // Melting/BoilingPoint in K
  // Density in kg/m³
  // BulkModulus in Pa
  // CohesionEnergy in J/mol
  strcpy(element[1].Symbol, "H");
  element[1].State = 0;
  element[1].Mass = 1.0080e-3;
  element[1].Radius = 120;
  element[1].IonizationEnergy = 13.598;
  element[1].MeltingPoint = 13.81;
  element[1].BoilingPoint = 20.28;
  element[1].Density = 0.00008988e3;
  strcpy(element[2].Symbol, "He");
  element[2].State = 0;
  element[2].Mass = 4.00260e-3;
  element[2].Radius = 140;
  element[2].IonizationEnergy = 24.587;
  element[2].MeltingPoint = 0.95;
  element[2].BoilingPoint = 4.22;
  element[2].Density = 0.0001785e3;
  strcpy(element[3].Symbol, "Li");
  element[3].State = 1;
  element[3].Mass = 7.0e-3;
  element[3].Radius = 182;
  element[3].IonizationEnergy = 5.392;
  element[3].MeltingPoint = 453.65;
  element[3].BoilingPoint = 1615;
  element[3].Density = 0.534e3;
  strcpy(element[4].Symbol, "Be");
  element[4].State = 1;
  element[4].Mass = 9.012183e-3;
  element[4].Radius = 153;
  element[4].IonizationEnergy = 9.323;
  element[4].MeltingPoint = 1560;
  element[4].BoilingPoint = 2744;
  element[4].Density = 1.85e3;
  strcpy(element[5].Symbol, "B");
  element[5].State = 1;
  element[5].Mass = 10.81e-3;
  element[5].Radius = 192;
  element[5].IonizationEnergy = 8.298;
  element[5].MeltingPoint = 2348;
  element[5].BoilingPoint = 4273;
  element[5].Density = 2.37e3;
  strcpy(element[6].Symbol, "C");
  element[6].State = 1;
  element[6].Mass = 12.011e-3;
  element[6].Radius = 170;
  element[6].IonizationEnergy = 11.260;
  element[6].MeltingPoint = 3823;
  element[6].BoilingPoint = 4098;
  element[6].Density = 2.2670e3;
  strcpy(element[7].Symbol, "N");
  element[7].State = 0;
  element[7].Mass = 14.007e-3;
  element[7].Radius = 155;
  element[7].IonizationEnergy = 14.534;
  element[7].MeltingPoint = 63.15;
  element[7].BoilingPoint = 77.36;
  element[7].Density = 0.0012506e3;
  strcpy(element[8].Symbol, "O");
  element[8].State = 0;
  element[8].Mass = 15.999e-3;
  element[8].Radius = 152;
  element[8].IonizationEnergy = 13.618;
  element[8].MeltingPoint = 54.36;
  element[8].BoilingPoint = 90.2;
  element[8].Density = 0.001429e3;
  strcpy(element[9].Symbol, "F");
  element[9].State = 0;
  element[9].Mass = 18.99840316e-3;
  element[9].Radius = 135;
  element[9].IonizationEnergy = 17.423;
  element[9].MeltingPoint = 53.53;
  element[9].BoilingPoint = 85.03;
  element[9].Density = 0.001696e3;
  strcpy(element[10].Symbol, "Ne");
  element[10].State = 0;
  element[10].Mass = 20.180e-3;
  element[10].Radius = 154;
  element[10].IonizationEnergy = 21.565;
  element[10].MeltingPoint = 24.56;
  element[10].BoilingPoint = 27.07;
  element[10].Density = 0.0008999e3;
  strcpy(element[11].Symbol, "Na");
  element[11].State = 1;
  element[11].Mass = 22.9897693e-3;
  element[11].Radius = 227;
  element[11].IonizationEnergy = 5.139;
  element[11].MeltingPoint = 370.95;
  element[11].BoilingPoint = 1156;
  element[11].Density = 0.97e3;
  strcpy(element[12].Symbol, "Mg");
  element[12].State = 1;
  element[12].Mass = 24.305e-3;
  element[12].Radius = 173;
  element[12].IonizationEnergy = 7.646;
  element[12].MeltingPoint = 923;
  element[12].BoilingPoint = 1363;
  element[12].Density = 1.74e3;
  element[12].Crystal = 3; // hcp
//element[12].Crystal = 1; // fcc
  element[12].BulkModulus = 35.4e9;
  element[12].CohesionEnergy = 1.51 / Jmol_to_eV;
  strcpy(element[13].Symbol, "Al");
  element[13].State = 1;
  element[13].Mass = 26.982e-3;
  element[13].Radius = 184;
  element[13].IonizationEnergy = 5.986;
  element[13].MeltingPoint = 933.437;
  element[13].BoilingPoint = 2792;
  element[13].Density = 2.70e3;
  element[13].Crystal = 1;
  element[13].BulkModulus = 0.48 * eVA3_to_Pa;
//element[13].BulkModulus = 77.0e9; // ml
  element[13].CohesionEnergy = 3.39 / Jmol_to_eV;
//element[13].CohesionEnergy = 322.3e3; // ml
  strcpy(element[14].Symbol, "Si");
  element[14].State = 1;
  element[14].Mass = 28.085e-3;
  element[14].Radius = 210;
  element[14].IonizationEnergy = 8.152;
  element[14].MeltingPoint = 1687;
  element[14].BoilingPoint = 3538;
  element[14].Density = 2.3296e3;
  strcpy(element[15].Symbol, "P");
  element[15].State = 1;
  element[15].Mass = 30.97376200e-3;
  element[15].Radius = 180;
  element[15].IonizationEnergy = 10.487;
  element[15].MeltingPoint = 317.3;
  element[15].BoilingPoint = 553.65;
  element[15].Density = 1.82e3;
  strcpy(element[16].Symbol, "S");
  element[16].State = 1;
  element[16].Mass = 32.07e-3;
  element[16].Radius = 180;
  element[16].IonizationEnergy = 10.360;
  element[16].MeltingPoint = 388.36;
  element[16].BoilingPoint = 717.75;
  element[16].Density = 2.067e3;
  strcpy(element[17].Symbol, "Cl");
  element[17].State = 0;
  element[17].Mass = 35.45e-3;
  element[17].Radius = 175;
  element[17].IonizationEnergy = 12.968;
  element[17].MeltingPoint = 171.65;
  element[17].BoilingPoint = 239.11;
  element[17].Density = 0.003214e3;
  strcpy(element[18].Symbol, "Ar");
  element[18].State = 0;
  element[18].Mass = 39.9e-3;
  element[18].Radius = 188;
  element[18].IonizationEnergy = 15.760;
  element[18].MeltingPoint = 83.8;
  element[18].BoilingPoint = 87.3;
  element[18].Density = 0.0017837e3;
  strcpy(element[19].Symbol, "K");
  element[19].State = 1;
  element[19].Mass = 39.0983e-3;
  element[19].Radius = 275;
  element[19].IonizationEnergy = 4.341;
  element[19].MeltingPoint = 336.53;
  element[19].BoilingPoint = 1032;
  element[19].Density = 0.89e3;
  strcpy(element[20].Symbol, "Ca");
  element[20].State = 1;
  element[20].Mass = 40.08e-3;
  element[20].Radius = 231;
  element[20].IonizationEnergy = 6.113;
  element[20].MeltingPoint = 1115;
  element[20].BoilingPoint = 1757;
  element[20].Density = 1.54e3;
  element[20].Crystal = 1;
  element[20].BulkModulus = 17e9;
  element[20].CohesionEnergy = 178e3;
  strcpy(element[21].Symbol, "Sc");
  element[21].State = 1;
  element[21].Mass = 44.95591e-3;
  element[21].Radius = 211;
  element[21].IonizationEnergy = 6.561;
  element[21].MeltingPoint = 1814;
  element[21].BoilingPoint = 3109;
  element[21].Density = 2.99e3;
  strcpy(element[22].Symbol, "Ti");
  element[22].State = 1;
  element[22].Mass = 47.867e-3;
  element[22].Radius = 187;
  element[22].IonizationEnergy = 6.828;
  element[22].MeltingPoint = 1941;
  element[22].BoilingPoint = 3560;
  element[22].Density = 4.5e3;
  strcpy(element[23].Symbol, "V");
  element[23].State = 1;
  element[23].Mass = 50.9415e-3;
  element[23].Radius = 179;
  element[23].IonizationEnergy = 6.746;
  element[23].MeltingPoint = 2183;
  element[23].BoilingPoint = 3680;
  element[23].Density = 6.0e3;
  strcpy(element[24].Symbol, "Cr");
  element[24].State = 1;
  element[24].Mass = 51.996e-3;
  element[24].Radius = 189;
  element[24].IonizationEnergy = 6.767;
  element[24].MeltingPoint = 2180;
  element[24].BoilingPoint = 2944;
  element[24].Density = 7.15e3;
  strcpy(element[25].Symbol, "Mn");
  element[25].State = 1;
  element[25].Mass = 54.93804e-3;
  element[25].Radius = 197;
  element[25].IonizationEnergy = 7.434;
  element[25].MeltingPoint = 1519;
  element[25].BoilingPoint = 2334;
  element[25].Density = 7.3e3;
  strcpy(element[26].Symbol, "Fe");
  element[26].State = 1;
  element[26].Mass = 55.84e-3;
  element[26].Radius = 194;
  element[26].IonizationEnergy = 7.902;
  element[26].MeltingPoint = 1811;
  element[26].BoilingPoint = 3134;
  element[26].Density = 7.874e3;
  strcpy(element[27].Symbol, "Co");
  element[27].State = 1;
  element[27].Mass = 58.93319e-3;
  element[27].Radius = 192;
  element[27].IonizationEnergy = 7.881;
  element[27].MeltingPoint = 1768;
  element[27].BoilingPoint = 3200;
  element[27].Density = 8.86e3;
  strcpy(element[28].Symbol, "Ni");
  element[28].State = 1;
  element[28].Mass = 58.693e-3;
  element[28].Radius = 163;
  element[28].IonizationEnergy = 7.640;
  element[28].MeltingPoint = 1728;
  element[28].BoilingPoint = 3186;
  element[28].Density = 8.912e3;
  element[28].Crystal = 1;
  element[28].BulkModulus = 187.5e9;
  element[28].CohesionEnergy = 428.4e3;
  strcpy(element[29].Symbol, "Cu");
  element[29].State = 1;
  element[29].Mass = 63.55e-3;
  element[29].Radius = 140;
  element[29].IonizationEnergy = 7.726;
  element[29].MeltingPoint = 1357.77;
  element[29].BoilingPoint = 2835;
  element[29].Density = 8.933e3;
  element[29].Crystal = 1;
  element[29].BulkModulus = 0.89 * eVA3_to_Pa;
  element[29].CohesionEnergy = 3.50 / Jmol_to_eV;
  strcpy(element[30].Symbol, "Zn");
  element[30].State = 1;
  element[30].Mass = 65.4e-3;
  element[30].Radius = 139;
  element[30].IonizationEnergy = 9.394;
  element[30].MeltingPoint = 692.68;
  element[30].BoilingPoint = 1180;
  element[30].Density = 7.134e3;
  strcpy(element[31].Symbol, "Ga");
  element[31].State = 1;
  element[31].Mass = 69.723e-3;
  element[31].Radius = 187;
  element[31].IonizationEnergy = 5.999;
  element[31].MeltingPoint = 302.91;
  element[31].BoilingPoint = 2477;
  element[31].Density = 5.91e3;
  strcpy(element[32].Symbol, "Ge");
  element[32].State = 1;
  element[32].Mass = 72.63e-3;
  element[32].Radius = 211;
  element[32].IonizationEnergy = 7.900;
  element[32].MeltingPoint = 1211.4;
  element[32].BoilingPoint = 3106;
  element[32].Density = 5.323e3;
  strcpy(element[33].Symbol, "As");
  element[33].State = 1;
  element[33].Mass = 74.92159e-3;
  element[33].Radius = 185;
  element[33].IonizationEnergy = 9.815;
  element[33].MeltingPoint = 1090;
  element[33].BoilingPoint = 887;
  element[33].Density = 5.776e3;
  strcpy(element[34].Symbol, "Se");
  element[34].State = 1;
  element[34].Mass = 78.97e-3;
  element[34].Radius = 190;
  element[34].IonizationEnergy = 9.752;
  element[34].MeltingPoint = 493.65;
  element[34].BoilingPoint = 958;
  element[34].Density = 4.809e3;
  strcpy(element[35].Symbol, "Br");
  element[35].State = 2;
  element[35].Mass = 79.90e-3;
  element[35].Radius = 183;
  element[35].IonizationEnergy = 11.814;
  element[35].MeltingPoint = 265.95;
  element[35].BoilingPoint = 331.95;
  element[35].Density = 3.11e3;
  strcpy(element[36].Symbol, "Kr");
  element[36].State = 0;
  element[36].Mass = 83.80e-3;
  element[36].Radius = 202;
  element[36].IonizationEnergy = 14.000;
  element[36].MeltingPoint = 115.79;
  element[36].BoilingPoint = 119.93;
  element[36].Density = 0.003733e3;
  strcpy(element[37].Symbol, "Rb");
  element[37].State = 1;
  element[37].Mass = 85.468e-3;
  element[37].Radius = 303;
  element[37].IonizationEnergy = 4.177;
  element[37].MeltingPoint = 312.46;
  element[37].BoilingPoint = 961;
  element[37].Density = 1.53e3;
  strcpy(element[38].Symbol, "Sr");
  element[38].State = 1;
  element[38].Mass = 87.62e-3;
  element[38].Radius = 249;
  element[38].IonizationEnergy = 5.695;
  element[38].MeltingPoint = 1050;
  element[38].BoilingPoint = 1655;
  element[38].Density = 2.64e3;
  strcpy(element[39].Symbol, "Y");
  element[39].State = 1;
  element[39].Mass = 88.90584e-3;
  element[39].Radius = 219;
  element[39].IonizationEnergy = 6.217;
  element[39].MeltingPoint = 1795;
  element[39].BoilingPoint = 3618;
  element[39].Density = 4.47e3;
  strcpy(element[40].Symbol, "Zr");
  element[40].State = 1;
  element[40].Mass = 91.22e-3;
  element[40].Radius = 186;
  element[40].IonizationEnergy = 6.634;
  element[40].MeltingPoint = 2128;
  element[40].BoilingPoint = 4682;
  element[40].Density = 6.52e3;
  element[40].Crystal = 3; // hcp
//element[40].Crystal = 1; // fcc
  element[40].BulkModulus = 97.2e9;               // chapter 5 of that thingy
  element[40].CohesionEnergy = 6.21 / Jmol_to_eV; // chapter 5 of that thingy
  strcpy(element[41].Symbol, "Nb");
  element[41].State = 1;
  element[41].Mass = 92.90637e-3;
  element[41].Radius = 207;
  element[41].IonizationEnergy = 6.759;
  element[41].MeltingPoint = 2750;
  element[41].BoilingPoint = 5017;
  element[41].Density = 8.57e3;
  strcpy(element[42].Symbol, "Mo");
  element[42].State = 1;
  element[42].Mass = 95.95e-3;
  element[42].Radius = 209;
  element[42].IonizationEnergy = 7.092;
  element[42].MeltingPoint = 2896;
  element[42].BoilingPoint = 4912;
  element[42].Density = 10.2e3;
  strcpy(element[43].Symbol, "Tc");
  element[43].State = 1;
  element[43].Mass = 96.90636e-3;
  element[43].Radius = 209;
  element[43].IonizationEnergy = 7.28;
  element[43].MeltingPoint = 2430;
  element[43].BoilingPoint = 4538;
  element[43].Density = 11e3;
  strcpy(element[44].Symbol, "Ru");
  element[44].State = 1;
  element[44].Mass = 101.1e-3;
  element[44].Radius = 207;
  element[44].IonizationEnergy = 7.361;
  element[44].MeltingPoint = 2607;
  element[44].BoilingPoint = 4423;
  element[44].Density = 12.1e3;
  strcpy(element[45].Symbol, "Rh");
  element[45].State = 1;
  element[45].Mass = 102.9055e-3;
  element[45].Radius = 195;
  element[45].IonizationEnergy = 7.459;
  element[45].MeltingPoint = 2237;
  element[45].BoilingPoint = 3968;
  element[45].Density = 12.4e3;
  element[45].Crystal = 1;
  element[45].BulkModulus = 1.68 * eVA3_to_Pa;
  element[45].CohesionEnergy = 5.75 / Jmol_to_eV;
  strcpy(element[46].Symbol, "Pd");
  element[46].State = 1;
  element[46].Mass = 106.42e-3;
  element[46].Radius = 202;
  element[46].IonizationEnergy = 8.337;
  element[46].MeltingPoint = 1828.05;
  element[46].BoilingPoint = 3236;
  element[46].Density = 12.0e3;
  element[46].Crystal = 1;
  element[46].BulkModulus = 1.22 * eVA3_to_Pa;
  element[46].CohesionEnergy = 3.94 / Jmol_to_eV;
  strcpy(element[47].Symbol, "Ag");
  element[47].State = 1;
  element[47].Mass = 107.868e-3;
  element[47].Radius = 172;
  element[47].IonizationEnergy = 7.576;
  element[47].MeltingPoint = 1234.93;
  element[47].BoilingPoint = 2435;
  element[47].Density = 10.501e3;
  element[47].Crystal = 1;
  element[47].BulkModulus = 0.68 * eVA3_to_Pa;
  element[47].CohesionEnergy = 2.96 / Jmol_to_eV;
  strcpy(element[48].Symbol, "Cd");
  element[48].State = 1;
  element[48].Mass = 112.41e-3;
  element[48].Radius = 158;
  element[48].IonizationEnergy = 8.994;
  element[48].MeltingPoint = 594.22;
  element[48].BoilingPoint = 1040;
  element[48].Density = 8.69e3;
  strcpy(element[49].Symbol, "In");
  element[49].State = 1;
  element[49].Mass = 114.818e-3;
  element[49].Radius = 193;
  element[49].IonizationEnergy = 5.786;
  element[49].MeltingPoint = 429.75;
  element[49].BoilingPoint = 2345;
  element[49].Density = 7.31e3;
  strcpy(element[50].Symbol, "Sn");
  element[50].State = 1;
  element[50].Mass = 118.71e-3;
  element[50].Radius = 217;
  element[50].IonizationEnergy = 7.344;
  element[50].MeltingPoint = 505.08;
  element[50].BoilingPoint = 2875;
  element[50].Density = 7.287e3;
  strcpy(element[51].Symbol, "Sb");
  element[51].State = 1;
  element[51].Mass = 121.760e-3;
  element[51].Radius = 206;
  element[51].IonizationEnergy = 8.64;
  element[51].MeltingPoint = 903.78;
  element[51].BoilingPoint = 1860;
  element[51].Density = 6.685e3;
  strcpy(element[52].Symbol, "Te");
  element[52].State = 1;
  element[52].Mass = 127.6e-3;
  element[52].Radius = 206;
  element[52].IonizationEnergy = 9.010;
  element[52].MeltingPoint = 722.66;
  element[52].BoilingPoint = 1261;
  element[52].Density = 6.232e3;
  strcpy(element[53].Symbol, "I");
  element[53].State = 1;
  element[53].Mass = 126.9045e-3;
  element[53].Radius = 198;
  element[53].IonizationEnergy = 10.451;
  element[53].MeltingPoint = 386.85;
  element[53].BoilingPoint = 457.55;
  element[53].Density = 4.93e3;
  strcpy(element[54].Symbol, "Xe");
  element[54].State = 0;
  element[54].Mass = 131.29e-3;
  element[54].Radius = 216;
  element[54].IonizationEnergy = 12.130;
  element[54].MeltingPoint = 161.36;
  element[54].BoilingPoint = 165.03;
  element[54].Density = 0.005887e3;
  strcpy(element[55].Symbol, "Cs");
  element[55].State = 1;
  element[55].Mass = 132.9054520e-3;
  element[55].Radius = 343;
  element[55].IonizationEnergy = 3.894;
  element[55].MeltingPoint = 301.59;
  element[55].BoilingPoint = 944;
  element[55].Density = 1.93e3;
  strcpy(element[56].Symbol, "Ba");
  element[56].State = 1;
  element[56].Mass = 137.33e-3;
  element[56].Radius = 268;
  element[56].IonizationEnergy = 5.212;
  element[56].MeltingPoint = 1000;
  element[56].BoilingPoint = 2170;
  element[56].Density = 3.62e3;
  strcpy(element[57].Symbol, "La");
  element[57].State = 1;
  element[57].Mass = 138.9055e-3;
  element[57].Radius = 240;
  element[57].IonizationEnergy = 5.577;
  element[57].MeltingPoint = 1191;
  element[57].BoilingPoint = 3737;
  element[57].Density = 6.15e3;
  strcpy(element[58].Symbol, "Ce");
  element[58].State = 1;
  element[58].Mass = 140.116e-3;
  element[58].Radius = 235;
  element[58].IonizationEnergy = 5.539;
  element[58].MeltingPoint = 1071;
  element[58].BoilingPoint = 3697;
  element[58].Density = 6.770e3;
  strcpy(element[59].Symbol, "Pr");
  element[59].State = 1;
  element[59].Mass = 140.90766e-3;
  element[59].Radius = 239;
  element[59].IonizationEnergy = 5.464;
  element[59].MeltingPoint = 1204;
  element[59].BoilingPoint = 3793;
  element[59].Density = 6.77e3;
  strcpy(element[60].Symbol, "Nd");
  element[60].State = 1;
  element[60].Mass = 144.24e-3;
  element[60].Radius = 229;
  element[60].IonizationEnergy = 5.525;
  element[60].MeltingPoint = 1294;
  element[60].BoilingPoint = 3347;
  element[60].Density = 7.01e3;
  strcpy(element[61].Symbol, "Pm");
  element[61].State = 1;
  element[61].Mass = 144.91276e-3;
  element[61].Radius = 236;
  element[61].IonizationEnergy = 5.55;
  element[61].MeltingPoint = 1315;
  element[61].BoilingPoint = 3273;
  element[61].Density = 7.26e3;
  strcpy(element[62].Symbol, "Sm");
  element[62].State = 1;
  element[62].Mass = 150.4e-3;
  element[62].Radius = 229;
  element[62].IonizationEnergy = 5.644;
  element[62].MeltingPoint = 1347;
  element[62].BoilingPoint = 2067;
  element[62].Density = 7.52e3;
  strcpy(element[63].Symbol, "Eu");
  element[63].State = 1;
  element[63].Mass = 151.964e-3;
  element[63].Radius = 233;
  element[63].IonizationEnergy = 5.670;
  element[63].MeltingPoint = 1095;
  element[63].BoilingPoint = 1802;
  element[63].Density = 5.24e3;
  strcpy(element[64].Symbol, "Gd");
  element[64].State = 1;
  element[64].Mass = 157.2e-3;
  element[64].Radius = 237;
  element[64].IonizationEnergy = 6.150;
  element[64].MeltingPoint = 1586;
  element[64].BoilingPoint = 3546;
  element[64].Density = 7.90e3;
  strcpy(element[65].Symbol, "Tb");
  element[65].State = 1;
  element[65].Mass = 158.92535e-3;
  element[65].Radius = 221;
  element[65].IonizationEnergy = 5.864;
  element[65].MeltingPoint = 1629;
  element[65].BoilingPoint = 3503;
  element[65].Density = 8.23e3;
  strcpy(element[66].Symbol, "Dy");
  element[66].State = 1;
  element[66].Mass = 162.500e-3;
  element[66].Radius = 229;
  element[66].IonizationEnergy = 5.939;
  element[66].MeltingPoint = 1685;
  element[66].BoilingPoint = 2840;
  element[66].Density = 8.55e3;
  strcpy(element[67].Symbol, "Ho");
  element[67].State = 1;
  element[67].Mass = 164.93033e-3;
  element[67].Radius = 216;
  element[67].IonizationEnergy = 6.022;
  element[67].MeltingPoint = 1747;
  element[67].BoilingPoint = 2973;
  element[67].Density = 8.80e3;
  strcpy(element[68].Symbol, "Er");
  element[68].State = 1;
  element[68].Mass = 167.26e-3;
  element[68].Radius = 235;
  element[68].IonizationEnergy = 6.108;
  element[68].MeltingPoint = 1802;
  element[68].BoilingPoint = 3141;
  element[68].Density = 9.07e3;
  strcpy(element[69].Symbol, "Tm");
  element[69].State = 1;
  element[69].Mass = 168.93422e-3;
  element[69].Radius = 227;
  element[69].IonizationEnergy = 6.184;
  element[69].MeltingPoint = 1818;
  element[69].BoilingPoint = 2223;
  element[69].Density = 9.32e3;
  strcpy(element[70].Symbol, "Yb");
  element[70].State = 1;
  element[70].Mass = 173.05e-3;
  element[70].Radius = 242;
  element[70].IonizationEnergy = 6.254;
  element[70].MeltingPoint = 1092;
  element[70].BoilingPoint = 1469;
  element[70].Density = 6.90e3;
  strcpy(element[71].Symbol, "Lu");
  element[71].State = 1;
  element[71].Mass = 174.9668e-3;
  element[71].Radius = 221;
  element[71].IonizationEnergy = 5.426;
  element[71].MeltingPoint = 1936;
  element[71].BoilingPoint = 3675;
  element[71].Density = 9.84e3;
  strcpy(element[72].Symbol, "Hf");
  element[72].State = 1;
  element[72].Mass = 178.49e-3;
  element[72].Radius = 212;
  element[72].IonizationEnergy = 6.825;
  element[72].MeltingPoint = 2506;
  element[72].BoilingPoint = 4876;
  element[72].Density = 13.3e3;
  strcpy(element[73].Symbol, "Ta");
  element[73].State = 1;
  element[73].Mass = 180.9479e-3;
  element[73].Radius = 217;
  element[73].IonizationEnergy = 7.89;
  element[73].MeltingPoint = 3290;
  element[73].BoilingPoint = 5731;
  element[73].Density = 16.4e3;
  strcpy(element[74].Symbol, "W");
  element[74].State = 1;
  element[74].Mass = 183.84e-3;
  element[74].Radius = 210;
  element[74].IonizationEnergy = 7.98;
  element[74].MeltingPoint = 3695;
  element[74].BoilingPoint = 5828;
  element[74].Density = 19.3e3;
  strcpy(element[75].Symbol, "Re");
  element[75].State = 1;
  element[75].Mass = 186.207e-3;
  element[75].Radius = 217;
  element[75].IonizationEnergy = 7.88;
  element[75].MeltingPoint = 3459;
  element[75].BoilingPoint = 5869;
  element[75].Density = 20.8e3;
  strcpy(element[76].Symbol, "Os");
  element[76].State = 1;
  element[76].Mass = 190.2e-3;
  element[76].Radius = 216;
  element[76].IonizationEnergy = 8.7;
  element[76].MeltingPoint = 3306;
  element[76].BoilingPoint = 5285;
  element[76].Density = 22.57e3;
  strcpy(element[77].Symbol, "Ir");
  element[77].State = 1;
  element[77].Mass = 192.22e-3;
  element[77].Radius = 202;
  element[77].IonizationEnergy = 9.1;
  element[77].MeltingPoint = 2719;
  element[77].BoilingPoint = 4701;
  element[77].Density = 22.42e3;
  element[77].Crystal = 1;
  element[77].BulkModulus = 2.31 * eVA3_to_Pa;
  element[77].CohesionEnergy = 6.93 / Jmol_to_eV;
  strcpy(element[78].Symbol, "Pt");
  element[78].State = 1;
  element[78].Mass = 195.08e-3;
  element[78].Radius = 209;
  element[78].IonizationEnergy = 9;
  element[78].MeltingPoint = 2041.55;
  element[78].BoilingPoint = 4098;
  element[78].Density = 21.46e3;
  element[78].Crystal = 1;
  element[78].BulkModulus = 1.80 * eVA3_to_Pa;
  element[78].CohesionEnergy = 5.86 / Jmol_to_eV;
  strcpy(element[79].Symbol, "Au");
  element[79].State = 1;
  element[79].Mass = 196.96657e-3;
  element[79].Radius = 166;
  element[79].IonizationEnergy = 9.226;
  element[79].MeltingPoint = 1337.33;
  element[79].BoilingPoint = 3129;
  element[79].Density = 19.282e3;
  element[79].Crystal = 1;
  element[79].BulkModulus = 1.03 * eVA3_to_Pa;
  element[79].CohesionEnergy = 3.78 / Jmol_to_eV;
  strcpy(element[80].Symbol, "Hg");
  element[80].State = 2;
  element[80].Mass = 200.59e-3;
  element[80].Radius = 209;
  element[80].IonizationEnergy = 10.438;
  element[80].MeltingPoint = 234.32;
  element[80].BoilingPoint = 629.88;
  element[80].Density = 13.5336e3;
  strcpy(element[81].Symbol, "Tl");
  element[81].State = 1;
  element[81].Mass = 204.383e-3;
  element[81].Radius = 196;
  element[81].IonizationEnergy = 6.108;
  element[81].MeltingPoint = 577;
  element[81].BoilingPoint = 1746;
  element[81].Density = 11.8e3;
  strcpy(element[82].Symbol, "Pb");
  element[82].State = 1;
  element[82].Mass = 207e-3;
  element[82].Radius = 202;
  element[82].IonizationEnergy = 7.417;
  element[82].MeltingPoint = 600.61;
  element[82].BoilingPoint = 2022;
  element[82].Density = 11.342e3;
  element[82].Crystal = 1;
  element[82].BulkModulus = 0.26 * eVA3_to_Pa;
  element[82].CohesionEnergy = 2.04 / Jmol_to_eV;
  strcpy(element[83].Symbol, "Bi");
  element[83].State = 1;
  element[83].Mass = 208.98040e-3;
  element[83].Radius = 207;
  element[83].IonizationEnergy = 7.289;
  element[83].MeltingPoint = 544.55;
  element[83].BoilingPoint = 1837;
  element[83].Density = 9.807e3;
  strcpy(element[84].Symbol, "Po");
  element[84].State = 1;
  element[84].Mass = 208.98243e-3;
  element[84].Radius = 197;
  element[84].IonizationEnergy = 8.417;
  element[84].MeltingPoint = 527;
  element[84].BoilingPoint = 1235;
  element[84].Density = 9.32e3;
  strcpy(element[85].Symbol, "At");
  element[85].State = 1;
  element[85].Mass = 209.98715e-3;
  element[85].Radius = 202;
  element[85].IonizationEnergy = 9.5;
  element[85].MeltingPoint = 575;
  element[85].Density = 7e3;
  strcpy(element[85].Symbol, "");
  strcpy(element[86].Symbol, "Rn");
  element[86].State = 0;
  element[86].Mass = 222.01758e-3;
  element[86].Radius = 220;
  element[86].IonizationEnergy = 10.745;
  element[86].MeltingPoint = 202;
  element[86].BoilingPoint = 211.45;
  element[86].Density = 0.00973e3;
  strcpy(element[87].Symbol, "Fr");
  element[87].State = 1;
  element[87].Mass = 223.01973e-3;
  element[87].Radius = 348;
  element[87].IonizationEnergy = 3.9;
  element[87].MeltingPoint = 300;
  strcpy(element[87].Symbol, "");
  strcpy(element[87].Symbol, "");
  strcpy(element[88].Symbol, "Ra");
  element[88].State = 1;
  element[88].Mass = 226.02541e-3;
  element[88].Radius = 283;
  element[88].IonizationEnergy = 5.279;
  element[88].MeltingPoint = 973;
  element[88].BoilingPoint = 1413;
  element[88].Density = 5e3;
  strcpy(element[89].Symbol, "Ac");
  element[89].State = 1;
  element[89].Mass = 227.02775e-3;
  element[89].Radius = 260;
  element[89].IonizationEnergy = 5.17;
  element[89].MeltingPoint = 1324;
  element[89].BoilingPoint = 3471;
  element[89].Density = 10.07e3;
  strcpy(element[90].Symbol, "Th");
  element[90].State = 1;
  element[90].Mass = 232.038e-3;
  element[90].Radius = 237;
  element[90].IonizationEnergy = 6.08;
  element[90].MeltingPoint = 2023;
  element[90].BoilingPoint = 5061;
  element[90].Density = 11.72e3;
  strcpy(element[91].Symbol, "Pa");
  element[91].State = 1;
  element[91].Mass = 231.03588e-3;
  element[91].Radius = 243;
  element[91].IonizationEnergy = 5.89;
  element[91].MeltingPoint = 1845;
  element[91].Density = 15.37e3;
  strcpy(element[92].Symbol, "U");
  element[92].State = 1;
  element[92].Mass = 238.0289e-3;
  element[92].Radius = 240;
  element[92].IonizationEnergy = 6.194;
  element[92].MeltingPoint = 1408;
  element[92].BoilingPoint = 4404;
  element[92].Density = 18.95e3;
  strcpy(element[93].Symbol, "Np");
  element[93].State = 1;
  element[93].Mass = 237.048172e-3;
  element[93].Radius = 221;
  element[93].IonizationEnergy = 6.266;
  element[93].MeltingPoint = 917;
  element[93].BoilingPoint = 4175;
  element[93].Density = 20.25e3;
  strcpy(element[94].Symbol, "Pu");
  element[94].State = 1;
  element[94].Mass = 244.06420e-3;
  element[94].Radius = 243;
  element[94].IonizationEnergy = 6.06;
  element[94].MeltingPoint = 913;
  element[94].BoilingPoint = 3501;
  element[94].Density = 19.84e3;
  strcpy(element[95].Symbol, "Am");
  element[95].State = 1;
  element[95].Mass = 243.061380e-3;
  element[95].Radius = 244;
  element[95].IonizationEnergy = 5.993;
  element[95].MeltingPoint = 1449;
  element[95].BoilingPoint = 2284;
  element[95].Density = 13.69e3;
  strcpy(element[96].Symbol, "Cm");
  element[96].State = 1;
  element[96].Mass = 247.07035e-3;
  element[96].Radius = 245;
  element[96].IonizationEnergy = 6.02;
  element[96].MeltingPoint = 1618;
  element[96].BoilingPoint = 3400;
  element[96].Density = 13.51e3;
  strcpy(element[97].Symbol, "Bk");
  element[97].State = 1;
  element[97].Mass = 247.07031e3;
  element[97].Radius = 244;
  element[97].IonizationEnergy = 6.23;
  element[97].MeltingPoint = 1323;
  element[97].Density = 14e3;
} //}}}

int main(int argc, char *argv[]) {
  // -h/--version options - print stuff and exit //{{{
  if (VersionOption(argc, argv)) {
    exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      Help(argv[0], false);
      exit(0);
    }
  }
  int req_args = 2; //}}}

  // check if correct number of arguments //{{{
  int count = 0;
  while ((count + 1) < argc && argv[count + 1][0] != '-') {
    count++;
  }
  // reverse bead type selection? ...do now to check correct number of arguments
  bool reverse = BoolOption(argc, argv, "--reverse");
  // possible to omit <bead name(s)> if '--reverse' is used
  if (count < (req_args - 1) || (count == (req_args - 1) && !reverse)) {
    ErrorArgNumber(count, req_args);
    Help(argv[0], true);
    exit(1);
  } //}}}

  // test if options are given correctly //{{{
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' && strcmp(argv[i], "-m") != 0 &&
        strcmp(argv[i], "-n") != 0 && strcmp(argv[i], "mu") != 0 &&
        strcmp(argv[i], "-B") != 0 && strcmp(argv[i], "-E") != 0 &&
        strcmp(argv[i], "-a") != 0 && strcmp(argv[i], "-mu") != 0 &&
        strcmp(argv[i], "-rho") != 0 && strcmp(argv[i], "--fcc") != 0 &&
        strcmp(argv[i], "--bcc") != 0 && strcmp(argv[i], "--hcp") != 0 &&
        strcmp(argv[i], "-i") != 0 && strcmp(argv[i], "-v") != 0 &&
        strcmp(argv[i], "--silent") != 0 && strcmp(argv[i], "--detailed") != 0 &&
        strcmp(argv[i], "--version") != 0 && strcmp(argv[i], "-h") != 0) {
      ErrorOption(argv[i]);
      Help(argv[0], true);
      exit(1);
    }
  } //}}}

  count = 0; // count mandatory arguments

  // <input> - input coordinate file //{{{
  char input_coor[LINE] = "", input_vsf[LINE] = "";
  snprintf(input_coor, LINE, "%s", argv[++count]);
  // test that <input> filename ends with '.vcf' or '.vtf'
  bool vtf;
  if (!InputCoor_old(&vtf, input_coor, input_vsf)) {
    Help(argv[0], true);
    exit(1);
  } //}}}

  // <element> - periodic table symbol for an element //{{{
  char species[3] = "\0";
  strncpy(species, argv[++count], 2);
  //}}}

  // options before reading system data
  bool silent, verbose, detailed;
  CommonOptions_old(argc, argv, input_vsf, LINE, &verbose, &silent, &detailed);

  // '-m' option //{{{
  int m_SC = 6; // no -m option
  if (IntegerOption(argc, argv, "-m", &m_SC)) {
    exit(1);
  } //}}}
  // '-n' option //{{{
  int n_SC = -1; // no -n option
  if (IntegerOption(argc, argv, "-n", &n_SC)) {
    exit(1);
  }
  if (n_SC != -1 && n_SC <= m_SC) {
    strcpy(ERROR_MSG, "in SC potentional, n > m!");
    PrintError();
    fprintf(stderr, "%s Supplied parameters: n=%s%d%s and m=%s%d%s\n",
            ErrRed(), ErrYellow(), n_SC,
            ErrRed(), ErrYellow(), m_SC, ErrColourReset());
    exit(1);
  } //}}}
  // options for element properties //{{{
  double Mw = -1, rho = -1, E_coh = -1, B = -1, a = -1;
  // atomic mass
  if (DoubleOption(argc, argv, "-mu", &Mw)) {
    exit(1);
  }
  // bulk modulus
  if (DoubleOption(argc, argv, "-B", &B)) {
    exit(1);
  }
  // cohesion energy
  if (DoubleOption(argc, argv, "-E", &E_coh)) {
    exit(1);
  }
  // lattice constant
  if (DoubleOption(argc, argv, "-a", &a)) {
    exit(1);
  }
  // density
  if (DoubleOption(argc, argv, "-rho", &rho)) {
    exit(1);
  }
  char lattice[4] = "\0";
  // lattice
  bool fcc = BoolOption(argc, argv, "--fcc");
  bool bcc = BoolOption(argc, argv, "--bcc");
  bool hcp = BoolOption(argc, argv, "--hcp");
  if ((fcc && bcc) || (fcc && hcp) || (bcc && hcp)) {
    strcpy(ERROR_MSG, "multiple lattices specified; \
preference order: fcc > bcc > hcp");
    PrintWarning();
  }
  if (fcc) {
    strcpy(lattice, "fcc");
  } else if (bcc) {
    strcpy(lattice, "bcc");
  } else if (hcp) {
    strcpy(lattice, "hcp");
  } //}}}

  // print command to stdout //{{{
  if (!silent) {
    PrintCommand(stdout, argc, argv);
  } //}}}

  // create the elements struct //{{{
  int n_elements = 118;
  ELEMENT *element = calloc(n_elements, sizeof *element);
  FillElements(element); //}}}

  // read information from vtf file(s) //{{{
  SYSTEM System = VtfReadStruct(input_vsf, detailed);
  VtfReadPBC(input_coor, &System.Box);
  if (!TriclinicCellData(&System.Box, 0)) {
    strcpy(ERROR_MSG, "wrong pbc data");
    PrintError();
    exit(1);
  } //}}}

  WarnChargedSystem(System, input_vsf, "\0", "\0");

  // read coordinates //{{{
  FILE *vcf = OpenFile(input_coor, "r");
  int count_vcf = 0,                         // count steps in the vcf file
      file_line_count = 0;                   // count lines in the vcf file
  char *stuff = calloc(LINE, sizeof *stuff); // array for the timestep preamble
  count_vcf++;
  if (!VtfReadTimestep(vcf, input_coor, &System, &file_line_count, stuff)) {
    count_vcf--;
  }
  fclose(vcf); //}}}

  // print information - verbose output //{{{
  if (verbose) {
    VerboseOutput(System);
  } //}}}

  COUNT *Count = &System.Count;

  double T_ref = 300; // reference temperature, K

  // identify the element //{{{
  int el = -1;
  for (int i = 0; i < 118; i++) {
    if (strcmp(species, element[i].Symbol) == 0) {
      el = i;
      break;
    }
  }
  if (el == -1) {
    strcpy(ERROR_MSG, "Unknown element!");
    printf("%sstring: %s%s%s\n", ErrRed(), ErrYellow(), species, ColourReset());
    PrintError();
    exit(1);
  } //}}}
  // error if missing data //{{{
  if (element[el].Mass == 0) {
    strcpy(ERROR_MSG, "Unspecified mass!");
    PrintError();
    printf("%selement %s%s%s\n", ErrRed(), ErrYellow(), species, ColourReset());
    exit(1);
  }
  if (element[el].Density == 0) {
    strcpy(ERROR_MSG, "Unspecified density!");
    PrintError();
    printf("%selement %s%s%s\n", ErrRed(), ErrYellow(), species, ColourReset());
    exit(1);
  }
  if (element[el].CohesionEnergy == 0) {
    strcpy(ERROR_MSG, "Unspecified cohesion energy!");
    PrintError();
    printf("%selement %s%s%s\n", ErrRed(), ErrYellow(), species, ColourReset());
    exit(1);
  }
  if (element[el].BulkModulus == 0) {
    strcpy(ERROR_MSG, "Unspecified bulk modulus!");
    PrintError();
    printf("%selement %s%s%s\n", ErrRed(), ErrYellow(), species, ColourReset());
    exit(1);
  } //}}}

  // element quantities from element array if not option-specified //{{{
  if (Mw == -1) {
    Mw = element[el].Mass;
  }
  if (rho == -1) {
    rho = element[el].Density;
  }
  if (E_coh == -1) {
    E_coh = element[el].CohesionEnergy;
  }
  if (B == -1) {
    B = element[el].BulkModulus;
  } //}}}
  // constants to recalculate stuff into different units
  double kT = k_B * T_ref, // to calculate in J
         eV = kT / e;      // to calculate in eV
  // system quantities
  double E = Count->BeadCoor / N_A * E_coh, // energy, J
         vol_a, // atomic volume, m^3
         vol; // volume of unit cell

  int n_unit = -1;
  if ((lattice[0] == '\0' && element[el].Crystal == 1) || fcc) {
    strcpy(lattice, "fcc");
    // 8 * 1/8 (corner atoms) +
    // 6 * 1/2 (face-centered atoms)
    n_unit = 4;
    vol_a = Mw / (N_A * rho);
    vol = n_unit * vol_a;
    if (a == -1) {
      a = pow(vol, 1.0 / 3) * 1e10; // lattice const, Å
    }
  } else if ((lattice[0] == '\0' && element[el].Crystal == 2) || bcc) {
    strcpy(lattice, "bcc");
    // 8 * 1/8 (corner atoms) +
    // 1 * 1   (body-centered atoms)
    n_unit = 2;
    vol_a = Mw / (N_A * rho);
    vol = n_unit * vol_a;
    if (a == -1) {
      a = pow(vol, 1.0 / 3) * 1e10; // lattice const, Å
    }
  } else if ((lattice[0] == '\0' && element[el].Crystal == 3) || hcp) {
    strcpy(lattice, "hcp");
    //  2 * 1/2 (face-centered atoms) +
    // 12 * 1/6 (corner atoms) +
    //  3 * 1   (body-centered atoms)
    n_unit = 6;
    vol_a = Mw / (N_A * rho);
    vol = n_unit * vol_a;
    /* volume of hcp cell
     * volume =
     * 3/2*sqrt(3)*a^2 (hexagon) *
     * sqrt(8/3)*a (height) =
     * 3*sqrt(2)a^3
     */
    if (a == -1) {
      a = pow(vol / (sqrt(2) * 3), 1.0 / 3) * 1e10; // lattice const, Å
    }
  }
  if (n_unit == -1) {
    strcpy(ERROR_MSG, "Wrong crystal structure! Should be fcc/bcc/hcp");
    PrintError();
    exit(1);
  }

//vol_a = vol / n_unit;

  // reduced quantities
  double vol_a_red = vol_a / vol;
  double E_red = E / (kT * Count->BeadCoor);
  double B_red = B * vol / kT;

  // unneeded reduced quantities //{{{
  // double vol_red = vol / a3;
  // double rho_red = Count->BeadCoor / vol_red;
  // double a_red = 1;
  // //}}}

  if (verbose) {
    printf("%sInput data for %s (%s):%s\n", Magenta(), species, lattice,
           ColourReset());
    printf("  Molar mass: %e kg/mol\n", Mw);
    printf("  Density: %e kg/m^3\n", rho);
    printf("  Cohesion energy: %e J/mol (%e eV)\n", E_coh, E_coh * Jmol_to_eV);
    printf("  Bulk modulus: %e Pa (%e eV/Å)\n\n", B, B / eVA3_to_Pa);
  }

  //// test prints //{{{
  // printf("%sCalculated quantities for %d atoms:%s\n",
  //        Magenta(), Count->BeadCoor, ColourReset());
  // printf("  Energy: %e J\n", U);
  // printf("  Volume: %e m^3\n", vol);
  // printf("  Lattice constant: %e m^3\n", a);
  // printf("  Atomic volume: %e m^3\n", vol_a);
  // printf("%sReduced quantities:%s\n", Magenta(), ColourReset());
  // printf("  Energy per atom: %e\n", U_red);
  // printf("  Volume: %e\n", vol_red);
  // printf("  Lattice constant: %e\n", a_red);
  // printf("  Atomic volume: %e\n", vol_a_red);
  // printf("  Density: %e\n", rho_red);
  // printf("  Bulk modulues: %e\n", B_red);
  //  //}}}

  // calculate n_SC as closest integer to 18*vol_a*B/(U*m) (in reduced units)
  double n_SC_dbl = -1;
  if (n_SC == -1) {
    n_SC_dbl = 18 * vol_a_red * B_red / (E_red * m_SC);
    n_SC = round(n_SC_dbl);
    if (n_SC <= m_SC) {
      printf("m isn't larger than n! (m=%d n=%d)\n", n_SC, m_SC);
      exit(0);
    }
  }
  // calculate how much the bulk modulus differs from input data
  double B_calc = n_SC * m_SC * E_red / (18 * vol_a_red);

  ToFractionalCoor(&System);
  // calculate lattice sums, S = sum_{i\neq j}(a/r_ij)^m_SC (and ^n_SC)
  count = 0;
  int values = 10; // TODO: make into option
  int n = Count->Bead / values;
  double avg_sum_m = 0, avg_sum_n = 0;
  for (int j = 0; j < Count->Bead; j += n) {
    count++;
    double sum_m = 0, sum_n = 0;
    int id1 = System.BeadCoor[j]; // atom to calculate distances from
    VECTOR *first = &System.Bead[id1].Position;
    for (int i = 0; i < Count->BeadCoor; i++) {
      int id2 = System.BeadCoor[i];
      if (id1 != id2) {
        VECTOR *pos = &System.Bead[id2].Position;
        VECTOR dist = Distance(*first, *pos, System.Box.Length);
        dist = FromFractional(dist, System.Box);
        double d = sqrt(SQR(dist.x) + SQR(dist.y) + SQR(dist.z));
        // printf("XXX %d dist: %lf (%lf %lf %lf) (%lf %lf %lf)\n",
        // id2, d, pos->x, pos->y, pos->z, dist.x, dist.y, dist.z);
        sum_m += 1 / pow(d, m_SC); // assumes reduced distance, a_red=1
        sum_n += 1 / pow(d, n_SC); //
      }
    }
    avg_sum_m += sum_m;
    avg_sum_n += sum_n;
  }
  avg_sum_m /= count;
  avg_sum_n /= count;

  //// test print sums //{{{
  // printf("%sLattice sums:%s\n", Magenta(), ColourReset());
  // printf("  sum_m (m=%d): %e\n", m_SC, sum_m);
  // printf("  sum_n (n=%d): %e\n", n_SC, sum_n);
  //  //}}}

  // calculate remaining parameters
  double c = n_SC * avg_sum_n / (m_SC * sqrt(avg_sum_m));
  double eps = 2 * m_SC * E_red / (avg_sum_n * (2 * n_SC - m_SC));

  printf("%sSutton-Chen parameters for %s (%s):%s\n", Magenta(), species,
         lattice, ColourReset());
  printf("   m: %d\n", m_SC);
  printf("   n: %d", n_SC);
  if (n_SC_dbl != -1) {
    printf(" (rounded from %lf)", n_SC_dbl);
  }
  putchar('\n');
  printf("   a: %.3f Å\n", a);
  printf("   c: %.3f\n", c);
  printf("   ε: %.3e K; %.3e eV\n", eps * T_ref, eps * eV);
  printf("   Ω: %.3f Å^3\n", vol_a*1e30);
  printf("   B: %.3e Pa; experiment %.3e Pa\n", B_calc * kT / vol, B);
  printf("  Mw: %.3f kg/m^3\n", Mw*1e3);

  // unneeded recalculated cohesion energy - same as input //{{{
  // double E_coh_calc = eps * sum_n * (2 * n_SC - m_SC) / (2 * m_SC);
  // printf("  coh E    %8.3f (%e J/mol; %e eV)\n",
  //       E_coh_calc, E_coh_calc*(k_B*T_ref*N_A), E_coh_calc*eV);
  // //}}}

  // free memory
  FreeSystem(&System);
  free(stuff);
  free(element);

  return 0;
}
