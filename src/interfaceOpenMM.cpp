! Copyright (c) 2021, Paolo Raiteri, Curtin University.
! All rights reserved.
! 
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License as published by the 
! Free Software Foundation; either version 3 of the License, or 
! (at your option) any later version.
!  
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! * Redistributions of source code must retain the above copyright notice, 
!   this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright notice, 
!   this list of conditions and the following disclaimer in the documentation 
!   and/or other materials provided with the distribution.
! * Neither the name of the <ORGANIZATION> nor the names of its contributors 
!   may be used to endorse or promote products derived from this software 
!   without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
#include "header.h"
#include <iomanip>
#include <string>
#include <sstream>

#ifdef GPTA_CUDA
extern "C" void registerAmoebaCudaKernelFactories();
#else
extern "C" void registerAmoebaReferenceKernelFactories();
#endif

extern "C" void openmm_init_()
{
#ifdef GPTA_CUDA
  registerAmoebaCudaKernelFactories();
  Platform::getPlatformByName("CUDA").setPropertyDefaultValue("Precision", "mixed");
#else
  registerAmoebaReferenceKernelFactories();
#endif
}

// Initialiase calculation
extern "C" void openmm_create_atoms_(int *n, double *mass)
{
  // Initialise particles
  numberOfAtoms = n[0];
  for (int i = 0; i < numberOfAtoms; i++)
    mySystem.addParticle(mass[i]);
}

extern "C" void openmm_setup_multipoles_(double *cutoff, double *ewaldPrecision, char *pType)
{
  amoebaMultipoleForce = new AmoebaMultipoleForce();
  mySystem.addForce(amoebaMultipoleForce);

  amoebaMultipoleForce->setNonbondedMethod(AmoebaMultipoleForce::PME);
  amoebaMultipoleForce->setEwaldErrorTolerance(*ewaldPrecision);

  string str;
  for (int i = 0; i < strlen(pType); i++)
    str.append(1, pType[i]);
  if (str == "mutual")
    amoebaMultipoleForce->setPolarizationType(AmoebaMultipoleForce::Mutual);
  else
    amoebaMultipoleForce->setPolarizationType(AmoebaMultipoleForce::Direct);
  amoebaMultipoleForce->setMutualInducedMaxIterations(1000);
  amoebaMultipoleForce->setMutualInducedTargetEpsilon(1e-6);

  // lame check that the input isn't already in nm
  if (*cutoff > 3.0)
    *cutoff /= 10;
  amoebaMultipoleForce->setCutoffDistance(*cutoff);

}

extern "C" void openmm_add_amoeba_multipoles_to_atom_(
    int *iType, double *multipole, double *pol, double *tho, int *axisType, int *refFrame, int *polarisationGroup)
{
  Charge[*iType] = multipole[0];
  Dipole[*iType][0] = multipole[1];
  Dipole[*iType][1] = multipole[2];
  Dipole[*iType][2] = multipole[3];

  Quadrupole[*iType][0] = multipole[4];
  Quadrupole[*iType][1] = multipole[5];
  Quadrupole[*iType][2] = multipole[6];
  Quadrupole[*iType][3] = multipole[7];
  Quadrupole[*iType][4] = multipole[8];
  Quadrupole[*iType][5] = multipole[9];
  Quadrupole[*iType][6] = multipole[10];
  Quadrupole[*iType][7] = multipole[11];
  Quadrupole[*iType][8] = multipole[12];

  Polarisation[*iType] = *pol;
  Thole[*iType] = *tho;

  int iatm = amoebaMultipoleForce->addMultipole(Charge[*iType],
                                                Dipole[*iType],
                                                Quadrupole[*iType],
                                                *axisType,
                                                refFrame[0] - 1, refFrame[1] - 1, refFrame[2] - 1,
                                                Thole[*iType],
                                                pow(Polarisation[*iType], 1.0 / 6.0), Polarisation[*iType]);

  int n = polarisationGroup[0];
  vector<int> polar(n);
  for (int i = 0; i < n; i++)
    polar[i] = polarisationGroup[i + 1] - 1;

  amoebaMultipoleForce->setCovalentMap(iatm, AmoebaMultipoleForce::PolarizationCovalent11, polar);
}

extern "C" void openmm_add_amoeba_vdw_(int *n, int *ivIndex, double *params,
                                       int *nShape, int *nCovalent12, int *lCovalent12,
                                       int *nCovalent13, int *lCovalent13)
{
  AmoebaVdwForce *amoebaVdwForce = new AmoebaVdwForce();

  std::string sigmaCombiningRule = std::string("CUBIC-MEAN");
  amoebaVdwForce->setSigmaCombiningRule(sigmaCombiningRule);

  std::string epsilonCombiningRule = std::string("HHG");
  amoebaVdwForce->setEpsilonCombiningRule(epsilonCombiningRule);

  amoebaVdwForce->setCutoff(1.);
  amoebaVdwForce->setNonbondedMethod(AmoebaVdwForce::CutoffPeriodic);

  int idx = 0;
  for (int i = 0; i < n[0]; i++)
  {
    int iv = ivIndex[i] - 1;
    double sigma = params[idx++] * 0.5; // radiussize DIAMETER
    double epsilon = params[idx++];
    double reduction = params[idx++];
    int idx = amoebaVdwForce->addParticle(iv, sigma, epsilon, reduction);
  }

  for (int i = 0; i < numberOfAtoms; i++)
  {
    std::vector<int> exclusions;
    exclusions.push_back(i);
    for (int j = 0; j < nCovalent12[i]; j++)
      exclusions.push_back(lCovalent12[i * nShape[0] + j] - 1);
    for (int j = 0; j < nCovalent13[i]; j++)
      exclusions.push_back(lCovalent13[i * nShape[0] + j] - 1);

    amoebaVdwForce->setParticleExclusions(i, exclusions);
  }

  // amoebaVdwForce->setUseDispersionCorrection(1);

  mySystem.addForce(amoebaVdwForce);
}

extern "C" void openmm_add_custom_non_bonded_(int *n, int *alist)
{

  string expression = string("A*exp(-r/B); A=0.5*(A1+A2); B=0.5*(B1+B2)");
  CustomNonbondedForce *buckingham = new CustomNonbondedForce(expression);

  buckingham->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
  buckingham->setCutoffDistance(1.);

  buckingham->addPerParticleParameter("A");
  buckingham->addPerParticleParameter("B");

  vector<double> parameters = {15437.6, .03};

  set<int> group1;
  set<int> group2;

  for (int i = 0; i < *n; i++)
    buckingham->addParticle(parameters);

  for (int i = 0; i < *n; i++)
    if (alist[i] == 1)
      group1.insert(i);
    else if (alist[i] == 2)
      group2.insert(i);

  buckingham->addInteractionGroup(group1, group2);

  mySystem.addForce(buckingham);
}

extern "C" void openmm_add_amoeba_bonds_(int *n, int *indices, double *params)
{
  AmoebaBondForce *amoebaBondForce = new AmoebaBondForce();
  amoebaBondForce->setAmoebaGlobalBondCubic(-25.5);
  amoebaBondForce->setAmoebaGlobalBondQuartic(379.3125);
  int idx = 0;
  int jdx = 0;
  for (int i = 0; i < n[0]; i++)
  {
    int iatm = indices[idx++] - 1;
    int jatm = indices[idx++] - 1;
    double distance = params[jdx++];
    double kappa = params[jdx++];
    amoebaBondForce->addBond(iatm, jatm, distance, kappa);
  }
  mySystem.addForce(amoebaBondForce);
}

extern "C" void openmm_add_amoeba_angles_(int *n, int *indices, double *params)
{
  AmoebaAngleForce *amoebaAngleForce = new AmoebaAngleForce();
  amoebaAngleForce->setAmoebaGlobalAngleCubic(-0.014);
  amoebaAngleForce->setAmoebaGlobalAngleQuartic(5.6e-05);
  amoebaAngleForce->setAmoebaGlobalAnglePentic(-7e-07);
  amoebaAngleForce->setAmoebaGlobalAngleSextic(2.2e-08);

  int idx = 0;
  int jdx = 0;
  for (int i = 0; i < n[0]; i++)
  {
    int iatm = indices[idx++] - 1;
    int jatm = indices[idx++] - 1;
    int katm = indices[idx++] - 1;
    double theta = params[jdx++];
    double kappa = params[jdx++];

    amoebaAngleForce->addAngle(iatm, jatm, katm, theta, kappa);
  }
  mySystem.addForce(amoebaAngleForce);
}

extern "C" void openmm_add_amoeba_torsions_(int *n, int *indices, double *params)
{
  PeriodicTorsionForce *periodicTorsionForce = new PeriodicTorsionForce();
  int idx = 0;
  int jdx = 0;
  for (int i = 0; i < n[0]; i++)
  {
    int iatm = indices[idx++] - 1;
    int jatm = indices[idx++] - 1;
    int katm = indices[idx++] - 1;
    int latm = indices[idx++] - 1;
    double phi1 = params[jdx++];
    double phi2 = params[jdx++];
    double phi3 = params[jdx++];
    double kappa1 = params[jdx++];
    double kappa2 = params[jdx++];
    double kappa3 = params[jdx++];

    if (kappa1 > 1e-6)
      periodicTorsionForce->addTorsion(iatm, jatm, katm, latm, 1, phi1, kappa1);
    if (kappa2 > 1e-6)
      periodicTorsionForce->addTorsion(iatm, jatm, katm, latm, 2, phi2, kappa2);
    if (kappa3 > 1e-6)
      periodicTorsionForce->addTorsion(iatm, jatm, katm, latm, 3, phi3, kappa3);
  }
  mySystem.addForce(periodicTorsionForce);
}

extern "C" void openmm_add_amoeba_out_of_plane_(int *n, int *list, double *kappa)
{
  AmoebaOutOfPlaneBendForce *amoebaOutOfPlaneBendForce = new AmoebaOutOfPlaneBendForce();
  amoebaOutOfPlaneBendForce->setAmoebaGlobalOutOfPlaneBendCubic(-0.1400000E-01);
  amoebaOutOfPlaneBendForce->setAmoebaGlobalOutOfPlaneBendQuartic(0.5600000E-04);
  amoebaOutOfPlaneBendForce->setAmoebaGlobalOutOfPlaneBendPentic(-0.7000000E-06);
  amoebaOutOfPlaneBendForce->setAmoebaGlobalOutOfPlaneBendSextic(0.2200000E-07);

  // double kOutOfPlaneBend = 0.30843384;
  // amoebaOutOfPlaneBendForce->addOutOfPlaneBend(0, 1, 2, 3, kOutOfPlaneBend);
  // The central atom should be the second in the list
  // Usual ambiguity for the order of the other three
  vector<int> id(4);
  int idx = 0;
  for (int i = 0; i < *n; i++)
  {
    id[1] = list[idx++] - 1;
    id[0] = list[idx++] - 1;
    id[2] = list[idx++] - 1;
    id[3] = list[idx++] - 1;
    amoebaOutOfPlaneBendForce->addOutOfPlaneBend(id[0], id[1], id[2], id[3], kappa[i]);
  }

  mySystem.addForce(amoebaOutOfPlaneBendForce);
}

extern "C" void openmm_add_amoeba_stretchbend_(int *n, int *indices, double *params)
{
  AmoebaStretchBendForce *amoebaStretchBendForce = new AmoebaStretchBendForce();
  int idx = 0;
  int jdx = 0;
  for (int i = 0; i < n[0]; i++)
  {
    int iatm = indices[idx++] - 1;
    int jatm = indices[idx++] - 1;
    int katm = indices[idx++] - 1;
    double abl = params[jdx++];
    double cbl = params[jdx++];
    double theta = params[jdx++];
    double k1 = params[jdx++];
    double k2 = params[jdx++];

    if (abs(k1) > 1e-6)
      amoebaStretchBendForce->addStretchBend(iatm, jatm, katm, abl, cbl, theta, k1, k2);
  }
  mySystem.addForce(amoebaStretchBendForce);
}

extern "C" void openmm_add_amoeba_ureybradley_(int *n, int *indices, double *params)
{
  HarmonicBondForce *harmonicBondForce = new HarmonicBondForce();
  int idx = 0;
  int jdx = 0;
  for (int i = 0; i < n[0]; i++)
  {
    int iatm = indices[idx++] - 1;
    int jatm = indices[idx++] - 1;
    int katm = indices[idx++] - 1;
    double distance = params[jdx++];
    double kappa = params[jdx++];

    if (abs(kappa) > 1e-6)
      harmonicBondForce->addBond(iatm, katm, distance, 2 * kappa);
  }
  mySystem.addForce(harmonicBondForce);
}

extern "C" void openmm_setup_exclusions_(int *nShape,
                                         int *nCovalent12, int *lCovalent12,
                                         int *nCovalent13, int *lCovalent13,
                                         int *nCovalent14, int *lCovalent14)
{
  // cout << "Setting 1-2 interactions\n";
  for (int i = 0; i < numberOfAtoms; i++)
  {
    vector<int> coval12(nCovalent12[i]);
    for (int j = 0; j < nCovalent12[i]; j++)
      coval12[j] = lCovalent12[i * nShape[0] + j] - 1;
    amoebaMultipoleForce->setCovalentMap(i, AmoebaMultipoleForce::Covalent12, coval12);
  }

  // cout << "Setting 1-3 interactions\n";
  for (int i = 0; i < numberOfAtoms; i++)
  {
    vector<int> coval13(nCovalent13[i]);
    for (int j = 0; j < nCovalent13[i]; j++)
      coval13[j] = lCovalent13[i * nShape[0] + j] - 1;
    amoebaMultipoleForce->setCovalentMap(i, AmoebaMultipoleForce::Covalent13, coval13);
  }

  // cout << "Setting 1-4 interactions\n";
  for (int i = 0; i < numberOfAtoms; i++)
  {
    vector<int> coval14(nCovalent14[i]);
    for (int j = 0; j < nCovalent14[i]; j++)
      coval14[j] = lCovalent14[i * nShape[0] + j] - 1;
    amoebaMultipoleForce->setCovalentMap(i, AmoebaMultipoleForce::Covalent14, coval14);
  }
}

extern "C" void openmm_compute_multipoles_(int *frameNumber, double *box, double *x, double *y, double *z, int *iType, char *dumpFile)
{
  mySystem.setDefaultPeriodicBoxVectors(
      Vec3(box[0], box[1], box[2]),
      Vec3(box[3], box[4], box[5]),
      Vec3(box[6], box[7], box[8]));

  LangevinIntegrator integrator(0.0, 0.1, 0.000001);
#ifdef GPTA_CUDA
  Context context(mySystem, integrator, Platform::getPlatformByName("CUDA"));
#else
  Context context(mySystem, integrator, Platform::getPlatformByName("Reference"));
#endif

  // Define PME parameters
  double actualAlpha;
  int actualSize[3];
  amoebaMultipoleForce->getPMEParametersInContext(context, actualAlpha, actualSize[0], actualSize[1], actualSize[2]);

  // Put atoms into context
  vector<Vec3> positions(numberOfAtoms);
  for (int i = 0; i < numberOfAtoms; i++)
    positions[i] = Vec3(x[i], y[i], z[i]);
  context.setPositions(positions);

  // define molecules
  molecules = context.getMolecules();

  // Compute energy
  myState = context.getState(State::Energy | State::Forces | State::Positions);

  // extract induced dipoles
  amoebaMultipoleForce->getInducedDipoles(context, inducedDipoles);

  // extract Total dipoles
  amoebaMultipoleForce->getTotalDipoles(context, totalDipole);

  string str;
  for (int i = 0; i < strlen(dumpFile); i++)
    str.append(1, dumpFile[i]);

  if (str != "NULL")
  {
    ofstream myfile;
    myfile.open(str);

    myfile << "# PME alpha " << actualAlpha << endl;
    myfile << "# PME Kx    " << actualSize[0] << endl;
    myfile << "# PME Ky    " << actualSize[1] << endl;
    myfile << "# PME Kx    " << actualSize[2] << endl;

    myfile << "# Number of molecules : " << molecules.size() << endl;

    double energy = myState.getPotentialEnergy();
    myfile << "# Electrostatic energy [kJ/mol] = " << energy << endl;

    for (int i = 0; i < numberOfAtoms; i++)
      positions[i] = myState.getPositions()[i];

    myfile << "# Positions | Atomic Polarisation | Induced Dipoles | Electric Field\n";
    for (int i = 0; i + 1 < numberOfAtoms; i++)
      if (Polarisation[iType[i]] > 0.0000)
      {
        myfile << setw(6) << i + 1 << " "
               << setw(15) << std::fixed << setprecision(7) << positions[i][0] * 10 << " "
               << setw(15) << std::fixed << setprecision(7) << positions[i][1] * 10 << " "
               << setw(15) << std::fixed << setprecision(7) << positions[i][2] * 10 << " "
               << setw(15) << std::fixed << setprecision(7) << Polarisation[iType[i]] << " "

               << setw(15) << std::fixed << setprecision(7) << 4.80321 * inducedDipoles[i][0] * 10 << " "
               << setw(15) << std::fixed << setprecision(7) << 4.80321 * inducedDipoles[i][1] * 10 << " "
               << setw(15) << std::fixed << setprecision(7) << 4.80321 * inducedDipoles[i][2] * 10 << " "

               << setw(15) << std::fixed << setprecision(7) << 0.1439 * inducedDipoles[i][0] / Polarisation[iType[i]] << " "
               << setw(15) << std::fixed << setprecision(7) << 0.1439 * inducedDipoles[i][1] / Polarisation[iType[i]] << " "
               << setw(15) << std::fixed << setprecision(7) << 0.1439 * inducedDipoles[i][2] / Polarisation[iType[i]] << " "
               << endl;
      }
    myfile.close();
  }
}

extern "C" void openmm_compute_energy_(double *energy)
{
  energy[0] = myState.getPotentialEnergy();
  // cout << " Electrostatic energy [kJ/mol] = " << energy[0] << endl;
}

extern "C" void openmm_compute_electric_field_(int *iType, double *eField, int *idir, double *z)
{
  int icart = *idir - 1;
  for (int i = 0; i < numberOfAtoms; i++)
    if (Polarisation[iType[i]] > 1e-8)
      eField[i] = inducedDipoles[i][icart] / Polarisation[iType[i]]; // units [e*nm / nm^3] <- z component only
    else
      eField[i] = 0.;

  // for (int i = 0; i < numberOfAtoms; i++)
  //   cout << z[i] << " " << eField[i] << endl;
}

extern "C" void openmm_compute_dipoles_(int *iType, double *moleculeDipole, double *x, double *y, double *z)
{
  int idx;
  int imol;
  int iatm = 0;
  Vec3 dij;
  Vec3 dipoleVector;
  
  for (imol = 0; imol < molecules.size(); imol++)
  {
    dipoleVector = totalDipole[iatm];
    for (int i = 1; i < molecules[imol].size(); i++)
    {
      dipoleVector += totalDipole[iatm + i];
      dij[0] = x[iatm + i] - x[iatm];
      dij[1] = y[iatm + i] - y[iatm];
      dij[2] = z[iatm + i] - z[iatm];
      dipoleVector[0] += Charge[iType[iatm + i]] * dij[0];
      dipoleVector[1] += Charge[iType[iatm + i]] * dij[1];
      dipoleVector[2] += Charge[iType[iatm + i]] * dij[2];
    }
//    dipoleVector *= 48.032047; // conversion from e*nm to Debye
    iatm += molecules[imol].size();
 
    idx = imol * 3;
    moleculeDipole[idx] = dipoleVector[0];
    moleculeDipole[idx + 1] = dipoleVector[1];
    moleculeDipole[idx + 2] = dipoleVector[2];
  }
}
