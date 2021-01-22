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
#include "openmm/Platform.h"
#include "openmm/Vec3.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/internal/AssertionUtilities.h"
#include "OpenMMAmoeba.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/CustomNonbondedForce.h"

// #include "AmoebaCudaKernelFactory.h"
// #include "AmoebaCudaKernels.h"
// #include "CudaPlatform.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string.h>

using namespace OpenMM;
using namespace std;

System mySystem;
State myState;
AmoebaMultipoleForce *amoebaMultipoleForce;

vector< vector<int> > molecules;
vector<Vec3> inducedDipoles;
vector<Vec3> totalDipole;

int numberOfAtoms;

double o_charge, o_polarisation;
vector<double> o_dipole(3), o_quadrupole(9, 0.0);

double h_charge, h_polarisation;
vector<double> h_dipole(3), h_quadrupole(9, 0.0);

// DCE multipoles
double cd_charge, cd_polarisation;
vector<double> cd_dipole(3, 0.), cd_quadrupole(9, 0.0);

double hd_charge, hd_polarisation;
vector<double> hd_dipole(3, 0.), hd_quadrupole(9, 0.0);

double cl_charge, cl_polarisation;
vector<double> cl_dipole(3, 0.), cl_quadrupole(9, 0.0);

vector<double> Charge(200, 0.);
vector< vector<double> > Dipole(200, vector<double>(3));
vector< vector<double> > Quadrupole(200, vector<double>(9));
vector<double> Polarisation(200);
vector<double> Thole(200);

void define_multipoles();
static void testMultipoleWaterPMEDirectPolarization();
