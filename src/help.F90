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
module moduleHelp
  
  implicit none

  public :: help, commandHelp
  private

contains

  subroutine help()
    use moduleVariables
    use moduleMessages
    implicit none
#ifdef GPTA_MPI
    integer :: idx
#endif

    call message(-10)

    call message(0,"GPTA HELP")
    call message(0,'For help about specific commands use the +help flag e.g.')
    call message(1,"gpta.x --gofr +help ")
    call message(0,'List of Available Commands in GPTA')
    
    call message(0,"%%%%% Global")
    call message(0," --nt          :: sets the number of openMP threads to use" )
    call message(0," --define      :: alters the value of some global parameters" )
    
    call message(0,"%%%%% Input/Output")
    call message(0," --i           :: reads the coordinates file(s)" )
    call message(0," --o           :: writes the coordinates file" )
    call message(0," --frame       :: selects one frame for processing" )
    call message(0," --frames      :: selects multiple frames for processing" )
    call message(0," --last        :: selects the last frame for processing" )
    call message(0," --skip        :: sets the stride for processing the frames" )
    call message(0," --log         :: sets the filename for the program output, it replaces the screen output" )

    call message(2)
    call message(0,"%%%%% Actions that change the atomic coordinates")
    call message(0," --pbc         :: applies PBC to the system")
    call message(0," --unwrap      :: removes the PBC wrapping in a trajectory")
    call message(0," --shift       :: translated the atoms")
    call message(0," --rescale     :: rescales the cell")
    call message(0," --fixcom      :: translated the centre of mass")
    call message(0," --repl        :: replicates the cell")
    call message(0," --mirror      :: adds the mirror image of the system")
    call message(0," --noclash     :: removes overlapping molecules")
    call message(0," --delete      :: deletes atoms/molecules in the system")
    call message(0," --set         :: changes the properties of the atoms")
    call message(0," --subs        :: replace molecules in the system")
    call message(0," --add         :: adds atoms/molecules to the system")
    call message(0," --surface     :: creates a surface of a crystal")
    call message(0," --cluster     :: extract spherical clusters from the coordinates")
    
    call message(2)
    call message(0,"%%%%% Actions that calculate properties ")
    call message(0," --top         :: calculates of the topology")
    call message(0," --extract     :: computes system-wide properties from the trajectory file")
    call message(0," --gofr        :: computes the radial pair distribution function")
    call message(0," --dmap1D      :: computes the 1D density map")
    call message(0," --dmap2D      :: computes the 2D density map")
    call message(0," --dmap3D      :: computes the 3D density map")
    call message(0," --solvation   :: computes the solvent density map around the solute")
    call message(0," --molprop     :: computes molecular properties")
    call message(0," --restime     :: computes the solvent residence time")
    call message(0," --xray        :: computes the x-ray powder diffraction spectrum")
    call message(0," --msd         :: computes the mean square displacement")

    call message(2)
    call message(0,"%%%%% Special actions ")
#ifdef GPTA_PLUMED
    call message(0, " --plumed      :: conputes properties of the system using PLUMED")
#endif

! openMM driver for AMOEBA
#ifdef GPTA_OPENMM
    call message(0, " --amoeba      :: computes the energy/induced dipole... using openMM")
#endif
    call message(1," --test        :: test routine for develpment")

#ifdef GPTA_MPI
    call MPI_Abort(MPI_COMM_WORLD, 9, idx )
#else
    stop
#endif

  end subroutine help

  subroutine commandHelp(cmd)
    use moduleMessages

    use moduleTemplateAction

    use moduleOpenMM
    use moduleAtomsAttributes
    use moduleDeleteAtoms
    use moduleMolecularTopology
    use moduleReplaceMolecule
    use moduleDumpCoordinates
    use moduleCreateSurface

    use moduleExtractSystemProperties
    use moduleRDF
    use moduleModifyCoordinates
    use moduleDensityProfile
    use moduleDensityMap2D
    use moduleDensityMap3D
    use moduleMolecularProperties
    use moduleAddAtoms
    use moduleResidenceTime
    use moduleXrayAction
    use moduleMeanSquareDisplacement
    use moduleExtractClustersAction
    use moduleSolvationShell

#ifdef GPTA_OPENMM
    use moduleAmoeba
#endif

#ifdef GPTA_PLUMED
    use modulePlumedInterface
#endif

    implicit none
    character(len=*), intent(in) :: cmd
    call message(-10)

    call message(0,"GPTA HELP for action "//trim(cmd))
    call message(0,"For a more detailed description of the command please refer to the manual.")
    call message(0,"..........................................................................................................")

    if (cmd == "--i" ) then
      call message(0,"This action reads the atomic coordinate from file(s).")
      call message(0,"It takes one or more filenames as argument; it can also be repeated.")
      call message(0,"There is also an alternative form, --iABC, where ABC is the typical file extension identifying the coordinates' type")
      call message(0,"These are the currently supported input coordinates file types")
      call message(0,"     xyz -> Extended XYZ")
      call message(0,"     pdb -> Protein Data Bank file")
      call message(0,"     dcd -> CHARMM DCD file ")
      call message(0,"    dcd2 -> CHARMM DCD file (safe/slow routine)")
      call message(0,"     arc -> TINKER ARC file (with connectivity)")
      call message(0,"    txyz -> TINKER XYZ input file")
      call message(0,"     gin -> GULP input file")
      call message(0,"     gau -> GAUSSIAN output file")
      call message(0,"     gro -> GROMOS input file")
      call message(0,"     lmp -> LAMMPS input file")
      call message(0,"  lmptrj -> LAMMPS custom trajectory file (partially supported)")
      call message(0,"     cif -> Crystallographic Information File ")
      call message(0,"     xtc -> GROMACS XTC trajectory file")
      call message(0,"     trr -> GROMACS TRR trajectory file")
      call message(0,"Examples:")
      call message(0,"  gpta.x --i coord.pdb trajectory.dcd")
      call message(0,"  gpta.x --ixyz coord.data --i traj.pdb")
    end if

    if (cmd == "--o" ) then
      call message(0,"This action write the atomic coordinates to file.")
      call message(0,"It takes one filenames as argument, it can be repeated.")
      call message(0,"There is also an alternative form, --oABC, where ABC is the typical file extension identifying the coordinates' type")
      call message(0,"These are the currently supported input coordinates file types")
      call message(0,"     xyz -> Extended XYZ")
      call message(0,"     pdb -> Protein Data Bank file")
      call message(0,"     pdb -> Protein Data Bank file (with CONECT section)")
      call message(0,"     dcd -> CHARMM DCD file ")
      call message(0,"     lmp -> LAMMPS input file")
      call message(0,"  lmptrj -> LAMMPS custom trajectory file (partially supported)")
      call message(0,"     psf -> CHARMM PSF file")
      call message(0,"     arc -> TINKER ARC file (with connectivity)")
      call message(0,"     gin -> GULP input file")
      call message(0,"    gin2 -> GULP input file (fractional coordinates)")
      call message(0,"     xtc -> GROMACS XTC trajectory file")
      call message(0,"Examples:")
      call message(0,"  gpta.x --i coord.pdb trajectory.dcd")
      call message(0,"  gpta.x --ixyz coord.data --i traj.pdb")
    end if

    if (cmd == "--test" ) then
      call message(0,"This action is used for development only and has custom flags.")
      call message(0,"Examples:")
      call message(0,"  gpta.x --i coord.pdb --test")
    end if

    if (cmd == "--top" ) then
      call message(0,"This action computes the system's topology and groups the atoms in molecules.")
      call message(0,"It is a prerequisity to many actions.")
      call message(0,"Examples:")
      call message(0,"  gpta.x --i coord.pdb --top")
    end if

    if (cmd == "--define"   ) call defineVariablesHelp()

    if (cmd == "--rescale"  ) call rescaleCellHelp()
    if (cmd == "--shift"    ) call shiftCoordinatesHelp()
    if (cmd == "--fixcom"   ) call shiftCOMHelp()
    if (cmd == "--pbc"      ) call applyPeriodicboundaryConditionsHelp()
    if (cmd == "--unwrap"   ) call unwrapCoordinatesHelp()
    if (cmd == "--repl"     ) call replicateCellHelp()
    if (cmd == "--mirror"   ) call mirrorCellHelp()
    if (cmd == "--noclash"  ) call removeOverlappingMoleculesHelp()
    
    if (cmd == "--delete"   ) call deleteAtomsHelp()
    if (cmd == "--set"      ) call setAtomAttributesHelp()
    if (cmd == "--subs"     ) call replaceMoleculesHelp()
    if (cmd == "--add"      ) call addAtomsHelp()
    if (cmd == "--surface"  ) call createSurfaceHelp()
    if (cmd == "--cluster"  ) call extractClustersHelp()

    if (cmd == "--extract"  ) call extractSystemPropertiesHelp()
    if (cmd == "--gofr"     ) call computeRadialPairDistributionHelp()
    if (cmd == "--dmap1D"   ) call computeDensityProfileHelp()
    if (cmd == "--dmap2D"   ) call computeDensityMap2DHelp()
    if (cmd == "--dmap3D"   ) call computeDensityMap3DHelp()
    if (cmd == "--solvation") call solvationShellHelp()
    if (cmd == "--molprop"  ) call computeMolecularPropertiesHelp()
    if (cmd == "--restime"  ) call computeResidenceTimeHelp()
    if (cmd == "--xray"     ) call computeXrayPowderHelp()
    if (cmd == "--msd"      ) call computeMSDHelp()

! openMM driver for AMOEBA
#ifdef GPTA_OPENMM
    if (cmd == "--amoeba" ) call amoebaHelp()
#endif

#ifdef GPTA_PLUMED
    ! if (cmd == "--plumed" ) call plumedInterfaceHelp()
#endif

    call message(-3)

  end subroutine

end module moduleHelp
