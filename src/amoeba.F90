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
module moduleAmoeba
  use moduleVariables
  use moduleFiles
  use moduleProperties

  implicit none

  public :: amoeba
  private

#ifdef GPTA_OPENMM 
  character(:), pointer :: actionCommand
  type(fileTypeDef), pointer :: outputFile
  character(len=STRLEN), pointer :: forcefieldFile
  integer, pointer :: tallyExecutions

  integer, pointer :: numberOfBins
  integer, pointer :: numberOfBinsZ

  real(8), allocatable, dimension(:,:) :: localPositions

  integer, pointer                 :: nProperties
  integer, pointer, dimension(:)   :: ID
#endif

contains

  subroutine initialiseAction(a)
#ifdef GPTA_OPENMM 
    use moduleStrings
    use moduleSystem
    use moduleOpenMM, only : polarisationType, cutoffOpenMM, ewaldPrecision
#endif
    implicit none
    type(actionTypeDef), target :: a

#ifdef GPTA_OPENMM       
    logical :: lflag

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0d0

    ! Local pointers
    actionCommand          => a % actionDetails
    tallyExecutions        => a % tallyExecutions
    outputFile             => a % outputFile
    forcefieldFile         => a % stringVariables(1)
    numberOfBins           => a % numberOfBins
    numberOfBinsZ          => a % integerVariables(1)
    
    nProperties            => a % integerVariables(2)
    ID(1:)                 => a % integerVariables(3:)

    call assignFlagValue(actionCommand,"+out", outputFile % fname, 'NULL')
    call assignFlagValue(actionCommand,"+f", forcefieldFile, 'water09.prm')
    
    call assignFlagValue(actionCommand,"+nbin", numberOfBins,100)
    call assignFlagValue(actionCommand,"+nbinZ", numberOfBinsZ,0)
    
    call assignFlagValue(actionCommand,"+rcut", cutoffOpenMM, 10.d0)
    call assignFlagValue(actionCommand,"+eps", ewaldPrecision, 1d-5)
    call assignFlagValue(actionCommand,"+direct", lflag, .false.)
    if (lflag) polarisationType = "direct"
    
    call openmm_init()
    
    ! allocate(efieldDistribution(numberOfBins) , source=0.d0)
    ! allocate(atomsdDistribution(numberOfBins) , source=0)
    
    ID = -1
    tallyExecutions = 0 
#else
    call message(-1,"Recompile with the interface to openMM")
#endif
  end subroutine initialiseAction

#ifdef GPTA_OPENMM 
  subroutine dumpScreenInfo(iFlag)
    use moduleMessages 
    use moduleOpenMM, only : polarisationType, cutoffOpenMM, ewaldPrecision
    implicit none
    integer, intent(in) :: iFlag

    call message(0,"Calling openMM ")
    call message(0,"...Forcefield file",str=forcefieldFile)
    call message(0,"...Real space cutoff",r=cutoffOpenMM)
    call message(0,"...PME precision",r=ewaldPrecision)
    call message(0,"...Polarisation type",str=polarisationType)
    
    if (iFlag == 0) then
      if (outputFile % fname == "NULL") outputFile % fname = "gpta_openmm.out"
      call message(0,"...Output file",str=outputFile % fname)
    else
      call message(0,"...Performing Full Analysis")
      call message(0,"......Number of bins for the distributions",i=numberOfBins)
      call message(0,"......Number of bins for the distributions in z",i=numberOfBinsZ)
    end if
    call message(2)
  end subroutine dumpScreenInfo
#endif


  subroutine amoeba(a)
    use moduleVariables
    use moduleSystem 
    use moduleOpenMM
    use moduleDistances
    use moduleMessages
    implicit none
    type(actionTypeDef), target :: a
  
#ifdef GPTA_OPENMM       
    real(8) :: energy
    character(len=100) :: str

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1
      
      if (firstAction) then
        
        if (numberOfMolecules == 0) call runInternalAction("topology","NULL")
        call dumpScreenInfo(0)

        ! Initilise OpenMM
        call initialiseOpenMM(forcefieldFile)

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      end if

      call computeAmoebaMutipoles(outputFile % fname)
      
      call openmm_compute_energy(energy)
      write(str,'("--> Frame ",i6," | Energy ",f12.2,1x,f12.5)') frame % nframe, energy, energy/96.485d0
      call message(0,str)

    end if

#endif    

  end subroutine amoeba

end module moduleAmoeba
