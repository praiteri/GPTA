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
module modulePlumedInterface
  use moduleVariables
  implicit none

  public :: plumedInterface
  private

#ifdef GPTA_PLUMED
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  character(len=STRLEN), pointer :: plumedInputFile
  character(len=STRLEN), pointer :: plumedOutputFile

  integer, pointer :: tallyExecutions

  real(8) :: energyUnits, lengthUnits, timeUnits
#endif

contains

  subroutine initialiseAction(a)
    use moduleStrings
    use moduleSystem
    implicit none
    type(actionTypeDef), target :: a

#ifdef GPTA_PLUMED
    a % actionInitialisation = .false.

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    plumedInputFile      => a % stringVariables(1)
    plumedOutputFile     => a % stringVariables(2)

    call assignFlagValue(actionCommand,"+f",plumedInputFile,'NULL')
    call assignFlagValue(actionCommand,"+out",plumedOutputFile,'plumed_gpta.out')

    call dumpScreenInfo()
    tallyExecutions = 0 
#else
    call message(-1,"Recompile with the interface to PLUMED 2.x")
#endif
  end subroutine initialiseAction

#ifdef GPTA_PLUMED
  subroutine dumpScreenInfo()
    use moduleMessages 
    implicit none
    call message(0,"Creating PLUMED interface")
  end subroutine dumpScreenInfo
#endif

  subroutine plumedInterface(a)
    use moduleVariables
    use moduleElements
    use moduleSystem
    implicit none
    type(actionTypeDef), target :: a

#ifdef GPTA_PLUMED

    integer :: plumedVariable
    real(8), pointer :: mass(:), force(:,:), virial(:,:)

    integer :: i
    real(8) :: etot

    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      if (firstAction) then
        energyUnits = 96.48530749925792d0
        lengthUnits = 0.1d0
        timeUnits = 1.0d0

        call plumed_f_installed(plumedVariable)
        if(plumedVariable<=0) call message(-1,"PLUMED 2 not available")

        if (plumedInputFile=="NULL") call message(-1,"Plumed input file must be specified with +f")
        
!        call message(1,"creating plumed from GPTA")
        call plumed_f_gcreate()
        call plumed_f_gcmd("setRealPrecision"//char(0),8)
        call plumed_f_gcmd("setMDEnergyUnits"//char(0),energyUnits)
        call plumed_f_gcmd("setMDLengthUnits"//char(0),lengthUnits)
        call plumed_f_gcmd("setMDTimeUnits"//char(0),timeUnits)
        call plumed_f_gcmd("setNatoms"//char(0),frame % natoms)
        call plumed_f_gcmd("setPlumedDat"//char(0),trim(plumedInputFile)//char(0))
        call plumed_f_gcmd("setLogFile"//char(0),trim(plumedOutputFile)//char(0))
        call plumed_f_gcmd("setMDEngine"//char(0),"driver");
        call plumed_f_gcmd("setTimestep"//char(0),1.0d0);
        call plumed_f_gcmd("init"//char(0),0);

        call dumpScreenInfo()

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      end if

      allocate(mass(frame % natoms))
      allocate(force(3,frame % natoms), source=0.d0)
      allocate(virial(3,3), source=0.d0)
      etot=0.d0
      do i=1,frame % natoms
        mass(i) = getElementMass(frame % lab(i))
      enddo
      CALL plumed_f_gcmd("setStep"//char(0),tallyExecutions)
      CALL plumed_f_gcmd("setMasses"//char(0),mass)
      CALL plumed_f_gcmd("setForces"//char(0),force)
      CALL plumed_f_gcmd("setPositions"//char(0),frame % pos)
      CALL plumed_f_gcmd("setBox"//char(0),frame % hmat)
      CALL plumed_f_gcmd("setVirial"//char(0),virial)
      CALL plumed_f_gcmd("setEnergy"//char(0),etot)
      CALL plumed_f_gcmd("calc"//char(0),0)
      deallocate(mass,force,virial)

    end if

#endif
  end subroutine plumedInterface

end module modulePlumedInterface
