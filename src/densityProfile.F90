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
module moduleDensityProfile
  use moduleVariables
  use moduleFiles
  use moduleMessages
  use moduleProperties

  implicit none

  public :: computeDensityProfile
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions

  integer, pointer :: numberOfBins
  integer, pointer :: iAxis
  integer, pointer :: ID

  real(8), pointer :: averageSize, averageVolume

  contains

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile

    numberOfBins         => a % numberOfBins
    iAxis                => a % integerVariables(1)
    ID                   => a % integerVariables(2)

    averageVolume        => a % doubleVariables(1)
    averageSize          => a % doubleVariables(2)

  end subroutine 

  subroutine initialiseAction(a)
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a

    integer :: i
    logical :: lflag(3)

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0d0

    call assignFlagValue(actionCommand,"+out",outputFile % fname,'dmap1D.out')

    call assignFlagValue(actionCommand,"+nbin",numberOfBins,100)

    call assignFlagValue(actionCommand,"+x",lflag(1),.false.)
    call assignFlagValue(actionCommand,"+y",lflag(2),.false.)
    call assignFlagValue(actionCommand,"+z",lflag(3),.false.)
    if (count(lflag)==0) call message(-1,"--dmap1D no direction specified +x/+y/+z")
    if (count(lflag)>1) call message(-1,"--dmap1D more than one direction specified +x/+y/+z")
    do i=1,3
      if (lflag(i)) iAxis = i
    enddo
    allocate(a % array1D(0:numberOfBins+1) , source=0.d0)

    tallyExecutions = 0 
    averageSize = 0.d0
    averageVolume = 0.d0

    call workData % initialise(ID, "dist1D", numberOfBins=[numberOfBins], lowerLimits=[0.d0], upperLimits=[1.d0])

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    use moduleMessages 
    implicit none
    call message(0,"Computing 1D density profile")
    call message(0,"...Output file",str=outputFile % fname)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleVariables
    use moduleSystem 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a

    integer :: iatm, nlocal
    real(8), allocatable, dimension(:), save :: fractionalCoord

    averageSize = averageSize + frame % hmat(iAxis,iAxis)
    averageVolume = averageVolume + frame % volume

    ! Local array with the coordinates of only the selected atoms
    nlocal = count(a % isSelected(:,1))
    allocate(fractionalCoord(nlocal))
    nlocal = 0
    do iatm=1,frame % natoms
      if (.not. a % isSelected(iatm,1)) cycle
      nlocal = nlocal + 1
      fractionalCoord(nlocal) = frame % frac(iAxis,iatm)
    enddo
    call workData % compute(ID, numberOfValues=nlocal, xValues=fractionalCoord)
    deallocate(fractionalCoord)

  end subroutine computeAction

  subroutine finaliseAction()
    implicit none
    ! type(actionTypeDef), target :: a
    real(8) :: dVolume

    call initialiseFile(outputFile,outputFile % fname)
    averageSize = averageSize / tallyExecutions
    dVolume = 1.d0/ (averageVolume / dble(numberOfBins) / tallyExecutions)
    write(outputFile % funit,'("# Position  | Density [atom/angstrom^3]")')
    call workData % dump(ID, outputFile % funit, upperLimits=[averageSize], normalisation=dVolume)
    close(outputFile % funit)

  end subroutine finaliseAction

  subroutine computeDensityProfile(a)
    use moduleSystem 
    implicit none
    type(actionTypeDef), target :: a

    call associatePointers(a)

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      if (firstAction) then
        call dumpScreenInfo()

        if (keepFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if

        call selectAtoms(1,actionCommand,a)
      
        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else

        if (a % updateAtomsSelection) call selectAtoms(1,actionCommand,a)

      end if
      
      call computeAction(a)

    end if

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) then
      call finaliseAction()
    end if

  end subroutine computeDensityProfile

end module moduleDensityProfile
