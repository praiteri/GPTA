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
module moduleMeanSquareDisplacement
  use moduleVariables
  use moduleFiles
  use moduleProperties
  use moduleMessages
  use moduleSystem
  use moduleStrings

  implicit none

  public :: msdAction
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions

  integer, pointer :: numberOfBins
  integer, pointer :: ID
  real(8), pointer :: timeInterval

contains

  subroutine msdAction(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: nsel

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if
    
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
        nsel = count(a % isSelected(:,1))
        allocate(a % localPositions(3,nsel))
        call createSelectionList(a,1)
        
        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else

        if (a % updateAtomsSelection) then
          call selectAtoms(1,actionCommand,a)
          call createSelectionList(a,1)
        end if

      end if

      if (tallyExecutions > numberOfBins+1) tallyExecutions = 1
      if (tallyExecutions == 1) then
        nsel = count(a % isSelected(:,1))
        call createReferenceCoordinatesFromList(nsel, a % idxSelection(:,1), a % localPositions)
      else
        call computeAction(a)
      end if
    end if

    if (endOfCoordinatesFiles) then
      call finaliseAction()
    end if 
  end subroutine msdAction

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    a % actionInitialisation = .false.

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile

    numberOfBins         => a % numberOfBins    
    ID                   => a % integerVariables(1)
    timeInterval         => a % doubleVariables(1)

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.

    call assignFlagValue(actionCommand,"+out",outputFile % fname,'msd.out')
    
    call assignFlagValue(actionCommand,"+nt",numberOfBins,1000)
    call assignFlagValue(actionCommand,"+dt",timeInterval,1.d0)
    
    call workData % initialise(ID, "avgDist ", numberOfBins=[numberOfBins], lowerLimits=[0.d0], upperLimits=[dble(numberOfBins)])

    tallyExecutions = 0
  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Computing mean square displacement")
    call message(0,"...Maximum number of frames",i=numberOfBins)
    call message(0,"...Output file",str=outputFile % fname)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: i, nsel, iatm
    real(8) :: msd, dij(3)
    real(8) :: time

    nsel = count(a % isSelected(:,1))
    msd = 0.d0
    do i=1,nsel
      iatm = a % idxSelection(i,1)
      dij = frame % pos(1:3,iatm) - a % localPositions(1:3,i)
      msd = msd + sum(dij*dij)
    end do
    msd = msd / dble(nsel)

    time = tallyExecutions - 1.01
    call workData % compute(ID, numberOfValues=1, xValues=[msd], yValues=[time])

  end subroutine computeAction

  subroutine finaliseAction()
    implicit none
    real(8), allocatable, dimension(:) :: msd
    integer, allocatable, dimension(:) :: tally
    real(8), allocatable, dimension(:) :: time
    integer :: i, idx, nbin, n
    real(8) :: alpha, beta, rfact

    ! call workData % extract(ID, msd)

    call workData % extract(ID, values=msd, tally=tally)

    if (me/=0) return

    nbin = size(tally)
    do i=1,nbin
      if (tally(i) == 0) then
        nbin = i-1
        exit
      end if
      msd(i) = msd(i) / tally(i)
    end do

    allocate(time(nbin))
    do i=1,nbin
      time(i) = i * timeInterval
    end do

    time = time * numberOfActiveProcesses

    idx = int(0.2*nbin)
    n = nbin - idx + 1

    ! y = alpha + beta*x
    call linearRegression(n,time(idx:nbin),msd(idx:nbin),alpha,beta,rfact)

    beta = beta / 6.d0
    beta = beta * 10.d0 ! conversion from angstrom^2/ps to 10^-5 cm^2/s

    call message(2)
    call message(0,"Self diffusion coefficient calculation")
    call message(0,"  Average D0 [10^-5 cm^2/s]",r=beta)
    call message(1,"  Average correlation factor",r=rfact)
    call initialiseFile(outputFile,outputFile % fname)
    write(outputFile % funit,'("# Time [ps] | MSD [A^3]")')
    do i=1,nbin
      write(outputFile % funit,'(f7.3,1x,f12.5)') time(i), msd(i)
    enddo
    close(outputFile % funit)

  end subroutine finaliseAction

  subroutine createReferenceCoordinatesFromList(natoms,isel,pos)
    implicit none
    integer, intent(in) :: natoms
    integer, intent(in), dimension(natoms) :: isel
    real(8), intent(out), dimension(:,:) :: pos

    integer :: iatm

    do iatm=1,natoms
      pos(1:3,iatm) = frame % pos(1:3,isel(iatm))
    end do

  end subroutine createReferenceCoordinatesFromList

end module moduleMeanSquareDisplacement
