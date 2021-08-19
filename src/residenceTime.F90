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
module moduleResidenceTime

  use moduleVariables
  use moduleFiles
  use moduleSystem
  use moduleMessages

  implicit none

  public :: computeResidenceTime, computeResidenceTimeHelp
  private
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  
  integer, pointer :: numberOfBins
  integer, pointer :: tallyExecutions
  real(8), pointer :: cutoffRadius
  integer, pointer, dimension(:,:,:) :: coordinationList

  integer, pointer :: numberOfCentres
  integer, pointer :: MAX_NEIGH
  integer, pointer :: threshold

contains

  subroutine computeResidenceTimeHelp()
    implicit none
    call message(0,"This action computes the solvent residence time around selected solute atoms.")
    call message(0,"For efficiency reasons it requires to provide the number of frames in the trajectory.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --restime +s Ca OW +rcut 3.2 +thres 1 +ntraj 1000")
  end subroutine computeResidenceTimeHelp

  subroutine computeResidenceTime(a)
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

        if (resetFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if

        ! create a list of the atoms' indices for each group
        call selectAtoms(2,actionCommand,a)
        call createSelectionList(a,2)

        numberOfCentres = count(a % isSelected(:,1))
        allocate(a % Iarray3D(MAX_NEIGH,numberOfCentres,numberOfBins) , source=0)
        call associatePointers(a)

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else

        if (a % updateAtomsSelection) call selectAtoms(1,actionCommand,a)

      end if

      call computeAction(a)

    end if

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) call finaliseAction()

  end subroutine computeResidenceTime

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile
    numberOfBins         => a % numberOfBins
    cutoffRadius         => a % doubleVariables(1)
    numberOfCentres      => a % integerVariables(1)
    MAX_NEIGH            => a % integerVariables(2)
    threshold            => a % integerVariables(3)
    coordinationList     => a % Iarray3D

  end subroutine associatePointers

  subroutine initialiseAction(a)
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .true.
    a % requiresNeighboursListDouble = .true.
    a % cutoffNeighboursList = 4.0d0

    call assignFlagValue(actionCommand,"+rcut",cutoffRadius,a % cutoffNeighboursList)
    call assignFlagValue(actionCommand,"+ntraj",numberOfBins,1000)
    call assignFlagValue(actionCommand,"+nmax",MAX_NEIGH,10)
    call assignFlagValue(actionCommand,"+thres",threshold,1)
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'restime.out')

    a % cutoffNeighboursList = cutoffRadius
    tallyExecutions = 0 

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    use moduleMessages 
    implicit none
    call message(0,"Residence time calculation setup")
    call message(0,"...Maximum number of neighbours",i=MAX_NEIGH)
    call message(0,"...Maximum trajectory length",i=numberOfBins)
    call message(0,"...Frames threshold",i=threshold)
    call message(0,"...Cutoff distance",r=cutoffRadius)
    call message(1,"...Output file",str=outputFile % fname)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleVariables
    use moduleSystem 
    implicit none
    type(actionTypeDef), target :: a

    integer :: iatm, jatm, ineigh
    integer, allocatable, dimension(:,:) :: localList
    integer :: icentre, nshell, nsel
    
    allocate(localList(MAX_NEIGH,numberOfCentres), source=0)
    
    icentre = 0
    nsel = numberOfCentres
    do icentre=1,nsel
      iatm = a % idxSelection(icentre,1)
      nshell = 0
      do ineigh = 1 , nneigh(iatm)
        if (rneigh(ineigh,iatm)>cutoffRadius) cycle
        jatm = lneigh(ineigh,iatm)
        if (a % isSelected(jatm,2)) then
          nshell = nshell + 1
          localList(nshell,icentre) = jatm
        end if
      end do
      if (nshell > MAX_NEIGH) call message(-1,"--restime | too many neighbours (+nmax)",i=nshell)
    end do
    if (frame % nframe > numberOfBins) call message(-1,"--restime | trajctory is too long (+ntraj)")
    coordinationList(:,:,frame % nframe) = localList

  end subroutine computeAction

  subroutine finaliseAction()
    implicit none
    integer :: idx, nf
    integer, allocatable, dimension(:,:) :: localList
    real(8), allocatable, dimension(:) :: residenceTime
    real(8), allocatable, dimension(:) :: exchangeProbability
#ifdef GPTA_MPI
    integer :: nn
#endif
    
    nf = frame % nframe

#ifdef GPTA_MPI
    call MPI_allreduce(MPI_IN_PLACE, nf, 1, MPI_INT, MPI_MAX, MPI_Working_Comm, ierr_mpi)
#endif

    allocate(residenceTime(0:nf), source=0.d0)
    allocate(exchangeProbability(0:nf), source=0.d0)
    allocate(localList(MAX_NEIGH,nf), source=0)

    do idx=1,numberOfCentres
      ! Coordination shell for one atoms in the list
      localList(1:MAX_NEIGH, 1:nf) = coordinationList(1:MAX_NEIGH, idx, 1:nf)

#ifdef GPTA_MPI
      nn = MAX_NEIGH * nf
      if (me == 0) then
        call MPI_reduce(MPI_IN_PLACE, localList, nn, MPI_INT, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)
      else 
        call MPI_reduce(localList, localList, nn, MPI_INT, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)
      end if
#endif

      call rest_original(MAX_NEIGH, nf, localList, threshold, residenceTime, exchangeProbability)
    end do

    if (me /= 0) return
    residenceTime = residenceTime / numberOfCentres
    residenceTime = residenceTime / residenceTime(0)

    exchangeProbability = exchangeProbability / numberOfCentres
    exchangeProbability = exchangeProbability / sum(exchangeProbability)

    call initialiseFile(outputFile,outputFile % fname)

    write(outputFile % funit,'("# Frames   | Survival func | Exchange probability")')
    do idx=0,nf
      write(outputFile % funit,'(i6,2(1x,f16.5))')idx, residenceTime(idx), exchangeProbability(idx)
    end do
  end subroutine finaliseAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rest_original(MAX_NEIGH,frames_read,local_list,thres,survivalFunction,distribution)
    implicit none
    integer, intent(in) :: frames_read, MAX_NEIGH

    integer, intent(in) :: thres
    integer, dimension(MAX_NEIGH,frames_read) :: local_list
    real(8), dimension(0:frames_read) :: survivalFunction
    real(8), dimension(0:frames_read) :: distribution

    integer, allocatable, dimension(:,:) :: coord

    integer :: i, j, n, n0

    ! coord(1,:) Current coordiantion shell
    ! coord(2,:) Time since the molecule has entered the coordination shell
    ! coord(3,:) Time since the molecule has left the coordination shell
    allocate(coord(3,MAX_NEIGH*2))
    coord=0

    n0 = count(local_list(:,1) > 0)

    ! Initial coordination shell
    coord(1,1:n0) = local_list(1:n0,1)

    do i=2,frames_read

      ! Number of atoms that were in the coordination shell in the previous step
      n0 = count(coord(1,:) > 0)
      ! Check if they are still alive
      do j=1,n0

        ! Yes
        if (any(local_list(:,i)==coord(1,j))) then
          if (coord(3,j)>0) then
            coord(2,j) = coord(2,j) + coord(3,j) ! Just passed out for a little while
            coord(3,j) = 0
          endif
          coord(2,j) = coord(2,j) + 1

          ! Not there anymore
        else
          coord(3,j) = coord(3,j) + 1
        endif

        ! It's been dead for too long - bury it
        ! -> add to survival function
        if ( coord(3,j) > thres) then
          if (coord(2,j) > thres) then
            survivalFunction(0:coord(2,j)) = survivalFunction(0:coord(2,j)) + 1.d0
            distribution(coord(2,j)) = distribution(coord(2,j)) + 1
          end if
          coord(:,j) = 0
        endif

      enddo

      ! Compact the list
      n=0
      do j=1,n0
        if (coord(1,j)>0) then
          n = n + 1
          coord(:,n) = coord(:,j)
        endif
      enddo

      ! Add new neighbours at the end of the list
      do j=1,count(local_list(:,i)>0)
        if (all(local_list(j,i)/=coord(1,:))) then
          n = n + 1
          coord(1,n) = local_list(j,i)
          coord(2:3,n) = 0
        endif
      enddo

      ! Clear the end of the array
      do j=n+1,MAX_NEIGH
        coord(:,j) = 0
      enddo

    enddo

    ! Add the molecules still in the list to the survival function
    do j=1,n
      survivalFunction(0:coord(2,j)) = survivalFunction(0:coord(2,j)) + 1.d0
    enddo

  end subroutine rest_original

end module moduleResidenceTime
