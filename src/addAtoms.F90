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
module moduleAddAtoms
  use moduleSystem
  use moduleVariables
  use moduleStrings
  use moduleRead
  use moduleFiles
  use moduleRandomNumbers
  use moduleDistances
  use moduleNeighbours
  use moduleElements
  use moduleMessages

  implicit none

  public :: addAtoms
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: inputFile
  integer, pointer :: numberOfNewMolecules
  real(8), pointer :: minimumDistance

  logical, pointer :: frameRead

  integer :: numberOfCreatedPositions = 0
  integer :: numberOfAddedMolecules = 0
  real(8), allocatable, dimension(:,:) :: newPositions

  real(8), pointer, dimension(:) :: localOrigin
  real(8), pointer, dimension(:) :: localCell

contains
  subroutine addAtoms(a)
    use moduleVariables
    use moduleSystem 
    implicit none
    type(actionTypeDef), target :: a

    call associatePointers(a)

    if (a % actionInitialisation) then
      a % actionInitialisation = .false.
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      if (firstAction) then

        if (frame % natoms > 0) then
          if (numberOfMolecules == 0) call runInternalAction("topology","NULL")
        end if

        call dumpScreenInfo(a)
        call extendFrame(frame, a % localFrame % natoms * numberOfNewMolecules)
        
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if

      call computeAction(a)

    end if

  end subroutine addAtoms

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    inputFile            => a % inputFile
    frameRead            => a % logicalVariables(1)
    numberOfNewMolecules => a % integerVariables(1)
    minimumDistance      => a % doubleVariables(1)

    localCell(1:9)       => a % doubleVariables(2:10)
    localOrigin(1:3)     => a % doubleVariables(11:13)

  end subroutine associatePointers

  subroutine initialiseAction(a)

    implicit none
    type(actionTypeDef), target :: a
    integer :: i, n
    real(8), dimension(3,3) :: hmat
    real(8), dimension(3) :: dij
    character(cp) :: atomName

    logical :: inputCell
    character(len=STRLEN) :: stringCell
    real(8), dimension(3,3) :: hnew

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0d0

    call assignFlagValue(actionCommand,"+cell",inputCell,.false.)
    if (inputCell) then
      call assignFlagValue(actionCommand,"+cell ", stringCell, "NONE")
      call readCellFreeFormat(stringCell, hnew)
      call hmat2cell(hnew, localCell, "DEG")
    else 
      call hmat2cell(frame % hmat, localCell, "DEG")
    end if
    
    ! Create a new cell if no coordinates had been read
    call getInverseCellMatrix(frame % hmat, frame % hinv, frame % volume)
    if (frame % volume < 1.1d0) then
      frame % hmat = hnew
      call getInverseCellMatrix(frame % hmat, frame % hinv, frame % volume)
      if (abs(frame % hmat(1,2)) + abs(frame % hmat(1,3)) + abs(frame % hmat(2,3)) .lt. 1.0d-6) then
        pbc_type = "ortho"
      else
        pbc_type = "tri"
      end if
      call hmat2cell(frame % hmat, frame % cell, "DEG")
      call initialisePBC(pbc_type)
      call initialiseNeighboursList()
    end if

    call assignFlagValue(actionCommand,"+origin",localOrigin,[0.d0,0.d0,0.d0])

    ! -> open file with coordinates
    call assignFlagValue(actionCommand,"+f",inputFile % fname,'NULL')
    
    ! -> or add an atom
    call assignFlagValue(actionCommand,"+atom",atomName,'NULL')
    
    if (inputFile % fname /= 'NULL') then
      call initialiseFile(inputFile, inputFile % fname)

      ! -> get number of atoms
      call getNumberOfAtoms(inputFile, n, hmat)
      rewind(inputFile % funit)

      call createSystemArrays(a % localFrame, n)

      call readCoordinates(frameRead, a % localFrame, inputFile)
      if (.not. frameRead) call message(-1,"--add | cannot read coordinates file")

    else if (atomName /= 'NULL') then
      call createSystemArrays(a % localFrame, 1)
      a % localFrame % lab(1) = atomName
      a % localFrame % pos = 0.d0
    
    else
      call message(-1,"--add | no atoms/molecule to add")
    endif

    ! translete com to the origin
    dij = a % localFrame % pos(:,1)
    do i=2,a % localFrame % natoms
      dij = dij + a % localFrame % pos(:,i)
    end do
    dij = dij / dble(a % localFrame % natoms)
    do i=1,a % localFrame % natoms
      a % localFrame % pos(1:3,i) = a % localFrame % pos(1:3,i) - dij(1:3)
    end do

    call assignFlagValue(actionCommand,"+n",numberOfNewMolecules,1)

    call assignFlagValue(actionCommand,"+rmin",minimumDistance,3.d0)

  end subroutine initialiseAction

  subroutine createRandomPositions()
    implicit none
    integer :: i, j, n
    real(8), dimension(3) :: rpos, dij
    real(8) :: rtmp, dmin
    real(8), allocatable, dimension(:,:) :: localPositions
    real(8), dimension(3,3) :: hmat

    dmin = minimumDistance**2
    call cell2hmat(localCell, hmat)
    if (numberOfCreatedPositions == 0) then
      allocate(newPositions(3,numberOfNewMolecules))
    else
      allocate(localPositions(3,numberOfCreatedPositions + numberOfNewMolecules))
      localPositions(:,1:numberOfCreatedPositions) = newPositions(:,1:numberOfCreatedPositions)
      call move_alloc(localPositions,newPositions)
    end if

    n = numberOfCreatedPositions
    do i=1,numberOfNewMolecules

      newpos : do 
        rpos(1:3) = localOrigin(1:3) + grnd() * hmat(:,1) + grnd() * hmat(:,2) + grnd() * hmat(:,3)

        do j=1,n
          dij = rpos(:) - newPositions(:,j)
          rtmp = computeDistanceSquaredPBC(dij)
          if (rtmp < dmin) cycle newpos
        end do

        n = n + 1
        newPositions(:,n) = rpos(:)
        exit
      end do newpos

    end do
    numberOfCreatedPositions = numberOfCreatedPositions + numberOfNewMolecules

  end subroutine createRandomPositions

  subroutine dumpScreenInfo(a)
    use moduleMessages     
    implicit none
    type(actionTypeDef), target :: a
    
    call message(0,"Add Atoms")
    call message(0,"...Size of the box to be filled",rv=localCell(1:3))
    call message(0,"...Angles of the box to be filled",rv=localCell(4:6))
    call message(0,"...Origin of the box to be filled",rv=localOrigin)
    if (inputFile % fname /= 'NULL') then
      call message(0,"...New moleculed from file",str=inputFile % fname)
      call message(0,"...Number of molcules to add",i=numberOfNewMolecules)
      call message(0,"...Minimum distance between the molecules",r=minimumDistance)
    else
      call message(0,"...New atom name",str=a % localFrame % lab(1))
      call message(0,"...Number of atoms to add",i=numberOfNewMolecules)
      call message(0,"...Minimum distance between the atoms",r=minimumDistance)
    end if

    call message(2)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: initialNumberOfAtoms
    integer :: i, n0, n
    real(8) :: dij(3), rtmp, rmax
    real(8) :: theta, phi, vec(3)
    real(8), allocatable, dimension(:,:) :: pos_tmp

    integer :: imol, idx, iatm, jatm
    logical, allocatable, dimension(:) :: ldelete
    logical :: newSystem = .false.

    call createRandomPositions()

    initialNumberOfAtoms = frame % natoms
    if (initialNumberOfAtoms == 0) newSystem = .true.

    frame % natoms = frame % natoms + (a % localFrame % natoms)*numberOfNewMolecules

    ! Array for the new molecules/atoms
    allocate(pos_tmp(3,a % localFrame % natoms), source=a % localFrame % pos)

    n = initialNumberOfAtoms
    do imol=1,numberOfNewMolecules
      
      ! Random rotation around the origin
      if (a % localFrame % natoms > 1)then
        theta  = grnd() * pi
        phi    = grnd() * twopi
        vec(1) = sin(theta) * cos(phi)
        vec(2) = sin(theta) * sin(phi)
        vec(3) = cos(theta)
        theta  = grnd() * twopi
        call rotateMolecule(vec(1:3), theta, a % localFrame % natoms, pos_tmp,-1)
      endif

      ! Add new molecules to the frame
      do i=1,a % localFrame % natoms
        frame % lab(n+i) = a % localFrame % lab(i)
        frame % pos(1:3,n+i) = pos_tmp(1:3,i) + newPositions(:,numberOfAddedMolecules + imol)
        frame % chg(n+i) = a % localFrame % chg(i) 
      end do

      n = n + a % localFrame % natoms
    enddo
    numberOfAddedMolecules = numberOfAddedMolecules + numberOfNewMolecules

    ! Removing initial molecules overlapping with the new ones
    if (newSystem) then
      call setUpNeigboursList()

    else
      allocate(ldelete(numberOfMolecules),source=.false.)
      rmax = (2.d0)**2
      m1 : do imol=1,numberOfMolecules
        do idx=1,listOfMolecules(imol) % numberOfAtoms
          iatm = listOfMolecules(imol) % listOfAtoms(idx)

          do jatm=initialNumberOfAtoms+1,frame % natoms
            dij = frame % pos(1:3,iatm) - frame % pos(1:3,jatm)
            rtmp = computeDistanceSquaredPBC(dij)
            if (rtmp < rmax) then
              ldelete(imol) = .true.
              cycle m1
            end if
          end do

        end do
      end do m1
    end if

    n0 = 0
    ! Remove the initial molecules that overlap with the solute
    do imol=1,numberOfMolecules
      if (ldelete(imol)) cycle
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        iatm = listOfMolecules(imol) % listOfAtoms(idx)
        n0 = n0 + 1
        frame % lab(n0) = frame % lab(iatm)
        frame % chg(n0) = frame % chg(iatm)
        frame % pos(:,n0) = frame % pos(:,iatm)
      end do
    end do

    ! Add the new molecules to the frame
    do iatm=initialNumberOfAtoms+1,frame % natoms
      n0 = n0 + 1
      frame % lab(n0) = frame % lab(iatm)
      frame % chg(n0) = frame % chg(iatm)
      frame % pos(:,n0) = frame % pos(:,iatm)
    end do
    frame % natoms = n0

    call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
    call updateNeighboursList(.true.)
    if (.not. newSystem) call runInternalAction("topology","+update")

  end subroutine computeAction


end module moduleAddAtoms

