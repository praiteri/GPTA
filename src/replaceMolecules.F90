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
module moduleReplaceMolecule
  use moduleVariables
  use moduleStrings
  use moduleFiles
  use moduleMessages 
  use moduleNeighbours
  use moduleRead
  use moduleAlignMolecules
  use moduleDistances
  use moduleProperties

  implicit none

  public :: replaceMolecules, replaceMoleculesHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  character(len=STRLEN), pointer :: targetMolecule
  integer, pointer :: numberOfLocalMolecules
  integer, pointer, dimension(:) :: localIndices

  type(fileTypeDef), pointer :: inputFile
  logical, pointer :: frameRead
  integer, pointer, dimension(:) :: currentIndices
  integer, pointer, dimension(:) :: referenceIndices

  logical, pointer :: computeOnlyRMSD
  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: ID

contains

  subroutine replaceMoleculesHelp()
    implicit none
    call message(0,"This action can be used to replace selected molecules with a different molecule read from a file.")
    call message(0,"The +is and +ir flags provide a mapping between the atoms in the system's and reference molecules.")
    call message(0,"Examples:")
    call message(0,"gpta.x --i co3.pdb --top --subs +mol M1 +is 1,2,3,4 +ir 1,4,3,2 +f hco3.pdb --o new.pdb")
    call message(0,"gpta.x --i slab.pdb --top --subs +mol M1 +all +f co3.pdb +rmsd rmsd.out")
  end subroutine replaceMoleculesHelp

  subroutine initialiseAction(a)
    use moduleStrings
    use moduleSystem
    implicit none
    type(actionTypeDef), target :: a
    integer :: i, n
    real(8), dimension(3,3) :: hmat

    integer, allocatable, dimension(:) :: idx
    logical :: lflag

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0d0

    ! molecule to replace
    call assignFlagValue(actionCommand,"+mol",targetMolecule,"NULL")

    ! file with the reference molecule
    call assignFlagValue(actionCommand,"+f",inputFile % fname,'NULL')
    if (inputFile % fname /= 'NULL') then
      call initialiseFile(inputFile, inputFile % fname)

      ! -> get number of atoms
      call getNumberOfAtoms(inputFile, n, hmat)
      rewind(inputFile % funit)

      call createSystemArrays(a % localFrame, n)

      call readCoordinates(frameRead, a % localFrame, inputFile)

      if (.not. frameRead) call message(-1,"--subs | cannot read coordinates file")

    else
      call message(-1,"--subs - missing reference structure file")

    end if

    call assignFlagValue(actionCommand,"+all",lflag,.false.)
    if (lflag) then
      n = a % localFrame % natoms
      allocate(idx(n))
      do i=1,n
        idx(i) = i
      enddo
      currentIndices(1:n) => a % integerVariables(101:100+n)
      currentIndices = idx
      referenceIndices(1:n) => a % integerVariables(151:150+n)
      referenceIndices = idx
      
    else
      ! building the matching arrays
      call assignFlagValue(actionCommand,"+is",idx)
      n = size(idx)
      currentIndices(1:n) => a % integerVariables(101:100+n)
      currentIndices = idx
      deallocate(idx)
      
      call assignFlagValue(actionCommand,"+ir",idx)
      n = size(idx)
      referenceIndices(1:n) => a % integerVariables(151:150+n)
      referenceIndices = idx
    end if
    
    call assignFlagValue(actionCommand,"+rmsd",computeOnlyRMSD,.false.)
    if (computeOnlyRMSD) then
      call assignFlagValue(actionCommand,"+rmsd",outputFile % fname,'rmsd.out')
      call initialiseFile(outputFile, outputFile % fname)
      call workData % initialise(ID, "dumpSeq", iounit=outputFile % funit)
    endif

    if ( size(referenceIndices) > 50 .or. size(currentIndices) > 50 ) call message(-1,"--subs - too many atoms in the molecule")

  end subroutine initialiseAction

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand           => a % actionDetails
    firstAction             => a % firstAction

    targetMolecule          => a % stringVariables(1)
    numberOfLocalMolecules  => a % integerVariables(1)
    localIndices            => a % localIndices

    inputFile               => a % inputFile
    frameRead               => a % logicalVariables(1)

    computeOnlyRMSD         => a % logicalVariables(2)
    ID                      => a % integerVariables(2)
    outputFile              => a % outputFile

  end subroutine associatePointers

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Replace Molecules Action")
    call message(0,"...Target molecule",str=targetMolecule)
    call message(0,"...New molecule read from",str=inputFile % fname)
    call message(0,"...Atoms used for the alignmen")
    call message(0,"......Current molecule",iv=currentIndices)
    call message(0,"......New molecule",iv=referenceIndices)
    call message(2)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleSystem
    implicit none
    type(actionTypeDef), target :: a

    integer :: nref, ncur, ilocal, imol
    real(8) :: rmsd, rvec(9), rotmat(3,3)
    real(8), allocatable, dimension(:,:) :: referencePositions
    real(8), allocatable, dimension(:,:) :: currentPositions
    real(8), allocatable, dimension(:) :: weights

    integer :: nidx0, nidx1, idx, jdx, iatm

    logical, allocatable, dimension(:) :: selectedMolecules
    
    integer :: newNumberOfAtoms
    real(8), allocatable, dimension(:,:) :: newPositions
    character(len=cp), allocatable, dimension(:) :: newLabels

    real(8), allocatable, dimension(:) :: localProperty

    nref = a % localFrame % natoms
    nidx0 = size(referenceIndices)

    imol = localIndices(1)
    ncur = listOfMolecules(imol) % numberOfAtoms
    nidx1 = size(currentIndices)

    if (nidx1 > nref) call message(-1,"--subs - wrong number of indices in +wr")
    if (nidx1 > ncur) call message(-1,"--subs - wrong number of indices in +ws")
    if (nidx1 /= nidx0) call message(-1,"--subs - +wr and +ws must have the same number of indices")

    ! allocating arrays for the alignment
    allocate(referencePositions(3,nref))
    allocate(currentPositions(3,nref))

    ! extracting the atoms for the reference structure to match
    allocate(weights(a % localFrame % natoms), source=0.d0)
    do idx=1,nidx0
      iatm = referenceIndices(idx)
      weights(iatm) = 1.d0
      referencePositions(1:3,idx) = a % localFrame % pos(1:3,iatm)
    end do

    ! shifting reference molecule and the subset of atoms used for alignement to the origin
    call CenterCoords(nidx0, referencePositions) 

    if (computeOnlyRMSD) then
      allocate(localProperty(numberOfLocalMolecules))
      do ilocal=1,numberOfLocalMolecules
        imol = localIndices(ilocal)
        do idx=1,nidx1
          jdx = currentIndices(idx)
          iatm = listOfMolecules(imol) % listOfAtoms(jdx)
          currentPositions(1:3,idx) = frame % pos(1:3,iatm)
        end do

        ! extracting the rotation matrix to rotate referencePositions onto currentPositions
        ! current positions are centered to the origin
        call Superimpose(nidx1, referencePositions, currentPositions, rmsd, rvec)
        localProperty(ilocal) = rmsd
      end do
      call workData % compute(ID, numberOfValues=numberOfLocalMolecules, xValues=localProperty)

    else

      ! using weighted function to get the same shift
      call CenterCoords(a % localFrame % natoms, a % localFrame % pos, weights) 

      ! allocating temporary arrays for the new positions
      allocate(selectedMolecules(numberOfMolecules), source=.false.)
      do ilocal=1,numberOfLocalMolecules
        imol = localIndices(ilocal)
        selectedMolecules(imol) = .true.
      end do
      newNumberOfAtoms = frame % natoms + count(selectedMolecules)*(nref-ncur)
      allocate(newPositions(3,newNumberOfAtoms))
      allocate(newLabels(newNumberOfAtoms))

      ! aligning and replacing the molecules
      newNumberOfAtoms = 0
      do imol=1,numberOfMolecules
      
        if (selectedMolecules(imol)) then
          ! extracting the atoms for the reference structure to match
          do idx=1,nidx1
            jdx = currentIndices(idx)
            iatm = listOfMolecules(imol) % listOfAtoms(jdx)
            currentPositions(1:3,idx) = frame % pos(1:3,iatm)
          end do

          ! extracting the rotation matrix to rotate referencePositions onto currentPositions
          ! current positions are centered to the origin
          call Superimpose(nidx1, referencePositions, currentPositions, rmsd, rvec)
          rotmat = reshape(rvec,[3,3])

          ! rotating the reference molecule + shifting its original centre of mass
          ! could lead to problems if old and new molecule are very different
          do idx=1,a % localFrame % natoms
            newPositions(1:3,newNumberOfAtoms+idx) = &
              matmul(rotmat , a % localFrame % pos(1:3,idx)) + listOfMolecules(imol) % centreOfMass
          end do
          do idx=1,a % localFrame % natoms
            newLabels(newNumberOfAtoms+idx) = a % localFrame % lab(idx)
          end do
          newNumberOfAtoms = newNumberOfAtoms + nref

        else

          ! ignoring this molecule, just copy it over the temporary arrays
          do idx=1,listOfMolecules(imol) % numberOfAtoms
            iatm = listOfMolecules(imol) % listOfAtoms(idx)
            newPositions(1:3,newNumberOfAtoms+idx) = frame % pos(1:3,iatm)
          end do
          do idx=1,listOfMolecules(imol) % numberOfAtoms
            iatm = listOfMolecules(imol) % listOfAtoms(idx)
            newLabels(newNumberOfAtoms+idx) = frame % lab(iatm)
          end do
          newNumberOfAtoms = newNumberOfAtoms + listOfMolecules(imol) % numberOfAtoms
        end if

      end do

      ! swapping allocations over
      frame % natoms = newNumberOfAtoms
      call move_alloc(newPositions, frame % pos)
      call move_alloc(newLabels, frame % lab)

      ! extending other frame arrays
      if (size(frame % chg) /= frame % natoms) then
        deallocate(frame % frac)
        deallocate(frame % chg)
        allocate(frame % frac(3,frame % natoms))
        allocate(frame % chg(frame % natoms), source=0.d0)
      end if

      ! recomputing topology
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
      call updateNeighboursList(.true.)
      call runInternalAction("topology","+update")
      
    end if ! RMSD only

  end subroutine computeAction

  subroutine replaceMolecules(a)
    use moduleVariables
    use moduleSystem 
    use moduleOpenMM
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a

    call associatePointers(a)
    
    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if
    
    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      
      if (firstAction) then
        call dumpScreenInfo()
        
        if (numberOfMolecules == 0) call runInternalAction("topology","NULL")
        
        ! get the indices of the molecules to replace
        call getLocalIndices(targetMolecule, a, numberOfLocalMolecules)
        localIndices => a % localIndices
        
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if

      call associatePointers(a)

      call computeAction(a)
    end if

  end subroutine replaceMolecules

end module moduleReplaceMolecule
