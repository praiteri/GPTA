! ! Copyright (c) 2021, Paolo Raiteri, Curtin University.
! ! All rights reserved.
! ! 
! ! This program is free software; you can redistribute it and/or modify it 
! ! under the terms of the GNU General Public License as published by the 
! ! Free Software Foundation; either version 3 of the License, or 
! ! (at your option) any later version.
! !  
! ! Redistribution and use in source and binary forms, with or without 
! ! modification, are permitted provided that the following conditions are met:
! ! 
! ! * Redistributions of source code must retain the above copyright notice, 
! !   this list of conditions and the following disclaimer.
! ! * Redistributions in binary form must reproduce the above copyright notice, 
! !   this list of conditions and the following disclaimer in the documentation 
! !   and/or other materials provided with the distribution.
! ! * Neither the name of the <ORGANIZATION> nor the names of its contributors 
! !   may be used to endorse or promote products derived from this software 
! !   without specific prior written permission.
! ! 
! ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! ! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! ! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! ! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! ! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! ! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! ! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! ! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! ! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! ! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ! 
module moduleSystem 
  use moduleVariables, only : frameTypeDef, fileTypeDef, moleculeTypeDef, MAXFILES, cp
  implicit none
 
  integer                                         :: numberInputFiles
  integer                                         :: currentInputFile
  type(fileTypeDef), dimension(MAXFILES), target  :: inputFileNames
  integer                                         :: numberOfFramesRead
  integer                                         :: numberOfFramesProcessed
  logical                                         :: frameReadSuccessfully
  logical                                         :: endOfCoordinatesFiles = .false.
  logical                                         :: lastFrameOnly = .false.
  integer                                         :: first_frame, last_frame, stride_frame, nlist_frames
  integer, allocatable, dimension(:)              :: list_frames
  type(frameTypeDef), target                      :: frame
  type(frameTypeDef), target                      :: referenceFrame

#ifdef GPTA_MPI
  logical, allocatable, dimension(:)              :: allFramesReadSuccessfully
  type(frameTypeDef), allocatable, dimension(:)   :: allFrames
#endif

  ! Global system properties  
  integer                                         :: inputNumberOfAtoms
  real(8)                                         :: totalMass
  real(8)                                         :: totalCharge
  integer                                         :: numberOfUniqueLabels
  character(len=cp), allocatable, dimension(:)    :: listOfUniqueLabels
  character(len=2), allocatable, dimension(:)     :: listOfElements

  logical                                         :: userDefinedCell = .false.
  real(8), dimension(3,3)                         :: userDefinedHMatrix
  character(len=5)                                :: pbc_type
  logical                                         :: inputCoordInNM = .false.
  logical                                         :: inputCoordInBohr = .false.
  logical                                         :: inputCoordInFractional = .false.
  logical                                         :: outputCoordInBohr = .false.

  ! Neighbours' list  
  logical                                         :: computeNeighboursList       = .false.
  logical                                         :: computeNeighboursListOnce   = .true.
  logical                                         :: computeNeighboursListDouble = .false.
  logical                                         :: computeNeighboursListVerlet = .false.
  real(8)                                         :: cutoffNeighboursListGlobal  = 3.d0
  integer                                         :: updatedListFrameNumber
  integer                                         :: neighboursListUpdates

  integer                                         :: neighmax = 100
  integer, allocatable, dimension(:)              :: nneigh
  integer, allocatable, dimension(:,:)            :: lneigh
  real(8), allocatable, dimension(:,:)            :: rneigh
  
  integer                                          :: numberOfMolecules
  type(moleculeTypeDef), allocatable, dimension(:) :: listOfMolecules
  integer                                          :: numberOfUniqueMoleculeTypes
  integer, allocatable, dimension(:)               :: uniqueMoleculesID
  type(moleculeTypeDef), allocatable, dimension(:) :: uniqueMoleculeTypes
  integer, allocatable, dimension(:)               :: atomToMoleculeIndex
  
  ! topology parameters
  integer, allocatable, dimension(:)               :: numberOfCovalentBondsPerAtom
  integer, allocatable, dimension(:,:)             :: listOfCovalentBondsPerAtom
  integer, allocatable, dimension(:)               :: numberOf_13_InteractionsPerAtom
  integer, allocatable, dimension(:,:)             :: listOf_13_InteractionsAtom
  integer, allocatable, dimension(:)               :: numberOf_14_InteractionsPerAtom
  integer, allocatable, dimension(:,:)             :: listOf_14_InteractionsAtom
         
  integer                                          :: numberOfUniqueBonds
  integer                                          :: numberOfUniqueAngles
  integer                                          :: numberOfUniqueTorsions
  integer                                          :: numberOfUniqueOutOfPlane
  integer, allocatable, dimension(:,:)             :: listOfUniqueBonds
  integer, allocatable, dimension(:,:)             :: listOfUniqueAngles
  integer, allocatable, dimension(:,:)             :: listOfUniqueTorsions
  integer, allocatable, dimension(:,:)             :: listOfUniqueOutOfPlane
           
  integer, parameter                               :: MAXBOND=10
  integer                                          :: sbond_n=0
  real(8), dimension(MAXBOND)                      :: sbond_d
  character(len=cp), dimension(2,MAXBOND)          :: sbond_l

  interface
    subroutine stuff(np, property, numberOfSelectedAtoms, selectionList)
      integer, intent(inout) :: np
      real(8), dimension(*), intent(out) :: property
      integer, intent(in), optional :: numberOfSelectedAtoms
      integer, dimension(:), intent(in), optional :: selectionList
    end subroutine stuff
  end interface

  interface
    subroutine dumpCoordinates(ftype, funit, natoms, pos, label, hmat, header)
      implicit none 
      character(len=*) :: ftype
      integer, intent(in) :: funit
      integer, intent(in) :: natoms
      real(8), dimension(3,natoms), intent(in) :: pos
      character(4), dimension(natoms), intent(in) :: label
      real(8), dimension(3,3), optional, intent(in) :: hmat
      logical, optional, intent(in) :: header
    end subroutine
  end interface

end module moduleSystem 

subroutine createSystemArrays(localFrame, numberOfAtoms)
  use moduleVariables, only : frameTypeDef
  use moduleMessages 
  implicit none
  type(frameTypeDef) :: localFrame
  integer, intent(in) :: numberOfAtoms
  
  localFrame % natoms = numberOfAtoms
  allocate(localFrame % lab(numberOfAtoms))
  allocate(localFrame % element(numberOfAtoms))
  allocate(localFrame % mass(numberOfAtoms))
  allocate(localFrame % pos(3,numberOfAtoms))
  allocate(localFrame % frac(3,numberOfAtoms))
  allocate(localFrame % chg(numberOfAtoms) , source=0.d0)
  
end subroutine createSystemArrays

subroutine deleteSystemArrays(localFrame)
  use moduleVariables, only : frameTypeDef
  use moduleMessages 
  implicit none
  type(frameTypeDef) :: localFrame
  
  localFrame % natoms = 0
  deallocate(localFrame % lab)
  deallocate(localFrame % element)
  deallocate(localFrame % mass)
  deallocate(localFrame % pos)
  deallocate(localFrame % frac)
  deallocate(localFrame % chg)
  
end subroutine deleteSystemArrays

subroutine extendFrame(n)
  use moduleSystem, only : frameTypeDef, frame, atomToMoleculeIndex
  integer, intent(in) :: n

  integer :: n0, nnew
  type(frameTypeDef) :: newFrame
  integer, allocatable, dimension(:) :: tmpIndices

  n0 = frame % natoms
  nnew = n0 + n

  if (n0 == 0) then
    if (allocated(frame % lab)) call deleteSystemArrays(frame)
    call createSystemArrays(frame, nnew)
    if (allocated(atomToMoleculeIndex)) deallocate(atomToMoleculeIndex)
    allocate(atomToMoleculeIndex(nnew))
    frame % natoms = 0
    return
  end if

  call createSystemArrays(newFrame, n0)

  newFrame % natoms         = frame % natoms 
  newFrame % lab (    1:n0) = frame % lab (    1:n0) 
  newFrame % pos (1:3,1:n0) = frame % pos (1:3,1:n0) 
  newFrame % frac(1:3,1:n0) = frame % frac(1:3,1:n0) 
  newFrame % chg (    1:n0) = frame % chg (    1:n0) 

  call deleteSystemArrays(frame)
  call createSystemArrays(frame, nnew)

  frame % natoms         = newFrame % natoms 
  frame % lab (    1:n0) = newFrame % lab (    1:n0) 
  frame % pos (1:3,1:n0) = newFrame % pos (1:3,1:n0) 
  frame % frac(1:3,1:n0) = newFrame % frac(1:3,1:n0) 
  frame % chg (    1:n0) = newFrame % chg (    1:n0) 

  allocate(tmpIndices(nnew))
  tmpIndices(1:n0) = atomToMoleculeIndex
  call move_alloc(tmpIndices, atomToMoleculeIndex)
  
end subroutine extendFrame

subroutine extendMolecules(n)
  use moduleSystem, only : moleculeTypeDef, listOfMolecules, numberOfMolecules!, atomToMoleculeIndex
  implicit none
  integer, intent(in) :: n
  type(moleculeTypeDef), allocatable, dimension(:) :: newMols
  integer :: n0, nnew
  integer :: imol

  n0 = numberOfMolecules
  nnew = n0 + n

  allocate(newMols(nnew))
  do imol=1,n0
    newMols(imol) % ID            = listOfMolecules(imol) % ID           
    newMols(imol) % resname       = listOfMolecules(imol) % resname      
    newMols(imol) % numberOfAtoms = listOfMolecules(imol) % numberOfAtoms
    newMols(imol) % centreOfMass  = listOfMolecules(imol) % centreOfMass 
    newMols(imol) % brokenBonds   = listOfMolecules(imol) % brokenBonds
!    allocate(newMols(imol) % listOfAtoms, source=listOfMolecules(imol) % listOfAtoms)
!    allocate(newMols(imol) % listOfLabels, source=listOfMolecules(imol) % listOfLabels)
    allocate(newMols(imol) % listOfAtoms(newMols(imol) % numberOfAtoms))
    allocate(newMols(imol) % listOfLabels(newMols(imol) % numberOfAtoms))
    newMols(imol) % listOfAtoms  = listOfMolecules(imol) % listOfAtoms
    newMols(imol) % listOfLabels = listOfMolecules(imol) % listOfLabels
  end do

  call move_alloc(newMols, listOfMolecules)
end subroutine extendMolecules

#ifdef GPTA_MPI
subroutine sendFrameToProcessor(f,idx)
  use moduleVariables
  implicit none
  type(frameTypeDef), intent(in) :: f 
  integer, intent(in) :: idx
  logical, save :: firstTimeIn(0:128) = .true.
  integer :: numberOfAtoms

  numberOfAtoms = f % natoms
  call MPI_Send(f % nframe, 1               , MPI_INT       , idx, 0, MPI_COMM_WORLD, ierr_mpi)
  call MPI_Send(f % hmat  , 9               , MPI_DOUBLE    , idx, 0, MPI_COMM_WORLD, ierr_mpi)
  call MPI_Send(f % pos   , numberOfAtoms*3 , MPI_DOUBLE    , idx, 0, MPI_COMM_WORLD, ierr_mpi)
  if (firstTimeIn(idx)) then
    firstTimeIn(idx) = .false.
    call MPI_Send(f % lab    , numberOfAtoms*cp, MPI_CHARACTER , idx, 0, MPI_COMM_WORLD, ierr_mpi)
    call MPI_Send(f % chg    , numberOfAtoms   , MPI_DOUBLE    , idx, 0, MPI_COMM_WORLD, ierr_mpi)
    call MPI_Send(f % element, numberOfAtoms*2 , MPI_CHARACTER , idx, 0, MPI_COMM_WORLD, ierr_mpi)
  end if
  return
end subroutine sendFrameToProcessor

subroutine receiveFrameFromMaster(f)
  use moduleVariables
  implicit none
  type(frameTypeDef), intent(inout) :: f 
  logical, save :: firstTimeIn = .true.
  character(cp), allocatable, dimension(:), save :: savedLabels
  real(8), allocatable, dimension(:), save :: savedCharges
  character(2), allocatable, dimension(:), save :: savedElements

  integer, save :: numberOfAtoms

  ! this is required becase atoms might be deleted by an action
  if (firstTimeIn) then
    numberOfAtoms = f % natoms
  else
    f % natoms = numberOfAtoms
  end if 

  call MPI_Recv(f % nframe, 1               , MPI_INT       , readingCPU, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)
  call MPI_Recv(f % hmat  , 9               , MPI_DOUBLE    , readingCPU, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)
  call MPI_Recv(f % pos   , numberOfAtoms*3 , MPI_DOUBLE    , readingCPU, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)
  if (firstTimeIn) then
    firstTimeIn = .false.
    call MPI_Recv(f % lab    , numberOfAtoms*cp, MPI_CHARACTER , readingCPU, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)
    call MPI_Recv(f % chg    , numberOfAtoms   , MPI_DOUBLE    , readingCPU, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)
    call MPI_Recv(f % element, numberOfAtoms*2 , MPI_CHARACTER , readingCPU, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)

    allocate(savedLabels(numberOfAtoms))
    allocate(savedCharges(numberOfAtoms))
    savedLabels = f % lab
    savedCharges = f % chg
    savedElements = f % element

  else

    f % lab = savedLabels
    f % chg = savedCharges
    f % element = savedElements

  end if

end subroutine receiveFrameFromMaster
#endif

