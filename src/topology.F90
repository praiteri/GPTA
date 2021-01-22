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
module moduleMolecularTopology
  use moduleVariables
  use moduleFiles
  use moduleElements 
  use moduleSystem
  use moduleStrings
  use m_refsor
  use moduleMessages 
  use moduleResizeArrays
  use moduleDistances
  use moduleNeighbours

  implicit none

  public :: computeTopology
  private

  character(:), pointer :: actionCommand
  integer, pointer :: tallyExecutions

  logical, pointer :: updateTopology
  logical, pointer :: updateTopologyOnce
  logical, pointer :: rebuildMolecules
  logical, pointer :: reorderAtoms
  logical, pointer :: userDefinedMolecules

  character(STRLEN), pointer, dimension(:) :: moleculesLabels

contains

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand          => a % actionDetails
    tallyExecutions        => a % tallyExecutions

    updateTopology         => a % logicalVariables(1)
    updateTopologyOnce     => a % logicalVariables(2)
    rebuildMolecules       => a % logicalVariables(3)
    reorderAtoms           => a % logicalVariables(4)
    
    userDefinedMolecules   => a % logicalVariables(5)
    moleculesLabels(1:)    => a % stringVariables(1:)

  end subroutine associatePointers

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    character(STRLEN), allocatable, dimension(:) :: tmpLabels
    character(STRLEN) :: flagString
    logical :: lflag

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0d0    
  
    call assignFlagValue(actionCommand,"+update",lflag,.false.)
    if (lflag) then
      updateTopologyOnce = .false.
    else
      updateTopologyOnce = .true.
    end if
    if (.not. updateTopologyOnce) a % requiresNeighboursListUpdates = .true.
    
    call assignFlagValue(actionCommand,"+ion",tmpLabels)
    if (allocated(tmpLabels)) then
      call defineIonicSpecies(tmpLabels)
      deallocate(tmpLabels)
    end if

    call assignFlagValue(actionCommand,"+rebuild",rebuildMolecules,.false.)
    
    call assignFlagValue(actionCommand,"+reorder",reorderAtoms,.false.)

    call assignFlagValue(actionCommand,"+def",userDefinedMolecules,.false.)
    if (userDefinedMolecules) then
      call extractFlag(actionCommand,"+def",flagString)
      call parse(flagString," ",moleculesLabels,numberOfUniqueMoleculeTypes)
    end if

    updateTopology = .true.
    tallyExecutions = 0

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    implicit none
    integer :: idx, jdx, kdx 
    integer :: iatm, n
    character(len=30) :: str, str1
    real(8) :: molecularCharge

    call message(0,"Topology information")
    call message(0,"...Number of molecules found",i=numberOfMolecules)
    call message(0,"...Number of unique molecules found",i=numberOfUniqueMoleculeTypes)
    do idx=1,numberOfUniqueMoleculeTypes
      kdx = uniqueMoleculesID(idx)
      write(str1,'(i0)')idx

      n = count(listOfMolecules(:) % ID == idx)
      if (n == 0) cycle
      
      call message(3)
      call message(0,"....Molecule type M"//trim(str1))

      call message(0,"......Number of molecules",i=n)
      call message(0,"......Number of atoms in molecule",i=listOfMolecules(kdx) % numberOfAtoms)

      str = trim(listOfMolecules(kdx) % listOfLabels(1))
      do jdx=2,listOfMolecules(kdx) % numberOfAtoms
        str = trim(str)//","//trim(listOfMolecules(kdx) % listOfLabels(jdx))
      enddo
      call message(0,"......Atoms types in molecule",str=str)

      call message(0,"......Number of atoms in molecule",i=listOfMolecules(kdx) % numberOfAtoms)

      
      iatm = listOfMolecules(kdx) % listOfAtoms(1)
      molecularCharge = frame % chg(iatm)
      do jdx=2,listOfMolecules(kdx) % numberOfAtoms
        iatm = listOfMolecules(kdx) % listOfAtoms(jdx)
        molecularCharge = molecularCharge + frame % chg(iatm) 
      enddo
      call message(0,"......Molecular charge",r=molecularCharge)
      if (frame % volume > 1e-6) &
        call message(0,"......Concentration [M]",r=1660.5d0*n/frame % volume)
    enddo
    call message(2)

  end subroutine dumpScreenInfo

  subroutine defineMoleculesFromDistances()
    implicit none
    ! type(actionTypeDef), target :: a

!    logical, external :: checkForBrokenMolecules

    integer :: iatm, jatm, ineigh, idx, jdx, itmp
    character(cp) :: l1 ,l2
    real(8) :: rmax

    logical, allocatable, dimension(:) :: atomsUsed
    type(moleculeTypeDef), allocatable, dimension(:) :: m

    allocate(numberOfCovalentBondsPerAtom(frame % natoms), source=0)
    allocate(listOfCovalentBondsPerAtom(neighmax,frame % natoms), source=0)
    allocate(atomToMoleculeIndex(frame % natoms), source=0)

    ! Compute list of all covalent bond for each atom
    numberOfCovalentBondsPerAtom = 0
    
    ! Double neighbours' list
    if (computeNeighboursListDouble) then
      do iatm=1,frame % natoms
        l1 = frame % lab(iatm)
        do ineigh=1,nneigh(iatm)
          jatm = lneigh(ineigh,iatm)
          l2 = frame % lab(jatm)
          rmax = getMaximumBondLength(l1,l2,distanceScaling)
          if (rneigh(ineigh,iatm) <= rmax) then
            numberOfCovalentBondsPerAtom(iatm) = numberOfCovalentBondsPerAtom(iatm) + 1
            listOfCovalentBondsPerAtom(numberOfCovalentBondsPerAtom(iatm),iatm) = jatm
          end if
        enddo
      enddo    
      ! Half neighbours' list
    else
      do iatm=1,frame % natoms
        l1 = frame % lab(iatm)
        do ineigh=1,nneigh(iatm)
          jatm = lneigh(ineigh,iatm)
          l2 = frame % lab(jatm)
          rmax = getMaximumBondLength(l1,l2,distanceScaling)
          if (rneigh(ineigh,iatm) <= rmax) then
            numberOfCovalentBondsPerAtom(iatm) = numberOfCovalentBondsPerAtom(iatm) + 1
            listOfCovalentBondsPerAtom(numberOfCovalentBondsPerAtom(iatm),iatm) = jatm
            numberOfCovalentBondsPerAtom(jatm) = numberOfCovalentBondsPerAtom(jatm) + 1
            listOfCovalentBondsPerAtom(numberOfCovalentBondsPerAtom(jatm),jatm) = iatm
          end if
        enddo
      enddo
    end if
    ! Search for molecules
    allocate(atomsUsed(frame % natoms), source=.false.)

    allocate(m(frame % natoms))
    do idx=1,frame % natoms
      allocate(m(idx) % listOfAtoms(100))
    enddo
    iatm = 0
    numberOfMolecules = 0
    do while (iatm<frame % natoms)
      iatm = iatm + 1
      if (atomsUsed(iatm)) cycle

      numberOfMolecules = numberOfMolecules + 1
      m(numberOfMolecules) % numberOfAtoms = 1
      m(numberOfMolecules) % listOfAtoms(1) = iatm

      idx = 0 
      do while (idx < m(numberOfMolecules) % numberOfAtoms)
        idx = idx + 1
        jatm = m(numberOfMolecules) % listOfAtoms(idx)
        call  addAtomsToMolecule(m(numberOfMolecules), numberOfCovalentBondsPerAtom(jatm), &
                                  listOfCovalentBondsPerAtom(1:numberOfCovalentBondsPerAtom(jatm),jatm))
      enddo
      ! Found a new molecule
      do idx=1,m(numberOfMolecules) % numberOfAtoms
        jatm = m(numberOfMolecules) % listOfAtoms(idx)
        atomsUsed(jatm) = .true.
        atomToMoleculeIndex(jatm) = numberOfMolecules
      enddo

    enddo

    ! Consolidate list of molecules into the global array
    allocate(listOfMolecules(numberOfMolecules))
    do idx=1,numberOfMolecules
      itmp = m(idx) % numberOfAtoms
      listOfMolecules(idx) % numberOfAtoms = itmp
      allocate(listOfMolecules(idx) % listOfAtoms(itmp))
      allocate(listOfMolecules(idx) % listOfLabels(itmp))
      call refsor(m(idx) % listOfAtoms(1:itmp))
      listOfMolecules(idx) % listOfAtoms = m(idx) % listOfAtoms(1:itmp)
      do jdx=1,itmp
        listOfMolecules(idx) % listOfLabels(jdx) = frame % lab(listOfMolecules(idx) % listOfAtoms(jdx))
      enddo
    enddo
    deallocate(m)

    call finaliseTopology()

  end subroutine defineMoleculesFromDistances

  subroutine finaliseTopology()
    implicit none
    integer :: idx, jdx, kdx, itmp, i, j, n
    character(len=30) :: str

    ! Search for unique molecules' types
    allocate(uniqueMoleculesID(numberOfMolecules)) 
    itmp = 1
    uniqueMoleculesID(itmp) = 1 ! Contains the index of the first molecule of each unique type
    listOfMolecules(itmp) % ID = 1 ! Contains the type ID of each molecule

    do idx=2,numberOfMolecules
      do jdx=1,itmp
        kdx = uniqueMoleculesID(jdx)
        if ( listOfMolecules(idx) % numberOfAtoms /= listOfMolecules(kdx) % numberOfAtoms ) cycle
        
        ! Found - most common case : atoms in the same order
        if ( all(listOfMolecules(idx) % listOfLabels == listOfMolecules(kdx) % listOfLabels) ) then
          listOfMolecules(idx) % ID = jdx
          exit
        end if
        
      enddo

      if (jdx>itmp) then
        itmp = itmp + 1
        uniqueMoleculesID(itmp) = idx
        listOfMolecules(idx) % ID = itmp
      end if
    enddo 
    numberOfUniqueMoleculeTypes = itmp

    ! Assign a name to the molecules
    do idx=1,numberOfMolecules
      write(str,'(i0)') listOfMolecules(idx) % ID
      listOfMolecules(idx) % resname = "M"//trim(str) 
    enddo

    ! Upper limit for bonds
    numberOfUniqueBonds = 0
    n = sum(numberOfCovalentBondsPerAtom)
    allocate(listOfUniqueBonds(2,n))

    do i=1,frame % natoms
      do j=1,numberOfCovalentBondsPerAtom(i)
        if (listOfCovalentBondsPerAtom(j,i) < i) cycle
        numberOfUniqueBonds = numberOfUniqueBonds + 1
        listOfUniqueBonds(1:2,numberOfUniqueBonds) = [i,listOfCovalentBondsPerAtom(j,i)]
      enddo
    enddo
    call resizeArray(listOfUniqueBonds,numberOfUniqueBonds)

  end subroutine finaliseTopology

  subroutine addAtomsToMolecule(mol,nat,lst) 
    implicit none
    type(moleculeTypeDef), intent(inout) :: mol
    integer, intent(in) :: nat
    integer, dimension(nat) :: lst
    integer :: iadd
    integer :: i, j, ncurr
    iadd = 0
    ncurr = mol % numberOfAtoms
    do i=1,nat
      j = lst(i)
      if (any(mol % listOfAtoms(1:ncurr) == j)) cycle
      mol % numberOfAtoms = mol % numberOfAtoms + 1
      iadd = iadd + 1
      mol % listOfAtoms(ncurr+iadd) = j
    enddo
    return
  end subroutine addAtomsToMolecule

  function checkPermutations(n,a,b,idx) result(found)
    use moduleVariables, only : cp
    implicit none
    integer, intent(in) :: n
    character(len=cp), dimension(n), intent(in) :: a, b
    integer, dimension(n), intent(inout) :: idx
    
    integer :: i
    logical :: found
    logical, external :: nextp
    ! integer, allocatable, dimension(:) :: idx
    character(cp), allocatable, dimension(:) :: c

    found = .false.

    ! this should not be necessary
    ! if ( all(a==b) ) then
    !   found = .true.
    !   return
    ! end if
    
    ! allocate(idx(n))
    allocate(c(n), source=b)
    do i=1,n
      idx(i) = i
    enddo
    do while (nextp(n,idx))
      do i=1,n
        c(i) = b(idx(i))
      enddo
      if ( all(a==c) ) then
        found = .true.
        return
      end if
    enddo

    return
  end function checkPermutations

  subroutine defineMoleculesFromLabels()
    implicit none
    integer :: imol, iatm, itmp, idx, jdx, i
    integer :: natoms, ntmp
    character(cp), dimension(50) :: labels
    integer, allocatable, dimension(:) :: midx
    character(cp), allocatable, dimension(:,:) :: mlab
    character(cp), allocatable, dimension(:) :: tmplab

    type(moleculeTypeDef), allocatable, dimension(:) :: m

    allocate(atomToMoleculeIndex(frame % natoms))

    natoms = 0
    do imol=1,numberOfUniqueMoleculeTypes
      call parse(moleculesLabels(imol),",",labels,ntmp)
      natoms = max(natoms,ntmp)
    end do
    allocate(mlab(natoms,numberOfUniqueMoleculeTypes))
    allocate(midx(numberOfUniqueMoleculeTypes))
    allocate(tmplab(natoms))
    do imol=1,numberOfUniqueMoleculeTypes
      call parse(moleculesLabels(imol),",",mlab(:,imol),midx(imol))
    end do

    numberOfMolecules = 0
    allocate(m(frame % natoms))
    do idx=1,frame % natoms
      allocate(m(idx) % listOfAtoms(50))
    enddo

    iatm = 0
    do while (iatm < frame % natoms)
      do imol=1,numberOfUniqueMoleculeTypes
        natoms = midx(imol)
        if (iatm+natoms > frame % natoms) cycle
        tmplab = frame % lab(iatm+1 : iatm+natoms )
        if ( all(tmplab(1:midx(imol)) == mlab(1:midx(imol),imol)) ) then
          numberOfMolecules = numberOfMolecules + 1
          m(numberOfMolecules) % numberOfAtoms = midx(imol)
          do i=1,midx(imol)
            m(numberOfMolecules) % listOfAtoms(i) = iatm+i
            atomToMoleculeIndex(iatm+i) = numberOfMolecules
          end do 
          iatm = iatm + midx(imol)
          exit
        end if
      end do
      if (imol > numberOfUniqueMoleculeTypes) call message(-1,"--top +def | No Molecule for atom ",i=iatm+1)
    end do

    ! Consolidate list of molecules into the global array
    allocate(listOfMolecules(numberOfMolecules))
    do idx=1,numberOfMolecules
      itmp = m(idx) % numberOfAtoms
      listOfMolecules(idx) % numberOfAtoms = itmp
      allocate(listOfMolecules(idx) % listOfAtoms(itmp))
      allocate(listOfMolecules(idx) % listOfLabels(itmp))
      call refsor(m(idx) % listOfAtoms(1:itmp))
      listOfMolecules(idx) % listOfAtoms = m(idx) % listOfAtoms(1:itmp)
      do jdx=1,itmp
        listOfMolecules(idx) % listOfLabels(jdx) = frame % lab(listOfMolecules(idx) % listOfAtoms(jdx))
      enddo
    enddo
    deallocate(m)

    call computeCovalentBondsFromMolecules()

    call finaliseTopology()

  end subroutine defineMoleculesFromLabels

  subroutine computeCovalentBondsFromMolecules()
    implicit none
    integer :: imol, idx, jdx, iatm, jatm
    character(cp) :: l1 ,l2
    real(8) :: rmax, dij(3), d2

    if (allocated(numberOfCovalentBondsPerAtom)) then
      deallocate(numberOfCovalentBondsPerAtom,listOfCovalentBondsPerAtom)
    end if
    
    allocate(numberOfCovalentBondsPerAtom(frame % natoms), source=0)
    allocate(listOfCovalentBondsPerAtom(neighmax,frame % natoms), source=0)
    do imol=1,numberOfMolecules
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        iatm = listOfMolecules(imol) % listOfAtoms(idx)
        l1 = listOfMolecules(imol) % listOfLabels(idx)
        do jdx=idx+1,listOfMolecules(imol) % numberOfAtoms
          jatm = listOfMolecules(imol) % listOfAtoms(jdx)
          l2 = listOfMolecules(imol) % listOfLabels(jdx)

          rmax = getMaximumBondLength(l1,l2,distanceScaling)
          dij = frame % pos(:,iatm) - frame % pos(:,jatm)
          d2 = computeDistanceSquaredPBC(dij)

          if (sqrt(d2) <= rmax) then
            numberOfCovalentBondsPerAtom(iatm) = numberOfCovalentBondsPerAtom(iatm) + 1
            listOfCovalentBondsPerAtom(numberOfCovalentBondsPerAtom(iatm),iatm) = jatm
            numberOfCovalentBondsPerAtom(jatm) = numberOfCovalentBondsPerAtom(jatm) + 1
            listOfCovalentBondsPerAtom(numberOfCovalentBondsPerAtom(jatm),jatm) = iatm
          end if

        end do
      end do
    end do

  end subroutine computeCovalentBondsFromMolecules

  subroutine computeTopology(a)
    implicit none
    type(actionTypeDef), target :: a
    logical :: brokenMolecules

    call associatePointers(a)

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      if (updateTopology) then
        if (allocated(atomToMoleculeIndex)) then
          deallocate(atomToMoleculeIndex)
          deallocate(listofmolecules)
          deallocate(uniqueMoleculesID)
          deallocate(numberOfCovalentBondsPerAtom)
          deallocate(listOfCovalentBondsPerAtom)
        end if
        
        if (allocated(listOfUniqueBonds)) then
          deallocate(listOfUniqueBonds)
        end if

        if (userDefinedMolecules) then
          call defineMoleculesFromLabels()
        else
          call defineMoleculesFromDistances()
          if (reorderAtoms) call reoderAtomsByMolecule()
        end if
        if (updateTopologyOnce) updateTopology = .false.
        call dumpScreenInfo()
        
        call checkUsedFlags(actionCommand)
      end if
      
      if (rebuildMolecules) then
        call reassembleAllMolecules()
      else 
        call checkForBrokenMolecules(brokenMolecules)
        if (brokenMolecules) call reassembleBrokenMolecules()
      end if

      call computeMoleculesCOM()
    end if

  end subroutine computeTopology

  subroutine reoderAtomsByMolecule()
    implicit none
    integer :: imol, jmol, idx, jdx, iatm, itmp
    real(8), allocatable, dimension(:,:) :: cartesianCoord
    real(8), allocatable, dimension(:) :: saveCharges
    character(cp), allocatable, dimension(:) :: saveLabels
    logical :: foundPerm
    logical, allocatable, dimension(:) :: removeMoleculeType
    integer, allocatable, dimension(:) :: newOrder
    integer, allocatable, dimension(:,:) :: newOrderAll
    integer, allocatable, dimension(:) :: newMoleculeType
    character(len=5) :: str

    allocate(cartesianCoord(3,frame % natoms), source=frame % pos(1:3,1:frame % natoms))
    allocate(saveLabels    (  frame % natoms), source=frame % lab(    1:frame % natoms))
    allocate(saveCharges   (  frame % natoms), source=frame % chg(    1:frame % natoms))
    ! cartesianCoord(1:3,1:frame % natoms) = frame % pos(1:3,1:frame % natoms)
    ! saveLabels    (    1:frame % natoms) = frame % lab(    1:frame % natoms)
    ! saveCharges   (    1:frame % natoms) = frame % chg(    1:frame % natoms)

    allocate(removeMoleculeType(numberOfUniqueMoleculeTypes), source=.false.)
    itmp = maxval(listOfMolecules(:) % numberOfAtoms)
    allocate(newOrder(itmp), source=0)
    allocate(newOrderAll(itmp,numberOfUniqueMoleculeTypes), source=0)
    
    allocate(newMoleculeType(numberOfUniqueMoleculeTypes))
    do idx=1,numberOfUniqueMoleculeTypes
      newMoleculeType(idx) = idx
    end do

    ! checking for reordering
    do idx=1,numberOfUniqueMoleculeTypes-1
        imol = uniqueMoleculesID(idx)

      ! Nothing to do for ions
      if (listOfMolecules(imol) % numberOfAtoms == 1) cycle
      
      do jdx=idx+1,numberOfUniqueMoleculeTypes
        
        ! cycle if it has already been removed
        if (removeMoleculeType(jdx)) cycle

        ! the molecule has a different number of atoms from the "Unique" reference
        jmol = uniqueMoleculesID(jdx)
        if ( listOfMolecules(imol) % numberOfAtoms /= listOfMolecules(jmol) % numberOfAtoms ) cycle

        ! check if the molecule has the same atoms types but in different order
        foundPerm = checkPermutations(listOfMolecules(imol) % numberOfAtoms, &
                                      listOfMolecules(imol) % listOfLabels,  &
                                      listOfMolecules(jmol) % listOfLabels,  &
                                      newOrder(1:listOfMolecules(imol) % numberOfAtoms) )

        if (foundPerm) then
          removeMoleculeType(jdx) = .true.
          newMoleculeType(jdx) = idx
          do itmp=1,listOfMolecules(jmol) % numberOfAtoms
            newOrderAll(itmp,jdx) = listOfMolecules(jmol) % listOfAtoms(newOrder(itmp))
          end do
        end if
        
      end  do 
    end  do
    
    if (count(removeMoleculeType) > 0) then
      do imol=1,numberOfMolecules
        idx = listOfMolecules(imol) % ID
        if (removeMoleculeType(idx)) then
          listOfMolecules(imol) % ID = newMoleculeType(idx)
          write(str,'(i0)') listOfMolecules(imol) % ID
          listOfMolecules(imol) % resname = "M"//trim(str) 
          listOfMolecules(imol) % listOfAtoms = newOrderAll(1:listOfMolecules(imol) % numberOfAtoms,idx)
        end if
      end do
    end if
    
    itmp = 0
    jdx = 0
    do while (jdx < numberOfUniqueMoleculeTypes)
      jdx = jdx + 1
      do imol=1,numberOfMolecules
        if (listOfMolecules(imol) % ID /= jdx) cycle
        do idx=1,listOfMolecules(imol) % numberOfAtoms
          iatm = listOfMolecules(imol) % listOfAtoms(idx)
          itmp = itmp + 1
          frame % pos(1:3,itmp) = cartesianCoord(1:3,iatm)
          frame % lab(    itmp) = saveLabels    (    iatm)
          frame % chg(    itmp) = saveCharges   (    iatm)
          atomToMoleculeIndex(itmp) = imol
          listOfMolecules(imol) % listOfAtoms(idx) = itmp
        end do
      end do
    end do
    call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
    call updateNeighboursList(.true.)
    call computeCovalentBondsFromMolecules()

  end subroutine reoderAtomsByMolecule

end module moduleMolecularTopology

subroutine computeConnectivity()
  use moduleVariables
  use moduleMessages 
  use moduleSystem 
  use moduleResizeArrays 
  implicit none

  integer :: n
  integer :: i, j, k
  integer :: idx, jdx, kdx, ldx
  integer :: iatm, jatm
  integer :: nShape(2)

  numberOfUniqueAngles = 0
  numberOfUniqueTorsions = 0
  numberOfUniqueOutOfPlane = 0
  
  ! Upper limit for angles
  n = numberOfUniqueBonds * 4
  allocate(listOfUniqueAngles(3,n))

  ! Loop over the atoms
  do i=1,frame % natoms
    ! Loop over the colvalent bonds of that atom
    do j=1,numberOfCovalentBondsPerAtom(i)-1
      jdx = listOfCovalentBondsPerAtom(j,i)
      ! Loop over the covalent bonds of the atoms bonded to i
      do k=j+1,numberOfCovalentBondsPerAtom(i)
        kdx = listOfCovalentBondsPerAtom(k,i)
        numberOfUniqueAngles = numberOfUniqueAngles + 1
        listOfUniqueAngles(1:3,numberOfUniqueAngles) = [jdx,i,kdx]
      enddo
    enddo
  enddo
  call resizeArray(listOfUniqueAngles,numberOfUniqueAngles)

  ! Skip if it's only ions in water
  if (any(listOfMolecules(:) % numberOfAtoms > 3)) then
    ! Upper limit for torsions
    n = numberOfUniqueBonds*4
    allocate(listOfUniqueTorsions(4,n))
    do i=1,numberOfUniqueAngles-1
      idx = listOfUniqueAngles(2,i)
      ! Skip if the molecule it too small
      if ( listOfMolecules(atomToMoleculeIndex(idx)) % numberOfAtoms < 3 ) cycle

      jdx = listOfUniqueAngles(1,i)
      kdx = listOfUniqueAngles(3,i)
      do j=i+1,numberOfUniqueAngles

        ! Cycle if the central atom is the same
        if (idx == listOfUniqueAngles(2,j))cycle

        ! Cycle if we are on different molecules
        ldx = listOfUniqueAngles(2,j)
        ! list of angles are assumed to be grouped by molecule
        if ( atomToMoleculeIndex(idx) /= atomToMoleculeIndex(ldx) ) exit

        if ( all([idx,jdx] == listOfUniqueAngles(1:2,j)) )then
          numberOfUniqueTorsions = numberOfUniqueTorsions + 1
          listOfUniqueTorsions(1:4,numberOfUniqueTorsions) = [kdx,listOfUniqueAngles(1:3,j)]

        else if (all([idx,jdx] == listOfUniqueAngles(3:2:-1,j)) ) then
          numberOfUniqueTorsions = numberOfUniqueTorsions + 1
          listOfUniqueTorsions(1:4,numberOfUniqueTorsions) = [listOfUniqueAngles(1:3,j),kdx]

        else if ( all([idx,kdx] == listOfUniqueAngles(1:2,j)) ) then
          numberOfUniqueTorsions = numberOfUniqueTorsions + 1
          listOfUniqueTorsions(1:4,numberOfUniqueTorsions) =  [jdx,listOfUniqueAngles(1:3,j)]

        else if (all([idx,kdx] == listOfUniqueAngles(3:2:-1,j)) ) then
          numberOfUniqueTorsions = numberOfUniqueTorsions + 1
          listOfUniqueTorsions(1:4,numberOfUniqueTorsions) =  [listOfUniqueAngles(1:3,j),jdx]
        end if
      enddo
    enddo
    call resizeArray(listOfUniqueTorsions,numberOfUniqueTorsions)

    ! Upper limit for improper torsions
    n = numberOfAtoms
    allocate(listOfUniqueOutOfPlane(4,n))
    do i=1,frame % natoms
      if (numberOfCovalentBondsPerAtom(i) == 3) then
        numberOfUniqueOutOfPlane = numberOfUniqueOutOfPlane + 1
        listOfUniqueOutOfPlane(1:4,numberOfUniqueOutOfPlane) = [i,listOfCovalentBondsPerAtom(1:3,i)]
      end if
    enddo
    call resizeArray(listOfUniqueOutOfPlane,numberOfUniqueOutOfPlane)
  end if

  ! Compute list of 1-n interactions
  nShape = 2*shape(listOfCovalentBondsPerAtom)
  allocate(numberOf_13_InteractionsPerAtom(frame % natoms), source=0)
  allocate(numberOf_14_InteractionsPerAtom(frame % natoms), source=0)
  allocate(listOf_13_InteractionsAtom(nShape(1),frame % natoms))
  allocate(listOf_14_InteractionsAtom(nShape(1),frame % natoms))

  do i=1,numberOfUniqueAngles
    iatm = listOfUniqueAngles(1,i)
    jatm = listOfUniqueAngles(3,i)
    numberOf_13_InteractionsPerAtom(iatm) = numberOf_13_InteractionsPerAtom(iatm) + 1
    numberOf_13_InteractionsPerAtom(jatm) = numberOf_13_InteractionsPerAtom(jatm) + 1
    listOf_13_InteractionsAtom(numberOf_13_InteractionsPerAtom(iatm),iatm) = jatm
    listOf_13_InteractionsAtom(numberOf_13_InteractionsPerAtom(jatm),jatm) = iatm
  enddo

  do i=1,numberOfUniqueTorsions
    iatm = listOfUniqueTorsions(1,i)
    jatm = listOfUniqueTorsions(4,i)
    numberOf_14_InteractionsPerAtom(iatm) = numberOf_14_InteractionsPerAtom(iatm) + 1
    numberOf_14_InteractionsPerAtom(jatm) = numberOf_14_InteractionsPerAtom(jatm) + 1
    listOf_14_InteractionsAtom(numberOf_14_InteractionsPerAtom(iatm),iatm) = jatm
    listOf_14_InteractionsAtom(numberOf_14_InteractionsPerAtom(jatm),jatm) = iatm
  enddo

  call message(0,"Connectivity summary")
  call message(0,"...Number of bonds",i=numberOfUniqueBonds)
  call message(0,"...Number of angles",i=numberOfUniqueAngles)
  call message(0,"...Number of dihedrals",i=numberOfUniqueTorsions)
  call message(1,"...Number of impropers",i=numberOfUniqueOutOfPlane)

end subroutine computeConnectivity
