!disclaimer
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
  logical, pointer :: checkForBroken
  logical, pointer :: reorderAtoms
  logical, pointer :: groupAtoms
  logical, pointer :: userDefinedMolecules
  logical, pointer :: verboseOutput

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
    checkForBroken         => a % logicalVariables(4)
    reorderAtoms           => a % logicalVariables(5)
    groupAtoms             => a % logicalVariables(6)
    verboseOutput          => a % logicalVariables(7)
    
    userDefinedMolecules   => a % logicalVariables(8)
    moleculesLabels(1:)    => a % stringVariables(1:)

  end subroutine associatePointers

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    character(STRLEN), allocatable, dimension(:) :: tmpLabels
    character(STRLEN) :: flagString
    logical :: lflag
    real(real64) :: rcut

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0_real64    
  
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

    call assignFlagValue(actionCommand,"+rebuild ",rebuildMolecules,.false.)
    
    call assignFlagValue(actionCommand,"+check ",checkForBroken,.false.)
    
    call assignFlagValue(actionCommand,"+reorder ",reorderAtoms,.false.)

    call assignFlagValue(actionCommand,"+group ",groupAtoms,.false.)

    call assignFlagValue(actionCommand,"+v ",verboseOutput,.false.)

    call assignFlagValue(actionCommand,"+def",userDefinedMolecules,.false.)
    if (userDefinedMolecules) then
      call extractFlag(actionCommand,"+def",flagString)
      call parse(flagString," ",moleculesLabels,numberOfUniqueMoleculeTypes)
    end if

    call assignFlagValue(actionCommand,"+rcut",rcut,a % cutoffNeighboursList)
    a % cutoffNeighboursList = rcut

    updateTopology = .true.
    tallyExecutions = 0

  end subroutine initialiseAction

  subroutine dumpScreenInfo(lflag)
    implicit none
    logical, intent(in) :: lflag
    integer :: idx, jdx, kdx 
    integer :: iatm, n
    character(len=30) :: str, str1
    real(real64) :: molecularCharge, molecularMass

    call message(0,"Topology information")
    call message(0,"...Number of molecules found",i=numberOfMolecules)
    call message(0,"...Number of unique molecules found",i=numberOfUniqueMoleculeTypes)
    do idx=1,numberOfUniqueMoleculeTypes
      kdx = uniqueMoleculesID(idx)
      write(str1,'(i0)')idx

      n = count(listOfMolecules(:) % ID == idx)
      if (n == 0) cycle
      
      call message(3)
      call message(0,"....Molecule type "//trim(listOfMolecules(kdx) % resname))

      call message(0,"......Number of molecules",i=n)
      call message(0,"......Number of atoms in molecule",i=listOfMolecules(kdx) % numberOfAtoms)

      str = trim(listOfMolecules(kdx) % listOfLabels(1))
      do jdx=2,listOfMolecules(kdx) % numberOfAtoms
        str = trim(str)//","//trim(listOfMolecules(kdx) % listOfLabels(jdx))
      enddo
      call message(0,"......Atoms types in molecule",str=str)

      call message(0,"......Number of atoms in molecule",i=listOfMolecules(kdx) % numberOfAtoms)

      molecularCharge = 0.0_real64
      molecularMass = 0.0_real64
      do jdx=1,listOfMolecules(kdx) % numberOfAtoms
        iatm = listOfMolecules(kdx) % listOfAtoms(jdx)
        molecularCharge = molecularCharge + frame % chg(iatm)
        molecularMass = molecularMass + getElementMass(frame % lab(iatm))
        enddo
      call message(0,"......Molecular charge",r=molecularCharge)
      if (frame % volume > 1e-6_real64) &
      call message(0,"......Concentration [M]",r=1660.5_real64*n/frame % volume)
      call message(0,"......Concentration [%w]",r=100*n*molecularMass/totalMass)
      call message(0,"......Concentration [%m]",r=100*dble(n)/numberOfMolecules)

    enddo
    if (lflag) then
      call message(3)
      call message(0,"Connectivity summary")
      call message(0,"  Number of bonds",i=numberOfUniqueBonds)
      if (numberOfUniqueBonds > 0) call writeUniqueCovalentTypes(1)
      call message(0,"  Number of angles",i=numberOfUniqueAngles)
      if (numberOfUniqueAngles > 0) call writeUniqueCovalentTypes(2)
      call message(0,"  Number of torsions",i=numberOfUniqueTorsions)
      if (numberOfUniqueTorsions > 0) call writeUniqueCovalentTypes(3)
      call message(0,"  Number of impropers",i=numberOfUniqueOutOfPlane)
      if (numberOfUniqueOutOfPlane > 0) call writeUniqueCovalentTypes(4)
    end if
    call message(2)

  end subroutine dumpScreenInfo

  subroutine defineMoleculesFromDistances()
    implicit none
    ! type(actionTypeDef), target :: a

    integer :: iatm, jatm, ineigh, idx, jdx, itmp
    character(cp) :: l1 ,l2
    real(real64) :: rmax

    logical, allocatable, dimension(:) :: atomsUsed
    type(moleculeTypeDef), allocatable, dimension(:) :: m

    if (allocated(numberOfCovalentBondsPerAtom)) then
      deallocate(numberOfCovalentBondsPerAtom)
      deallocate(listOfCovalentBondsPerAtom)
      deallocate(atomToMoleculeIndex)
    end if    

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
      allocate(m(idx) % listOfAtoms(200))
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
    integer :: isWater, indices(4)
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

    ! identity water molecules
    isWater=-1
    do idx=1,numberOfUniqueMoleculeTypes
      itmp = uniqueMoleculesID(idx)
      if (listOfMolecules(itmp) % numberOfAtoms == 3) then
        do jdx=1,listOfMolecules(itmp) % numberOfAtoms
          indices(jdx) = getAtomicNumber(listOfMolecules(itmp) % listOfLabels(jdx))
        end do
        if (count(indices(1:3)==8)==1 .and. count(indices(1:3)==1)==2) isWater = idx
      else if (listOfMolecules(itmp) % numberOfAtoms == 4) then
        do jdx=1,listOfMolecules(itmp) % numberOfAtoms
          indices(jdx) = getAtomicNumber(listOfMolecules(itmp) % listOfLabels(jdx))
        end do
        if (count(indices(1:4)==8)==1 .and. count(indices(1:4)==1)==2 .and. count(indices(1:4)==0)==1 ) isWater = idx
      end if
    end do

    ! Assign a name to the molecules
    do idx=1,numberOfMolecules
      if (isWater==listOfMolecules(idx) % ID) then
        listOfMolecules(idx) % resname = "HOH"
      else
        write(str,'(i0)') listOfMolecules(idx) % ID
        listOfMolecules(idx) % resname = "M"//trim(str) 
      end if
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
      if (imol > numberOfUniqueMoleculeTypes) then
        write(0,*)tmplab(1:midx(imol-1))
        write(0,*)mlab(1:midx(imol-1),imol-1)
        write(0,*)tmplab(1:midx(imol-1)) == mlab(1:midx(imol-1),imol-1)        
        call message(1,"--top +def | No Molecule for atom ",i=iatm+1)
      end if
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
    real(real64) :: rmax, dij(3), d2

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

#ifdef DEBUG
    write(0,*) " --> Entering computeTopology <-- "
#endif

    call associatePointers(a)

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      if (updateTopology) then
        if (allocated(atomToMoleculeIndex         ))  deallocate(atomToMoleculeIndex)
        if (allocated(listofmolecules             ))  deallocate(listofmolecules)
        if (allocated(uniqueMoleculesID           ))  deallocate(uniqueMoleculesID)
        if (allocated(numberOfCovalentBondsPerAtom))  deallocate(numberOfCovalentBondsPerAtom)
        if (allocated(listOfCovalentBondsPerAtom  ))  deallocate(listOfCovalentBondsPerAtom)
        if (allocated(listOfUniqueBonds           ))  deallocate(listOfUniqueBonds)

        
        if (allocated(numberOfCovalentBondsPerAtom)) deallocate(numberOfCovalentBondsPerAtom)
        if (allocated(listOfCovalentBondsPerAtom  )) deallocate(listOfCovalentBondsPerAtom)
        if (allocated(atomToMoleculeIndex         )) deallocate(atomToMoleculeIndex)
    
        if (userDefinedMolecules) then
          call defineMoleculesFromLabels()
        else
          call defineMoleculesFromDistances()
          if (reorderAtoms) call reorderAtomsByMolecule()
          if (groupAtoms) call groupAtomsByMolecule()
        end if
        if (updateTopologyOnce) updateTopology = .false.
        
        call computeConnectivity() 
        call dumpScreenInfo(verboseOutput)

        call checkUsedFlags(actionCommand)
      end if
      
      if (rebuildMolecules) call reassembleAllMolecules()

      if (checkForBroken) then
        call checkForBrokenMolecules(brokenMolecules)
        if (brokenMolecules) call reassembleBrokenMolecules()
      end if

      call computeMoleculesCOM(0)

    end if

  end subroutine computeTopology

  subroutine reorderAtomsByMolecule()
    implicit none
    integer :: imol, jmol, idx, jdx, iatm, itmp
    real(real64), allocatable, dimension(:,:) :: cartesianCoord
    real(real64), allocatable, dimension(:) :: saveCharges
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
                                      newOrderAll(1:listOfMolecules(imol) % numberOfAtoms,jdx))

        if (foundPerm) then
          removeMoleculeType(jdx) = .true.
          newMoleculeType(jdx) = idx
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
          do itmp=1,listOfMolecules(imol) % numberOfAtoms
            newOrder(itmp) = listOfMolecules(imol) % listOfAtoms(newOrderAll(itmp,idx))
          end do

          listOfMolecules(imol) % listOfAtoms = newOrder(1:listOfMolecules(imol) % numberOfAtoms)
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

  end subroutine reorderAtomsByMolecule

  subroutine groupAtomsByMolecule()
    implicit none
    integer :: imol, jmol, idx, jdx, iatm, itmp, jatm
    real(real64), allocatable, dimension(:,:) :: cartesianCoord
    real(real64), allocatable, dimension(:) :: saveCharges
    character(cp), allocatable, dimension(:) :: saveLabels
    logical, allocatable, dimension(:) :: lused
    character(len=5) :: str

    allocate(cartesianCoord(3,frame % natoms), source=frame % pos(1:3,1:frame % natoms))
    allocate(saveLabels    (  frame % natoms), source=frame % lab(    1:frame % natoms))
    allocate(saveCharges   (  frame % natoms), source=frame % chg(    1:frame % natoms))
    allocate(lused         (  frame % natoms), source=.false.)

    jatm = 0
    do imol=1,numberOfMolecules
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        iatm = listOfMolecules(imol) % listOfAtoms(idx)
        jatm = jatm + 1
        frame % pos(1:3,jatm) = cartesianCoord(1:3,iatm)
        frame % lab(    jatm) = saveLabels    (    iatm)
        frame % chg(    jatm) = saveCharges   (    iatm)
        lused(iatm) = .true.
      end do
    end do

    call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
    call updateNeighboursList(.true.)
    call runInternalAction("topology","+update")
    ! call computeCovalentBondsFromMolecules()

  end subroutine groupAtomsByMolecule

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
  if (allocated(listOfUniqueAngles)) deallocate(listOfUniqueAngles)
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
    n = numberOfUniqueBonds*6
    if (allocated(listOfUniqueTorsions)) deallocate(listOfUniqueTorsions)
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
    if (allocated(listOfUniqueOutOfPlane)) deallocate(listOfUniqueOutOfPlane)

    n = frame % natoms
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
  if (allocated(numberOf_13_InteractionsPerAtom)) then
    deallocate(numberOf_13_InteractionsPerAtom)
    deallocate(numberOf_14_InteractionsPerAtom)
    deallocate(listOf_13_InteractionsAtom)
    deallocate(listOf_14_InteractionsAtom)
  end if

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

end subroutine computeConnectivity

subroutine writeUniqueCovalentTypes(ntype)
  use moduleVariables
  use moduleMessages 
  use moduleSystem 

  integer, intent(in) :: ntype

  integer :: i, j, nlist
  character(cp), allocatable, dimension(:,:) :: list
  character(cp), dimension(4) :: labels
  character(len=100) :: str

  if (ntype == 1) then
    allocate(list(2,numberOfUniqueBonds))
    nlist = 1
    labels(1:2) = [frame%lab(listOfUniqueBonds(1,nlist)), \
                   frame%lab(listOfUniqueBonds(2,nlist))]
    list(1:2,nlist) = labels(1:2)
    a : do i=2,numberOfUniqueBonds
      labels(1:2) = [frame%lab(listOfUniqueBonds(1,i)), \
                     frame%lab(listOfUniqueBonds(2,i))]

      do j=1,nlist
        if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) ) cycle a
        if ( (list(1,j) == labels(2)) .and. (list(2,j) == labels(1)) ) cycle a
      end do
      nlist = nlist + 1 
      labels(1:2) = [frame%lab(listOfUniqueBonds(1,i)), \
                     frame%lab(listOfUniqueBonds(2,i))]
      list(1:2,nlist) = labels(1:2)
    end do a

    call message(0,"    Number of unique bond types",i=nlist)
    do i=1,nlist
      str = list(1,i)//" - "//list(2,i)
      call message(0,"      "//trim(str))
    end do

  else if (ntype == 2) then
    allocate(list(3,numberOfUniqueAngles))

    nlist = 1
    labels(1:3) = [frame%lab(listOfUniqueAngles(1,nlist)), \
                   frame%lab(listOfUniqueAngles(2,nlist)), \
                   frame%lab(listOfUniqueAngles(3,nlist))]
    list(1:3,nlist) = labels(1:3)
    b : do i=2,numberOfUniqueAngles
      labels(1:3) = [frame%lab(listOfUniqueAngles(1,i)), \
                     frame%lab(listOfUniqueAngles(2,i)), \
                     frame%lab(listOfUniqueAngles(3,i))]

      do j=1,nlist
        if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(3)) ) cycle b
        if ( (list(1,j) == labels(3)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(1)) ) cycle b
      end do
      nlist = nlist + 1 
      labels(1:3) = [frame%lab(listOfUniqueAngles(1,i)), \
                     frame%lab(listOfUniqueAngles(2,i)), \
                     frame%lab(listOfUniqueAngles(3,i))]
      list(1:3,nlist) = labels(1:3)
    end do b

    call message(0,"    Number of unique angle types",i=nlist)
    do i=1,nlist
      str = list(1,i)//" - "//list(2,i)//" - "//list(3,i)
      call message(0,"      "//trim(str))
    end do

  else if (ntype == 3) then
    allocate(list(4,numberOfUniqueTorsions))

    nlist = 1
    labels(1:4) = [frame%lab(listOfUniqueTorsions(1,nlist)), \
                   frame%lab(listOfUniqueTorsions(2,nlist)), \
                   frame%lab(listOfUniqueTorsions(3,nlist)), \
                   frame%lab(listOfUniqueTorsions(4,nlist))]
    list(1:4,nlist) = labels(1:4)
    c : do i=2,numberOfUniqueTorsions
      labels(1:4) = [frame%lab(listOfUniqueTorsions(1,i)), \
                     frame%lab(listOfUniqueTorsions(2,i)), \
                     frame%lab(listOfUniqueTorsions(3,i)), \
                     frame%lab(listOfUniqueTorsions(4,i))]

      do j=1,nlist
        if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(3)) .and. (list(4,j) == labels(4)) ) cycle c
        if ( (list(1,j) == labels(4)) .and. (list(2,j) == labels(3)) .and. (list(3,j) == labels(2)) .and. (list(4,j) == labels(1)) ) cycle c
      end do
      nlist = nlist + 1 
      labels(1:4) = [frame%lab(listOfUniqueTorsions(1,i)), \
                     frame%lab(listOfUniqueTorsions(2,i)), \
                     frame%lab(listOfUniqueTorsions(3,i)), \
                     frame%lab(listOfUniqueTorsions(4,i))]
      list(1:4,nlist) = labels(1:4)
    end do c

    call message(0,"    Number of unique torsion types",i=nlist)
    do i=1,nlist
      str = list(1,i)//" - "//list(2,i)//" - "//list(3,i)//" - "//list(4,i)
      call message(0,"      "//trim(str))
    end do

  else if (ntype == 4) then
    allocate(list(4,numberOfUniqueOutOfPlane))

    nlist = 1
    labels(1:4) = [frame%lab(listOfUniqueOutOfPlane(1,nlist)), \
                   frame%lab(listOfUniqueOutOfPlane(2,nlist)), \
                   frame%lab(listOfUniqueOutOfPlane(3,nlist)), \
                   frame%lab(listOfUniqueOutOfPlane(4,nlist))]
    list(1:4,nlist) = labels(1:4)
    d : do i=2,numberOfUniqueOutOfPlane
      labels(1:4) = [frame%lab(listOfUniqueOutOfPlane(1,i)), \
                     frame%lab(listOfUniqueOutOfPlane(2,i)), \
                     frame%lab(listOfUniqueOutOfPlane(3,i)), \
                     frame%lab(listOfUniqueOutOfPlane(4,i))]

      do j=1,nlist
        if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(3)) .and. (list(4,j) == labels(4)) ) cycle d
        if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(4)) .and. (list(4,j) == labels(3)) ) cycle d
        if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(3)) .and. (list(3,j) == labels(2)) .and. (list(4,j) == labels(4)) ) cycle d
        if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(3)) .and. (list(3,j) == labels(4)) .and. (list(4,j) == labels(2)) ) cycle d
        if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(4)) .and. (list(3,j) == labels(2)) .and. (list(4,j) == labels(3)) ) cycle d
        if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(4)) .and. (list(3,j) == labels(3)) .and. (list(4,j) == labels(2)) ) cycle d
      end do
      nlist = nlist + 1 
      labels(1:4) = [frame%lab(listOfUniqueOutOfPlane(1,i)), \
                     frame%lab(listOfUniqueOutOfPlane(2,i)), \
                     frame%lab(listOfUniqueOutOfPlane(3,i)), \
                     frame%lab(listOfUniqueOutOfPlane(4,i))]
      list(1:4,nlist) = labels(1:4)
    end do d

    call message(0,"    Number of unique improper types",i=nlist)
    do i=1,nlist
      str = list(1,i)//" - "//list(2,i)//" - "//list(3,i)//" - "//list(4,i)
      call message(0,"      "//trim(str))
    end do

  end if
  
end subroutine writeUniqueCovalentTypes
