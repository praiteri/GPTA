module addDummyModule

  use moduleVariables
  use moduleSystem
  use moduleMessages

  public :: addDummyParticles

  private
  
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  integer, pointer :: numberOfSites
  integer, pointer :: numberOfNewMolecules

  type dummyStructure
    character(fp) :: resname
    integer :: numberOfAtoms
    integer, allocatable, dimension(:) :: listOfAtoms
    character(cp) :: label
    real(real64) :: distance
  end type
  type(dummyStructure) :: dummy

  contains

  subroutine addDummyParticles(a)
    implicit none
    type(actionTypeDef), target :: a

    call associatePointers(a)

    if (a % actionInitialisation) call initialiseAction(a)

    if (frameReadSuccessfully) then
      if (firstAction) then

        call extractMoleculesIndices(a)

        call extendFrame(numberOfNewMolecules)
        
        block ! add molecule to global array
          integer :: i, iatm
          integer :: currentMolecules
                    
          currentMolecules = numberOfMolecules
          call extendMolecules(numberOfNewMolecules)
          iatm = frame % natoms
          do i=currentMolecules+1, currentMolecules+numberOfNewMolecules
            listOfMolecules(i) % ID = numberOfUniqueMoleculeTypes + 1
            listOfMolecules(i) % numberOfAtoms = 1
            listOfMolecules(i) % resname = "COM"

            allocate(listOfMolecules(i) % listOfAtoms(1))
            allocate(listOfMolecules(i) % listOfLabels(1) , source="XC  ")
            iatm = iatm + 1
            listOfMolecules(i) % listOfAtoms(1) = iatm

            atomToMoleculeIndex(iatm) = i
          end do
        end block

        call dumpScreenInfo()
        
        call checkUsedFlags(actionCommand)
        firstAction = .false.

      end if
    
      call computeAction(a)
    end if

  end subroutine addDummyParticles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand        => a % actionDetails
    firstAction          => a % firstAction

    numberOfSites        => a % integerVariables(1)
    numberOfNewMolecules => a % integerVariables(2)

  end subroutine associatePointers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialiseAction(a)
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a

    character(len=STRLEN), allocatable, dimension(:) :: moleculeID
    character(len=STRLEN), allocatable, dimension(:) :: atomsIndices

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.

    call assignFlagValue(actionCommand,"+mol",moleculeID)
    numberOfSites = size(moleculeID)

    if (numberOfSites > 1) then
      call message(-1,"--add dummy - cannot add more than one dummy site at a time")
    end if
    read(moleculeID(1),*) dummy % resname 

    call assignFlagValue(actionCommand,"+i ",atomsIndices)
    dummy % numberOfAtoms = size(atomsIndices)
    if (dummy % numberOfAtoms /= 3) then
      call message(-1,"--add dummy - dummy site requires three atoms")
    end if

    allocate(dummy % listOfAtoms(dummy % numberOfAtoms))

    block
      integer :: i
      do i=1,dummy % numberOfAtoms
        read(atomsIndices(i),*) dummy % listOfAtoms(i)
      end do
    end block

    call assignFlagValue(actionCommand,"+dist ",dummy % distance,0.0_real64)
    
    call assignFlagValue(actionCommand,"+l ",dummy % label,"MW")

  end subroutine initialiseAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: idx, imol, n0
    integer :: i, j, k
    integer :: ii, jj ,kk
    real(real64), dimension(3) :: vec

    ii = dummy % listOfAtoms(1)
    jj = dummy % listOfAtoms(2)
    kk = dummy % listOfAtoms(3)

    n0 = frame % natoms
    do idx=1,numberOfNewMolecules
      imol = a % localIndices(idx)
      
      i = listOfMolecules(imol) % listOfAtoms(ii)
      j = listOfMolecules(imol) % listOfAtoms(jj)
      k = listOfMolecules(imol) % listOfAtoms(kk)

      vec(1:3) = frame % pos(1:3,j) + frame % pos(1:3,k) - 2.0_real64 * frame % pos(1:3,i)
      vec = vec / sqrt(sum(vec*vec))
      
      frame % pos(1:3 , n0+idx) = frame % pos(1:3,i) + vec(1:3) * dummy % distance
      frame % lab(n0+idx) = dummy % label
    end do

    frame % natoms = frame % natoms + numberOfNewMolecules

  end subroutine computeAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Add dummy particle to molecules")

    call message(0,"...Molecule's name",str=dummy % resname)
    call message(0,"...Number of molecules found",i=numberOfMolecules)
    call message(0,"...Indices of the atoms",iv=dummy % listOfAtoms)
    call message(0,"...Distance from the first atom",r=dummy % distance)

    call message(2)

  end subroutine dumpScreenInfo

  subroutine extractMoleculesIndices(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: imol, idx

    numberOfNewMolecules = 0
    do imol=1,numberOfMolecules
      if (listOfMolecules(imol) % resname == dummy % resname) numberOfNewMolecules = numberOfNewMolecules + 1
    end do
    allocate(a % localIndices(numberOfNewMolecules))

    idx = 0
    do imol=1,numberOfMolecules
      if (listOfMolecules(imol) % resname == dummy % resname) then
        idx = idx + 1
        a % localIndices(idx) = imol
      end if
    end do

  end subroutine

end module addDummyModule
