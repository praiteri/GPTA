module addCentresModule

  use moduleVariables
  use moduleSystem
  use moduleMessages

  public :: addCentresParticles

  private
  
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  integer, pointer :: numberOfSites
  integer, pointer :: numberOfNewMolecules

  type centresStructure
    character(fp) :: resname
    integer :: numberOfAtoms = 0
    integer, allocatable, dimension(:) :: listOfAtoms
    character(cp) :: label
  end type
  type(centresStructure) :: centres

  contains

  subroutine addCentresParticles(a)
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

  end subroutine addCentresParticles

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
      call message(-1,"--add centres - cannot add more than one centres site at a time")
    end if
    read(moleculeID(1),*) centres % resname 

    call assignFlagValue(actionCommand,"+i ",atomsIndices)
    if (.not. allocated(atomsIndices)) then
      centres % numberOfAtoms = 0
    else
      centres % numberOfAtoms = size(atomsIndices)
    end if

    if (centres % numberOfAtoms > 0) then
      allocate(centres % listOfAtoms(centres % numberOfAtoms))
      
      block
        integer :: i
        do i=1,centres % numberOfAtoms
          read(atomsIndices(i),*) centres % listOfAtoms(i)
        end do
      end block
    end if

    call assignFlagValue(actionCommand,"+l ",centres % label,"XC")

  end subroutine initialiseAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: idx, jdx, kdx, imol, n0, iatm
    real(real64), dimension(3) :: xcom

    n0 = frame % natoms

    if (centres % numberOfAtoms == 0) then
      do idx=1,numberOfNewMolecules
        imol = a % localIndices(idx)
        frame % pos(1:3 , n0+idx) = listOfMolecules(imol) % centreOfMass
        frame % lab(n0+idx) = centres % label
      end do
      
    else
      
      do idx=1,numberOfNewMolecules
        imol = a % localIndices(idx)
        xcom = 0.0_real64
        do jdx=1,centres % numberOfAtoms
          kdx = centres % listOfAtoms(jdx)
          iatm = listOfMolecules(imol) % listOfAtoms(kdx)
          xcom = xcom + frame % pos(1:3,iatm)
        end do
        frame % pos(1:3 , n0+idx) = xcom / centres % numberOfAtoms
        frame % lab(n0+idx) = centres % label
      end do

    end if

    frame % natoms = frame % natoms + numberOfNewMolecules

  end subroutine computeAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dumpScreenInfo()
    implicit none

    call message(0,"Add centres of mass to molecules")

    call message(0,"...Molecule's name",str=centres % resname)
    call message(0,"...Number of molecules found",i=numberOfMolecules)
    ! if (centres % numberOfAtoms > 0) \
      ! call message(0,"...Indices of the atoms",iv=centres % listOfAtoms)

    call message(2)

  end subroutine dumpScreenInfo

  subroutine extractMoleculesIndices(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: imol, idx

    numberOfNewMolecules = 0
    do imol=1,numberOfMolecules
      if (listOfMolecules(imol) % resname == centres % resname) numberOfNewMolecules = numberOfNewMolecules + 1
    end do
    allocate(a % localIndices(numberOfNewMolecules))

    idx = 0
    do imol=1,numberOfMolecules
      if (listOfMolecules(imol) % resname == centres % resname) then
        idx = idx + 1
        a % localIndices(idx) = imol
      end if
    end do

  end subroutine
  
end module addCentresModule
