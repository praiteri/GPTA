module addFileModule

  use moduleVariables
  use moduleSystem
  use moduleMessages

  public :: addSystemFromFile

  private
  
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  integer, pointer :: numberOfFiles
  character(len=STRLEN), pointer, dimension(:) :: fileNames

  type newMolecule
    integer :: natoms
    character(cp), allocatable, dimension(:) :: lab
    real(real64), allocatable, dimension(:,:) :: pos
    real(real64), allocatable, dimension(:) :: chg
    real(real64) :: centre(3)
  end type 
  type(newMolecule), allocatable, dimension(:) :: localMolecules

contains

  subroutine addSystemFromFile(a)
    use moduleFiles
    use moduleRead
    implicit none
    type(actionTypeDef), target :: a

    type(fileTypeDef), allocatable, dimension(:) :: inputFile
    logical :: frameRead
    integer :: n, iFile
    real(real64) :: hmat(3,3)

    call associatePointers(a)

    if (a % actionInitialisation) call initialiseAction(a)

    if (frameReadSuccessfully) then
      ! if (firstAction) then

        allocate(localMolecules(numberOfFiles))
        allocate(inputFile(numberOfFiles))
        
        do iFile=1,numberOfFiles
          call initialiseFile(inputFile(iFile), fileNames(iFile))
        
          ! -> get number of atoms
          call getNumberOfAtoms(inputFile(iFile), n, hmat)
          rewind(inputFile(iFile) % funit)
        
          call createSystemArrays(a % localFrame, n)
          call readCoordinates(frameRead, a % localFrame, inputFile(iFile))

          localMolecules(iFile) % natoms = a % localFrame % natoms

          allocate(localMolecules(iFile) % lab(  localMolecules(iFile) % natoms))
          allocate(localMolecules(iFile) % pos(3,localMolecules(iFile) % natoms))
          allocate(localMolecules(iFile) % chg(  localMolecules(iFile) % natoms))
          localMolecules(iFile) % lab = a % localFrame % lab
          localMolecules(iFile) % pos = a % localFrame % pos
          localMolecules(iFile) % chg = a % localFrame % chg

          call deleteSystemArrays(a % localFrame)
        end do
        
        ! extend coordinates arrays
        n = sum(localMolecules(1:numberOfFiles) % natoms)
        call extendFrame(n) 

        !extend molecules arrays
        call extendMolecules(n)

        call dumpScreenInfo()

        call checkUsedFlags(actionCommand)
        firstAction = .false.

        ! end if
      call computeAction(a)
        
      deallocate(localMolecules)
      deallocate(inputFile)

    end if

  end subroutine addSystemFromFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand            => a % actionDetails
    firstAction              => a % firstAction

    numberOfFiles            => a % integerVariables(1)
    fileNames(1:)            => a % stringVariables(1:)

  end subroutine associatePointers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialiseAction(a)
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a

    character(len=STRLEN), allocatable, dimension(:) :: localString

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.
    
    call assignFlagValue(actionCommand,"+f",localString)
    numberOfFiles = size(localString)

    fileNames(1:numberOfFiles) = localString(1:numberOfFiles)
    deallocate(localString)

  end subroutine initialiseAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: iFile, m, n, i
    integer :: currentAtoms
    integer :: newAtoms

    n = size(frame % lab)
    m = frame % natoms

    call createSystemArrays(a % localFrame, n)
    
    a % localFrame % lab(1:m)   = frame % lab(1:m)  
    a % localFrame % pos(:,1:m) = frame % pos(:,1:m)
    a % localFrame % chg(1:m)   = frame % chg(1:m)  
    
    newAtoms = n - m

    ! add the new atoms to the frame at the end
    currentAtoms = 0
    do iFile=1,numberOfFiles
      do i=1,localMolecules(iFile) % natoms
        currentAtoms = currentAtoms + 1
        frame % pos(1:3,currentAtoms) = localMolecules(iFile) % pos(1:3,i) 
        frame % lab(currentAtoms) = localMolecules(iFile) % lab(i) 
        frame % chg(currentAtoms) = localMolecules(iFile) % chg(i) 
      end do
    end do

    ! newAtoms = currentAtoms
    do i=1,frame % natoms
      currentAtoms = currentAtoms + 1
      frame % pos(1:3,currentAtoms) = a % localFrame % pos(1:3,i) 
      frame % lab(currentAtoms) = a % localFrame % lab(i) 
      frame % chg(currentAtoms) = a % localFrame % chg(i) 
    end do

    frame % natoms = currentAtoms

  end subroutine computeAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dumpScreenInfo()
    implicit none
    integer :: i, j
    character(len=STRLEN) :: str
 
    call message(0,"Add atoms from file(s)")
    do i=1,numberOfFiles
      call message(0,"...File name",str=fileNames(i))
    end do

    call message(2)

  end subroutine dumpScreenInfo

end module addFileModule
