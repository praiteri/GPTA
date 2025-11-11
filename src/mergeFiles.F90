!disclaimer
module moduleMergeFiles
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleFiles
  use moduleMessages
  use moduleRead
  use moduleDistances
  use moduleNeighbours

  implicit none

  public :: MergeFiles, MergeFilesHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions
  
  integer, pointer :: numberOfFiles
  character(len=STRLEN), pointer, dimension(:) :: fileNames
  character(len=STRLEN), pointer:: stack
  real(real64), pointer :: shift

contains

subroutine MergeFilesHelp()
  use moduleMessages
  implicit none
  call message(0,"This action merges coordinates from different files.")
  call message(0,"Examples:")
  call message(0,"  gpta.x --i coord.pdb --merge +f f1.pdb,f2.pdb")
end subroutine MergeFilesHelp

subroutine MergeFiles(a)
  implicit none
  type(actionTypeDef), target :: a

  call associatePointers(a)
  if (a % actionInitialisation) then
    call initialiseAction(a)
    return
  end if

  if (frameReadSuccessfully) then
    tallyExecutions = tallyExecutions + 1

    ! Atoms' selection
    if (firstAction) then
      ! dump info about the action on the screen
      call dumpScreenInfo()

      ! Throw a warning for unused flags
      call checkUsedFlags(actionCommand)
      firstAction = .false.

    end if

    call computeAction(a)
  end if

  if (endOfCoordinatesFiles) then
    call finaliseAction()
  end if 

end subroutine MergeFiles

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions

    numberOfFiles        => a % integerVariables(1)
    stack                => a % stringVariables(1)
    fileNames(1:)        => a % stringVariables(2:)
    shift                => a % doubleVariables(1)

  end subroutine associatePointers

  subroutine initialiseAction(a)

    implicit none
    type(actionTypeDef), target :: a
    character(len=STRLEN), allocatable, dimension(:) :: localString

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.2d0

    tallyExecutions = 0

    call assignFlagValue(actionCommand,"+f",localString)
    numberOfFiles = size(localString)

    if (numberOfFiles == 0) then
      call message(0,"Please provide a list of files to be processed")
      call message(0,"Example: gpta.x --i file1 --merge +f file2 file3")
      call message(-1,"No input files provided")
    end if
    fileNames(1:numberOfFiles) = localString(1:numberOfFiles)
    deallocate(localString)

    call assignFlagValue(actionCommand,"+stack ", stack, "NONE")
    call assignFlagValue(actionCommand,"+shift ", shift, 0.d0)

    call message(0,"Number of files to be merged",i=numberOfFiles)
    block
      integer :: i
      call message(0,"  File(s) to be merged",str=fileNames(1))
      do i=2,numberOfFiles
        call message(1,"  ",str=fileNames(i))
      end do
    end block

  end subroutine initialiseAction

  subroutine finaliseAction()
    implicit none
  end subroutine finaliseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Test interface Merging coordinates")
  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: i, iFile, n
    logical :: frameRead
    type(fileTypeDef), allocatable, dimension(:) :: inputFile
    type(frameTypeDef), dimension(:), allocatable :: newCoordinates

    integer :: initialAtoms, currentAtoms
    real(real64) :: hmat(3,3)
    integer :: newNumberOfAtoms, idir
    real(real64), dimension(:,:), allocatable :: newPositions
    character(len=STRLEN), dimension(:), allocatable :: newLabels
    real(real64) :: vshift(3)
    real(real64) :: maxHeight, minHeight, delta

    allocate(inputFile(numberOfFiles))
    allocate(newCoordinates(numberOfFiles))

    newNumberOfAtoms = frame % natoms
    do iFile=1,numberOfFiles
      call initialiseFile(inputFile(iFile), fileNames(iFile))
    
      ! -> get number of atoms
      call getNumberOfAtoms(inputFile(iFile), n, hmat)
      rewind(inputFile(iFile) % funit)
    
      call createSystemArrays(newCoordinates(iFile), n)
      call readCoordinates(frameRead, newCoordinates(iFile), inputFile(iFile))

      if (.not. frameRead) then
        call message(-1,"Error reading file",str=fileNames(iFile))
      end if

      call message(0,"  File",str=fileNames(iFile))
      call message(0,"    Number of atoms",i=n)

      newNumberOfAtoms = newNumberOfAtoms + n
    end do

    initialAtoms = frame % natoms
    currentAtoms = frame % natoms

    ! -> merge coordinates 
    allocate(newPositions(3,newNumberOfAtoms))
    allocate(newLabels(newNumberOfAtoms))

    hmat = 0.d0
    do i=1,3
      hmat(i,i) = frame % hmat(i,i)
    end do
    if (stack /= "NONE") then
      select case (trim(stack))
        case ("x")
          idir = 1
        case ("y")
          idir = 2
        case ("z")
          idir = 3
      end select
      hmat(idir,idir) = hmat(idir,idir) + shift
    end if

    currentAtoms = initialAtoms
    do i=1,frame % natoms
      newPositions(:,i) = frame % pos(:,i) 
      newLabels(i) = frame % lab(i)
    end do

    vshift = [0.d0,0.d0,0.d0]
    do iFile=1,numberOfFiles
      n = newCoordinates(iFile) % natoms
      if (stack /= "NONE") then
        select case (trim(stack))
          case ("x")
            idir = 1
          case ("y")
            idir = 2
          case ("z")
            idir = 3
        end select
        maxHeight = maxval(newPositions(idir,1:currentAtoms))
        minHeight = minval(newCoordinates(iFile) % pos(idir,:))
        delta = shift + maxHeight - minHeight
        vshift(idir) = delta
      end if
      do i=1,n
        newPositions(:,currentAtoms+i) = newCoordinates(iFile) % pos(:,i) + vshift
      end do
      newLabels(currentAtoms+1:currentAtoms+n) = newCoordinates(iFile) % lab
      currentAtoms = currentAtoms + size(newCoordinates(iFile) % pos,2)

      if (stack /= "NONE") then
        select case (trim(stack))
          case ("x")
            hmat(1,1) = hmat(1,1) + newCoordinates(iFile) % hmat(1,1) + shift*2
          case ("y")
            hmat(2,2) = hmat(2,2) + newCoordinates(iFile) % hmat(2,2) + shift*2
          case ("z")
            hmat(3,3) = hmat(3,3) + newCoordinates(iFile) % hmat(3,3) + shift*2
        end select
      end if
  
    end do

    call deleteSystemArrays(frame)
    call createSystemArrays(frame, newNumberOfAtoms)
    frame % natoms = newNumberOfAtoms
    frame % hmat = hmat 
    frame % pos = newPositions
    frame % lab = newLabels

    call hmat2cell (frame % hmat, frame % cell, "DEG")
    call getInverseCellMatrix (frame % hmat, frame % hinv, frame % volume)

    call message(0,"New total number of atoms",i=frame % natoms)
    if (frame % volume > 1.1_real64) then
      call message(0,"New cell")
      call message(0,"...Cell vector A",rv=frame % hmat(1:3,1))
      call message(0,"...Cell vector B",rv=frame % hmat(1:3,2))
      call message(0,"...Cell vector C",rv=frame % hmat(1:3,3))
      call message(0,"...Cell lengths", rv=frame % cell(1:3))
      call message(0,"...Cell angles",  rv=frame % cell(4:6))
      call message(1,"...Volume",       r =frame % volume)
    end if

    call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

    call setUpNeigboursList()
    call initialiseNeighboursList()
    call updateNeighboursList(.true.)

  end subroutine computeAction

end module moduleMergeFiles
