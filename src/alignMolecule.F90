!disclaimer
module moduleAlignMolecule
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleFiles
  use moduleSuperimposeMolecules
  use moduleDistances
  use moduleMessages
  
  implicit none

  public :: alignMolecule, alignMoleculeHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions
  
  integer, pointer :: moleculeID
  integer, pointer :: numberOfSwaps

  character(len=STRLEN), pointer, dimension(:) :: listOfSwaps
  character(len=STRLEN), pointer, dimension(:) :: listOfEquiv

  character(len=STRLEN), pointer :: fileName
  character(len=STRLEN) :: moleculeName

  integer, pointer, dimension(:) :: atomIndices
  real(real64), pointer, dimension(:,:) :: referencePositions
  real(real64), pointer, dimension(:,:) :: currentPositions

contains

  subroutine alignMoleculeHelp()
    implicit none
    call message(0,"This action rotates the atomic coordinates to align the selected molecule to their initial position.")
    call message(0,"It works best if before molecule is initially at the centre of the box.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --top --align +id M2")
  end subroutine alignMoleculeHelp

  subroutine alignMolecule(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: nsel
    character(len=STRLEN) :: string

    call associatePointers(a)

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      ! Atoms' selection
      if (firstAction) then
        if (fileName == "NULL") then
          call selectCentralMolecule(a)
        else
          call selectCentralMoleculeFromFile(a,fileName)
        end if

        ! dump info about the action on the screen
        call dumpScreenInfo()

        ! select two groups of atoms
        call getNumberOfAtomSelections(actionCommand,nsel)
        if (nsel == 0) then
          string = trim(actionCommand)//" +s all"
          call selectAtoms(1,string,a)
        else
          call selectAtoms(1,actionCommand,a)
        end if

        ! create a list of the atoms' indices for each group
        call createSelectionList(a,1)

        ! Throw a warning for unused flags
        call checkUsedFlags(actionCommand)
        firstAction = .false.          
      end if

      call associatePositionsPointers(a)

      call computeAction(a)
    end if

  end subroutine alignMolecule

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile

    moleculeID           => a % integerVariables(1)
    numberOfSwaps        => a % integerVariables(6)

    fileName             => a % stringVariables(1)
    listOfSwaps(1:)      => a % stringVariables(2:)
    listOfEquiv(1:)      => a % stringVariables(2:)

  end subroutine associatePointers

  subroutine associatePositionsPointers(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: n1, n2, n3
    n1 = size(a % localIndices)
    n2 = n1 * 2
    n3 = n1 * 3

    atomIndices        => a % localIndices(1:n1)

    referencePositions => a % localPositions(1:3,   1:n1)
    currentPositions   => a % localPositions(1:3,n1+1:n2)

  end subroutine associatePositionsPointers

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    logical :: lflag
    character(STRLEN) :: flagString

    a % actionInitialisation = .false.

    ! get output file name from the command line, if present
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'align.pdb')
    call initialiseFile(outputFile,outputFile % fname)

    ! select molecule
    call assignFlagValue(actionCommand,"+id",moleculeName,"NULL")
    if (moleculeName == "NULL") call message(-1,"--align - the solute molecule must be chosen with the +id flag")

    call assignFlagValue(actionCommand,"+ref",fileName,"NULL")

    call assignFlagValue(actionCommand,"+swap",lflag,.false.)
    if (lflag) then
      call extractFlag(actionCommand,"+swap",flagString)
      call parse(flagString," ",listOfSwaps,numberOfSwaps)
    else
      numberOfSwaps = 0
    end if

    tallyExecutions = 0

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    implicit none
    integer :: i
    character(len=STRLEN) :: str
    call message(0,"Aligning molecules")
    call message(0,"...Solute molecule ID",str=moleculeName)
    call message(0,"...Output file",str=outputFile % fname)
    if (numberOfSwaps > 0) then
      str = listOfSwaps(1)
      do i=2,numberOfSwaps
        str = trim(str) // " - " // trim(listOfSwaps(i))
      end do
      call message(0,"...Checking atoms' permutations for rotation",str=str)
    end if
  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: i, iatm, nRef, nSolvent

    real(real64), save :: shift0(3)
    real(real64), dimension(3) :: shift, dij
    real(real64), allocatable, dimension(:,:) :: tmpArray
    real(real64), allocatable, dimension(:,:) :: systemPositions
    
    real(real64), dimension(3,3) :: rotationalMatrix
    
    integer :: irot, nRotations
    integer, allocatable, dimension(:,:) :: localIndices

    nRef = size(atomIndices)
    nSolvent = count(a % isSelected(:,1))

    allocate(localIndices(nRef,9))
    call rotateIndices(nRotations,localIndices)

    do irot=1,nRotations
    
      ! Extract coordinated of current molecule ...
      do i=1,nRef
        iatm = localIndices(i,irot)
        currentPositions(1:3,i) =  frame % pos(1:3,iatm)
      end do
      call CenterCoords(shift, nRef, currentPositions)
      if (tallyExecutions == 1) shift0 = shift

      ! ... and solvent atoms
      allocate(systemPositions(3,frame % natoms))
      do i=1,frame % natoms
        systemPositions(1:3,i) = frame % pos(1:3,i) - shift(1:3)
      enddo
      
      allocate(tmpArray(3,frame % natoms))
      call cartesianToFractionalNINT(frame % natoms, systemPositions, tmpArray)
      call fractionalToCartesian(frame % natoms, tmpArray, systemPositions)
      if (numberOfMolecules > 0) call reassembleAllMolecules2(frame % natoms , systemPositions)
      
      ! Align current molcule to reference
      call findBestOverlap(nRef, referencePositions, currentPositions, rotationalMatrix)
  
      ! Rotate all the solvent atoms
      do i=1,nSolvent
        iatm = a % idxSelection(i,1)
        dij(1:3) = systemPositions(1:3,iatm)
        tmpArray(1:3,i) = matmul(rotationalMatrix , dij(1:3)) + shift0(1:3)
      end do

      block
        character(len=cp), allocatable, dimension(:) :: tmpLab
        integer :: xx
        
        if ( outputFile % ftype == "dcd" ) then
          if (tallyExecutions == 1) then
            call dumpDCD(outputFile % funit, nSolvent, tmpArray(1:3,1:nSolvent), frame%hmat, .true.)
          else 
            call dumpDCD(outputFile % funit, nSolvent, tmpArray(1:3,1:nSolvent), frame%hmat, .false.)
          end if
        else
          allocate(tmpLab(nSolvent))
          do xx=1,nSolvent
            iatm = a % idxSelection(xx,1)
            tmpLab(xx) = frame % lab(iatm)
          enddo
          call dumpPDB(outputFile % funit, nSolvent, tmpArray(1:3,1:nSolvent), tmpLab(1:nSolvent), frame%hmat)
        end if
      end block

    end do

  end subroutine computeAction

  subroutine rotateIndices(nRot, localIndex)
    integer, intent(out) :: nRot
    integer, dimension(:,:), intent(inout) :: localIndex

      nRot = 1
      localIndex(:,nRot) = atomIndices

  end subroutine rotateIndices

  subroutine findBestOverlap(nat, referencePositions, currentPositions, rotmat)
    integer, intent(in) :: nat
    real(real64), intent(inout), dimension(3,nat) :: referencePositions, currentPositions
    real(real64), intent(out), dimension(3,3) :: rotmat
    
    real(real64) :: rmsd, rvec(9)
    integer :: i, j, k, ntot, n1, n2, n, ntmp
    integer, allocatable, dimension(:) :: order
    logical, external :: nextp
    integer, allocatable, dimension(:) :: numberOfIndices
    integer, allocatable, dimension(:) :: listOfIndices
    
    real(real64), allocatable, dimension(:,:) :: currentPositions1
    real(real64), allocatable, dimension(:,:) :: currentPositions2
    real(real64) :: rmsd2, rvec2(9)

    character(len=200), dimension(100) :: tokens

    allocate(currentPositions1(3,nat))
    allocate(currentPositions2(3,nat))
    currentPositions1 = currentPositions
    currentPositions2 = currentPositions

    ! current positions are centered to the origin
    call Superimpose(nat, referencePositions, currentPositions, rmsd, rvec)

    ! check if swapped atoms give a lower rmsd
    if ( numberOfSwaps > 0) then
     allocate(numberOfIndices(numberOfSwaps))
      
      do i=1,numberOfSwaps
        call parse(listOfSwaps(i),",",tokens,numberOfIndices(i))
      end do
      ntot = sum(numberOfIndices)
      allocate(listOfIndices(ntot))
      ntot = 0
      do i=1,numberOfSwaps
        call parse(listOfSwaps(i),",",tokens,ntmp)
        do j=1,ntmp
          read(tokens(j),*) listOfIndices(ntot+j) 
        end do
        ntot = ntot + ntmp
      end do
      allocate(order(ntot))
      do i=1,ntot
        order(i) = i
      enddo
    
      main : do while (nextp(ntot,order)) 
        n=0
        do k=1,numberOfSwaps
          n1 = n+1
          n2 = n+numberOfIndices(k)
          if (any(order(n1:n2) < n1) .or. any(order(n1:n2) > n2)) cycle main
          n = n + numberOfIndices(k)
        end do
      
        do i=1,ntot
          j = listOfIndices(i)
          k = listOfIndices(order(i))
          currentPositions2(:,j) = currentPositions1(:,k)
        end do

        call Superimpose(nat, referencePositions, currentPositions2, rmsd2, rvec2)
        if (rmsd2 < rmsd) then
          rmsd = rmsd2
          rvec = rvec2
          currentPositions = currentPositions2
        end if
        
      end do main
    end if

    rotmat = transpose(reshape(rvec,[3,3]))

  end subroutine

  subroutine selectCentralMolecule(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: imol, n, iatm, idx
    real(real64), dimension(3) :: shift
    
    if (moleculeName == "NULL") then
      call message(-1,"Align - +id flag is required")
    else
      do imol=1,numberOfMolecules
        if (listOfMolecules(imol) % resname == moleculeName) then
          moleculeID = imol
          exit
        end if
      end do
    end if

    n = listOfMolecules(imol) % numberOfAtoms
    allocate(a % localLabels(n))
    allocate(a % localIndices(n))
    allocate(a % localPositions(3,n*3), source=0.0_real64)    

    a % localLabels(:) = listOfMolecules(imol) % listOfLabels(:)
    a % localIndices(:) = listOfMolecules(imol) % listOfAtoms(:)
    
    ! this can be replaced if an external frame is required
    do idx=1,n
      iatm = listOfMolecules(imol) % listOfAtoms(idx)
      a % localPositions(1:3,idx) = frame % pos(1:3,iatm)
    end do

    call CenterCoords(shift, n, a % localPositions) 

  end subroutine selectCentralMolecule

  subroutine selectCentralMoleculeFromFile(a,f)
    use moduleFiles
    use moduleRead

    implicit none
    type(actionTypeDef), target :: a
    character(len=STRLEN) :: f

    type(fileTypeDef) :: inputFile

    integer :: imol, n
    logical :: frameRead
    real(real64), dimension(3,3) :: hmat
    
    ! Initialise coordinates file
    call initialiseFile(inputFile, f)
    
    ! -> get number of atoms
    call getNumberOfAtoms(inputFile, n, hmat)
    rewind(inputFile % funit)
  
    call createSystemArrays(a % localFrame, n)
    call readCoordinates(frameRead, a % localFrame, inputFile)

    allocate(a % localLabels(n))
    a % localLabels = a % localFrame % lab

    allocate(a % localPositions(3,n*2))    
    a % localPositions(1:3,1:n) = a % localFrame % pos(1:3,1:n)
    
    do imol=1,numberOfMolecules
      if (listOfMolecules(imol) % resname == moleculeName) then
        moleculeID = imol
        exit
      end if
    end do

    a % localIndices = listOfMolecules(moleculeID) % listOfAtoms
    
    call deleteSystemArrays(a % localFrame)

  end subroutine selectCentralMoleculeFromFile

end module moduleAlignMolecule