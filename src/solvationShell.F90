!disclaimer
module moduleSolvationShell
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleFiles
  use moduleSuperimposeMolecules
  use moduleDistances
  use moduleMessages
  
  implicit none

  public :: solvationShell, solvationShellHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions
  
  integer, pointer :: moleculeID
  integer, pointer :: smoothingType
  integer, pointer :: numberOfSwaps
  integer, pointer :: numberOfAtomIndices

  integer, pointer, dimension(:) :: numberOfBins
  real(real64), pointer, dimension(:) :: densityBox
  character(len=STRLEN), pointer, dimension(:) :: listOfSwaps
  character(len=STRLEN), pointer :: selectedAtomIndices

  character(len=STRLEN), pointer :: fileName
  character(len=STRLEN) :: moleculeName

  integer, pointer, dimension(:) :: atomIndices
  real(real64), pointer, dimension(:,:) :: referencePositions
  real(real64), pointer, dimension(:,:) :: currentPositions
  real(real64), pointer, dimension(:,:) :: averagePositions

contains

  subroutine solvationShellHelp()
    implicit none
    call message(0,"This action computes the 3D density map for the selected atoms around a solute molecule.")
    call message(0,"The solvent density is computed only inside the largest ellipsoid the fits inside the box around the solvent.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --top --solvation +s O2 +id M2  +cell 10,12,14 +nbin 50,50,50")
    call message(0,"  gpta.x --i coord.pdb --top --solvation +s O2 +id M2  +cell 10,12,14 +nbin 50,50,50 + indices 1,2,5")
    call message(0,"  gpta.x --i coord.pdb --top --solvation +id M2 +s O +cell 9 +nbin 20,20,20 +smooth +ref ace.xyz")
    call message(0,"  gpta.x --i coord.pdb --top --solvation +id M2 +s O +cell 9 +nbin 20,20,20 +ref ace.xyz +swap 2,3,4 6,7")
    ! call message(0,"  gpta.x --i coord.pdb --top --solvation +id M2 +s O +cell 9 +nbin 20,20,20 +ref ace.xyz +equiv 2,3,4 6,7")
  end subroutine solvationShellHelp

  subroutine solvationShell(a)
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
        if (fileName == "NULL") then
          call selectCentralMolecule(a)
        else
          call selectCentralMoleculeFromFile(a,fileName)
        end if

        ! dump info about the action on the screen
        call dumpScreenInfo()

        ! select two groups of atoms
        call selectAtoms(1,actionCommand,a)
        if (count(a % isSelected(:,1)) == 0) call message(-1,"--solvation - no atoms selected in the first group")

        ! create a list of the atoms' indices for each group
        call createSelectionList(a,1)

        ! Throw a warning for unused flags
        call checkUsedFlags(actionCommand)
        firstAction = .false.          
      end if

      call associatePositionsPointers(a)

      call computeAction(a)
    end if

    if (endOfCoordinatesFiles) then
      if (smoothingType == 1) call smoothDistribution3D(numberOfBins(1),numberOfBins(2),numberOfBins(3),a % array3D)
      call finaliseAction(a)
    end if 

  end subroutine solvationShell

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand       => a % actionDetails
    firstAction         => a % firstAction
    tallyExecutions     => a % tallyExecutions
    outputFile          => a % outputFile

    moleculeID          => a % integerVariables(1)
    numberOfBins(1:3)   => a % integerVariables(2:4)    
    smoothingType       => a % integerVariables(5)
    numberOfSwaps       => a % integerVariables(6)
    numberOfAtomIndices => a % integerVariables(7)

    densityBox(1:9)     => a % doubleVariables(1:9)

    fileName            => a % stringVariables(1)
    selectedAtomIndices => a % stringVariables(2)
    listOfSwaps(1:)     => a % stringVariables(3:)

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
    averagePositions   => a % localPositions(1:3,n2+1:n3)

  end subroutine associatePositionsPointers

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    logical :: lflag
    character(STRLEN) :: flagString

    a % actionInitialisation = .false.

    ! get output file name from the command line, if present
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'solvation.cube')
    
    ! get number of bins for the distribution from the command line, if present
    call assignFlagValue(actionCommand,"+nbin",numberOfBins,[50,50,50])
    ! make sure nuberOfBins is an even number
    numberOfBins = 2 * int(numberOfBins / 2)
    allocate(a % array3D(numberOfBins(1),numberOfBins(2),numberOfBins(3)) , source=0.0_real64)

    ! select molecule
    call assignFlagValue(actionCommand,"+id",moleculeName,"NULL")
    if (moleculeName == "NULL") call message(-1,"--solvation - the solute molecule must be chosen with the +id flag")

    call assignFlagValue(actionCommand,"+ref",fileName,"NULL")

    ! set size of the 3D volume
    block
      real(real64), dimension(3,3) :: htmp
      real(real64) :: rtmp
      character(len=STRLEN) :: stringCell
  
      call assignFlagValue(actionCommand,"+cell ", stringCell, "NONE")
      if (stringCell == "NONE") then
        stringCell = "10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0"
      end if
      call readCellFreeFormat(stringCell, htmp)
      densityBox = reshape(htmp,[9])
      rtmp = htmp(1,1) + htmp(2,2) + htmp(3,3)
      if ( sum(abs(densityBox)) - rtmp > 1e-6_real64 ) then
        call message(-1,"--solvation | currently only orthorhombic boxes can be used")
      end if
    end block

    ! smoothing
    call assignFlagValue(actionCommand,"+smooth",lflag,.false.)
    if (lflag) then
      smoothingType = 1
    else
      smoothingType = 0
    end if

    call assignFlagValue(actionCommand,"+swap",lflag,.false.)
    if (lflag) then
      call extractFlag(actionCommand,"+swap",flagString)
      call parse(flagString," ",listOfSwaps,numberOfSwaps)
    else
      numberOfSwaps = 0
    end if

    call assignFlagValue(actionCommand,"+indices",lflag,.false.)
    if (lflag) then
      call extractFlag(actionCommand,"+indices",selectedAtomIndices)
    else
      selectedAtomIndices = "NONE"
    end if

    write(0,*)selectedAtomIndices,"AA"

  tallyExecutions = 0

  end subroutine initialiseAction

  subroutine finaliseAction(a)
    use moduleElements

    implicit none
    type(actionTypeDef), target, intent(inout) :: a
    
    integer :: funit
    real(real64) :: hmat(3,3), hinv(3,3), dvol, dh(3,3), origin(3)
    integer :: i, ix, iy, n

    call initialiseFile(outputFile,outputFile % fname)
  
    ! write cube file
    funit = outputFile % funit

    hmat = reshape(densityBox,[3,3])
    call getInverseCellMatrix(hmat,hinv,dvol)
    dvol = dvol / product(numberOfBins)
    a % array3D = a % array3D / dvol / tallyExecutions

    do i=1,3
      dh(:,i) = hmat(:,i) / numberOfBins(i)
    end do

    origin(1) = -densityBox(1) / 2.0_real64
    origin(2) = -densityBox(5) / 2.0_real64
    origin(3) = -densityBox(9) / 2.0_real64

    n = size(a % localIndices)

    write(funit,'("CUBE file generate by GPTA")')
    write(funit,'("Density reported in atoms/angstom^3")')
    write(funit,'(i6,3f13.5)')n, origin / rbohr

    ! Cell
    write(funit,'(i6,3f13.5)')numberOfBins(1),dh(:,1) / rbohr
    write(funit,'(i6,3f13.5)')numberOfBins(2),dh(:,2) / rbohr
    write(funit,'(i6,3f13.5)')numberOfBins(3),dh(:,3) / rbohr

    ! average molecule
    do i=1,n
      ix = getAtomicNumber(a % localLabels(i))
      write(funit,'(i6,4f13.5)')ix, 0.0_real64, averagePositions(1:3,i) / rbohr / tallyExecutions
    end do

    do ix=1,numberOfBins(1)
      do iy=1,numberOfBins(2)
        write(funit,'(6e13.5)') a % array3D(ix,iy,1:numberOfBins(3))
      enddo
    enddo
    call flush(funit)

    close(outputFile % funit)

  end subroutine finaliseAction

  subroutine dumpScreenInfo()
    implicit none
    integer :: i
    character(len=STRLEN) :: str
    call message(0,"Computing solvation shell")
    call message(0,"...Solute molecule ID",str=moleculeName)
    call message(0,"...Output file",str=outputFile % fname)
    call message(0,"...Number of bins",iv=numberOfBins)
    call message(0,"...Region size A",rv=densityBox(1:3))
    call message(0,"...Region size B",rv=densityBox(4:6))
    call message(0,"...Region size C",rv=densityBox(7:9))
    if (smoothingType == 1) call message(0,"...Density smoothing with tirnagular kernel")
    if (numberOfSwaps > 0) then
      str = listOfSwaps(1)
      do i=2,numberOfSwaps
        str = trim(str) // " - " // trim(listOfSwaps(i))
      end do
      call message(0,"...Checking atoms' permutations for rotation",str=str)
    end if
    if (selectedAtomIndices /= "NONE") then
      call message(0,"...Selecting atoms from molecule",str=selectedAtomIndices)
    end if
  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: i, idx, iatm, nRef, nSolvent

    real(real64), dimension(3) :: shift, dij
    real(real64), allocatable, dimension(:,:) :: tmpArray
    real(real64), allocatable, dimension(:,:) :: solventPositions
    
    real(real64), dimension(3,3) :: rotationalMatrix
    
    real(real64) :: pinv(3,3), dp(3), boundaries(3)
    integer :: ipos(3)

    nRef = size(atomIndices)
    nSolvent = count(a % isSelected(:,1))

    ! Reference molecule's com is the origin
    ! do i=1,nRef
    !   write(0,*)i,referencePositions(1:3,i)
    ! end do
  
    ! Extract coordinated of current molecule ...
    do i=1,nRef
      iatm = atomIndices(i)
      currentPositions(1:3,i) =  frame % pos(1:3,iatm)
    end do
    call CenterCoords(shift, nRef, currentPositions)
    ! Reference molecule's com is the origin
    ! do i=1,nRef
    !   write(0,*)i,currentPositions(1:3,i)
    ! end do

    ! ... and solvent atoms
    allocate(solventPositions(3,nSolvent))
    do i=1,nSolvent
      iatm = a % idxSelection(i,1)
      solventPositions(1:3,i) = frame % pos(1:3,iatm) - shift(1:3)
    enddo
    allocate(tmpArray(3,nSolvent))
    call cartesianToFractionalNINT(nSolvent, solventPositions, tmpArray)
    call fractionalToCartesian(nSolvent, tmpArray, solventPositions)
    deallocate(tmpArray)
    
    ! Compute number of permutations to enhance statistics <-- hard task, need to preseve chirality -->

    ! Align current molcule to reference
    ! -- currentPositions are updated inside
    ! *=- check for rotations (+swap) is done inside -=*
    call findBestOverlap(nRef, referencePositions, currentPositions, rotationalMatrix)
    ! write(0,*)rotationalMatrix
    ! Average positions
    do i=1,nRef
      dij = matmul(rotationalMatrix, currentPositions(1:3,i))
      averagePositions(1:3,i) = averagePositions(1:3,i) + dij
    end do

    ! Rotate all the solvent atoms
    do i=1,nSolvent
      dij(1:3) = solventPositions(1:3,i)
      solventPositions(1:3,i) = matmul(rotationalMatrix , dij(1:3))
    end do

    allocate(tmpArray(3,nSolvent))
    call cartesianToFractionalNINT(nSolvent, solventPositions, tmpArray)
    call move_alloc(tmpArray, solventPositions)

    ! Tally the density map
    call cartesianToFractionalNINT(3,densityBox/2.0_real64,pinv)
    do idx=1,3
      boundaries(idx) = abs(pinv(idx,idx))
      dp(idx) = 2.0_real64 * pinv(idx,idx) / numberOfBins(idx)
    end do

    do idx=1,nSolvent
      ! consider only the points inside an ellypsoid that fits into the box
      if ( sum( (solventPositions(:,idx) / boundaries(:))**2 ) > 1 ) cycle
      ipos(1:3) = int( (solventPositions(1:3,idx) + boundaries(1:3)) / dp(1:3)) + 1
      a % array3D(ipos(1),ipos(2),ipos(3)) = a % array3D(ipos(1),ipos(2),ipos(3))  + 1.0_real64
    end do    

    deallocate(solventPositions)

  end subroutine computeAction

  ! subroutine rotateIndices(nRot, localIndex)
  !   integer, intent(out) :: nRot
  !   integer, dimension(:,:), intent(inout) :: localIndex

  !   integer :: i, j, nat, ntot

  !   integer, allocatable, dimension(:) :: numberOfIndices
  !   integer, allocatable, dimension(:) :: listOfIndices
  !   character(len=200), dimension(100) :: tokens

  !   nat = size(atomIndices)

  !   if (numberOfAtomIndices == 0) then
  !     nRot = 1
  !     localIndex(:,nRot) = atomIndices
  !   else

  !     allocate(numberOfIndices(numberOfAtomIndices))
  !     do i=1,numberOfAtomIndices
  !       call parse(selectedAtomIndices(i),",",tokens,numberOfIndices(i))
  !     end do
  !     ntot = sum(numberOfIndices)
  !     allocate(listOfIndices(ntot))
  !     ntot = 0
  !     do i=1,numberOfAtomIndices
  !       call parse(selectedAtomIndices(i),",",tokens,nat)
  !       do j=1,nat
  !         read(tokens(j),*) listOfIndices(ntot+j)
  !       end do
  !       ntot = ntot + nat
        
  !       ! do j=1,nat
  !       !   write(0,*)trim(tokens(j))
  !       ! end do
        
  !       ! ... do something smart here ...

  !     end do

  !     nRot = 9
  !     localIndex(:,1) = [1,2,3,4,5]
  !     localIndex(:,2) = [1,2,4,5,3]
  !     localIndex(:,3) = [1,2,5,3,4]
  !     localIndex(:,4) = [1,4,3,5,2]
  !     localIndex(:,5) = [1,5,3,2,4]
  !     localIndex(:,6) = [1,3,5,4,2]
  !     localIndex(:,7) = [1,5,2,4,3]
  !     localIndex(:,8) = [1,3,4,2,5]
  !     localIndex(:,9) = [1,4,2,3,5]

  !   end if

  ! end subroutine rotateIndices

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

    integer :: imol, iatm, idx, jdx
    integer, allocatable, dimension(:) :: itmp
    character(len=200), dimension(100) :: tokens

    real(real64), dimension(3) :: shift
    
    if (moleculeName == "NULL") then
      call message(-1,"Solvation - +id flag is required")
    else
      do imol=1,numberOfMolecules
        if (listOfMolecules(imol) % resname == moleculeName) then
          moleculeID = imol
          exit
        end if
      end do
    end if

    if (selectedAtomIndices == "NONE") then
      numberOfAtomIndices = listOfMolecules(imol) % numberOfAtoms
      allocate(itmp(numberOfAtomIndices))
      do idx=1,numberOfAtomIndices
        itmp(idx) = idx
      end do
    else
      call parse(selectedAtomIndices,",",tokens,numberOfAtomIndices)
      allocate(itmp(numberOfAtomIndices))
      do idx=1,numberOfAtomIndices
        read(tokens(idx),*) itmp(idx)
      end do
    end if

    allocate(a % localLabels(numberOfAtomIndices))
    allocate(a % localIndices(numberOfAtomIndices))
    allocate(a % localPositions(3,numberOfAtomIndices*3), source=0.0_real64)    

    ! this can be replaced if an external frame is required
    do jdx=1,numberOfAtomIndices
      idx = itmp(jdx)
      iatm = listOfMolecules(imol) % listOfAtoms(idx)
      a % localIndices(jdx) = iatm
      a % localLabels(jdx) = frame % lab(iatm)
      a % localPositions(1:3,idx) = frame % pos(1:3,iatm)
    end do

    ! write(0,*)numberOfAtomIndices
    ! write(0,*)a % localIndices
    ! write(0,*)a % localLabels
    ! do idx=1,numberOfAtomIndices
    !   write(0,*)idx, a % localPositions(1:3,idx)
    ! end do

    call CenterCoords(shift, numberOfAtomIndices, a % localPositions) 

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

end module moduleSolvationShell

subroutine smoothDistribution3D(nx,ny,nz,array3D)
  use moduleVariables, only: real64
  implicit none
  integer, intent(in) :: nx, ny, nz
  real(real64), dimension(nx,ny,nz) :: array3D
  integer :: ix, iy, iz
  real(real64) :: val
  real(real64), allocatable, dimension(:,:,:) :: smooth
  
  integer :: i, n
  real(real64), allocatable, dimension(:) :: w
  
  allocate(smooth(nx,ny,nz), source=0.0_real64)
  n = 1
  allocate(w(-n:n))

  do i=-n,n
    w(i) = abs(real(i,8))/real(n+1,8)
  end do
  
  do ix=2,nx-1
    do iy=2,ny-1
      do iz=2,nz-1
        val = array3D(ix,iy,iz)
        smooth(ix,iy,iz) = val
        do i=-n,n
          smooth(ix+i,iy,iz) = smooth(ix+i,iy,iz) + val * w(i)
        end do
        do i=-n,n
          smooth(ix,iy+i,iz) = smooth(ix,iy+i,iz) + val * w(i)
        end do
        do i=-n,n
          smooth(ix,iy,iz+i) = smooth(ix,iy,iz+i) + val * w(i)
        end do
      enddo
    enddo
  enddo

  array3D = smooth / sum(w)

end subroutine smoothDistribution3D
