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
  integer, pointer :: numberOfEquiv

  integer, pointer, dimension(:) :: numberOfBins
  real(8), pointer, dimension(:) :: densityBox
  character(len=STRLEN), pointer, dimension(:) :: listOfSwaps
  character(len=STRLEN), pointer, dimension(:) :: listOfEquiv

  character(len=STRLEN), pointer :: fileName
  character(len=STRLEN) :: moleculeName

  integer, pointer, dimension(:) :: atomIndices
  real(8), pointer, dimension(:,:) :: referencePositions
  real(8), pointer, dimension(:,:) :: currentPositions
  real(8), pointer, dimension(:,:) :: averagePositions

contains

  subroutine solvationShellHelp()
    implicit none
    call message(0,"This action computes the 3D density map for the selected atoms around a solute molecule.")
    call message(0,"The solvent density is computed only inside the largest ellipsoid the fits inside the box around the solvent.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --top --solvation +s O2 +id M2  +cell 10,12,14  +nbin 50,50,50")
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
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile

    moleculeID           => a % integerVariables(1)
    numberOfBins(1:3)    => a % integerVariables(2:4)    
    smoothingType        => a % integerVariables(5)
    numberOfSwaps        => a % integerVariables(6)
    numberOfEquiv        => a % integerVariables(7)

    densityBox(1:9)      => a % doubleVariables(1:9)

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
    allocate(a % array3D(numberOfBins(1),numberOfBins(2),numberOfBins(3)) , source=0.d0)

    ! select molecule
    call assignFlagValue(actionCommand,"+id",moleculeName,"NULL")
    if (moleculeName == "NULL") call message(-1,"--solvation - the solute molecule must be chosen with the +id flag")

    call assignFlagValue(actionCommand,"+ref",fileName,"NULL")

    ! set size of the 3D volume
    block
      real(8), dimension(3,3) :: htmp
      real(8) :: rtmp
      character(len=STRLEN) :: stringCell
  
      call assignFlagValue(actionCommand,"+cell ", stringCell, "NONE")
      if (stringCell == "NONE") then
        stringCell = "10.d0, 0.d0, 0.d0, 0.d0, 10.d0, 0.d0, 0.d0, 0.d0, 10.d0"
      end if
      call readCellFreeFormat(stringCell, htmp)
      densityBox = reshape(htmp,[9])
      rtmp = htmp(1,1) + htmp(2,2) + htmp(3,3)
      if ( sum(abs(densityBox)) - rtmp > 1e-6 ) then
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

    call assignFlagValue(actionCommand,"+equi",lflag,.false.)
    if (lflag) then
      call extractFlag(actionCommand,"+equi",flagString)
      call parse(flagString," ",listOfEquiv,numberOfEquiv)
    else
      numberOfEquiv = 0
    end if

    tallyExecutions = 0

  end subroutine initialiseAction

  subroutine finaliseAction(a)
    use moduleElements

    implicit none
    type(actionTypeDef), target, intent(inout) :: a
    
    integer :: funit
    real(8) :: hmat(3,3), hinv(3,3), dvol, dh(3,3), origin(3)
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

    origin(1) = -densityBox(1) / 2.d0
    origin(2) = -densityBox(5) / 2.d0
    origin(3) = -densityBox(9) / 2.d0

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
      write(funit,'(i6,4f13.5)')ix, 0.d0, averagePositions(1:3,i) / rbohr / tallyExecutions
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
    if (numberOfEquiv > 0) then
      str = listOfEquiv(1)
      do i=2,numberOfEquiv
        str = trim(str) // " - " // trim(listOfEquiv(i))
      end do
      call message(0,"...Checking atoms' permutations for equivalent atoms",str=str)
    end if
  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: i, idx, iatm, nRef, nSolvent

    real(8), dimension(3) :: shift, dij
    real(8), allocatable, dimension(:,:) :: tmpArray
    real(8), allocatable, dimension(:,:) :: solventPositions
    
    real(8), dimension(3,3) :: rotationalMatrix
    
    real(8) :: pinv(3,3), dp(3), boundaries(3)
    integer :: ipos(3)

    integer :: irot, nRotations
    integer, allocatable, dimension(:,:) :: localIndices

    nRef = size(atomIndices)
    nSolvent = count(a % isSelected(:,1))

    ! Reference molecule's com is the origin
    ! do i=1,nRef
    !   write(0,*)i,referencePositions(1:3,i)
    ! end do
  
    allocate(localIndices(nRef,9))
    call rotateIndices(nRotations,localIndices)

    do irot=1,nRotations
    
      ! Extract coordinated of current molecule ...
      do i=1,nRef
        iatm = localIndices(i,irot)
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
        averagePositions(1:3,i) = averagePositions(1:3,i) + dij / nRotations
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
      call cartesianToFractionalNINT(3,densityBox/2.d0,pinv)
      do idx=1,3
        boundaries(idx) = abs(pinv(idx,idx))
        dp(idx) = 2.d0 * pinv(idx,idx) / numberOfBins(idx)
      end do
  
      do idx=1,nSolvent
        ! consider only the points inside an ellypsoid that fits into the box
        if ( sum( (solventPositions(:,idx) / boundaries(:))**2 ) > 1 ) cycle
        ipos(1:3) = int( (solventPositions(1:3,idx) + boundaries(1:3)) / dp(1:3)) + 1
        a % array3D(ipos(1),ipos(2),ipos(3)) = a % array3D(ipos(1),ipos(2),ipos(3))  + 1.d0/nRotations
      end do    

      deallocate(solventPositions)

    end do
    ! stop
  end subroutine computeAction

  subroutine rotateIndices(nRot, localIndex)
    integer, intent(out) :: nRot
    integer, dimension(:,:), intent(inout) :: localIndex

    integer :: i, j, nat, ntot

    integer, allocatable, dimension(:) :: numberOfIndices
    integer, allocatable, dimension(:) :: listOfIndices
    character(len=200), dimension(100) :: tokens

    nat = size(atomIndices)

    if (numberOfEquiv == 0) then
      nRot = 1
      localIndex(:,nRot) = atomIndices
    else

      allocate(numberOfIndices(numberOfEquiv))
      do i=1,numberOfEquiv
        call parse(listOfEquiv(i),",",tokens,numberOfIndices(i))
      end do
      ntot = sum(numberOfIndices)
      allocate(listOfIndices(ntot))
      ntot = 0
      do i=1,numberOfEquiv
        call parse(listOfEquiv(i),",",tokens,nat)
        do j=1,nat
          read(tokens(j),*) listOfIndices(ntot+j)
        end do
        ntot = ntot + nat
        
        ! do j=1,nat
        !   write(0,*)trim(tokens(j))
        ! end do
        
        ! ... do something smart here ...

      end do

      nRot = 9
      localIndex(:,1) = [1,2,3,4,5]
      localIndex(:,2) = [1,2,4,5,3]
      localIndex(:,3) = [1,2,5,3,4]
      localIndex(:,4) = [1,4,3,5,2]
      localIndex(:,5) = [1,5,3,2,4]
      localIndex(:,6) = [1,3,5,4,2]
      localIndex(:,7) = [1,5,2,4,3]
      localIndex(:,8) = [1,3,4,2,5]
      localIndex(:,9) = [1,4,2,3,5]

    end if

  end subroutine rotateIndices

  subroutine findBestOverlap(nat, referencePositions, currentPositions, rotmat)
    integer, intent(in) :: nat
    real(8), intent(inout), dimension(3,nat) :: referencePositions, currentPositions
    real(8), intent(out), dimension(3,3) :: rotmat
    
    real(8) :: rmsd, rvec(9)
    integer :: i, j, k, ntot, n1, n2, n, ntmp
    integer, allocatable, dimension(:) :: order
    logical, external :: nextp
    integer, allocatable, dimension(:) :: numberOfIndices
    integer, allocatable, dimension(:) :: listOfIndices
    
    real(8), allocatable, dimension(:,:) :: currentPositions1
    real(8), allocatable, dimension(:,:) :: currentPositions2
    real(8) :: rmsd2, rvec2(9)

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
      ! write(0,*)listOfIndices
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
    real(8), dimension(3) :: shift
    
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

    n = listOfMolecules(imol) % numberOfAtoms
    allocate(a % localLabels(n))
    allocate(a % localIndices(n))
    allocate(a % localPositions(3,n*3), source=0.d0)    

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
    real(8), dimension(3,3) :: hmat
    
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
  implicit none
  integer, intent(in) :: nx, ny, nz
  real(8), dimension(nx,ny,nz) :: array3D
  integer :: ix, iy, iz
  real(8) :: val
  real(8), allocatable, dimension(:,:,:) :: smooth
  
  integer :: i, n
  real(8), allocatable, dimension(:) :: w
  
  allocate(smooth(nx,ny,nz), source=0.d0)
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
