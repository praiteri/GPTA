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
module moduleSolvationShell
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleFiles
  use moduleAlignMolecules
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
  
  integer, pointer, dimension(:) :: numberOfBins
  real(8), pointer, dimension(:) :: densityBox

  character(len=STRLEN), pointer :: fileName
  character(len=STRLEN) :: moleculeName

contains

  subroutine solvationShellHelp()
    implicit none
    call message(0,"This action computes the 3D density map for the selected atoms around a solute molecule.")
    call message(0,"The solvent density is computed only inside the largest ellipsoid the fits inside the box around the solvent.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --top --solvation +s O2 +id M2  +cell 10,12,14  +nbin 50,50,50")
    call message(0,"  gpta.x --i coord.pdb --top --solvation +id M2 +s O +cell 9 +nbin 20,20,20 +smooth +ref ace.xyz")
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

        block
          integer :: n
          n = size(a % localIndices)
          call CenterCoords(n, a % localPositions) 
        end block

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

    densityBox(1:9)      => a % doubleVariables(1:9)

    fileName             => a % stringVariables(1)

  end subroutine associatePointers

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    logical :: lflag

    a % actionInitialisation = .false.

    ! get output file name from the command line, if present
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'solvation.cube')
    
    ! get number of bins for the distribution from the command line, if present
    call assignFlagValue(actionCommand,"+nbin",numberOfBins,[50,50,50])
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
        htmp = frame % hmat
      else
        call readCellFreeFormat(stringCell, htmp)
      end if
      densityBox = reshape(htmp,[9])
      rtmp = htmp(1,1) + htmp(2,2) + htmp(3,3)
      if ( sum(abs(densityBox)) - rtmp > 1e-6 ) then
        call message(-1,"--solvation | only orthorhombic boxes can be used")
      end if
    end block

    ! smoothing
    call assignFlagValue(actionCommand,"+smooth",lflag,.false.)
    if (lflag) then
      smoothingType = 1
    else
      smoothingType = 0
    end if
    tallyExecutions = 0

  end subroutine initialiseAction

  subroutine finaliseAction(a)
    use moduleElements

    implicit none
    type(actionTypeDef), target, intent(inout) :: a
    
    integer :: funit
    real(8) :: hmat(3,3), hinv(3,3), dvol, dh(3,3), origin(3)
    integer :: i, ix, iy, iz, n

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
    do i=n+1,2*n
      ix = getAtomicNumber(a % localLabels(i-n))
      write(funit,'(i6,4f13.5)')ix, 0.d0, a % localPositions(1:3,i) / rbohr / tallyExecutions
    end do

    i=0
    do ix=1,numberOfBins(1)
      do iy=1,numberOfBins(2)
        do iz=1,numberOfBins(3)
          i=i+1
          write(funit,'(e13.5)',advance='no')a % array3D(ix,iy,iz)
          if(mod(i,6)==0)write(funit,*)
        enddo
        if (mod(numberOfBins(3),6)/=0) write(funit,*)
      enddo
    enddo
    call flush(funit)

    close(outputFile % funit)

  end subroutine finaliseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Computing solvation shell")
    call message(0,"...Output file",str=outputFile % fname)
    call message(0,"...Number of bins",iv=numberOfBins)
    call message(0,"...Region size A",rv=densityBox(1:3))
    call message(0,"...Region size B",rv=densityBox(4:6))
    call message(0,"...Region size C",rv=densityBox(7:9))

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: nsel, idx, iatm, n
    real(8), allocatable, dimension(:,:) :: currentPositions
    real(8), allocatable, dimension(:,:) :: solventPositionsCart
    real(8), allocatable, dimension(:,:) :: solventPositionsFrac
    real(8), dimension(3) :: xcom
    real(8), dimension(3) :: boundaries, dp
    integer, dimension(3) :: nb
    integer :: ipos(3)
    real(8) :: rmsd, rotmat(3,3), dij(3), rvec(9), rtmp

    n = size(a % localIndices)

    ! centre of mass of the solute
    xcom = 0.d0
    allocate(currentPositions(3,n))
    do idx=1,n
      iatm = a % localIndices(idx)
      currentPositions(1:3,idx) =  frame % pos(1:3,iatm) 
    end do
    xcom = sum(currentPositions,dim=2) / n
    do idx=1,n
      currentPositions(1:3,idx) = currentPositions(1:3,idx) - xcom(1:3)
    end do

    ! current positions are centered to the origin
    call Superimpose(n, a % localPositions, currentPositions, rmsd, rvec)
    rotmat = transpose(reshape(rvec,[3,3]))

    ! Average positions
    do idx=1,n
      dij = matmul(rotmat, currentPositions(1:3,idx))
      a % localPositions(1:3,idx+n) = a % localPositions(1:3,idx+n) + dij(1:3)
    end do      

    nsel = count(a % isSelected(:,1))
    allocate(solventPositionsCart(3,nsel))
    allocate(solventPositionsFrac(3,nsel))

    ! loop over over the first group of atoms
    do idx=1,nsel
      iatm = a % idxSelection(idx,1)
      dij(1:3) = frame % pos(1:3,iatm) - xcom(1:3)
      solventPositionsCart(1:3,idx) = matmul(rotmat , dij(1:3))
    end do

    call cartesianToFractionalNINT(nsel, solventPositionsCart, solventPositionsFrac)

    ! requires an even number of bins
    nb = (int((numberOfBins-1) / 2.d0) + 1) *2
    block
      real(8), dimension(3,3) :: pinv
      call cartesianToFractionalNINT(3,densityBox/2.d0,pinv)
      do idx=1,3
        boundaries(idx) = abs(pinv(idx,idx))
        dp(idx) = 2.d0 * pinv(idx,idx) / nb(idx)
      end do
    end block

    do idx=1,nsel
      ! consider only the points inside an ellypsoid that fits into the box
      rtmp = sum( (solventPositionsFrac(:,idx) / boundaries(:))**2 )
      if (rtmp > 1) cycle
      ipos(1:3) = int( (solventPositionsFrac(1:3,idx) + boundaries(1:3)) / dp(1:3)) + 1
      a % array3D(ipos(1),ipos(2),ipos(3)) = a % array3D(ipos(1),ipos(2),ipos(3))  + 1
    end do 

  end subroutine computeAction

  subroutine selectCentralMolecule(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: imol, n, iatm, idx
    
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
    allocate(a % localPositions(3,n*2))    
    allocate(a % localIndices(n))

    a % localLabels(:) = listOfMolecules(imol) % listOfLabels(:)
    a % localIndices(:) = listOfMolecules(imol) % listOfAtoms(:)
    
    ! this can be replaced if an external frame is required
    do idx=1,n
      iatm = listOfMolecules(imol) % listOfAtoms(idx)
      a % localPositions(1:3,idx) = frame % pos(1:3,iatm)
    end do

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

    array3D = smooth

  end subroutine smoothDistribution3D
