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
module moduleDensityMap3D
  use moduleVariables
  use moduleStrings
  use moduleMessages 
  use moduleSystem 
  use moduleFiles
  use moduleDistances

  public :: computeDensityMap3D, computeDensityMap3DHelp
  private

  character(:), pointer :: actionCommand

  integer, pointer :: tallyExecutions
  type(fileTypeDef), pointer :: outputFile
  
  integer, pointer, dimension(:) :: numberOfBins
  integer, pointer :: numberOfLocalAtoms
  real(8), pointer :: origin(:)
  real(8), pointer, dimension(:) :: densityBox
  
  real(8), pointer, dimension(:,:,:) :: dmap

contains

  subroutine computeDensityMap3DHelp()
    implicit none
    call message(0,"This action computes the 3D density map for the selected atoms in a portion of the system.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --dmap3D +out dmap3D.cube +s O2 +origin 12,34,56 +cell 10")
    call message(0,"  gpta.x --i coord.pdb --dmap3D +out dmap3D.cube +s O2 +origin 12,34,56 +cell 10,10,2")
    call message(0,"  gpta.x --i coord.pdb --dmap3D +out dmap3D.cube +s O2 +cell 10,10,10,90,90,90")
    call message(0,"  gpta.x --i coord.pdb --dmap3D +out dmap3D.cube +s O2 +cell 10,0,0,0,10,0,0,0,10")
  end subroutine computeDensityMap3DHelp

  subroutine computeDensityMap3D(a)

    implicit none
    type(actionTypeDef), target :: a

    integer, dimension(3) :: itmp
    real(8), allocatable, dimension(:,:) :: cartesianCoord
    real(8), allocatable, dimension(:,:) :: fractionalCoord

    integer :: iatm
    real(8), dimension(3) :: p, pinv, dp
    real(8), dimension(3,3) :: htmp
    character(len=STRLEN) :: stringCell

    ! Associate variables
    actionCommand          => a % actionDetails
    tallyExecutions        => a % tallyExecutions

    outputFile             => a % outputFile

    numberOfBins(1:3)      => a % integerVariables(1:3)
    numberOfLocalAtoms     => a % integerVariables(4)

    origin(1:3)            => a % doubleVariables(1:3)
    densityBox(1:9)        => a % doubleVariables(4:12)

    if (a % actionInitialisation ) then
      a % actionInitialisation  = .false.
      a % requiresNeighboursList = .false.

      call assignFlagValue(actionCommand,"+out",outputFile % fname,'dmap3D.cube')
      
      call assignFlagValue(actionCommand,"+nbin",itmp,[50,50,50])
      allocate(a % array3D(itmp(1),itmp(2),itmp(3)) , source=0.d0)
      numberOfBins(1:3) = itmp
      
      call assignFlagValue(actionCommand,"+origin",origin,[0.d0,0.d0,0.d0])
      
      call assignFlagValue(actionCommand,"+cell ", stringCell, "NONE")
      if (stringCell == "NONE") then
        htmp = frame % hmat
      else
        call readCellFreeFormat(stringCell, htmp)
      end if
      densityBox = reshape(htmp,[9])

      tallyExecutions = 0
      return
    end if

    ! Associate distribution array
    dmap => a % array3D

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      ! Select atoms
      if (a % updateAtomsSelection) then
        call dumpSetupInfo()
        call selectAtoms(1,actionCommand,a)
        call createSelectionList(a,1)
        a % updateAtomsSelection=.false.

        numberOfLocalAtoms = count(a % isSelected(:,1))
      end if
      allocate(cartesianCoord(3,numberOfLocalAtoms))
      allocate(fractionalCoord(3,numberOfLocalAtoms))
      
      ! Local array with the coordinates of only the selected atoms
      do iatm=1,numberOfLocalAtoms
        cartesianCoord(1:3,iatm) = frame % pos(1:3,a % idxSelection(iatm,1)) - origin(1:3)
      enddo
      call cartesianToFractional(numberOfLocalAtoms, cartesianCoord, fractionalCoord)

      p(1:3) = densityBox(1:3) + densityBox(4:6) + densityBox(7:9)
      call cartesianToFractional(1,p,pinv)
      where (pinv < 1.d-3) pinv = 1.d0
      dp = pinv / real(numberOfBins,8)

      do iatm=1,numberOfLocalAtoms
        itmp(1:3) = int(fractionalCoord(1:3,iatm) / dp(1:3)) + 1
        if (any(itmp <= 0)) cycle
        if (any(itmp > numberOfBins)) cycle

        dmap(itmp(1),itmp(2),itmp(3)) = dmap(itmp(1),itmp(2),itmp(3)) + 1
      end do

    end if

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) then
      dmap = dmap / tallyExecutions
      call dumpDensityCube()
    end if

  end subroutine ComputeDensityMap3D

  subroutine dumpSetupInfo()
    implicit none
    call message(0,"Compute 3D density map")
    call message(0,"...Output file",str=outputFile % fname)
    call message(0,"...Number of bins",iv=numberOfBins)
    call message(0,"...Origin",rv=origin)
    call message(0,"...Region size A",rv=densityBox(1:3))
    call message(0,"...Region size B",rv=densityBox(4:6))
    call message(0,"...Region size C",rv=densityBox(7:9))
  end subroutine dumpSetupInfo

  subroutine dumpDensityCube()
    implicit none
    real(8), dimension(3,3) :: hmat, hinv, dh
    real(8) :: dvol
    integer :: i, ix, iy, iz
    integer :: funit=123

    call initialiseFile(outputFile, outputFile % fname)
    funit = outputFile % funit

    hmat(1:3,1) = densityBox(1:3) 
    hmat(1:3,2) = densityBox(4:6) 
    hmat(1:3,3) = densityBox(7:9)
    call getInverseCellMatrix(hmat,hinv,dvol)
    dvol = dvol / product(numberOfBins)
    dmap = dmap / dvol

    do i=1,3
      dh(:,i) = hmat(:,i) / numberOfBins(i)
    end do

    write(funit,'("CUBE file generate by GPTA")')
    write(funit,'("Density reported in atoms/angstom^3")')
    write(funit,'(i6,3f13.5)')8, origin / rbohr
    write(funit,'(i6,3f13.5)')numberOfBins(1),dh(:,1) / rbohr
    write(funit,'(i6,3f13.5)')numberOfBins(2),dh(:,2) / rbohr
    write(funit,'(i6,3f13.5)')numberOfBins(3),dh(:,3) / rbohr

    do ix=0,1
      do iy=0,1
        do iz=0,1
          write(funit,'(i6,4f13.5)')1, 1.d0, (origin + ix*hmat(:,1) + iy*hmat(:,2) + iz*hmat(:,3)) / rbohr
        end do
      end do
    end do

    i=0
    do ix=1,numberOfBins(1)
      do iy=1,numberOfBins(2)
        do iz=1,numberOfBins(3)
          i=i+1
          write(funit,'(e13.5)',advance='no')dmap(ix,iy,iz)
          if(mod(i,6)==0)write(funit,*)
        enddo
        if (mod(numberOfBins(3),6)/=0) write(funit,*)
      enddo
    enddo
    call flush(funit)

  end subroutine dumpDensityCube

end module moduleDensityMap3D

