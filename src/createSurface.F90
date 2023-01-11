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
module moduleCreateSurface
  use moduleVariables
  use moduleFiles
  use moduleMessages
  use moduleSystem
  use moduleStrings
  use m_mrgrnk
  
  implicit none

  public :: createSurface, createSurfaceHelp
  private

  character(:), pointer :: actionCommand
  type(fileTypeDef), pointer :: outputFile
  logical, pointer :: firstAction

  integer, pointer, dimension(:) :: imiller

contains

  subroutine createSurfaceHelp()
    implicit none
    call message(0,"This action creates a new system unit cell oriented with the chosen Miller index parallel to the z axis.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --surface +hkl 1,0,4 +out surface.pdb")
  end subroutine createSurfaceHelp


  subroutine initialiseAction(a)

    implicit none
    type(actionTypeDef), target :: a

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0d0

    ! Local pointers
    actionCommand          => a % actionDetails
    firstAction            => a % firstAction

    outputFile             => a % outputFile
    imiller(1:3)           => a % integerVariables(1:3)

    call assignFlagValue(actionCommand,"+out",outputFile % fname,'surface.pdb')
    call initialiseFile(outputFile, outputFile % fname)

    call assignFlagValue(actionCommand,"+hkl",imiller,[0,0,0])

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Create Surface")
    call message(1,"...Miller indices",iv=imiller)

  end subroutine dumpScreenInfo

  subroutine computeAction()
    implicit none
    real(8), dimension(3) :: dij, sij, rdist, xcom
    real(8) :: rtmp
    integer :: nmols_new, natoms_new
    
    integer :: h, k, l, kmax, ii, kk, jj
    integer :: iatm, imol, jmol, idx
    integer, allocatable, dimension(:) :: order
    real(8), dimension(3,3) :: hmat_new, hinv_new
    real(8), dimension(3,3) :: hmat_fin, hinv_fin
    real(8), dimension(3,3) :: hnew

    type(moleculeTypeDef), allocatable, dimension(:) :: molecule_new
    real(8), allocatable, dimension(:,:) :: pos_new
    real(8), allocatable, dimension(:) :: chg_new
    character(cp), allocatable, dimension(:) :: lab_new
    integer, allocatable, dimension(:) :: atm2mol_new
    real(8), allocatable, dimension(:,:) :: zmol
   
    real(8), pointer, dimension(:,:) :: pos_fin
    character(cp), pointer, dimension(:) :: lab_fin
    integer, pointer, dimension(:) :: atm2mol_fin
    real(8), dimension(3) :: shift, shift0, dipole, rcom

    integer :: nmax

    ! Determine the vector normal to the surface from the reciprocal space vectors
    call defineSurfaceVectors(imiller, frame % hmat, hmat_new)

    ! make sure that the three axes form a right-handed cell
    call getInverseCellMatrix(hmat_new,hinv_new,rtmp)
    if (rtmp<0.0) then
      hmat_new(1:3,2)= -hmat_new(1:3,2)
      call getInverseCellMatrix(hmat_new,hinv_new,rtmp)
    endif

    if ( abs((rtmp/frame % volume)-nint((rtmp/frame % volume))) > 1.0d-4 ) &
      call message(-1,"Incommensurate volume for the new cell",r=rtmp/frame % volume)
  
    ! Look for molecules inside the new cell
    iatm=0
    nmols_new=0
    nmax = 100
    allocate(molecule_new(numberOfMolecules*nmax))
    allocate(zmol(3,numberOfMolecules*nmax))
    natoms_new=0
    allocate(pos_new(3,  frame % natoms*nmax))
    allocate(lab_new(    frame % natoms*nmax))
    allocate(chg_new(    frame % natoms*nmax))
    allocate(atm2mol_new(frame % natoms*nmax))

    do imol=1,numberOfMolecules
      xcom=0.0d0
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        xcom = xcom + frame % pos(1:3,listOfMolecules(imol)%listOfAtoms(idx))
      enddo
      xcom=xcom/listOfMolecules(imol) % numberOfAtoms

      kmax=10
      do h=-kmax,kmax
        do k=-kmax,kmax
          do l=-kmax,kmax
            dij = xcom + h * frame % hmat(:,1) + k * frame % hmat(:,2) + l * frame % hmat(:,3)

            sij(1) = hinv_new(1,1)*dij(1) + hinv_new(1,2)*dij(2) + hinv_new(1,3)*dij(3)
            sij(2) = hinv_new(2,1)*dij(1) + hinv_new(2,2)*dij(2) + hinv_new(2,3)*dij(3)
            sij(3) = hinv_new(3,1)*dij(1) + hinv_new(3,2)*dij(2) + hinv_new(3,3)*dij(3)
            if (sij(1)>=0.0d0 .and. sij(1)<=1.0d0 .and. &
                sij(2)>=0.0d0 .and. sij(2)<=1.0d0 .and. &
                sij(3)>=0.0d0 .and. sij(3)<=1.0d0) then

              ! Remove molecule it it overlaps with something already in the cell
              do jmol=1,nmols_new
                rdist(1:3) = dij(1:3) - zmol(1:3,jmol)
                sij(1) = hinv_new(1,1)*rdist(1) + hinv_new(1,2)*rdist(2) + hinv_new(1,3)*rdist(3)
                sij(2) = hinv_new(2,1)*rdist(1) + hinv_new(2,2)*rdist(2) + hinv_new(2,3)*rdist(3)
                sij(3) = hinv_new(3,1)*rdist(1) + hinv_new(3,2)*rdist(2) + hinv_new(3,3)*rdist(3)
                sij(1:3) = sij(1:3)-nint(sij(1:3))
                rdist(1) = hmat_new(1,1)*sij(1) + hmat_new(1,2)*sij(2) + hmat_new(1,3)*sij(3)
                rdist(2) = hmat_new(2,1)*sij(1) + hmat_new(2,2)*sij(2) + hmat_new(2,3)*sij(3)
                rdist(3) = hmat_new(3,1)*sij(1) + hmat_new(3,2)*sij(2) + hmat_new(3,3)*sij(3)
                rtmp = sqrt(sum(rdist*rdist))
                if (rtmp < 0.1d0) exit
              enddo
              if (jmol<=nmols_new) cycle

              ! New molecule found
              nmols_new = nmols_new+1
              natoms_new = natoms_new+listOfMolecules(imol) % numberOfAtoms

              ! Add the new molecule's info
              molecule_new(nmols_new) % ID            = listOfMolecules(imol) % ID
              molecule_new(nmols_new) % numberOfAtoms = listOfMolecules(imol) % numberOfAtoms
              allocate(molecule_new(nmols_new) % listOfAtoms(listOfMolecules(imol) % numberOfAtoms))
              zmol(1:3,nmols_new) = dij(1:3)

              ! Add the new atoms' info
              do idx=1,listOfMolecules(imol) % numberOfAtoms
                iatm = iatm + 1
                molecule_new(nmols_new) % listOfAtoms(idx) = iatm
                ii = listOfMolecules(imol) % listOfAtoms(idx)
                lab_new(iatm) = frame % lab(ii)
                pos_new(1:3,iatm) = frame % pos(1:3,ii) + h * frame % hmat(:,1) + k * frame % hmat(:,2) + l * frame % hmat(:,3)
                chg_new(iatm) = frame % chg(ii)
                atm2mol_new(iatm) = nmols_new
              enddo

            endif
          enddo
        enddo
      enddo
    enddo

! Rotate cell, the atoms' positions and molecules' com
    hmat_fin = hmat_new
    call makeUpperTriangularCell(hmat_new,pos_new,natoms_new)
    call makeUpperTriangularCell(hmat_fin,zmol,nmols_new)
    call straightenCell(hmat_fin,"z",hnew) ; hmat_fin=hnew
    call getInverseCellMatrix(hmat_fin,hinv_fin,rtmp)

    allocate(order(nmols_new))
    call mrgrnk(zmol(3,1:nmols_new) , order(1:nmols_new))
    ! call rnkpar(zmol(3,1:nmols_new) , order , nmols_new)

    allocate(pos_fin(3,natoms_new))
    allocate(lab_fin(natoms_new))
    allocate(atm2mol_fin(natoms_new))

    rtmp=0.0d0
    shift0=-1.0d0
    numberOfMolecules=0

    kk=0
    do ii=1,nmols_new
      idx=order(ii)
      ! Generate all the possible surface shifts
      ! from the centres of masses of the molecules
      shift(1:3) = (/0.0d0 , 0.0d0 , zmol(3,order(ii))/)

      ! New positions
      do iatm=1,natoms_new
        pos_fin(1:3,iatm) = pos_new(1:3,iatm) - shift(1:3) -0.0001 + hmat_fin(:,3)
      enddo

      ! Apply PBC to leave the right molecule on the surface
      do imol=1,nmols_new
        rcom=0.0d0
        do jj=1,molecule_new(imol) % numberOfAtoms
          iatm=molecule_new(imol) % listOfAtoms(jj)
          rcom(1:3)=rcom(1:3)+pos_fin(1:3,iatm)
        enddo
        rcom(1:3)=rcom(1:3)/molecule_new(imol) % numberOfAtoms
        sij(1) = hinv_fin(1,1)*rcom(1) + hinv_fin(1,2)*rcom(2) + hinv_fin(1,3)*rcom(3)
        sij(2) = hinv_fin(2,1)*rcom(1) + hinv_fin(2,2)*rcom(2) + hinv_fin(2,3)*rcom(3)
        sij(3) = hinv_fin(3,1)*rcom(1) + hinv_fin(3,2)*rcom(2) + hinv_fin(3,3)*rcom(3)
        sij(1:3) = sij(1:3)-floor(sij(1:3))
        dij(1) = hmat_fin(1,1)*sij(1) + hmat_fin(1,2)*sij(2) + hmat_fin(1,3)*sij(3)
        dij(2) = hmat_fin(2,1)*sij(1) + hmat_fin(2,2)*sij(2) + hmat_fin(2,3)*sij(3)
        dij(3) = hmat_fin(3,1)*sij(1) + hmat_fin(3,2)*sij(2) + hmat_fin(3,3)*sij(3)

        dij(1:3) = dij(1:3) - rcom(1:3)

        do jj=1,molecule_new(imol) % numberOfAtoms
          iatm=molecule_new(imol) % listOfAtoms(jj)
          pos_fin(1:3,iatm) = pos_fin(1:3,iatm) + dij(1:3)
        enddo
      enddo

      ! Calculate the cell dipole moment
      dipole=0.0d0
      do iatm=1,natoms_new
        dipole = dipole + chg_new(iatm)*(pos_fin(1:3,iatm)-(/5.,5.,5./))
      enddo
      where(dipole<1e-5) dipole=0.0d0

      kk=kk+1

      ! Convet to Debye when writing  
      ! write(str,'("  Cell dipole (Debye) for configuration ",i0)')kk
      ! call message(0,str,rv=dipole)

      call dumpCoordinates(outputFile % ftype, outputFile % funit, natoms_new, pos_fin, lab_new, hmat_fin)
    enddo
  
  end subroutine computeAction

  subroutine createSurface(a)
    implicit none
    type(actionTypeDef), target :: a

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then

      if (firstAction) then
        if (numberOfMolecules == 0) call runInternalAction("topology","NULL")

        call dumpScreenInfo()
        call computeAction()
  
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if

    end if

  end subroutine createSurface

end module moduleCreateSurface
