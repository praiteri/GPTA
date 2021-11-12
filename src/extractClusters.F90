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
module moduleExtractClustersAction
  use moduleVariables
  use moduleFiles
  use moduleMessages
  use moduleSystem
  use moduleStrings
  use moduleDistances
  
  implicit none

  public :: extractClusters, extractClustersHelp
  private

  character(:), pointer :: actionCommand
  type(fileTypeDef), pointer :: outputFile
  logical, pointer :: firstAction

  integer, pointer :: numberOfBins
  integer, pointer :: ID
  real(8), pointer :: clusterRadiusSQ

contains

  subroutine extractClustersHelp()
    implicit none
    call message(0,"This action extracts a spherical cluster of radius R around the selected atom.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --cluster +i 12 +r 6.0 +out clusters.xyz")
  end subroutine extractClustersHelp

  subroutine extractClusters(a)
    implicit none
    type(actionTypeDef), target :: a

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then

      if (firstAction) then
        call dumpScreenInfo()
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if
      
      call computeAction()
    end if

  end subroutine extractClusters

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    a % actionInitialisation = .false.

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    outputFile           => a % outputFile

    numberOfBins         => a % numberOfBins    
    ID                   => a % integerVariables(1)
    clusterRadiusSQ      => a % doubleVariables(1)

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .true.
    a % requiresNeighboursListDouble = .true.
    a % cutoffNeighboursList = 3.2d0

    call assignFlagValue(actionCommand,"+out",outputFile % fname,'cluster')

    call assignFlagValue(actionCommand,"+i",ID,0)
    if (ID==0) call message(-1,"cluster | central atoms must be specified with +i")

    call assignFlagValue(actionCommand,"+r",clusterRadiusSQ,0.d0)
    if (clusterRadiusSQ < 1d-3) call message(-1,"cluster | cluster radius must be specified with +r")
    clusterRadiusSQ = clusterRadiusSQ**2

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Extract Clusters Action")
    call message(0,"...Central atom",i=ID)
    call message(1,"...Radius",r=sqrt(clusterRadiusSQ))

  end subroutine dumpScreenInfo

  subroutine computeAction()
    implicit none
    integer :: idx
    integer :: iatm, jatm
    real(8) :: dij(3), dist
    logical, allocatable, dimension(:) :: lmol
    integer :: imol, nsel
    character(cp), allocatable, dimension(:) :: llab
    real(8), allocatable, dimension(:,:) :: lpos
    character(len=6) :: str
    type(fileTypeDef) :: localfile
    real(8), dimension(3) :: shift

    nsel = 0

    if (numberOfMolecules > 0) then
      allocate(lmol(numberOfMolecules), source=.false.)
      do iatm=1,frame % natoms
        dij(1:3) = frame % pos(1:3,iatm) - frame % pos(1:3,ID)
        dist = computeDistanceSquaredPBC(dij)
        if (dist < clusterRadiusSQ) then
          imol = atomToMoleculeIndex(iatm)
          if (lmol(imol)) cycle
          lmol(imol) = .true.
          nsel = nsel + listOfMolecules(imol) % numberOfAtoms
        end if
      end do

      allocate(llab(nsel))
      allocate(lpos(3,nsel))

      idx = 0
      do imol=1,numberOfMolecules
        if (lmol(imol)) then
          do iatm=1,listOfMolecules(imol) % numberOfAtoms
            idx = idx + 1
            jatm = listOfMolecules(imol) % listOfAtoms(iatm)
            lpos(:,idx) =  frame % pos(:,jatm)
            llab(idx) = frame % lab(jatm)
          end do
        end if
      end do

    else
      allocate(lmol(frame % natoms), source=.false.)
      do iatm=1,frame % natoms
        dij(1:3) = frame % pos(1:3,iatm) - frame % pos(1:3,ID)
        dist = computeDistanceSquaredPBC(dij)
        if (dist < clusterRadiusSQ) then
          lmol(iatm) = .true.
          nsel = nsel + 1
        end if
      end do
      allocate(llab(nsel))
      allocate(lpos(3,nsel))

      idx = 0
      do iatm=1,frame % natoms
        if (lmol(iatm)) then
          idx = idx + 1
          lpos(:,idx) =  frame % pos(:,iatm)
          llab(idx) = frame % lab(iatm)
        end if
      end do

    end if

    ! translate the cluster to the centre of the cell
    shift = -frame % pos(1:3,ID)
    shift = shift + 0.5d0 * (frame % hmat(:,1) + frame % hmat(:,2) + frame % hmat(:,3))

    do iatm=1,nsel
      lpos(1:3,iatm) = lpos(1:3,iatm) + shift(1:3)
    end do
    call wrapCoordinates(nsel,lpos)
    do iatm=1,nsel
      lpos(1:3,iatm) = lpos(1:3,iatm) - shift(1:3)
    end do

    write(str,'(i0)') frame % nframe
    localFile % fname = trim(outputFile % fname)//"."//trim(str)//".xyz"
    call initialiseFile(localFile, localFile % fname)
    ! call dumpXYZ(localFile % funit, nsel, lpos, llab)

    call dumpCoordinates("xyz", localFile % funit, nsel, lpos, llab)
 
    close(localFile % funit)

    deallocate(llab)
    deallocate(lpos)
    
  end subroutine computeAction

end module moduleExtractClustersAction
