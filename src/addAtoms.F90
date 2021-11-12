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
module moduleAddAtoms
  use moduleActions 
  use moduleMessages
  
  implicit none

  public :: addAtoms, addAtomsHelp
  private

  type :: localWorkProcedure
    procedure (command), pointer, nopass :: work => null ()
  end type localWorkProcedure
  type(localWorkProcedure) :: addVirtualSite(100)

  integer :: lastFrameCalled = -11
  integer :: nCalled

contains

  subroutine addAtomsHelp()
    implicit none
    call message(0,"This action adds new particles to a frame.")
    call message(0,"It can be used to add solute/solvent atoms or molecules, or to add dummy particles to the system.")
    call message(0,"The new species are added randomly to the cell.")
    call message(0,"This action is a wrapper for different actions (solute, solvent, centre, dummy), which work in slightly different ways.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --rescale +cell 10")
    call message(0,"  gpta.x --i coord.pdb --top --add solute +atom Na,Cl +n 2,2 +rmin 4.0")
    call message(0,"  gpta.x --i coord.pdb --top --add solvent +f benzene.pdb +n 5 +rmin 5.0 --o new.pdb")
    call message(0,"  gpta.x --add solvent +f h2o.pdb +box 12. +n 57 +rmin 2 --top --o water.pdb")
    call message(0,"  gpta.x --i coord.pdb --top --add centre +f M1")
    call message(0,"  gpta.x --i benzene.pdb --top --add centre +f M1 +i 1,3,5,7,9,11")
    call message(0,"  gpta.x --i water.pdb --top --add dummy +i 1,2,3 +dist 0.1577 --o tip4p_ice.pdb")
  end subroutine addAtomsHelp
  
  subroutine addAtoms(a)
    use moduleStrings
    use moduleSystem
    use moduleDistances
    use moduleNeighbours

    use addDummyModule
    use addMoleculesModule
    use addElementsModule
    use addCentresModule

    implicit none
    type(actionTypeDef), target :: a
    character(len=10) :: option
    integer :: removeOverlap ! (0=no, 1=initial, 2=added)
    integer :: originalNumberOfAtoms
    integer :: originalNumberOfMolecules
    logical :: updateTopology

    if (frame % nframe > lastFrameCalled) then
      lastFrameCalled = frame % nframe
      nCalled = 0
    end if
    nCalled = nCalled + 1

    read(a % actionDetails,*) option

    removeOverlap = 0
    select case (option)
      case default
        call message(-1,"--add - nothing to do")

      case("box")
        block
          integer :: idx
          character(len=100) :: stringCell
          idx = index(a % actionDetails,"box") + 3
          read(a % actionDetails(idx:),'(a100)')stringCell
          call readCellFreeFormat(stringCell, frame % hmat)
        end block
        if (abs(frame % hmat(1,2)) + abs(frame % hmat(1,3)) + abs(frame % hmat(2,3)) .lt. 1.0d-6) then
          pbc_type = "ortho"
        else
          pbc_type = "tri"
        end if
        call initialisePBC(pbc_type)
        ! Initialise neighbours list
        call hmat2cell (frame % hmat, frame % cell, "DEG")
        call getInverseCellMatrix (frame % hmat, frame % hinv, frame % volume)
        call initialiseNeighboursList()
        return

      case("solute" , "solvent")
        updateTopology = .false.
        if (index(a % actionDetails,"+atom") > 0) then
          addVirtualSite(nCalled) % work => addElements

        else if (index(a % actionDetails,"+f") > 0) then
          addVirtualSite(nCalled) % work => addMoleculesFromFile

        else
          call message(-1,"--add solute/solvent requires either +f or +atom")

        end if

        if (option == "solute") then
          removeOverlap = 1
        else if (option == "solvent") then
          removeOverlap = 2
        end if

      case("dummy")
        updateTopology = .true.
        addVirtualSite(nCalled) % work => addDummyParticles

      case("centre")
        updateTopology = .true.
        addVirtualSite(nCalled) % work => addCentresParticles

    end select
    
    originalNumberOfAtoms = frame % natoms
    originalNumberOfMolecules = numberOfMolecules
    
    if (frameReadSuccessfully) then
      call addVirtualSite(nCalled) % work(a)

      if (originalNumberOfAtoms == 0) then
        call setUpNeigboursList()
        call initialiseNeighboursList()
      end if
  
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
      call updateNeighboursList(.true.)
      !      if (updateTopology) call runInternalAction("topology","+update +reorder")
      
      !call computeMoleculesCOM(originalNumberOfMolecules)
      if (removeOverlap /= 0) call removeClashes(frame,originalNumberOfMolecules,removeOverlap)
    end if

  end subroutine addAtoms

  !!!!!!!!!!!!!
  subroutine removeClashes(f,originalNumberOfMolecules,iType)
    use moduleSystem
    use moduleDistances
    implicit none
    type(frameTypeDef), intent(inout) :: f
    integer, intent(in) :: originalNumberOfMolecules
    integer, intent(in) :: iType

    integer :: imol, jmol, idx, jdx, iatm, jatm
    real(8) :: rmax, distanceCOM, dij(3), rtmp
    logical, allocatable, dimension(:) :: ldelete
    integer :: n0
    real(8), allocatable, dimension(:,:) :: xx

    integer :: nMol1, nMol2, nMol3, nMol4

    allocate(ldelete(numberOfMolecules),source=.false.)
    n0 = maxval(listOfMolecules(1:numberOfMolecules) % numberOfAtoms)
    rmax = (1.5d0)**2

    allocate(xx(3,n0))
    
    if (iType == 1) then
      nMol1 = 1
      nMol2 = originalNumberOfMolecules
      nMol3 = originalNumberOfMolecules + 1
      nMol4 = numberOfMolecules
    else if (iType == 2) then
      nMol1 = originalNumberOfMolecules + 1
      nMol2 = numberOfMolecules
      nMol3 = 1
      nMol4 = originalNumberOfMolecules
    end if

    toDelete : do imol=nMol1,nMol2
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        iatm = listOfMolecules(imol) % listOfAtoms(idx)
        xx(1:3,idx) = f % pos(1:3,iatm)
      end do

      do jmol=nMol3,nMol4
        dij(1:3) = listOfMolecules(jmol) % centreOfMass - listOfMolecules(imol) % centreOfMass
        distanceCOM = computeDistanceSquaredPBC(dij)
        if (distanceCOM < 100.d0) then
          do idx=1,listOfMolecules(imol) % numberOfAtoms
            do jdx=1,listOfMolecules(jmol) % numberOfAtoms
              jatm = listOfMolecules(jmol) % listOfAtoms(jdx)
              dij(1:3) = xx(1:3,idx) - f % pos(1:3,jatm)
              rtmp = computeDistanceSquaredPBC(dij)
              if (rtmp < rmax) then
                ldelete(imol) = .true.
                cycle toDelete
              end if
            end do
          end do
        end if
      end do
    end do toDelete

    n0 = 0
    ! Remove molecules 
    do imol=1,numberOfMolecules
      if (ldelete(imol)) cycle
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        iatm = listOfMolecules(imol) % listOfAtoms(idx)
        n0 = n0 + 1
        f % lab(n0) = f % lab(iatm)
        f % chg(n0) = f % chg(iatm)
        f % pos(:,n0) = f % pos(:,iatm)
      end do
    end do
    f % natoms = n0

  end subroutine removeClashes

end module moduleAddAtoms
