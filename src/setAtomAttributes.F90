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
module moduleAtomsAttributes
  use moduleMessages 
  
contains

  subroutine setAtomAttributesHelp()
    implicit none
    call message(0,"This action overrides some of the basic properties of the atoms in the system that can be used by other actions.")
    call message(0,"Examples:")
    call message(0,"gpta.x --i coord.pdb --top --set +s O2,H2 +l OW,HW")
    call message(0,"gpta.x --i coord.pdb --top --set +i 1:300:3 +l OW")
    call message(0,"gpta.x --i coord.pdb --top --set +mol M1 +l OW,HW,HW")
    call message(0,"gpta.x --i coord.pdb --top --set +s O2,H2 +tq 1,2")
  end subroutine setAtomAttributesHelp

  subroutine setAtomAttributes(a)
    use moduleVariables
    use moduleStrings
    use moduleSystem 
    use moduleElements
    use moduleOpenMM, only : qType

    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    logical, pointer :: actionInitialisation
    logical, pointer :: firstAction

    logical, pointer :: labelSelection
    logical, pointer :: indexSelection
    logical, pointer :: moleculeSelection
    logical, pointer :: renameAtoms
    logical, pointer :: rechargeAtoms
    logical, pointer :: tinkerAtoms

    integer :: iatm, natoms, idx, imol, itmp, jdx
    integer :: indices(3)
    integer, allocatable, dimension(:) :: atomFlag
    character(cp), allocatable, dimension(:) :: oldLabels, moleculeLabels
    character(len=40), allocatable, dimension(:) :: newLabels
    character(len=40) :: ctmp 
    character(len=STRLEN) :: str

    actionCommand        => a % actionDetails

    ! Associate variables
    actionInitialisation => a % actionInitialisation
    firstAction          => a % firstAction

    labelSelection       => a % logicalVariables(1)
    indexSelection       => a % logicalVariables(2)
    moleculeSelection    => a % logicalVariables(3)
    renameAtoms          => a % logicalVariables(5)
    rechargeAtoms        => a % logicalVariables(6)
    tinkerAtoms          => a % logicalVariables(7)

    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.

      call assignFlagValue(actionCommand,"+l",renameAtoms,.false.)
      call assignFlagValue(actionCommand,"+q",rechargeAtoms,.false.)
      call assignFlagValue(actionCommand,"+tq",tinkerAtoms,.false.)

      call assignFlagValue(actionCommand,"+s",labelSelection,.false.)
      call assignFlagValue(actionCommand,"+i",indexSelection,.false.)
      call assignFlagValue(actionCommand,"+mol",moleculeSelection,.false.)
      if ( count([labelSelection,indexSelection,moleculeSelection]) /= 1) call message(-1,"--set | mixed selection not available")

      call checkUsedFlags(actionCommand)
      return
    end if

    if (frameReadSuccessfully) then

      if (firstAction) then
        allocate(atomFlag(frame % natoms), source=0)

        if (labelSelection) then
          call assignFlagValue(actionCommand,"+s",oldLabels)
          if (oldLabels(1) == "all") then
            deallocate(oldLabels)
            allocate(oldLabels(numberOfUniqueLabels), source=listOfUniqueLabels)
          end if
          do iatm=1,frame % natoms
            do idx=1,size(oldLabels)
              if (frame % lab(iatm) == oldLabels(idx)) atomFlag(iatm) = idx
            enddo
          enddo

        else if (moleculeSelection) then
          call assignFlagValue(actionCommand,"+mol",moleculeLabels)
          idx = 0
          do itmp=1,size(moleculeLabels)
            do imol=1,numberOfMolecules
              if (listOfMolecules(imol) % resname == moleculeLabels(itmp)) then
                natoms = listOfMolecules(imol) % numberOfAtoms
                do iatm=1,natoms
                  jdx = listOfMolecules(imol) % listOfAtoms(iatm)
                  atomFlag(jdx) = idx + iatm
                end do
              end if
            end do
            idx = idx + natoms
          end do

        else if (indexSelection) then
          call assignFlagValue(actionCommand,"+i",newLabels)
          do idx=1,size(newLabels)
            indices = extractSelectionIndices(newLabels(idx))
            do jdx=indices(1),indices(2),indices(3)
              atomFlag(jdx) = 1
            end do
          end do

        else 
          call message(-1,"--set | Unknown selection")
        end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! call message(0,"Modify atoms' Attributes")
        if (labelSelection) then
          str = trim(oldLabels(1))
          do idx=2,size(oldLabels)
            str = trim(str)//" "//trim(oldLabels(idx))
          end do
          call message(0,"...Atoms selected by label",str=str)

        else if (indexSelection) then
          str = trim(newLabels(1))
          do idx=2,size(newLabels)
            str = trim(str)//" "//trim(newLabels(idx))
          end do
          call message(0,"...Atoms selected by index",str=str)

        else if (moleculeSelection) then
          str = trim(moleculeLabels(1))
          do idx=2,size(moleculeLabels)
            str = trim(str)//" "//trim(moleculeLabels(idx))
          end do
          call message(0,"...Atoms selected by molecule",str=str)

        end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (renameAtoms) then
          allocate(a % localLabels(frame % natoms))
          call assignFlagValue(actionCommand,"+l",newLabels)
          if (index(newLabels(1),"elem") > 0) then
            deallocate(newLabels)
            allocate(newLabels(numberOfUniqueLabels))
            newLabels = " "
            do idx=1,numberOfUniqueLabels
              newLabels(idx)(1:2) = listOfElements(idx)
            end do
          else
            if (moleculeSelection) then
              if (size(newLabels) /= maxval(atomFlag)) &
                call message(-1,"--set | +mol and +l have different number of elements",iv=[size(oldLabels),size(newLabels)])
            else
              if (size(newLabels) /= size(oldLabels)) then
                if (size(newLabels) /= 1) call message(-1,"--set | +s and +l have different number of elements",iv=[size(oldLabels),size(newLabels)])
                ctmp = newLabels(1)
                deallocate(newLabels)
                allocate(newLabels(size(oldLabels)), source=ctmp)
              end if
            end if
          end if

          do iatm=1,frame % natoms
            if (atomFlag(iatm)==0) then
              a % localLabels(iatm) = frame % lab(iatm)
            else
              a % localLabels(iatm) = newLabels(atomFlag(iatm))(1:4)
            end if
          enddo
          str = trim(newLabels(1))
          do idx=2,size(newLabels)
            str = trim(str)//" "//trim(newLabels(idx))
          enddo
          call message(0,"...Renaming selected atoms as",str=str)
        end if

        if (rechargeAtoms) then
          allocate(a % localCharges(frame % natoms))
          call assignFlagValue(actionCommand,"+q",newLabels)
          if (size(newLabels) /= maxval(atomFlag)) then
            if (size(newLabels) /= 1) call message(-1,"select +s and +q have different number of elements")
            ctmp = newLabels(1)
            deallocate(newLabels)
            allocate(newLabels(size(oldLabels)), source=ctmp)
          end if

          do iatm=1,frame % natoms
            if (atomFlag(iatm)==0) then
              a % localCharges(iatm) = frame % chg(iatm)
            else
              read(newLabels(atomFlag(iatm)),*) a % localCharges(iatm) 
            end if
          enddo

          str = trim(newLabels(1))
          do idx=2,size(newLabels)
            str = trim(str)//" "//trim(newLabels(idx))
          enddo
          call message(0,"...Setting atomic charges to",str=str)
        end if

        if (tinkerAtoms) then
          allocate(qType(frame % natoms))
          allocate(a % localIndices(frame % natoms))
          call assignFlagValue(actionCommand,"+tq",newLabels)
          if (size(newLabels) /= maxval(atomFlag)) then
            if (size(newLabels) /= 1) call message(-1,"select +s and +tq have different number of elements",iv=[size(oldLabels),size(newLabels),maxval(atomFlag)])
            ctmp = newLabels(1)
            deallocate(newLabels)
            allocate(newLabels(size(oldLabels)), source=ctmp)
          end if

          do iatm=1,frame % natoms
            if (atomFlag(iatm)==0) then
              call message(-1,"All atoms must have a TINKER type")
            else
              read(newLabels(atomFlag(iatm)),*) a % localIndices(iatm) 
            end if
          enddo

          str = trim(newLabels(1))
          do idx=2,size(newLabels)
            str = trim(str)//" "//trim(newLabels(idx))
          enddo
          call message(0,"...Setting TINKER types to",str=str)

        end if

        call message(2)

        call checkUsedFlags(actionCommand)
        firstAction = .false.
      
      end if

      if (renameAtoms) frame % lab = a % localLabels

      if (rechargeAtoms) frame % chg = a % localCharges

      if (tinkerAtoms) qType = a % localIndices

    end if

    if (endOfCoordinatesFiles) return

  end subroutine setAtomAttributes

end module moduleAtomsAttributes
