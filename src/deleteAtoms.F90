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
module moduleDeleteAtoms
  
contains

  subroutine deleteAtoms(a)
    use moduleVariables
    use moduleStrings
    use moduleSystem 
    use moduleMessages 

    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    logical, pointer :: actionInitialisation
    logical, pointer :: firstAction

    logical, pointer :: labelSelection
    logical, pointer :: indexSelection
    logical, pointer :: moleculeSelection

    integer :: iatm, natoms

    actionCommand        => a % actionDetails

    ! Associate variables
    actionInitialisation => a % actionInitialisation
    firstAction          => a % firstAction

    labelSelection       => a % logicalVariables(1)
    indexSelection       => a % logicalVariables(2)
    moleculeSelection    => a % logicalVariables(3)

    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.

      call assignFlagValue(actionCommand,"+s",labelSelection,.false.)
      call assignFlagValue(actionCommand,"+i",indexSelection,.false.)
      call assignFlagValue(actionCommand,"+mol",moleculeSelection,.false.)
      if ( count([labelSelection,indexSelection,moleculeSelection]) /= 1) call message(-1,"--set | mixed selection not available")

      call checkUsedFlags(actionCommand)
      return
    end if

    if (frameReadSuccessfully) then

      ! Atoms' selection
      if (firstAction) then
        call message(0,"Deleting Atoms")

        if (keepFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if
        call selectAtoms(1,actionCommand,a)

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      ! Repeat selection for reactive trajectories
      else
        if (a % updateAtomsSelection) call selectAtoms(1,actionCommand,a)
      end if

      ! Delete atoms
      natoms = 0
      do iatm=1,frame % natoms
        if (a % isSelected(iatm, 1)) cycle
        natoms = natoms + 1
        frame % pos(1:3,natoms) = frame % pos(1:3,iatm)
        frame % lab(    natoms) = frame % lab(    iatm)
        frame % chg(    natoms) = frame % chg(    iatm)
      enddo
      frame % natoms = natoms

    end if

    if (endOfCoordinatesFiles) return

  end subroutine deleteAtoms

end module moduleDeleteAtoms
