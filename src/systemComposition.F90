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
subroutine systemComposition(localFrame)
  use moduleVariables, only : fileTypeDef
  use moduleSystem , only : cp, frameTypeDef, numberOfUniqueLabels, listOfUniqueLabels, listOfElements
  use moduleMessages 
  use moduleElements, only : getElement
  implicit none

  type(frameTypeDef), intent(in) :: localFrame

  integer :: i, iatm
  integer :: nspecies
  integer, dimension(100) ::  ncount
  character(cp), dimension(100) :: species
  character(2), dimension(100) :: element
  character(5)  :: str


  call message(0,"System composition")
  call message(0,'...Number of atoms found',i=localFrame % natoms)

  ncount = 0
  nspecies = 1
  species(1) = localFrame % lab(1)
  element(1) = localFrame % element(1)

  ncount(1) = 1
  do iatm=2,localFrame % natoms
    if (all(species(1:nspecies) /= localFrame % lab(iatm))) then
      nspecies = nspecies + 1
      species(nspecies) = localFrame % lab(iatm)
      element(nspecies) = localFrame % element(iatm)

    end if
    where (species(1:nspecies) == localFrame % lab(iatm)) ncount(1:nspecies) = ncount(1:nspecies) + 1
  enddo
  do i=1,nspecies
    element(i) = getElement(species(i))
    write(str,'("[",a2,"] ")') adjustl(element(i))
    call message(0,"......Number of unique species "//str//trim(species(i)),i=ncount(i))
  enddo

  call message(3)

  numberOfUniqueLabels = nspecies
  allocate(listOfUniqueLabels(numberOfUniqueLabels))
  allocate(listOfElements(numberOfUniqueLabels))
  listOfUniqueLabels = species(1:numberOfUniqueLabels)
  listOfElements = element(1:numberOfUniqueLabels)

end subroutine systemComposition
