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
subroutine selectAtoms(nsel,cmd,a)
  use moduleVariables
  use moduleSystem 
  use moduleStrings
  use moduleMessages 
  implicit none
  integer, intent(in) :: nsel
  character(len=*), intent(inout) :: cmd
  type(actionTypeDef), intent(inout) :: a

  integer :: numberOfFlags, numberOfWords
  character(len=STRLEN), dimension(100) :: listOfWords
  character(len=STRLEN), dimension(100) :: listOfFlags
  
  character(len=STRLEN) :: selectionString
  character(len=STRLEN) :: flagString
  integer :: selectionID
  
  integer :: ntmp
  integer :: idx, jdx, indices(3)
  integer :: i, j, itmp, imol, iatm
  character(len=1) :: chr
  
  character(cp), dimension(100) :: localLabels
  character(len=10), dimension(100) :: localString

  selectionString = ''

  if (.not. allocated(a % isSelected)) then
    allocate(a % isSelected(frame % natoms , nsel), source=.false.)
  else
    a % isSelected =.false.
  end if 

  selectionID = 0
  call parse2(cmd,"+",listOfFlags,numberOfFlags)
  do idx=1,numberOfFlags

    ! select all atoms
    if (index(listOfFlags(idx),"all") > 0 .or. index(listOfFlags(idx),"ALL") > 0) then
        selectionID = selectionID + 1
        a % isSelected(: , selectionID) = .true.
        cycle
    end if

    ! select atoms by label
    if ( index(listOfFlags(idx),"+s ") > 0) then
      selectionString = trim(selectionString) // " " // trim(listOfFlags(idx))
      
      call extractFlag(listOfFlags(idx),"+s",flagString)
      call parse(flagString," ",listOfWords,numberOfWords)
      
      do jdx=1,numberOfWords
        selectionID = selectionID + 1
        if (selectionID > nsel) call message(-1,"Too many selection flags or blocks",str=selectionString)

        call parse(listOfWords(jdx),",",localLabels,ntmp)    
        do i=1,frame % natoms
          if (any(localLabels(1:ntmp)==frame % lab(i))) a % isSelected(i , selectionID) = .true.
        end do
      end do

    ! select atoms by index
    else if ( index(listOfFlags(idx),"+i ") > 0) then
      selectionString = trim(selectionString) // " " // trim(listOfFlags(idx))

      call extractFlag(listOfFlags(idx),"+i ",flagString)
      call parse(flagString," ",listOfWords,numberOfWords)

      do jdx=1,numberOfWords
        selectionID = selectionID + 1
        if (selectionID > nsel) call message(-1,"Too many selection flags or blocks",str=selectionString)

        call parse(trim(listOfWords(jdx)),",",localString,ntmp) 
        if (ntmp > 0) then
          do j=1,ntmp
            indices = extractSelectionIndices(localString(j))
            do i=indices(1),indices(2),indices(3)
              a % isSelected(i , selectionID) = .true.
            end do
          end do

        else
          call message(-1,"Cannot parse selection by index",str=listOfWords(jdx))
        end if
      end do

    ! select atoms by molecule
    else if ( index(listOfFlags(idx),"+mol ") > 0) then
      selectionString = trim(selectionString) // " " // trim(listOfFlags(idx))
  
      if (numberOfMolecules == 0) call message(-1,"Cannot perform selection by molecules, use --top")

      call extractFlag(listOfFlags(idx),"+mol ",flagString)
      call parse(flagString," ",listOfWords,numberOfWords)

      do jdx=1,numberOfWords
        selectionID = selectionID + 1
        if (selectionID > nsel) call message(-1,"Too many selection flags or blocks",str=selectionString)

        call parse(listOfWords(jdx),",",localLabels,ntmp)    

        do itmp=1,ntmp
          do imol=1,numberOfMolecules
            if (listOfMolecules(imol) % resname == localLabels(itmp)) then
              do iatm=1,listOfMolecules(imol) % numberOfAtoms
                i = listOfMolecules(imol) % listOfAtoms(iatm)
                a % isSelected(i , selectionID) = .true.
              end do
            end if
          end do
        end do

      end do

    else
      cycle
    end if
  end do

  if (selectionID < nsel) call message(-1,"Not enough selection flags or blocks",str=selectionString)

  ! mark the selection flags as used - not very elegant
  do 
    idx = maxval( [index(cmd,"+s") , index(cmd,"+i") , index(cmd,"+mol")] )
    if (idx == 0) exit
    cmd(idx:idx) = '%'
  end do

  if (a % firstAction) then
    call message(0,"...Atoms selection for "//trim(a % name))
    call message(0,"......Selection command",str=selectionString)
    do i=1,nsel
      write(chr,'(i1)') i
      call message(0,"......Atoms selected in group "//chr,i=count(a % isSelected(:,i)))
    end do
    call message(2)
  end if
  
end subroutine selectAtoms

