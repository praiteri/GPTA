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
subroutine parseActionLine()
  use moduleActions, only : numberOfActions, actionType, actionDetails
  use moduleMessages
  use moduleHelp

  implicit none

  integer :: iarg
  integer :: nargs
  character(len=200) :: word

  ! Number of command line arguments
  nargs=iargc()

  ! Count the commands
  numberOfActions=0
  iarg=0
  do while(iarg<nargs)
    iarg = iarg + 1
    call getarg(iarg,word)
    if (word(1:2)=="--") then
      numberOfActions=numberOfActions+1
      actionType(numberOfActions) = trim(word)
      actionDetails(numberOfActions) = ""
    else 
      if (word(1:1)=="-1") then
        if ( ((ichar(word(2:2)) >= 48) .and. (ichar(word(2:2)) <= 57)) .or. (ichar(word(2:2)) == 46) )then
          actionDetails(numberOfActions) = trim(actionDetails(numberOfActions))//" "//trim(word)
        else
          write(0,*)ichar(word(2:2)),word(2:2)
          call message(-1,"Unknown command",str=word)
        end if
      else
        actionDetails(numberOfActions) = trim(actionDetails(numberOfActions))//" "//trim(word)
      end if  
    end if
  enddo

  if (any(actionType == "--help")) call help()

  do iarg=1,numberOfActions
    if ( index(actionDetails(iarg),"+help") > 0) call commandHelp(actionType(iarg))
  end do

end subroutine parseActionLine
