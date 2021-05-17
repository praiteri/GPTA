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
subroutine defineVariables(cmd)
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleNeighbours 
  use moduleElements
  use moduleMessages
  use moduleNeighbours

  implicit none

  character(len=STRLEN), intent(in) :: cmd

  integer :: iword, i
  integer :: numberOfWords, nw
  character(len=STRLEN), dimension(100) :: listOfWords, words

  character(len=STRLEN) :: wstr
  ! character(len=STRLEN), external :: lowercase

  call parse(cmd, " ", listOfWords, numberOfWords)

  iword = 0
  do while (iword < numberOfWords)
    iword = iword + 1
    i = len_trim(listOfWords(iword)) + 1
    call lowercase(listOfWords(iword),wstr)
    wstr(i:) = ""
    select case(wstr)
      case default
        call message(-1,"--define - unknown variable",str=wstr)

      case ("seed")
        iword = iword + 1
        read(listOfWords(iword),*) randomNumberSeed

      case("rscale")
        iword = iword + 1
        read(listOfWords(iword),*) distanceScaling
        
      case("safedist")
        safeDistances = .true.

      case("verlet")
        forceVerletList = .true.

      case("rcov")

        iword = iword + 1
        call parse(listOfWords(iword),"=",words,nw)
        if (nw /= 2) then
          call message(-1,"--define - cannot read covalent radius",str=listOfWords(iword))
        else
          do i=1,nelement
            if (words(1) == atom(i) % lab) then
              read(words(2),*) atom(i) % rcov
              exit
            end if
          end do
        end if

      case ("bond")
        iword = iword + 1
        call parse(listOfWords(iword),"=",words,nw)
        if (nw /= 2) call message(-1,"Cannot parse special bond",str=listOfWords(iword))
        sbond_n = sbond_n + 1
        read(words(2),*) sbond_d(sbond_n)
        
        wstr = words(1)
        call parse(wstr,"-",words,nw)
        if (nw /= 2) then
          call parse(wstr,",",words,nw)
          if (nw /= 2) call message(-1,"Cannot parse special bond",str=wstr)
        end if
        read(words(1),*) sbond_l(1,sbond_n)
        read(words(2),*) sbond_l(2,sbond_n)

      case ("reax")
        keepFrameLabels = .false.
        keepFrameCharges = .false.
        
    end select
  end do

end subroutine defineVariables

subroutine dumpVariablesInfo()
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleNeighbours 
  use moduleMessages
  
  integer :: i

  call message(0,"Global variables settings")
  call message(0,"...Random number seed",i=randomNumberSeed)
  call message(0,"...Scaling factor for covalent bonds",r=distanceScaling)
  if (safeDistances) call message(0,"...Compute distances with loop over neighbour cells")
  if (sbond_n>0) then
    do i=1,sbond_n
      call message(0,"...Special bond length for "//trim(sbond_l(1,i))//"-"//trim(sbond_l(2,i)),r=sbond_d(i))
    end do
  end if
  call message(2)

end subroutine dumpVariablesInfo
