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
subroutine getNumberOfAtomsGromacs(iounit,natoms,hmat)
  use moduleStrings
  implicit none
  integer, intent(in) :: iounit
  integer, intent(out) :: natoms
  real(8), dimension(3,3), intent(inout) :: hmat

  integer :: ios, i, nw
  character(len=50) :: words(10)
  character(len=200) :: line
 
  natoms = 0

  read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
  read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
  read (line,*) natoms
  do i=1,natoms
    read(iounit,*)
  end do
  read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
  call parse(line," ",words,nw)
  if (nw == 6) then
    read(line,*) hmat(1,1),hmat(2,2),hmat(3,3)
  else if (nw == 9) then
    read(line,*)hmat(1,1),hmat(2,2),hmat(3,3),hmat(2,1),hmat(3,1),hmat(1,2),hmat(3,2),hmat(1,3),hmat(2,3)
  else
    hmat(1,1) = 1.d0
    hmat(2,2) = 1.d0
    hmat(3,3) = 1.d0
  end if

end subroutine getNumberOfAtomsGromacs

subroutine readCoordinatesGromacs(iounit,n,pos,lab,chg,hmat,go)
  use moduleStrings
  use moduleElements
  use moduleMessages
  implicit none
  integer, intent(in) :: iounit
  integer, intent(in) :: n
  real(8), dimension(3,n), intent(out) :: pos
  real(8), dimension(n), intent(out) :: chg
  character(len=4), dimension(n), intent(out) :: lab
  real(8), dimension(3,3), intent(inout) :: hmat
  logical, intent(out) :: go

  integer :: ios, i, nw, natoms, itmp, id
  character(len=50) :: words(10), str
  character(len=200) :: line
  real(8), dimension(3) :: ptmp

  go = .false.
  chg = 0.d0 

  read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
  read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
  read (line,*) natoms
  if (n/=natoms) call message(-1,"Wrong number of atoms in the GROMACS input file")
  do i=1,natoms
    read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
    read(line,'(i5,2a5,i5,3f8.3)')itmp, str, str, itmp, ptmp
    str = adjustl(str)
    id = getAtomicNumber(str)
    lab(itmp) = trim(atom(id)%lab)
    pos(1:3,itmp) = 10.d0*ptmp
  end do
  read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
  call parse(line," ",words,nw)
  if (nw == 6) then
    read(line,*) hmat(1,1),hmat(2,2),hmat(3,3)
  else if (nw == 9) then
    read(line,*)hmat(1,1),hmat(2,2),hmat(3,3),hmat(2,1),hmat(3,1),hmat(1,2),hmat(3,2),hmat(1,3),hmat(2,3)
  else
    hmat(1,1) = 1.d0
    hmat(2,2) = 1.d0
    hmat(3,3) = 1.d0
  end if
  
  hmat = hmat * 10.d0
  go = .true.

end subroutine readCoordinatesGromacs
