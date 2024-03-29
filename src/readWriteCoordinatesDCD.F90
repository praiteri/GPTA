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
subroutine readCoordinatesDCD(uinp,first,natoms,pos,h,go)
  use moduleMessages 
  implicit none

  integer, intent(in) :: uinp 
  logical, intent(inout) :: first
  integer, intent(in) :: natoms
  real(8), dimension(3,natoms), intent(out) :: pos
  real(8), intent(out) :: h(3,3)
  logical, intent(out) :: go
!
! DCD header
  real(4)            :: delta
  integer(4)         :: nframes
  integer(4)         :: nstart
  integer(4)         :: nsanc
  integer(4)         :: nstep
  integer(4)         :: i4(4)
  integer(4)         :: namin
  integer(4), save   :: icell
  integer(4)         :: i9(9)
  integer(4)         :: ntitle
  character(len=4)   :: car4
  character(len=80)  :: car(10)
!
! DCD variables
  integer :: j
  integer(4) :: nn
  real(8), dimension(6) :: cell
  real(8), dimension(6) :: cell_tmp
  real(4), allocatable, dimension(:), save  :: x, y, z

  logical, save :: actionInitialisation=.true.
  character(len=20) :: str

  integer, save :: ncount

  if (actionInitialisation)then
    allocate(x(natoms),y(natoms),z(natoms))
    actionInitialisation=.false.
    ncount=0
  end if

  go=.true.
  ncount = ncount + 1
  if(first)then
    first = .false.
    read(uinp,end=112,err=112) car4, nframes, nstart, nsanc, nstep, i4, namin, delta, icell, i9
    read(uinp,end=111,err=111) ntitle, (car(j),j=1,ntitle)
    read(uinp,end=111,err=111) nn

    write(str,'(i0," /= ",i0)')nn,natoms
    if (nn/=natoms) call message(-1,"Wrong number of atoms in the DCD file",str=str)
  end if

  if(icell==1)then
    read(uinp,end=112,err=112)cell(1),cell(6),cell(2),cell(5),cell(4),cell(3)

    if ( all(cell(4:6) >= -1.d0) .and.  all(cell(4:6) <= 1.d0) ) then
      cell(4) = dacos(cell(4))
      cell(5) = dacos(cell(5))
      cell(6) = dacos(cell(6))
    end if
    
    cell_tmp=real(cell,8)
    call cell2hmat(cell_tmp,h)
  else
    h=0.d0
    h(1,1) = 1.d0
    h(2,2) = 1.d0
    h(3,3) = 1.d0
  end if

  read(uinp,end=111,err=111)x
  read(uinp,end=111,err=111)y
  read(uinp,end=111,err=111)z
  pos(1,:) = real(x,8)
  pos(2,:) = real(y,8)
  pos(3,:) = real(z,8)

  return

111 write(0,*) "DCD error reading frame ",ncount
112 go=.false.
  return

end subroutine readCoordinatesDCD

subroutine writeCoordinatesDCD(uout,write_header)
  use moduleMessages 
  use moduleSystem
  implicit none
  integer, intent(in)    :: uout
  logical, intent(inout) :: write_header

  integer :: j
  real*8 :: cell(6)
  real*4, allocatable, dimension (:)  :: x!, y, z
!
! DCD header
  real(4)            :: delta
  integer(4)         :: nframes
  integer(4)         :: nstart
  integer(4)         :: nsanc
  integer(4)         :: nstep
  integer(4)         :: i4(4)
  integer(4)         :: namin
  integer(4), save   :: icell
  integer(4)         :: i9(9)
  integer(4)         :: ntitle
  character(len=4)   :: car4
  character(len=80)  :: car(10)

  call hmat2cell(frame % hmat,cell,"COS")

  if(all(cell(1:3)>1.d0))then
    icell=1                 !! flag to indicate whether the cell is written in the dcd
  else
    icell=0
  end if

  ! DCD HEADER TAKEN FROM NAMD2-6
  if(write_header)then
    write_header = .false.
    car4="CORD"
    nframes=1                 !! frames in the output dcd (it'll give a worning if wrong)
    nstart=1                  !! first timestep in output
    nsanc=1                   !! dcd frequency
    nstep=1                   !! timestep of the last dcd
    i4=(/0,0,0,0/)            !! four useless integers :)
    namin=0                   !! endian type
    delta=0.01                !! timestep size
    i9=(/0,0,0,0,0,0,0,0,24/) !! nine useless integers :)
    ntitle=1                  !! number of dcd title lines
    car(1:ntitle)=(/"REMARKS DCD FILE CREATED BY GPTA"/)

    write(uout) car4, nframes, nstart, nsanc, nstep, i4, namin, delta, icell, i9
    write(uout) ntitle, (car(j),j=1,ntitle)
    write(uout) frame % natoms
  end if

  if(icell==1)then
    write(uout)cell(1),cell(6),cell(2),cell(5),cell(4),cell(3)
  end if

!  allocate(x(numberOfAtoms))
!  allocate(y(numberOfAtoms))
!  allocate(z(numberOfAtoms))
!
!  do j=1,numberOfAtoms
!    x(j)=real(frame % pos(1,j),4)
!    y(j)=real(frame % pos(2,j),4)
!    z(j)=real(frame % pos(3,j),4)
!  enddo
!  write(uout)x(1:numberOfAtoms)
!  write(uout)y(1:numberOfAtoms)
!  write(uout)z(1:numberOfAtoms)

  block
    integer :: numberOfAtoms
    numberOfAtoms = frame % natoms
    allocate(x(numberOfAtoms))
    do j=1,3
      x = real(frame % pos(j,1:numberOfAtoms),4)
      write(uout)x(1:numberOfAtoms)
    enddo
  end block

  return
end subroutine writeCoordinatesDCD
