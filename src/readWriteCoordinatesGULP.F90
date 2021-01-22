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
subroutine getNumberOfAtomsGULP(uinp,natoms,hmat)
  use moduleVariables
  use moduleMessages 
  use moduleSymmetry , only : rot, sgnops, sgncen, sgnops_tot, cen, get_spg_name
  use moduleStrings
  implicit none
  integer, intent(in)  :: uinp
  integer, intent(out) :: natoms
  real(8), dimension(3,3), intent(out) :: hmat

  character(len=100)  :: line
  integer :: nasym
  real(8), dimension(3) :: rtmp
  real(8) :: rqq, cell(6)
  character(cp) :: ctmp

  character(len=20)  :: str
  integer :: i, j, k, l, ierr
  logical :: lerr

  character(len=1), external :: uppercase

! GSAS routine
  CHARACTER*20  :: SPG
  INTEGER*4     :: LAUE,SGINV,SGLATT,SGUNIQ,IERR_LOCAL
  ! REAL*4        :: SGMTRX(24,3,3),SGTRNS(24,3),SGGEN(24)
  REAL*4        :: RT(5,4,25)
  INTEGER*4     :: JRT(3,5,24)
  INTEGER*4     :: SGPOL

  spg = "P 1"
  nasym = 0
  do
    ! read(uinp,'(a100)',end=546)line
    call readline(uinp,line,IERR_LOCAL)
    if (IERR_LOCAL/=0) goto 546

    ! Remove leading tab
    if ( index(line,achar(9)) > 0 ) line(index(line,achar(9)):index(line,achar(9))) = " "

    if (len_trim(line)==0) cycle
    if ( line(1:1) == "#") cycle

    read(line,*)str
    if (str(1:4)=='cell') then
      ! read(uinp,'(a100)',end=547)line
      call readline(uinp,line,ierr)
      if (ierr/=0) goto 546

      read(line,*)cell
      cell(4) = cell(4) * pi / 180.0d0
      cell(5) = cell(5) * pi / 180.0d0
      cell(6) = cell(6) * pi / 180.0d0
      call cell2hmat(cell,hmat)
    end if

    if (str(1:4)=='vect') then
      read(uinp,*,end=546)hmat(:,1)
      read(uinp,*,end=546)hmat(:,2)
      read(uinp,*,end=546)hmat(:,3)
    end if

    if (str(1:4)=='svec') then
      read(uinp,*,end=546)hmat(1:2,1)
      read(uinp,*,end=546)hmat(1:2,2)
      hmat(1,1)=abs(hmat(1,1))
      hmat(2,2)=abs(hmat(2,2))
    end if

    if (str(1:4)=='frac' .or. str(1:4)=='sfra' .or. str(1:4)=='cart') then
      do
        call readline(uinp,line,IERR_LOCAL)
        if (IERR_LOCAL/=0) goto 546
        call parseGulpCoordinates(line,ctmp,rtmp,lerr,rqq)
        if (lerr) exit
        nasym = nasym + 1
      enddo
      natoms = nasym
    end if
    read(line,*)str

    if (str(1:4) == "spac") then
      do
        call readline(uinp,line,IERR_LOCAL)
        if (IERR_LOCAL /= 0) call message(-1,"unknown space group "//trim(spg))
        read(line,'(a16)')str
        str=trim(adjustl(str))
        if ( ichar(str(1:1)) > ichar("0") .and. ichar(str(1:1)) <= ichar("9") ) then
          read(str,*)i
          spg = get_spg_name(i)
        else
          do i=1,len_trim(str)
            spg(i:i)=uppercase(str(i:i))
          enddo
        end if
        
        exit
      enddo
    end if
  enddo
  546 continue
  
  sgnops_tot = 0
  call sgroupnp(spg,laue,sguniq,sginv,sglatt,sgnops,sgpol,jrt,cen,sgncen,rt,IERR_LOCAL)
  do l=0,sginv
    do k=1,sgnops
      sgnops_tot = sgnops_tot + 1
      do j=1,3
        do i=1,3
          rot(i,j,sgnops_tot) = real(jrt(i,j,k),8) * (1-2*l)
        end do
        rot(j,4,sgnops_tot) = real(jrt(j,4,k),8)/12.d0
      end do
    end do
  enddo

  call message(0,"Parsing GULP input file")
  call message(0,"...Space group name",str=spg)
  call message(0,"......Number of representative symmetry operations",i=sgnops)
  call message(0,"......Space group has inversion centre (0=no; 1=yes)",i=sginv)
  call message(0,"......Total number of symmetry operations",i=sgnops_tot)
  if (ldebug) then
    call message(2)
    do k=1,sgnops_tot
      do i=1,3
        call message(0,"    Rotation matrix",rv=rot(1:3,i,k))
      enddo
      call message(0,"    Translation vector",rv=rot(1:3,4,k))
      if (k<sgnops_tot) call message(3)
    enddo
    call message(3)
    do l=1,sgncen
      call message(0,"  Centering vector(s)",rv=real(cen(1:3,l),8))
    enddo
  end if
  call message(2)

  natoms = nasym * sgnops_tot * sgncen

end subroutine getNumberOfAtomsGULP

subroutine readCoordinatesGULP(uinp,natoms,pos,label,hmat,chg,go)
  use moduleVariables
  use moduleMessages 
  use moduleSymmetry , only : rot, sgnops_tot, sgncen, cen, get_spg_name
  use moduleStrings, only : readline

  implicit none
  integer, intent(in) :: uinp
  integer, intent(inout) :: natoms
  real(8), dimension(3,natoms), intent(out)  :: pos
  character(cp), dimension(natoms), intent(out)  :: label
  real(8), dimension(3,3), intent(out) :: hmat
  real(8), dimension(natoms), intent(out) :: chg
  logical, intent(out) :: go

  real(8) :: cell(6)
  integer :: iatm, nasym
  character(len=100)  :: line
  character(len=16)  :: str(10)
  integer :: icoord = 0
  real(8) :: sij(3), ptmp(3)
  real(8), allocatable, dimension(:,:) :: dij
  character(cp), allocatable, dimension(:)  :: label1
  real(8), allocatable, dimension(:) :: chg1
  real(8) :: rqq
  real(8), dimension(3) :: rtmp
  character(cp) :: ctmp
  logical :: lerr
  integer :: i, j, k, l, idx, ierr
  logical, allocatable, dimension(:) :: lremove
  real(8) :: hinv(3,3), volume, rdist

  go=.false.
  nasym = 0

  allocate(label1(natoms))
  label=''
  label1=''
  allocate(chg1(natoms))
  chg1=0.0d0

  hmat=0.0d0
  pos=0.0d0
  allocate(dij(3,natoms))

  do
    call readline(uinp,line,ierr)
    if (ierr<0) goto 547
    if (ierr>0) goto 548

    read(line,*)str(1)
    if (str(1)(1:4)=='cell') then
      call readline(uinp,line,ierr)
      if (ierr/=0) goto 547

      read(line,*)cell
      cell(4) = cell(4) * pi / 180.0d0
      cell(5) = cell(5) * pi / 180.0d0
      cell(6) = cell(6) * pi / 180.0d0
      call cell2hmat(cell,hmat)

      call getInverseCellMatrix(hmat,hinv,volume)

    end if

    if (str(1)(1:4)=='vect') then
      read(uinp,*,end=547)hmat(:,1)
      read(uinp,*,end=547)hmat(:,2)
      read(uinp,*,end=547)hmat(:,3)
    end if

    if (str(1)(1:4)=='svec') then
      read(uinp,*,end=547)hmat(1:2,1)
      read(uinp,*,end=547)hmat(1:2,2)
      hmat(1,1)=abs(hmat(1,1))
      hmat(2,2)=abs(hmat(2,2))
    end if

    if (str(1)(1:4)=='frac' .or. str(1)(1:4)=='sfra' .or. str(1)(1:4)=='cart') then
      icoord=0
      if (str(1)(1:4)=='frac') icoord = 1
      if (str(1)(1:4)=='sfra') icoord = 2
      do
        call readline(uinp,line,ierr)
        if (ierr/=0) goto 547
        call parseGulpCoordinates(line,ctmp,rtmp,lerr,rqq)
        if (lerr) exit
        go=.true.
        nasym=nasym+1
        dij(1:3,nasym) = rtmp
        label1(nasym) = ctmp
        chg1(nasym) = rqq
      enddo
    end if

  enddo

547 continue

  iatm = 0
  do idx=1,nasym
    if (icoord == 0) then
      sij(1) = hinv(1,1)*dij(1,idx) + hinv(1,2)*dij(2,idx) + hinv(1,3)*dij(3,idx)
      sij(2) = hinv(2,1)*dij(1,idx) + hinv(2,2)*dij(2,idx) + hinv(2,3)*dij(3,idx)
      sij(3) = hinv(3,1)*dij(1,idx) + hinv(3,2)*dij(2,idx) + hinv(3,3)*dij(3,idx)
    else if (icoord == 1) then
      sij(1:3) = dij(1:3,idx)
    else
      write(0,*) "Wrong icoord in readGulp"
    end if

    do l=1,sgncen
      do k=1,sgnops_tot
        iatm=iatm+1
        do i=1,3
          ptmp(i) = rot(i,4,k)
          do j=1,3
            ptmp(i) = ptmp(i) + rot(i,j,k)*sij(j)
          enddo
        enddo
        ptmp(1:3) = ptmp(1:3) + real(cen(1:3,l),8)
        ptmp = ptmp - floor(ptmp)
        pos(1,iatm) = hmat(1,1)*ptmp(1) + hmat(1,2)*ptmp(2) + hmat(1,3)*ptmp(3)
        pos(2,iatm) = hmat(2,1)*ptmp(1) + hmat(2,2)*ptmp(2) + hmat(2,3)*ptmp(3)
        pos(3,iatm) = hmat(3,1)*ptmp(1) + hmat(3,2)*ptmp(2) + hmat(3,3)*ptmp(3)

        label(iatm) = label1(idx)
        chg(iatm) = chg1(idx)
      enddo
    enddo
  enddo

  allocate(lremove(natoms))
  lremove=.false.
  do i=1,natoms
    if (lremove(i)) cycle
    do j=i+1,natoms
      if (lremove(j)) cycle
      rtmp(1:3) = pos(1:3,j) - pos(1:3,i)
      rdist = distance(rtmp,hmat)
      if ( rdist < 1d-3 ) lremove(j) = .true.
    enddo
  enddo
  j=0
  do i=1,natoms
    if (lremove(i)) cycle
    j=j+1
    pos(1:3,j) = pos(1:3,i)
    label(  j) = label(  i)
  enddo
  natoms = j

  deallocate(dij)
  deallocate(label1)
  deallocate(chg1)

  return

548 go=.false.
  return

  contains

  function distance(dij,hmat) result(rdist)
    implicit none
    real(8), intent(in) :: dij(3), hmat(3,3)
    real(8) :: rdist
    real(8) :: sij(3), rtmp
    integer :: i, j, k
    rdist = sum(dij*dij)
    do i=-1,1
      do j=-1,1
        do k=-1,1
          if (i==0 .and. j==0 .and. k==0) cycle
          sij = dij + i*hmat(:,1) + j*hmat(:,2) + k*hmat(:,3)
          rtmp=sum(sij*sij)
          if (rtmp<rdist) rdist = rtmp
        enddo
      enddo
    enddo
    rdist = sqrt(rdist)

  end function distance
  
end subroutine readCoordinatesGULP

subroutine parseGulpCoordinates(line,lab,rtmp,lerr,chg)
  use moduleVariables, only : cp
  implicit none

  character(len=100), intent(inout) :: line
  character(cp) :: lab
  real(8), dimension(3), intent(out) :: rtmp
  real(8), intent(out) :: chg
  logical, intent(out) :: lerr

  integer :: i, iw, nw, ia, ib, ic, ix, iy, ipos
  logical :: noword
  integer, dimension(2,20) :: iword

  chg=0.0d0

  lerr=.false.

  ! do while( index(line,achar(9)) > 0)
  !   i = index(line,achar(9))
  !   line(i:i) = " "
  ! enddo

  nw=0
  noword=.true.
  iword=0
  do i=1,len_trim(line)
    if (noword) then
      if ( line(i:i) == " " ) then
        cycle
      else
        noword=.false.
        nw = nw + 1
        iword(1,nw) = i
      end if
    else
      if ( line(i:i) == " " ) then
        noword=.true.
        iword(2,nw) = i-1
      end if
    end if
  enddo
  if (iword(2,nw) == 0) iword(2,nw) = len_trim(line)

! read the label
  iw=1
  ia=iword(1,iw)
  ib=iword(2,iw)
  read(line( ia:ib ),*) lab
  if ( ib-ia < 2 ) then
    read(line( ia:ib ),*) lab
  else
! a ->  97 ;  z -> 122 ;  A ->  65 ;  Z ->  90
    if ( ichar(line( ia+2:ia+2 )) >= 97 .and. ichar(line( ia+2:ia+2 )) <= 122 .or. &
         ichar(line( ia+2:ia+2 )) >= 65 .and. ichar(line( ia+2:ia+2 )) <=  90 ) then
      lerr=.true.
      return
    else
      read(line( ia:ib ),*) lab
    end if
  end if

! read the coordinates
  ic=0
  do while (ic<3)
    iw=iw+1
    ia=iword(1,iw)
    ib=iword(2,iw)
    if ( line( ia:ia+2 ) == "cor" ) cycle
    if ( line( ia:ia+2 ) == "she" ) then
      lab = trim(lab)//"S"
      cycle
    end if
    ic=ic+1
    
    ipos=0
    do i=1,ib-ia
      if (line( ia+i:ia+i ) == "/" ) then
        ipos=ia+i
        exit
      end if
    enddo
    
    if ( ipos>0 ) then
      read(line( ia:ipos-1 ),*) ix
      read(line( ipos+1:ib ),*) iy
      rtmp(ic) = real(ix,8)/real(iy,8)
    else
      read(line( ia:ib ),*) rtmp(ic)
    end if
  enddo

  if (iw>=nw) return
  iw=iw+1
  ia=iword(1,iw)
  ib=iword(2,iw)
  read(line( ia:ib ),*) chg

  return
end subroutine parseGulpCoordinates

subroutine writeGulpCoordinates(uout)
  use moduleSystem 
  implicit none
  integer, intent(in) :: uout

  integer :: iatom, imol

  write(uout,'("vectors")')
  write(uout,'(3f20.6)') frame % hmat(:,1)
  write(uout,'(3f20.6)') frame % hmat(:,2)
  write(uout,'(3f20.6)') frame % hmat(:,3)

  write(uout,'("cartesian")')

  do iatom=1,frame % natoms
    write(uout,'(a8,3f20.6)') adjustl(trim(frame % lab(iatom))),frame % pos(1:3,iatom)
  enddo

  if (numberOfMolecules < numberOfAtoms) then
    write(uout,*)
    write(uout,'("# molecular connectivity")')
    do imol=1,numberOfMolecules
      write(uout,'("molatom ")',advance='no')
      do iatom=1,listOfMolecules(imol) % numberOfAtoms
        if (mod(iatom,6)==0) write(uout,'(" &"/,"        ")',advance='no')
        write(uout,'(i6)',advance='no')listOfMolecules(imol) % listOfAtoms(iatom)
      enddo
      write(uout,*)
    enddo
  end if

  return
end subroutine writeGulpCoordinates

subroutine writeGulpCoordinatesFractional(uout)
  use moduleSystem 
  implicit none
  integer, intent(in) :: uout

  integer :: iatom, imol

  write(uout,'("vectors")')
  write(uout,'(3f20.6)') frame % hmat(:,1)
  write(uout,'(3f20.6)') frame % hmat(:,2)
  write(uout,'(3f20.6)') frame % hmat(:,3)

  write(uout,'("fractional")')

  do iatom=1,frame % natoms
    write(uout,'(a8,3f20.6)') adjustl(trim(frame % lab(iatom))),frame % frac(1:3,iatom)
  enddo

  if (numberOfMolecules < numberOfAtoms) then
    write(uout,*)
    write(uout,'("# molecular connectivity")')
    do imol=1,numberOfMolecules
      write(uout,'("molatom ")',advance='no')
      do iatom=1,listOfMolecules(imol) % numberOfAtoms
        if (mod(iatom,6)==0) write(uout,'(" &"/,"        ")',advance='no')
        write(uout,'(i6)',advance='no')listOfMolecules(imol) % listOfAtoms(iatom)
      enddo
      write(uout,*)
    enddo
  end if

  return
end subroutine writeGulpCoordinatesFractional


