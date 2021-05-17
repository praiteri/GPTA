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
subroutine debugTiming_old(idx)
  use moduleVariables
  use moduleMessages
  implicit none
  integer, intent(in) :: idx
  logical, save :: firstTimeIn = .true.
  integer, save :: numberOfSplits
  real(8), save, allocatable, dimension(:) :: timeTally
  real(8), save, allocatable, dimension(:) :: startClock
  logical, save, allocatable, dimension(:) :: splitOn
  
  integer :: i
  character(len=2) :: str

  if (firstTimeIn) then
    firstTimeIn = .false.
    numberOfSplits = 20
    allocate(startClock(numberOfSplits))
    allocate(splitOn(numberOfSplits), source=.false.)
    allocate(timeTally(numberOfSplits), source=-1.d0)
  end if

  if (idx<0) then
    if (all(timeTally < 0.d0)) return
    call message(0,"---> Partial Debug Timings <---")
    do i=1,numberOfSplits
      if (timeTally(i) < 0.d0) exit
      write(str,'(i0)')i
      if (splitOn(i)) then
        call message(0,"...Unended split time "//trim(str))
      else
        call message(0,"...Timing for block "//trim(str),r=timeTally(i))
      end if
    end do
    call message(2)
    return
  end if
  
  if (splitOn(idx)) then
    timeTally(idx) = timeTally(idx) + timing() - startClock(idx)
    splitOn(idx) = .false.
  else
    splitOn(idx) = .true.
    startClock(idx) = timing()
    if (timeTally(idx)<0.d0) timeTally(idx)=0.d0
  end if

end subroutine debugTiming_old

subroutine debugTiming(idx)
  use moduleVariables, only : timing
  use moduleMessages
  implicit none
  integer, intent(in) :: idx

  logical, save :: firstTimeIn = .true.
  integer, save :: numberOfSplits
  real(8), save, allocatable, dimension(:) :: timeTally
  real(8), save :: startClock, stopClock
  integer, save :: splitID = 0
  integer :: i
  character(len=2) :: str

  if (firstTimeIn) then
    firstTimeIn = .false.
    numberOfSplits = 20
    allocate(timeTally(numberOfSplits), source=-1.d-8)
  end if

  if (idx < 0) then
    do i=1,numberOfSplits
      if (timeTally(i) < 0.d0) exit
      write(str,'(i0)')i
      call message(0,"...Timing for block "//trim(str),r=timeTally(i))
    end do
    call message(2)

  else if (idx == 0) then
    splitID = 1
    startClock = timing()

  else
    stopClock = timing()
    timeTally(splitID) = timeTally(splitID)  + stopClock - startClock
    splitID = splitID + 1
    startClock = timing()
  
  end if


end subroutine debugTiming

function uppercase(strin) result(strout)
  implicit none
  character,intent(in) :: strin
  character            :: strout
  if (ichar(strin).ge.ichar('a') .and. ichar(strin).le.ichar('z')) then
    strout = char( ichar('A')+ichar(strin)-ichar("a"))
  else
    strout = strin
  end if
  return
end function uppercase

! function uppercase(str) result(ucstr)

! ! convert string to upper case

! character (len=*):: str
! character (len=len_trim(str)):: ucstr

! ilen=len_trim(str)
! ioffset=iachar('A')-iachar('a')
! iquote=0
! ucstr=str
! do i=1,ilen
!   iav=iachar(str(i:i))
!   if(iquote==0 .and. (iav==34 .or.iav==39)) then
!     iquote=1
!     iqc=iav
!     cycle
!   end if
!   if(iquote==1 .and. iav==iqc) then
!     iquote=0
!     cycle
!   end if
!   if (iquote==1) cycle
!   if(iav >= iachar('a') .and. iav <= iachar('z')) then
!     ucstr(i:i)=achar(iav+ioffset)
!   else
!     ucstr(i:i)=str(i:i)
!   end if
! end do
! return

! end function uppercase

subroutine lowercase(str,lcstr)

! convert string to lower case

character (len=*):: str
character (len=len_trim(str)):: lcstr

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')
iquote=0
lcstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('A') .and. iav <= iachar('Z')) then
    lcstr(i:i)=achar(iav-ioffset)
  else
    lcstr(i:i)=str(i:i)
  end if
end do

end subroutine lowercase

subroutine cell2hmat(cell,h)
  use moduleVariables, only : pi, pih
  implicit none
  real(8) :: cell(6)
  real(8) :: h(3,3)

  if (any(cell(4:6)>pi)) then
    cell(4:6) = cell(4:6)*pi/180.0d0
  end if

  if (cell(2)<1.0d-6) cell(2) = cell(1)
  if (cell(3)<1.0d-6) cell(3) = cell(2)
  if (cell(4)<1.0d-6) cell(4) = pih
  if (cell(5)<1.0d-6) cell(5) = pih
  if (cell(6)<1.0d-6) cell(6) = pih

!!p cell  :: a, b, c, alpha, beta, gamma in radians
  h = 0d0
  h(1,1) = cell(1)
  h(1,2) = cell(2) * cos(cell(6))
  h(2,2) = cell(2) * sin(cell(6))
  h(1,3) = cell(3) * cos(cell(5))
  h(2,3) = (cell(3) * cell(2) * cos(cell(4)) - h(1,2) * h(1,3)) / h(2,2)
  h(3,3) = sqrt(cell(3)*cell(3) - h(1,3)*h(1,3) - h(2,3)*h(2,3))
!
  if ( abs(h(1,2)) < 1.d-4 ) h(1,2) = 0.0d0
  if ( abs(h(1,3)) < 1.d-4 ) h(1,3) = 0.0d0
  if ( abs(h(2,3)) < 1.d-4 ) h(2,3) = 0.0d0
!
  return
end subroutine cell2hmat

subroutine hmat2cell(hmat,cell,angle)
  use moduleVariables, only : pi, pih
  implicit none
  real(8) :: hmat(3,3)
  real(8) :: cell(6)
  character(len=3)  :: angle

!!p a, b, c
  cell(1) = sqrt(dot_product(hmat(:,1),hmat(:,1)))
  cell(2) = sqrt(dot_product(hmat(:,2),hmat(:,2)))
  cell(3) = sqrt(dot_product(hmat(:,3),hmat(:,3)))
!!p alpha, beta, gamma in degrees
  if(angle=="DEG")then
    cell(4) = acos(dot_product(hmat(:,2),hmat(:,3))/cell(2)/cell(3))/pi*180d0
    cell(5) = acos(dot_product(hmat(:,1),hmat(:,3))/cell(1)/cell(3))/pi*180d0
    cell(6) = acos(dot_product(hmat(:,1),hmat(:,2))/cell(1)/cell(2))/pi*180d0
!!p alpha, beta, gamma in radians
  else if(angle=="RAD")then
    cell(4) = acos(dot_product(hmat(:,2),hmat(:,3))/cell(2)/cell(3))
    cell(5) = acos(dot_product(hmat(:,1),hmat(:,3))/cell(1)/cell(3))
    cell(6) = acos(dot_product(hmat(:,1),hmat(:,2))/cell(1)/cell(2))
!!p cosines of alpha, beta, gamma
  else if(angle=="COS")then
    cell(4) = dot_product(hmat(:,2),hmat(:,3))/cell(2)/cell(3)
    cell(5) = dot_product(hmat(:,1),hmat(:,3))/cell(1)/cell(3)
    cell(6) = dot_product(hmat(:,1),hmat(:,2))/cell(1)/cell(2)
  end if
  return
end subroutine hmat2cell

subroutine getInverseCellMatrix (hmat,hmati,deth)
  implicit none
  real(8), dimension(3,3)  :: hmat, hmati
  real(8) :: deth
  real(8) :: odet

  deth = &
       hmat(1,1) * ( hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2) ) + &
       hmat(1,2) * ( hmat(2,3)*hmat(3,1)-hmat(2,1)*hmat(3,3) ) + &
       hmat(1,3) * ( hmat(2,1)*hmat(3,2)-hmat(2,2)*hmat(3,1) )
  if ( deth < 1.d-4 ) return
  odet = 1.d0 / deth
  hmati(1,1) = (hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2))*odet
  hmati(2,2) = (hmat(1,1)*hmat(3,3)-hmat(1,3)*hmat(3,1))*odet
  hmati(3,3) = (hmat(1,1)*hmat(2,2)-hmat(1,2)*hmat(2,1))*odet
  hmati(1,2) = (hmat(1,3)*hmat(3,2)-hmat(1,2)*hmat(3,3))*odet
  hmati(2,1) = (hmat(3,1)*hmat(2,3)-hmat(2,1)*hmat(3,3))*odet
  hmati(1,3) = (hmat(1,2)*hmat(2,3)-hmat(1,3)*hmat(2,2))*odet
  hmati(3,1) = (hmat(2,1)*hmat(3,2)-hmat(3,1)*hmat(2,2))*odet
  hmati(2,3) = (hmat(1,3)*hmat(2,1)-hmat(2,3)*hmat(1,1))*odet
  hmati(3,2) = (hmat(3,1)*hmat(1,2)-hmat(3,2)*hmat(1,1))*odet

  return
end subroutine getInverseCellMatrix
 
subroutine makeUpperTriangularCell(hmat,pos,nn)
  implicit none
  real(8), intent(inout) :: hmat(3,3)
  real(8), dimension(3,nn), intent(inout) :: pos
  real(8) :: hnew(3,3), rotchi(3,3)
  real(8) :: u(3), v(3), rotmat(3,3)
  integer :: nn
  integer :: i

  u = hmat(:,1)
  v = hmat(:,2)
  call rotateCell(u,v,rotmat)
  do i=1,3
    hnew(:,i)=matmul(rotmat, hmat(:,i))
  enddo

  ! One more rotation if the b vector points towards -y
  if ( dot_product(hnew(:,2) , [0.d0,1.d0,0.d0]) < 0.d0 ) then
    rotchi(1,1) = 1.0d0
    rotchi(1,2) = 0.0d0
    rotchi(1,3) = 0.0d0
    rotchi(2,1) = 0.0d0
    rotchi(2,2) = -1.d0
    rotchi(2,3) = 0.d0
    rotchi(3,1) = 0.0d0
    rotchi(3,2) = 0.d0
    rotchi(3,3) = -1.d0
  
    rotmat  = matmul(rotchi,rotmat)
    do i=1,3
      hnew(:,i)=matmul(rotmat, hmat(:,i))
    enddo
  end if

  do i=1,nn
    u(1:3)=pos(1:3,i)
    pos(1:3,i)=matmul(rotmat,u)
  enddo
  hmat = hnew

  return
end subroutine makeUpperTriangularCell

subroutine rotateCell(u,v,rot)
  implicit none
  real(8)  :: u(3), v(3)
  real(8)  :: normu, u1(3)
  real(8)  :: cospsi, costheta, coschi, sinchi, sintheta, sinpsi
  real(8)  :: rottheta(3,3), rotpsi(3,3), rot1(3,3), rotchi(3,3)
  real(8)  :: vnew(3), rot(3,3)
  integer :: i

  rot=0.0d0
  do i=1,3
    rot(i,i) = 1.0d0
  enddo

! Align a with x
  if(abs(u(2))>1.d-6 .or. abs(u(3))>1.d-6)then
    normu = sqrt(dot_product(u,u))

    costheta = u(2)/sqrt(u(2)**2+u(3)**2)
    sintheta = u(3)/sqrt(u(2)**2+u(3)**2)

    rottheta(1,1) = 1.0d0
    rottheta(1,2) = 0.0d0
    rottheta(1,3) = 0.0d0
    rottheta(2,1) = 0.0d0
    rottheta(2,2) = costheta
    rottheta(2,3) = sintheta
    rottheta(3,1) = 0.0d0
    rottheta(3,2) = -sintheta
    rottheta(3,3) = costheta

    u1 = MATMUL(rottheta, u )

    cospsi = u1(1)/normu
    sinpsi = u1(2)/normu

    rotpsi(1,1) = cospsi
    rotpsi(1,2) = sinpsi
    rotpsi(1,3) = 0.0d0
    rotpsi(2,1) = -sinpsi
    rotpsi(2,2) = cospsi
    rotpsi(2,3) = 0.0d0
    rotpsi(3,1) = 0.0d0
    rotpsi(3,2) = 0.0d0
    rotpsi(3,3) = 1.0d0

    rot1 = matmul(rotpsi,rottheta)
  else
    rot1=0.0d0
    do i=1,3
      rot1(i,i)=1.0d0
    enddo
  end if
  vnew = matmul(rot1, v)

!!p rotate around x in order to place b on the xy plane

  if( abs(vnew(3)) > 1.d-6)then
    coschi = vnew(2)/sqrt(vnew(2)**2+vnew(3)**2)
    sinchi = vnew(3)/sqrt(vnew(2)**2+vnew(3)**2)
  else 
    coschi = 1.d0
    sinchi = 0.d0
  end if

  rotchi(1,1) = 1.0d0
  rotchi(1,2) = 0.0d0
  rotchi(1,3) = 0.0d0
  rotchi(2,1) = 0.0d0
  rotchi(2,2) = coschi
  rotchi(2,3) = sinchi
  rotchi(3,1) = 0.0d0
  rotchi(3,2) = -sinchi
  rotchi(3,3) = coschi
  
  rot  = matmul(rotchi,rot1)

  return
end subroutine rotateCell

subroutine straightenCell(hmat,dir,hnew)
  use moduleMessages
  implicit none
  real(8), intent(in) :: hmat(3,3)
  character(len=*) :: dir
  real(8), intent(out), target :: hnew(3,3)

  integer :: nn, ia, ib
  real(8), dimension(3) :: sij, dij
  real(8), pointer, dimension(:) :: vec1, vec2, vec3
  real(8) :: len0, len1

  hnew = hmat
  if (dir == "x") then
    vec1 => hnew(:,3)
    vec2 => hnew(:,2)
    vec3 => hnew(:,1)
  elseif (dir == "y") then
    vec1 => hnew(:,1)
    vec2 => hnew(:,3)
    vec3 => hnew(:,2)
  elseif (dir == "z") then
    vec1 => hnew(:,1)
    vec2 => hnew(:,2)
    vec3 => hnew(:,3)
  else
    call message(-1,"Unknown axis to straighten")
  endif

  nn=10
  sij(1:3) = vec3(1:3)
  len0=sum(sij(1:3)*sij(1:3))

  do ia=-nn,nn
    dij(1:3) = vec3(1:3) + real(ia,8) * vec1(1:3)
    len1 = sum(dij(1:3)*dij(1:3))
    if (len1<len0) then
      sij(1:3) = dij(1:3)
      len0 = len1
    endif
  enddo

  do ib=-nn,nn
    dij(1:3) = vec3(1:3) + real(ib,8) * vec2(1:3)
    len1 = sum(dij(1:3)*dij(1:3))
    if (len1<len0) then
      sij(1:3) = dij(1:3)
      len0 = len1
    endif
  enddo
  vec3(1:3) = sij(1:3)

  return
end subroutine straightenCell

function cross_product(a,b) result(c)
  implicit none
  real(8), dimension (3), intent(in)  :: a
  real(8), dimension (3), intent(in)  :: b
  real(8), dimension (3) :: c

  c(1) =  a(2)*b(3) - a(3)*b(2)
  c(2) = -a(1)*b(3) + a(3)*b(1)
  c(3) =  a(1)*b(2) - a(2)*b(1)

  return
end function cross_product

subroutine readCellFreeFormat(string, hmat)
  use moduleStrings
  use moduleMessages
  implicit none
  character(len=*), intent(in) :: string
  real(8), dimension(3,3), intent(out) :: hmat

  integer :: numberOfWords
  character(len=STRLEN), dimension(100) :: listOfWords
  real(8), dimension(6) :: cell

  call parse(string,",",listOfWords,numberOfWords)
  write(0,*)"#"//trim(string)//"#"
  write(0,*)numberOfWords
  if (numberOfWords == 1) then
    read(listOfWords(1),*) hmat(1,1)
    hmat(2,2) = hmat(1,1)
    hmat(3,3) = hmat(1,1)

  else if (numberOfWords == 3) then
    read(listOfWords(1),*) hmat(1,1)
    read(listOfWords(2),*) hmat(2,2)
    read(listOfWords(3),*) hmat(3,3)

  else if (numberOfWords == 6) then
    read(listOfWords(1),*) cell(1)
    read(listOfWords(2),*) cell(2)
    read(listOfWords(3),*) cell(3)
    read(listOfWords(4),*) cell(4)
    read(listOfWords(5),*) cell(5)
    read(listOfWords(6),*) cell(6)
    call cell2hmat(cell,hmat)

  else if (numberOfWords == 9) then
    read(listOfWords(1),*) hmat(1,1)
    read(listOfWords(2),*) hmat(2,1)
    read(listOfWords(3),*) hmat(3,1)
    read(listOfWords(4),*) hmat(1,2)
    read(listOfWords(5),*) hmat(2,2)
    read(listOfWords(6),*) hmat(3,2)
    read(listOfWords(7),*) hmat(1,3)
    read(listOfWords(8),*) hmat(2,3)
    read(listOfWords(9),*) hmat(3,3)

  else
    call message(-1,"Cannot read cell")
  endif

end subroutine readCellFreeFormat

! writes all permutations of the indices of a
! a must start as a sequence (1,2,...,n)
  function nextp ( n, a )
    logical :: nextp
    integer,intent(in) :: n
    integer,dimension(n),intent(inout) :: a

    integer i,j,k,t

    i = n-1
 10 if ( a(i) .lt. a(i+1) ) goto 20
    i = i-1
    if ( i .eq. 0 ) goto 20
    goto 10
 20 j = i+1
    k = n
 30 t = a(j)
    a(j) = a(k)
    a(k) = t
    j = j+1
    k = k-1
    if ( j .lt. k ) goto 30
    j = i
    if (j .ne. 0 ) goto 40

    nextp = .false.

    return

 40 j = j+1
    if ( a(j) .lt. a(i) ) goto 40
    t = a(i)
    a(i) = a(j)
    a(j) = t

    nextp = .true.

    return
  end function
  
subroutine  rotateMolecule(vec,theta,natoms,pos,pivot)
  implicit none
  real(8), intent(in) :: vec(3), theta
  integer, intent(in) :: natoms
  real(8), dimension(3,natoms), intent(inout) :: pos
  integer, intent(in) :: pivot

  integer :: i
  real(8) :: v(3), v2(3), costheta, sintheta, u(3), rcom(3)
  real(8), dimension (3,3) :: rotmat
  real(8) :: norm

  norm=sqrt(sum(vec*vec))
  v=vec/norm
  v2=v*v
  costheta=cos(theta)
  sintheta=sin(theta)

  rotmat(1,1)=v2(1)+(1.-v2(1))*costheta
  rotmat(1,2)= v(1)*v(2)*(1-costheta)-v(3)*sintheta
  rotmat(1,3)= v(1)*v(3)*(1-costheta)+v(2)*sintheta

  rotmat(2,1)= v(1)*v(2)*(1-costheta)+v(3)*sintheta
  rotmat(2,2)=v2(2)+(1.-v2(2))*costheta
  rotmat(2,3)= v(2)*v(3)*(1-costheta)-v(1)*sintheta

  rotmat(3,1)= v(1)*v(3)*(1-costheta)-v(2)*sintheta
  rotmat(3,2)= v(2)*v(3)*(1-costheta)+v(1)*sintheta
  rotmat(3,3)=v2(3)+(1.-v2(3))*costheta

  ! rotates around one atom
  if (pivot >= 1 .and. pivot <= natoms) then
    rcom = pos(1:3,pivot)

  ! rotates around the centre of mass
  elseif (pivot == 0) then
    rcom=0.0d0
    do i=1,natoms
      rcom=rcom+pos(1:3,i)
    enddo
    rcom=rcom/real(natoms,8)

  ! rotates around the origin
  else
    rcom = 0.d0
  endif

  do i=1,natoms
    u=pos(1:3,i)-rcom(1:3)
    pos(1:3,i)=matmul(rotmat,u(1:3))+rcom(1:3)
  enddo

end subroutine rotateMolecule

subroutine linearRegression(n,x,y,alpha,beta,rfact)
  ! y = alpha + beta*x
  implicit none
  integer, intent(in) :: n
  real(8), dimension(n), intent(in) :: x
  real(8), dimension(n), intent(in) :: y
  real(8), intent(out) :: alpha, beta, rfact

  real(8) :: xavg, yavg, xyavg, x2avg, y2avg
  real(8) :: ssxx, ssyy, ssxy
  integer :: i

  xavg  = 0.0d0
  yavg  = 0.0d0
  xyavg = 0.0d0
  x2avg = 0.0d0
  y2avg = 0.0d0
  do i=1,n
    xavg  = xavg  + x(i)
    yavg  = yavg  + y(i)
    xyavg = xyavg + x(i)*y(i)
    x2avg = x2avg + x(i)**2
    y2avg = y2avg + y(i)**2
  enddo

  xavg  = xavg  / real(n,8)
  yavg  = yavg  / real(n,8)

  ssxx = x2avg - real(n,8) * xavg**2
  ssyy = y2avg - real(n,8) * yavg**2
  ssxy = xyavg - real(n,8) * xavg*yavg

  beta = ssxy/ssxx
  alpha = yavg - xavg*(ssxy/ssxx)
  rfact = ssxy**2 / ssxx / ssyy

end subroutine linearRegression

subroutine createSelectionList(a,n)
  use moduleVariables, only : actionTypeDef
  use moduleSystem, only : frame 
  use moduleMessages
  implicit none
  type(actionTypeDef), target :: a
  integer, intent(in) :: n
  integer :: ntmp, i, idx
  if (.not. allocated(a % isSelected)) call message(-1,"Cannot creat selection list",str=a % actionDetails)

  ntmp = count(a % isSelected(:,1))
  do idx=1,n
    ntmp = max(ntmp , count(a % isSelected(:,idx)))
  end do

  if (.not. allocated(a % idxSelection)) allocate(a % idxSelection(ntmp,n))
  
  do idx=1,n
    ntmp = 0
    do i=1,frame % natoms
      if (.not. a % isSelected(i, idx)) cycle
      ntmp = ntmp + 1
      a % idxSelection(ntmp, idx) = i
    end do
  end do
end subroutine createSelectionList

subroutine defineSurfaceVectors(imiller, hmat, hsurf)
  ! use m_rnkpar
  use m_mrgrnk
  implicit none
  integer, dimension(3), intent(in) :: imiller
  real(8), dimension(3,3), intent(in) :: hmat
  real(8), dimension(3,3), intent(out) :: hsurf
  
  real(8), dimension(3,3) :: hinv
  integer :: kmax, maxsurf, itmp
  integer :: h, k, l, ii, i
  real(8), dimension(3) :: vec0, vec1, vec2
  real(8) :: nrm1, nrm2, rtmp1, d0
  real(8), allocatable, dimension(:,:) :: hkl
  real(8), allocatable, dimension(:) :: dist
  integer, allocatable, dimension(:) :: order

  ! Vector normal to the surface with length equal to to the spacing between the planes
  call getInverseCellMatrix(hmat,hinv,rtmp1)
  vec0(1:3) =             imiller(1) * hinv(1,:)
  vec0(1:3) = vec0(1:3) + imiller(2) * hinv(2,:)
  vec0(1:3) = vec0(1:3) + imiller(3) * hinv(3,:)
  vec0 = vec0 / sqrt(sum(vec0**2))

  kmax = 4
100 continue
  maxsurf=(kmax*2+1)**3
  allocate(dist(maxsurf))
  allocate(hkl(3,maxsurf))
  allocate(order(maxsurf))
 
  itmp=0
  do l=kmax,-kmax,-1
    do k=kmax,-kmax,-1
      do h=kmax,-kmax,-1
        vec1(1:3) = real(h,8) * hmat(1:3,1) + real(k,8) * hmat(1:3,2) + real(l,8) * hmat(1:3,3)
        nrm1=sqrt(sum(vec1*vec1))
        if (nrm1<0.1) cycle
        rtmp1 = dot_product(vec0,vec1) / nrm1
        if (abs(rtmp1)>0.001d0) cycle
        itmp=itmp+1
        dist(itmp) = nrm1
        hkl(1:3,itmp) = [h,k,l]
      enddo
    enddo
  enddo
 
  ! Extend the search in the k-space if not enough vecors are found
  if (itmp<4) then
    deallocate(dist,hkl,order)
    kmax=int(kmax*1.5)
    goto 100
  endif

  ! Rank the vecors and start from the shortest
  maxsurf=itmp
  call mrgrnk(dist(1:maxsurf) , order(1:maxsurf))
  ! call rnkpar(dist(1:maxsurf) , order(1:maxsurf) , maxsurf)

  ii=order(1)
  d0 = dist(ii)
  vec1(1:3) = hkl(1,ii) * hmat(:,1) + hkl(2,ii) * hmat(:,2) + hkl(3,ii) * hmat(:,3)
  nrm1=sqrt(sum(vec1*vec1))
  do itmp=2,maxsurf
    ii=order(itmp)
    if ( (dist(ii) - d0) < 0.01) cycle

    vec2(1:3) = hkl(1,ii) * hmat(:,1) + hkl(2,ii) * hmat(:,2) + hkl(3,ii) * hmat(:,3)
    nrm2=sqrt(sum(vec2*vec2))
    rtmp1 = dot_product(vec1,vec2) / nrm1 / nrm2

    ! Limit the space to get the more orthorhombic possible - maybe useless/dangerous
    if (rtmp1 > 0.7 .or. rtmp1 < -0.7) cycle

    hsurf(1:3,1) = vec1
    hsurf(1:3,2) = vec2
    ! hsurf(1:3,3) = vec0 ! vector normal to the surface

    ! This is the new c vector for a 3D periodic cell
    do i=1,3
      if (imiller(i)/=0) hsurf(1:3,3) = hsurf(1:3,3) + hmat(1:3,i)
    enddo
    exit
  enddo
  deallocate(dist,hkl,order)

end subroutine defineSurfaceVectors

logical function integer2logical(i)
  implicit none
  integer, intent(in) :: i
  if (i>0) then
    integer2logical = .true.
  else
    integer2logical = .false.
  end if
end function integer2logical
