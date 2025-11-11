!disclaimer
subroutine debugTiming_old(idx)
  use moduleVariables
  use moduleMessages
  implicit none
  integer, intent(in) :: idx
  logical, save :: firstTimeIn = .true.
  integer, save :: numberOfSplits
  real(real64), save, allocatable, dimension(:) :: timeTally
  real(real64), save, allocatable, dimension(:) :: startClock
  logical, save, allocatable, dimension(:) :: splitOn
  
  integer :: i
  character(len=2) :: str

  if (firstTimeIn) then
    firstTimeIn = .false.
    numberOfSplits = 20
    allocate(startClock(numberOfSplits))
    allocate(splitOn(numberOfSplits), source=.false.)
    allocate(timeTally(numberOfSplits), source=-1.0_real64)
  end if

  if (idx<0) then
    if (all(timeTally < 0.0_real64)) return
    call message(0,"---> Partial Debug Timings <---")
    do i=1,numberOfSplits
      if (timeTally(i) < 0.0_real64) exit
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
    if (timeTally(idx)<0.0_real64) timeTally(idx)=0.0_real64
  end if

end subroutine debugTiming_old

subroutine debugTiming(idx)
  use moduleVariables, only : timing, real64
  use moduleMessages
  implicit none
  integer, intent(in) :: idx

  logical, save :: firstTimeIn = .true.
  integer, save :: numberOfSplits
  real(real64), save, allocatable, dimension(:) :: timeTally
  real(real64), save :: startClock, stopClock
  integer, save :: splitID = 0
  integer :: i
  character(len=2) :: str

  if (firstTimeIn) then
    firstTimeIn = .false.
    numberOfSplits = 20
    allocate(timeTally(numberOfSplits), source=-1.0e-8_real64)
  end if

  if (idx < 0) then
    do i=1,numberOfSplits
      if (timeTally(i) < 0.0_real64) exit
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
do i=ilen+1,len(lcstr)
  lcstr(i:i)=" "
end do

end subroutine lowercase

subroutine cell2hmat(cell,h)
  use moduleVariables, only : pi, pih, real64
  implicit none
  real(real64) :: cellInput(6)
  real(real64) :: cell(6)
  real(real64) :: h(3,3)

  cellInput = cell

  if (any(cell(4:6)>pi)) then
    cell(4:6) = cell(4:6)*pi/180.0_real64
  end if

  if (cell(2)<1.0e-6_real64) cell(2) = cell(1)
  if (cell(3)<1.0e-6_real64) cell(3) = cell(2)
  if (cell(4)<1.0e-6_real64) cell(4) = pih
  if (cell(5)<1.0e-6_real64) cell(5) = pih
  if (cell(6)<1.0e-6_real64) cell(6) = pih

!!p cell  :: a, b, c, alpha, beta, gamma in radians
  h = 0_real64
  h(1,1) = cell(1)
  h(1,2) = cell(2) * cos(cell(6))
  h(2,2) = cell(2) * sin(cell(6))
  h(1,3) = cell(3) * cos(cell(5))
  h(2,3) = (cell(3) * cell(2) * cos(cell(4)) - h(1,2) * h(1,3)) / h(2,2)
  h(3,3) = sqrt(cell(3)*cell(3) - h(1,3)*h(1,3) - h(2,3)*h(2,3))
!
  if ( abs(h(1,2)) < 1.0e-4_real64 ) h(1,2) = 0.0_real64
  if ( abs(h(1,3)) < 1.0e-4_real64 ) h(1,3) = 0.0_real64
  if ( abs(h(2,3)) < 1.0e-4_real64 ) h(2,3) = 0.0_real64
!
  cell = cellInput
  return
end subroutine cell2hmat

subroutine hmat2cell(hmat,cell,angle)
  use moduleVariables, only : pi, pih, real64
  implicit none
  real(real64) :: hmat(3,3)
  real(real64) :: cell(6)
  character(len=3)  :: angle

!!p a, b, c
  cell(1) = sqrt(dot_product(hmat(:,1),hmat(:,1)))
  cell(2) = sqrt(dot_product(hmat(:,2),hmat(:,2)))
  cell(3) = sqrt(dot_product(hmat(:,3),hmat(:,3)))
!!p alpha, beta, gamma in degrees
  if(angle=="DEG")then
    cell(4) = acos(dot_product(hmat(:,2),hmat(:,3))/cell(2)/cell(3))/pi*180_real64
    cell(5) = acos(dot_product(hmat(:,1),hmat(:,3))/cell(1)/cell(3))/pi*180_real64
    cell(6) = acos(dot_product(hmat(:,1),hmat(:,2))/cell(1)/cell(2))/pi*180_real64
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
  use moduleVariables, only: real64
  implicit none
  real(real64), dimension(3,3)  :: hmat, hmati
  real(real64) :: deth
  real(real64) :: odet

  deth = &
       hmat(1,1) * ( hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2) ) + &
       hmat(1,2) * ( hmat(2,3)*hmat(3,1)-hmat(2,1)*hmat(3,3) ) + &
       hmat(1,3) * ( hmat(2,1)*hmat(3,2)-hmat(2,2)*hmat(3,1) )
  if ( deth < 1.0e-4_real64 ) return
  odet = 1.0_real64 / deth
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
  use moduleVariables, only: real64
  implicit none
  integer :: nn
  real(real64), intent(inout) :: hmat(3,3)
  real(real64), dimension(3,nn), intent(inout) :: pos
  real(real64) :: hnew(3,3), rotchi(3,3)
  real(real64) :: u(3), v(3), rotmat(3,3)
  integer :: i

  u = hmat(:,1)
  v = hmat(:,2)
  call rotateCell(u,v,rotmat)
  do i=1,3
    hnew(:,i)=matmul(rotmat, hmat(:,i))
  enddo

  ! One more rotation if the b vector points towards -y
  if ( dot_product(hnew(:,2) , [0.0_real64,1.0_real64,0.0_real64]) < 0.0_real64 ) then
    rotchi(1,1) = 1.0_real64
    rotchi(1,2) = 0.0_real64
    rotchi(1,3) = 0.0_real64
    rotchi(2,1) = 0.0_real64
    rotchi(2,2) = -1.0_real64
    rotchi(2,3) = 0.0_real64
    rotchi(3,1) = 0.0_real64
    rotchi(3,2) = 0.0_real64
    rotchi(3,3) = -1.0_real64
  
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
  use moduleVariables, only: real64
  implicit none
  real(real64)  :: u(3), v(3)
  real(real64)  :: normu, u1(3)
  real(real64)  :: cospsi, costheta, coschi, sinchi, sintheta, sinpsi
  real(real64)  :: rottheta(3,3), rotpsi(3,3), rot1(3,3), rotchi(3,3)
  real(real64)  :: vnew(3), rot(3,3)
  integer :: i

  rot=0.0_real64
  do i=1,3
    rot(i,i) = 1.0_real64
  enddo

! Align a with x
  if(abs(u(2))>1.0e-6_real64 .or. abs(u(3))>1.0e-6_real64)then
    normu = sqrt(dot_product(u,u))

    costheta = u(2)/sqrt(u(2)**2+u(3)**2)
    sintheta = u(3)/sqrt(u(2)**2+u(3)**2)

    rottheta(1,1) = 1.0_real64
    rottheta(1,2) = 0.0_real64
    rottheta(1,3) = 0.0_real64
    rottheta(2,1) = 0.0_real64
    rottheta(2,2) = costheta
    rottheta(2,3) = sintheta
    rottheta(3,1) = 0.0_real64
    rottheta(3,2) = -sintheta
    rottheta(3,3) = costheta

    u1 = MATMUL(rottheta, u )

    cospsi = u1(1)/normu
    sinpsi = u1(2)/normu

    rotpsi(1,1) = cospsi
    rotpsi(1,2) = sinpsi
    rotpsi(1,3) = 0.0_real64
    rotpsi(2,1) = -sinpsi
    rotpsi(2,2) = cospsi
    rotpsi(2,3) = 0.0_real64
    rotpsi(3,1) = 0.0_real64
    rotpsi(3,2) = 0.0_real64
    rotpsi(3,3) = 1.0_real64

    rot1 = matmul(rotpsi,rottheta)
  else
    rot1=0.0_real64
    do i=1,3
      rot1(i,i)=1.0_real64
    enddo
  end if
  vnew = matmul(rot1, v)

!!p rotate around x in order to place b on the xy plane

  if( abs(vnew(3)) > 1.0e-6_real64)then
    coschi = vnew(2)/sqrt(vnew(2)**2+vnew(3)**2)
    sinchi = vnew(3)/sqrt(vnew(2)**2+vnew(3)**2)
  else 
    coschi = 1.0_real64
    sinchi = 0.0_real64
  end if

  rotchi(1,1) = 1.0_real64
  rotchi(1,2) = 0.0_real64
  rotchi(1,3) = 0.0_real64
  rotchi(2,1) = 0.0_real64
  rotchi(2,2) = coschi
  rotchi(2,3) = sinchi
  rotchi(3,1) = 0.0_real64
  rotchi(3,2) = -sinchi
  rotchi(3,3) = coschi
  
  rot  = matmul(rotchi,rot1)

  return
end subroutine rotateCell

subroutine straightenCell(hmat,dir,hnew)
  use moduleVariables, only : cp, real64
  use moduleMessages
  implicit none
  real(real64), intent(in) :: hmat(3,3)
  character(len=*) :: dir
  real(real64), intent(out), target :: hnew(3,3)

  integer :: nn, ia, ib
  real(real64), dimension(3) :: sij, dij
  real(real64), pointer, dimension(:) :: vec1, vec2, vec3
  real(real64) :: len0, len1

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
  use moduleVariables, only : cp, real64
  implicit none
  real(real64), dimension (3), intent(in)  :: a
  real(real64), dimension (3), intent(in)  :: b
  real(real64), dimension (3) :: c

  c(1) =  a(2)*b(3) - a(3)*b(2)
  c(2) = -a(1)*b(3) + a(3)*b(1)
  c(3) =  a(1)*b(2) - a(2)*b(1)

  return
end function cross_product

subroutine readCellFreeFormat(string, hmat)
  use moduleVariables, only: real64
  use moduleStrings
  use moduleMessages
  implicit none
  character(len=*), intent(in) :: string
  real(real64), dimension(3,3), intent(out) :: hmat

  integer :: numberOfWords
  character(len=STRLEN), dimension(100) :: listOfWords
  real(real64), dimension(6) :: cell

  hmat = 0.0_real64
  call parse(string,",",listOfWords,numberOfWords)
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
  
subroutine rotateMolecule(vec,theta,natoms,pos,pivot)
  use moduleVariables, only: real64
  implicit none
  real(real64), intent(in) :: vec(3), theta
  integer, intent(in) :: natoms
  real(real64), dimension(3,natoms), intent(inout) :: pos
  integer, intent(in) :: pivot

  integer :: i
  real(real64) :: v(3), v2(3), costheta, sintheta, u(3), rcom(3)
  real(real64), dimension (3,3) :: rotmat
  real(real64) :: norm

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
    rcom=0.0_real64
    do i=1,natoms
      rcom=rcom+pos(1:3,i)
    enddo
    rcom=rcom/real(natoms,8)

  ! rotates around the origin
  else
    rcom = 0.0_real64
  endif

  do i=1,natoms
    u=pos(1:3,i)-rcom(1:3)
    pos(1:3,i)=matmul(rotmat,u(1:3))+rcom(1:3)
  enddo

end subroutine rotateMolecule

subroutine linearRegression(n,x,y,alpha,beta,rfact)
  use moduleVariables, only: real64
  ! y = alpha + beta*x
  implicit none
  integer, intent(in) :: n
  real(real64), dimension(n), intent(in) :: x
  real(real64), dimension(n), intent(in) :: y
  real(real64), intent(out) :: alpha, beta, rfact

  real(real64) :: xavg, yavg, xyavg, x2avg, y2avg
  real(real64) :: ssxx, ssyy, ssxy
  integer :: i

  xavg  = 0.0_real64
  yavg  = 0.0_real64
  xyavg = 0.0_real64
  x2avg = 0.0_real64
  y2avg = 0.0_real64
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
  use moduleVariables, only : actionTypeDef, real64
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

  if (.not. allocated(a % idxSelection)) then
    allocate(a % idxSelection(ntmp,n))
    allocate(a % idxToSelection(frame%natoms,n))
  end if
  
  do idx=1,n
    ntmp = 0
    do i=1,frame % natoms
      if (.not. a % isSelected(i, idx)) cycle
      ntmp = ntmp + 1
      a % idxSelection(ntmp, idx) = i
      a % idxToSelection(i, idx) = ntmp
    end do
  end do
end subroutine createSelectionList

subroutine createInvertedSelectionList(a,n)
  use moduleVariables, only : actionTypeDef, real64
  use moduleSystem, only : frame 
  use moduleMessages
  implicit none
  type(actionTypeDef), target :: a
  integer, intent(in) :: n
  integer :: ntmp, i, idx
  if (.not. allocated(a % isSelected)) call message(-1,"Cannot creat selection list",str=a % actionDetails)

  ntmp = frame % natoms - count(a % isSelected(:,1))
  do idx=2,n
    ntmp = max(ntmp , frame % natoms - count(a % isSelected(:,idx)))
  end do

  if (.not. allocated(a % idxSelection)) allocate(a % idxSelection(ntmp,n))
  
  do idx=1,n
    ntmp = 0
    do i=1,frame % natoms
      if (a % isSelected(i, idx)) cycle
      ntmp = ntmp + 1
      a % idxSelection(ntmp, idx) = i
    end do
  end do
end subroutine createInvertedSelectionList

subroutine defineSurfaceVectors(imiller, hmat, hsurf)
  use moduleVariables, only: real64
  ! use m_rnkpar
  use m_mrgrnk
  implicit none
  integer, dimension(3), intent(in) :: imiller
  real(real64), dimension(3,3), intent(in) :: hmat
  real(real64), dimension(3,3), intent(out) :: hsurf
  
  real(real64), dimension(3,3) :: hinv
  integer :: kmax, maxsurf, itmp
  integer :: h, k, l, ii, i
  real(real64), dimension(3) :: vec0, vec1, vec2
  real(real64) :: nrm1, nrm2, rtmp1, dist0
  real(real64), allocatable, dimension(:,:) :: hkl
  real(real64), allocatable, dimension(:) :: dist
  integer, allocatable, dimension(:) :: order

  ! Vector normal to the surface with length equal to to the spacing between the planes
  call getInverseCellMatrix(hmat,hinv,rtmp1)
  vec0(1:3) =             imiller(1) * hinv(1,:)
  vec0(1:3) = vec0(1:3) + imiller(2) * hinv(2,:)
  vec0(1:3) = vec0(1:3) + imiller(3) * hinv(3,:)
  vec0 = vec0 / sqrt(sum(vec0**2))
  kmax = 10
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
        if (nrm1<0.01) cycle
        rtmp1 = dot_product(vec0,vec1) / nrm1
        if (abs(rtmp1)>0.001_real64) cycle
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

  ! do itmp=1,maxsurf
  !   ii=order(itmp)
  !   write(0,*)itmp,dist(ii),hkl(:,ii)
  ! end do

  ii=order(1)
  dist0 = dist(ii)
  vec1(1:3) = hkl(1,ii) * hmat(:,1) + hkl(2,ii) * hmat(:,2) + hkl(3,ii) * hmat(:,3)
  nrm1=sqrt(sum(vec1*vec1))
  ! do itmp=maxsurf,2,-1
  do itmp=2,maxsurf
    ii=order(itmp)
    if ( (dist(ii) - dist0) < 0.01) cycle
    
    vec2(1:3) = hkl(1,ii) * hmat(:,1) + hkl(2,ii) * hmat(:,2) + hkl(3,ii) * hmat(:,3)
    nrm2=sqrt(sum(vec2*vec2))
    rtmp1 = dot_product(vec1,vec2) / nrm1 / nrm2
    
    ! Limit the space to get the more orthorhombic possible - maybe useless/dangerous
    if (rtmp1 > 0.2 .or. rtmp1 < -0.2) cycle

    hsurf(1:3,1) = vec1
    hsurf(1:3,2) = vec2
    ! hsurf(1:3,3) = vec0 ! vector normal to the surface

    ! This is the new c vector for a 3D periodic cell
    do i=1,3
      if (imiller(i)/=0) hsurf(1:3,3) = hsurf(1:3,3) + hmat(1:3,i) * sign(1,imiller(i))
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

subroutine compute_histogram(data, n_data, min_val, max_val, n_bins, bins, hist)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    
    ! Arguments
    real(real64), intent(in) :: data(n_data)        ! Input data array
    integer, intent(in) :: n_data           ! Number of data points
    real(real64), intent(in) :: min_val             ! Minimum value for binning
    real(real64), intent(in) :: max_val             ! Maximum value for binning
    integer, intent(in) :: n_bins           ! Number of bins
    real(real64), intent(OUT) :: bins(n_bins)       ! Bin edges
    real(real64), intent(OUT) :: hist(n_bins)    ! Histogram counts
    
    ! Local variables
    integer :: i, bin_index
    real :: bin_width
    
    ! Initialize histogram counts to zero
    hist = 0
    
    ! Calculate bin width
    bin_width = (max_val - min_val) / n_bins
    
    ! Calculate bin edges
    do i = 1, n_bins
        bins(i) = min_val + (i-1) * bin_width
    end do
    
    ! Compute histogram
    do i = 1, n_data
        if (data(i) >= min_val .AND. data(i) < max_val) then
            bin_index = INT((data(i) - min_val) / bin_width) + 1
            if (bin_index <= n_bins) then
                hist(bin_index) = hist(bin_index) + 1
            end if
        end if
    end do
end subroutine compute_histogram

subroutine get_time_seed(seed)
    implicit none
    integer, intent(out) :: seed
    
    ! Local variables
    integer :: time_values(8)
    integer :: milliseconds
    
    ! Get current date and time
    call date_and_time(values=time_values)
    
    ! time_values contains:
    ! (1) year, (2) month, (3) day, (4) time difference with UTC in minutes
    ! (5) hour, (6) minute, (7) second, (8) milliseconds
    
    ! Create seed from time components
    ! Use seconds, minutes, hours, and milliseconds for maximum variability
    seed = time_values(8) + &           ! milliseconds (0-999)
           time_values(7) * 1000 + &    ! seconds 
           time_values(6) * 6000        ! minutes 
    
    ! Ensure seed is positive (some RNGs require positive seeds)
    seed = abs(seed)
    
    ! Avoid zero seed (some RNGs fail with zero)
    if (seed == 0) seed = 1
    
end subroutine get_time_seed

subroutine get_filename(path, filename)
    implicit none
    character(len=*), intent(in) :: path
    character(len=256), intent(out) :: filename
    integer :: i, last_slash

    ! Find the last '/' in the path
    last_slash = 0
    do i = 1, len_trim(path)
        if (path(i:i) == '/') then
            last_slash = i
        end if
    end do

    ! Extract filename (everything after the last '/')
    if (last_slash > 0) then
        filename = path(last_slash+1:len_trim(path))
    else
        filename = path  ! No '/' found, entire string is filename
    end if

end subroutine get_filename