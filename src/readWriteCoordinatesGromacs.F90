!disclaimer
subroutine getNumberOfAtomsGromacs(iounit,natoms,hmat)
  use moduleVariables, only: real64
  use moduleStrings
  implicit none
  integer, intent(in) :: iounit
  integer, intent(out) :: natoms
  real(real64), dimension(3,3), intent(inout) :: hmat
  real(real64) :: cell(6)
  real(real64), parameter :: pi = 3.1415926535898_real64
  
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
  if (nw == 3) then
    read(line,*) hmat(1,1),hmat(2,2),hmat(3,3)

  else if (nw == 6) then
    read(line,*) cell
    cell(4) = cell(4) * pi / 180.0_real64
    cell(5) = cell(5) * pi / 180.0_real64
    cell(6) = cell(6) * pi / 180.0_real64
    call cell2hmat(cell,hmat)
  
  else if (nw == 9) then
    read(line,*)hmat(1,1),hmat(2,2),hmat(3,3),hmat(2,1),hmat(3,1),hmat(1,2),hmat(3,2),hmat(1,3),hmat(2,3)
  else
    hmat(1,1) = 1.0_real64
    hmat(2,2) = 1.0_real64
    hmat(3,3) = 1.0_real64
  end if
  hmat = hmat*10
  
end subroutine getNumberOfAtomsGromacs

subroutine readCoordinatesGromacs(iounit,n,pos,lab,chg,hmat,go)
  use moduleStrings
  use moduleElements
  use moduleMessages
  implicit none
  integer, intent(in) :: iounit
  integer, intent(in) :: n
  real(real64), dimension(3,n), intent(out) :: pos
  real(real64), dimension(n), intent(out) :: chg
  character(len=4), dimension(n), intent(out) :: lab
  real(real64), dimension(3,3), intent(inout) :: hmat
  logical, intent(out) :: go

  integer :: ios, i, nw, natoms, itmp, id
  character(len=50) :: words(10), str
  character(len=200) :: line
  real(real64), dimension(3) :: ptmp

  go = .false.
  chg = 0.0_real64 

  read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
  read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
  read (line,*) natoms
  if (n/=natoms) call message(-1,"Wrong number of atoms in the GROMACS input file")
  do i=1,natoms
    read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
    read(line,*) str, str, itmp, ptmp
    ! read(line,'(i5,2a5,i5,3f8.3)')itmp, str, str, itmp, ptmp
    str = adjustl(str)
    id = getAtomicNumber(str)
    lab(itmp) = trim(atom(id)%lab)
    pos(1:3,itmp) = 10.0_real64*ptmp
  end do
  read(iounit,'(a200)',iostat=ios) line ; if (ios/=0) return
  call parse(line," ",words,nw)
  if (nw == 3) then
    read(line,*) hmat(1,1),hmat(2,2),hmat(3,3)
  else if (nw == 9) then
    read(line,*)hmat(1,1),hmat(2,2),hmat(3,3),hmat(2,1),hmat(3,1),hmat(1,2),hmat(3,2),hmat(1,3),hmat(2,3)
  else
    hmat(1,1) = 1.0_real64
    hmat(2,2) = 1.0_real64
    hmat(3,3) = 1.0_real64
  end if
  
  hmat = hmat * 10.0_real64
  go = .true.

end subroutine readCoordinatesGromacs
