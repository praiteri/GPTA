!disclaimer
subroutine dumpXYZ(iounit,nat,pos,lab,hmat)
  use moduleVariables, only : cp, real64
  implicit none
  integer, intent(in) :: iounit
  integer, intent(in) :: nat
  real(real64), dimension(3,nat), intent(in) :: pos
  character(cp), dimension(nat), optional, intent(in) :: lab
  real(real64), dimension(3,3), optional, intent(in) :: hmat
  integer :: i
  
  character(len=200) :: line

  write(iounit,'(i0)')nat

  if (present(hmat)) then
    write(line,'(9(f17.8,1x))') hmat
  else
    line = ""
  end if
  write(iounit,'(a)') trim(line)

  if (present(lab)) then
    do i=1,nat
      write(iounit,'(a4,1x,4(f12.5,1x))')lab(i), pos(1:3,i)
    enddo
  else
    do i=1,nat
      write(iounit,'(a4,1x,4(f12.5,1x))')"Ar", pos(1:3,i)
    enddo
  end if
  call flush(iounit)

end subroutine

subroutine dumpPDB(uout,natoms,pos,label,hmat)
  use moduleVariables, only : cp, real64
  use moduleElements, only : getAtomicNumber, atom

  implicit none
  integer, intent(in) :: uout
  integer, intent(in) :: natoms
  real(real64), dimension(3,natoms), intent(in) :: pos
  character(cp), dimension(natoms), intent(in) :: label
  real(real64), intent(in) :: hmat(3,3)

  integer :: i,nn
  character(len=80)  :: line
  real(real64) :: cell(6)
  real(real64)  :: occ, beta, tmp
  character(len=5) :: chr
  real(real64), allocatable, dimension(:,:), target :: pos_tmp
  real(real64) :: hnew(3,3), hmati(3,3), volume
  integer :: id

  allocate(pos_tmp(3,natoms), source=pos)

  call getInverseCellMatrix (hmat,hmati,volume)
  if (volume>1e-6_real64) then
    tmp=abs(hmat(2,1)) + abs(hmat(3,1)) + abs(hmat(3,2))
    if (tmp.gt.1.0e-9_real64) then
      hnew = hmat
      call makeUpperTriangularCell(hnew,pos_tmp,natoms)
      call hmat2cell(hnew,cell,"DEG")
    else
      call hmat2cell(hmat,cell,"DEG")
    endif
    if(volume>1.0_real64)write(uout,'(a6,3f9.3,3f7.2)')"CRYST1",cell
  endif

  occ=1.0_real64
  beta=0.0_real64

  do i=1,natoms
    id=getAtomicNumber(label(i))
    line=" "

    write(line( 1: 6),'(a6   )')"ATOM  "                    ! a literal "ATOM  " (note two trailing spaces).
    if (natoms<=99999) then
      write(line( 7:11),'(i5   )')i                         ! atom serial number, e.g. "   86".
    else
      write(line( 7:11),'(z5   )')i                         ! NON-STANDARD atom serial number, e.g. "   86".
    endif
    write(line(12:12),'(a1   )')" "                         ! space
    write(line(13:16),'(a4   )')adjustl(label(i))           ! Atom role name, e.g. " CG1;".
    write(line(17:17),'(a1   )')" "                         ! atom variant, officially called the "alternate location indicator". This is usually " "

    write(line(18:22),'(a4,1x)')"UNK"                       ! amino acid abbreviation, e.g. "ARG".

    if (nn<=99999) then
      write(chr,'(i5)')i
    else
      write(chr,'(z5)')i
    endif
    write(line(23:27),'(a5   )')adjustl(chr)                ! residue sequence number (I4) and insertion code (A1), e.g. "  11 " or " 256C".

    write(line(28:30),'(a3   )')"   "                       ! NON STANDARD one spaces
    write(line(31:54),'(3f8.3)')pos_tmp(:,i)                ! atom coordinates (X, Y, Z)
    write(line(55:60),'(f6.2 )')occ                         ! atom occupancy, usually "  1.00".
    write(line(61:66),'(f6.2 )')beta                        ! B value or temperature factor.
    write(line(67:72),'(a6   )')"      "                    ! Atom's species
    write(line(73:76),'(a4   )')"    "                      ! segment identifier, left-justified. [format version 2.0 and later.]
    write(line(77:78),'(a2   )')adjustr(trim(atom(id)%lab)) ! element symbol, right-justified. [format version 2.0 and later.]
    write(line(79:80),'(a2   )')"  "                        ! charge on the atom. [format version 2.0 and later.]
    write(uout,'(a80)')line
  enddo
  write(uout,'("END")')

  deallocate(pos_tmp)

end subroutine dumpPDB

subroutine dumpGULP(uout, natoms, pos, lab, hmat)
  use moduleVariables, only : cp, real64
  implicit none
  integer, intent(in) :: uout
  integer, intent(in) :: natoms
  real(real64), dimension(3,natoms), intent(in) :: pos
  character(cp), dimension(natoms) :: lab
  real(real64), dimension(3,3), intent(in) :: hmat

  integer :: iatom

  write(uout,'("vectors")')
  write(uout,'(3f22.8)') hmat(:,1)
  write(uout,'(3f22.8)') hmat(:,2)
  write(uout,'(3f22.8)') hmat(:,3)

  write(uout,'("cartesian")')

  do iatom=1,natoms
    write(uout,'(a8,3f22.8)') adjustl(trim(lab(iatom))), pos(1:3,iatom)
  enddo

end subroutine dumpGULP

subroutine dumpDCD(uout, natoms, pos, hmat, write_header)
  use moduleVariables, only: real64
  use moduleMessages 
  use moduleSystem
  implicit none
  integer, intent(in)    :: uout
  integer, intent(in) :: natoms
  real(real64), dimension(3,natoms), intent(in) :: pos
  real(real64), dimension(3,3), intent(in) :: hmat
  logical, intent(in) :: write_header

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

  call hmat2cell(hmat,cell,"COS")

  if(all(cell(1:3)>1.0_real64))then
    icell=1                 !! flag to indicate whether the cell is written in the dcd
  else
    icell=0
  end if

  ! DCD HEADER TAKEN FROM NAMD2-6
  if(write_header)then
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
    write(uout) natoms
  end if

  if(icell==1)then
    write(uout)cell(1),cell(6),cell(2),cell(5),cell(4),cell(3)
  end if

  block
    integer :: numberOfAtoms
    numberOfAtoms = natoms
    allocate(x(numberOfAtoms))
    do j=1,3
      x = real(pos(j,1:numberOfAtoms),4)
      write(uout)x(1:numberOfAtoms)
    enddo
  end block

  return
end subroutine dumpDCD

subroutine dumpCoordinates(ftype, funit, natoms, pos, label, hmat, header)
  use moduleVariables
  use moduleMessages
  implicit none 
  character(len=*) :: ftype
  integer, intent(in) :: funit
  integer, intent(in) :: natoms
  real(real64), dimension(3,natoms), intent(in) :: pos
  character(4), dimension(natoms), intent(in) :: label
  real(real64), dimension(3,3), optional, intent(in) :: hmat
  logical, optional, intent(in) :: header

  select case (ftype)
    case default
      call message(-1,"Unknown coordinates type",str=ftype)
    case ("xyz" , "XYZ")

      if (present(hmat)) then
        call dumpXYZ(funit, natoms, pos, label, hmat)
      else
        call dumpXYZ(funit, natoms, pos, label)
      end if

    case ("pdb" , "PDB")
      call dumpPDB(funit, natoms, pos, label, hmat)
      
    case ("gin")
      call dumpGULP(funit, natoms, pos, label, hmat)

    case ("dcd")
      call dumpDCD(funit, natoms, pos, hmat, header)

  end select

  call flush(funit)

end subroutine dumpCoordinates
