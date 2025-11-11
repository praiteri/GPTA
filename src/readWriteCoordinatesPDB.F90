!disclaimer
subroutine getNumberOfAtomsPDB(uinp,natoms,hmat)
  use moduleVariables, only: real64
  implicit none
  integer, intent(in)  :: uinp
  integer, intent(out) :: natoms
  real(real64), dimension(3,3), intent(out) :: hmat
  real(real64), parameter :: pi = 3.1415926535898_real64
  character(len=80)  :: line
  character(len=10)  :: str
  real(real64), dimension(6) :: cell

  natoms=0
  do
    read(uinp,'(a80)',end=546)line
    read(line(1:6),*)str
    
    if(str == "CRYST1") then
      read(line,*)str,cell
      cell(4) = cell(4) * pi / 180.0_real64
      cell(5) = cell(5) * pi / 180.0_real64
      cell(6) = cell(6) * pi / 180.0_real64
      call cell2hmat(cell,hmat)
      
    else if(str=="ATOM" .or. str=="HETATM")then
      natoms=natoms+1

    else if(str=="END" .or. str=="ENDMDL")then
      exit
    else
      cycle
    end if
  enddo
546 continue

end subroutine getNumberOfAtomsPDB

subroutine readCoordinatesPDB(uinp,natoms,pos,label,charge,hmat,element,go)
  use moduleVariables, only: real64
  use moduleMessages 
  implicit none
  integer, intent(in) :: uinp
  integer, intent(in) :: natoms
  real(real64), dimension(3,natoms), intent(inout) :: pos
  character(*), dimension(natoms), intent(inout) :: label
  real(real64), dimension(natoms), intent(inout) :: charge
  real(real64), dimension(3,3), intent(inout) :: hmat
  character(*), dimension(natoms), intent(inout) :: element
  logical, intent(inout) :: go

  real(real64), parameter :: pi = 3.1415926535898_real64

  integer :: iatom
  real(real64) :: cell(6)
  character(len=100) :: line
  character(len=10)  :: str
  integer :: ios

  iatom = 0
  go=.true.
  do
    read(uinp,'(a100)',iostat=ios)line
    if (ios/=0) then
      if (is_iostat_end(ios) .and. iatom==natoms) return 
      go=.false.
      return
    end if

    ! if (len_trim(line)==0) cycle

    read(line(1:6),*)str
    if(str=="ATOM" .or. str=="HETATM")then
      iatom=iatom+1
      if (iatom>natoms) then
        call message(-1,"###### WARNING ###### too many atoms in the PDB")
        go=.false.
        return
      end if
! Labels and positions
      read(line(13:16),*)label(iatom) 
      read(line(31:54),'(3f8.3)')pos(:,iatom)
      ! read(line(77:78),*)element(iatom)

! charges - beyond character 80
      charge(iatom)=0.0_real64
      if (len_trim(line)>80) then
        read(line(81:100),*,iostat=ios) charge(iatom) 
      end if
 
    else if(str == "CRYST1") then
      read(line,*)str,cell
      cell(4) = cell(4) * pi / 180.0_real64
      cell(5) = cell(5) * pi / 180.0_real64
      cell(6) = cell(6) * pi / 180.0_real64
      call cell2hmat(cell,hmat)

    else if(str=="END" .or. str=="ENDMDL")then
      exit
    else
      cycle
    end if
  enddo

end subroutine readCoordinatesPDB

subroutine writeCoordinatesPDB(uout,lpdb2)
  use moduleVariables, only: real64
  use moduleSystem 
  use moduleElements 
  implicit none
  integer, intent(in) :: uout
  logical, intent(in) :: lpdb2

  character(len=100) :: line
  integer :: iatm
  integer :: id, jd, itmp
  character(len=4) :: fmtAtom
  character(len=5) :: chr

  write(uout,'(a18,i6)')"REMARK --- frame:  ",frame % nframe
  if(frame % volume>1.0_real64)write(uout,'(a6,3f9.3,3f7.2)')"CRYST1",frame % cell

  if (frame % natoms<=99999) then
    fmtAtom = '(i5)'
  else
    fmtAtom = '(z5)'
  end if

  do iatm=1,frame % natoms
    id=getAtomicNumber(frame % lab(iatm))
    if (numberOfMolecules>0) then
      jd=atomToMoleculeIndex(iatm)
    else
      jd=iatm
    end if
    line=" "
 
    write(line( 1: 6),'(a6   )')"ATOM  "                           ! a literal "ATOM  " (note two trailing spaces).
    write(line( 7:11),fmtAtom) iatm                                ! atom serial number, e.g. "   86". 
    write(line(12:12),'(a1   )')" "                                ! space
    write(line(13:16),'(a4   )')adjustl(frame % lab(iatm))         ! Atom role name, e.g. " CG1;". 
    write(line(17:17),'(a1   )')" "                                ! atom variant, officially called the "alternate location indicator". This is usually " " 
 
    if (numberOfMolecules>0) then
      if (jd<1 .or. jd>numberOfMolecules) then
        write(line(18:21),'(a3,2x)')"UNK"                            ! amino acid abbreviation, e.g. "ARG". 
      else
        write(line(18:22),'(a4,1x)')listOfMolecules(jd) % resname    ! amino acid abbreviation, e.g. "ARG". 
      end if
    else
      write(line(18:22),'(a3,2x)')"UNK"                            ! amino acid abbreviation, e.g. "ARG". 
    end if
 
    write(chr,fmtAtom)jd
    write(line(23:27),'(a5   )')adjustl(chr)                       ! residue sequence number (I4) and insertion code (A1), e.g. "  11 " or " 256C". 
 
    write(line(28:30),'(a3   )')"   "                              ! NON STANDARD one spaces
    write(line(31:54),'(3f8.3)')frame % pos(:,iatm)                ! atom coordinates (X, Y, Z)
    write(line(55:60),'(f6.2 )')1.0_real64                               ! atom occupancy, usually "  1.00".
    write(line(61:66),'(f6.2 )')0.0_real64                               ! B value or temperature factor.
    write(line(67:72),'(a6   )')"      "                           ! Atom's species
    write(line(73:76),'(a4   )')"    "                             ! segment identifier, left-justified. [format version 2.0 and later.]
    write(line(77:78),'(a2   )')adjustr(trim(atom(id) % lab))      ! element symbol, right-justified. [format version 2.0 and later.]
    write(line(79:80),'(a2   )')"  "                               ! charge on the atom. [format version 2.0 and later.]
    if (.not.lpdb2) write(line(81:100),'(f12.8)')frame % chg(iatm) ! charges 
    write(uout,'(a100)')line
  enddo

  if (numberOfMolecules>0) then

    block
      integer, allocatable, dimension(:) :: nBonds
      integer, allocatable, dimension(:,:) :: newList
      allocate(nBonds, source=numberOfCovalentBondsPerAtom)
      allocate(newList, source=listOfCovalentBondsPerAtom)

      do iatm=1,frame % natoms
        nBonds(iatm) = 0
        if (frame%lab(iatm) == "M" .or. frame%lab(iatm) == "MW") cycle
        do itmp=1,numberOfCovalentBondsPerAtom(iatm)
          id = listOfCovalentBondsPerAtom(itmp,iatm)
          if (frame%lab(id) == "M" .or. frame%lab(id) == "MW") cycle
          nBonds(iatm) = nBonds(iatm) + 1
          newList(nBonds(iatm),iatm) = id
        end do
      end do

      write(uout,'("ENDMDL")')
      do iatm=1,frame % natoms
        itmp = nBonds(iatm)
        if (itmp>4) then
          write(uout,'("CONECT",5i5)') iatm , newList(1:4,iatm)
          write(uout,'("CONECT",5i5)') iatm , newList(5:itmp,iatm)
        else if (itmp>0) then
          write(uout,'("CONECT",5i5)') iatm , newList(1:itmp,iatm)
        end if
        enddo
    end block

  else
    write(uout,'("END")')
  end if

end subroutine writeCoordinatesPDB
