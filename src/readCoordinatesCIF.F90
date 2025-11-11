module moduleCIF
  use moduleSymmetry
  implicit none

  private

  public :: getNumberOfAtomsCIF, readCoordinatesCIF

  CHARACTER*20  :: SPG
  INTEGER*4     :: LAUE,SGINV,SGLATT,SGUNIQ,IERR_LOCAL
  REAL*4        :: RT(5,4,25)
  INTEGER*4     :: JRT(3,5,24)
  INTEGER*4     :: SGPOL

  integer, parameter :: cp = 4
  integer, parameter :: lineLength=200
  character(len=20) :: fmtLine

  integer :: numberOfAtomsCIF

  character(len=20) :: spaceGroup_HM = ""
  character(len=20) :: spaceGroup_HM_alt = ""
  character(len=20) :: spaceGroup_Hall = ""
  integer :: spaceGroupNumber = 0

  integer, parameter :: maxNumberOfAtoms = 1000
  real(real64), dimension(6) :: cell
  real(real64), dimension(3,maxNumberOfAtoms) :: localPositions
  character(len=cp), dimension(maxNumberOfAtoms) :: localLabels
  integer :: ios

contains

  subroutine readLine(uinp,ios,line)
    integer, intent(in) :: uinp
    integer, intent(out) :: ios
    character(len=lineLength), intent(out) :: line
    read(uinp,fmtLine,iostat=ios) line
    line = adjustl(line)
  end subroutine readLine

  subroutine getNumberOfAtomsCIF(uinp,natoms,hmat)
    use moduleMessages
    implicit none

    integer, intent(in) :: uinp
    integer, intent(out) :: natoms
    real(real64), dimension(3,3), intent(out) :: hmat

    character(len=lineLength) :: line
    character(len=lineLength) :: line2

    integer :: nw
    character(len=lineLength) :: words(100)
    character(len=lineLength) :: field(100)

    integer :: iField
    natoms = 0

    write(fmtLine,('("(a",i0,")")')) lineLength

    ! read a new line
    call readLine(uinp,ios,line)

    mainLoop: do
      ! empty line
      ! nothing to do read new line
      if (len_trim(line) == 0) then
        call readLine(uinp,ios,line)
        if (ios/=0) exit
        cycle mainLoop
      end if

      ! comment
      ! nothing to do read new line
      if (line(1:1) == "#") then
        call readLine(uinp,ios,line)
        if (ios/=0) exit
        cycle mainLoop
      end if

      ! loop_ (more than one item in the next lines)
      if (line(1:5) == "loop_") then
        nw = 0
        do
          call readLine(uinp,ios,line)
          line = adjustl(line)
          if (line(1:1) /= "_") exit

          nw = nw + 1
          read(line,*) field(nw)
        end do
        do
          read(line,*,iostat=ios)words(1:nw)
          if (ios/=0) then
            call readLine(uinp,ios,line2)
            line = trim(line)//" "//trim(line2)
            read(line,*,iostat=ios)words(1:nw)
          end if

          do iField=1,nw
            call extractField(field(iField),words(iField))
          end do
          call readLine(uinp,ios,line)
          if (ios/=0) exit mainLoop
          if (len_trim(line) == 0) exit
          if (line(1:1) == "_") cycle mainLoop
          if (line(1:5) == "loop_") cycle mainLoop
        end do
      else if (line(1:1) == "_") then
        nw = 1
        read(line,*,iostat=ios)field(nw), words(nw)
        if (ios == 0) then
          call extractField(field(nw),words(nw))
        else
          call readFieldFromMultipleLines(uinp)
        end if
      
      else

        ! nothing to do read new line
        call readLine(uinp,ios,line)
        if (ios/=0) exit
        cycle mainLoop
      
      end if

      ! read a new line
      call readLine(uinp,ios,line)
      if (ios/=0) exit

    end do mainLoop
   
    if (spaceGroupNumber > 0) then
      spg = get_spg_name(spaceGroupNumber)
      
    else if (len_trim(spaceGroup_HM) > 0) then
      spg = trim(spaceGroup_HM)
      
    else if (len_trim(spaceGroup_HM_alt) > 0) then
      spg = trim(spaceGroup_HM_alt)
      
    else if (len_trim(spaceGroup_Hall) > 0) then
      spg = trim(spaceGroup_Hall)
      
    else
      spg = "P 1"
      
    end if
    call message(0,"Parsing CIF coordinates file")
    call message(0,"......Space group name",str=spg)
    call message(0,"......Space group number",i=spaceGroupNumber)
    call message(0,"......Space group name (H-M)",str=spaceGroup_HM)
    call message(0,"......Space group name (H-M alt)",str=spaceGroup_HM_alt)
    call message(0,"......Space group name (Hall)",str=spaceGroup_Hall)

    sgnops_tot = 0
    call sgroupnp(spg,laue,sguniq,sginv,sglatt,sgnops,sgpol,jrt,cen,sgncen,rt,IERR_LOCAL)

    call message(0,"......Number of representative symmetry operations",i=sgnops)
    call message(0,"......Space group has inversion centre (0=no; 1=yes)",i=sginv)

    block
      integer :: i, j, k, l
      do l=0,sginv
        do k=1,sgnops
          sgnops_tot = sgnops_tot + 1
          do j=1,3
            do i=1,3
              rot(i,j,sgnops_tot) = real(jrt(i,j,k),8) * (1-2*l)
            end do
            rot(j,4,sgnops_tot) = real(jrt(j,4,k),8)/12.0_real64
          end do
        end do
      enddo
    end block

    natoms = numberOfAtomsCIF * sgnops_tot * sgncen
      
    call message(0,"......Total number of symmetry operations",i=sgnops_tot)
    call message(1,"......Number of atoms in the asymmetric unit",i=numberOfAtomsCIF)
      
    call cell2hmat(cell,hmat)

  end subroutine getNumberOfAtomsCIF

  subroutine readCoordinatesCIF(natoms,pos,label,hmat,go)
    integer, intent(inout) :: natoms
    real(real64), dimension(3,natoms), intent(out)  :: pos
    character(cp), dimension(natoms), intent(out)  :: label
    real(real64), dimension(3,3), intent(out) :: hmat
    logical, intent(out) :: go

    integer :: i, j, k, l, idx
    integer :: iatm
    logical, allocatable, dimension(:) :: lremove
    real(real64) :: hinv(3,3), volume, rdist
    real(real64), dimension(3) :: sij, ptmp, rtmp
    real(real64), allocatable, dimension(:,:) :: dij

    logical, save :: firstTimeIn = .true.

    if (firstTimeIn) then
      firstTimeIn = .false.
      go = .true.
    else
      go = .false.
      return
    end if

    call cell2hmat(cell,hmat)
    call getInverseCellMatrix(hmat,hinv,volume)

    label=''
    allocate(dij(3,natoms))
  
    iatm = 0
    do idx=1,numberOfAtomsCIF
      
      ! sij contains the factional coordinates
      sij(1:3) = localPositions(:,idx)

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
          label(iatm) = localLabels(idx)
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
        if ( rdist < 1.0e-3_real64 ) lremove(j) = .true.
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
  
  end subroutine readCoordinatesCIF

  subroutine extractField(keyword,stringValue) 
    implicit none
    character(len=*), intent(in) :: keyword
    character(len=*), intent(in) :: stringValue

    if (index(keyword,"space_group_name_H-M "      ) > 0) read (stringValue,'(a20)') spaceGroup_HM
    if (index(keyword,"space_group_name_H-M_alt"   ) > 0) read (stringValue,'(a20)') spaceGroup_HM_alt
    if (index(keyword,"space_group_name_Hall"      ) > 0) read (stringValue,'(a20)') spaceGroup_Hall
    if (index(keyword,"space_group_IT_number"      ) > 0) read (stringValue,*) spaceGroupNumber
    if (index(keyword,"_symmetry_Int_Tables_number") > 0) read (stringValue,*) spaceGroupNumber

    select case(keyword)
      case default
        return

      case ("_cell_length_a")
        cell(1) = readFloat(stringValue)
      case ("_cell_length_b")
        cell(2) = readFloat(stringValue)
      case ("_cell_length_c")
        cell(3) = readFloat(stringValue)

      case ("_cell_angle_alpha")
        cell(4) = readFloat(stringValue)
      case ("_cell_angle_beta")
        cell(5) = readFloat(stringValue)
      case ("_cell_angle_gamma")
        cell(6) = readFloat(stringValue)

      case("_atom_site_label")
        numberOfAtomsCIF = numberOfAtomsCIF + 1
        read(stringValue,*) localLabels(numberOfAtomsCIF)
        
      case("_atom_site_fract_x")
        localPositions(1,numberOfAtomsCIF) = readFloat(stringValue)
      case("_atom_site_fract_y")
        localPositions(2,numberOfAtomsCIF) = readFloat(stringValue)
      case("_atom_site_fract_z")
        localPositions(3,numberOfAtomsCIF) = readFloat(stringValue)

    end select

  end subroutine extractField

  function readFloat(str) result(val)
    implicit none
    character(len=*) :: str
    real(real64) :: val
    integer :: i

    do i=1,len_trim(str)
      if (str(i:i) == "(") exit
    end do
    read(str(1:i-1),*) val

  end function readFloat

  subroutine readFieldFromMultipleLines(io)
    implicit none
    integer, intent(in) :: io
    character(len=100) :: localLine

    read(io,fmtLine,iostat=ios) localLine
    if (localLine(1:1) == ";") then
      do
        read(io,fmtLine,iostat=ios) localLine
        if (localLine(1:1) == ";") return
      end do
    else
      return
    end if

  end subroutine readFieldFromMultipleLines

  function distance(dij,hmat) result(rdist)
    implicit none
    real(real64), intent(in) :: dij(3), hmat(3,3)
    real(real64) :: rdist
    real(real64) :: sij(3), rtmp
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
  
end module moduleCIF
