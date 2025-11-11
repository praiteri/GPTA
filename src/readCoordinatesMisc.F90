!disclaimer
  subroutine getNumberOfAtomsGaussian(iounit,natoms)
    integer, intent(in) :: iounit
    integer, intent(out) :: natoms
    integer :: ierr
    character(len=100) :: line

    main : do
      read(iounit,'(a100)',iostat=ierr)line
      if (ierr /= 0) exit main

      if (len_trim(line) == 0) cycle

      if (index(line,"Standard orientation") > 0) then
        do i=1,4
          read(iounit,'(a100)',iostat=ierr)line
          if (ierr /= 0) exit main
        end do
      
        natoms = 0
        do
          read(iounit,'(a100)',iostat=ierr)line
          if (ierr /= 0) exit main
          if (index(line,"--") > 0) exit main
          natoms = natoms + 1
        end do
      
      end if

    end do main

  end subroutine getNumberOfAtomsGaussian

  subroutine readCoordinatesGaussian(iounit,n,pos,lab,go)
    use moduleVariables, only: real64
    use moduleElements
    implicit none
    integer, intent(in) :: iounit
    integer, intent(in) :: n
    real(real64), dimension(3,n), intent(out) :: pos
    character(len=4), dimension(n), intent(out) :: lab
    logical, intent(out) :: go

    integer :: iatm, i, ierr, nw
    character(len=100) :: line

!                              Standard orientation:
!  ---------------------------------------------------------------------
!  Center     Atomic      Atomic             Coordinates (Angstroms)
!  Number     Number       Type             X           Y           Z
!  ---------------------------------------------------------------------
!       1          6           0        1.415358    0.307944    1.371401
!       ...
!      60          1           0       -5.696168   -0.945555   -1.360240
! ---------------------------------------------------------------------
!  Rotational constants (GHZ):           0.1348428           0.1085367           0.0968257

    go = .true.
    main : do
      read(iounit,'(a100)',iostat=ierr)line
      if (ierr /= 0) exit main

      if (len_trim(line) == 0) cycle

      if (index(line,"Standard orientation") > 0) then
        do i=1,4
          read(iounit,'(a100)',iostat=ierr)line
          if (ierr /= 0) exit main
        end do
      
        do
          read(iounit,'(a100)',iostat=ierr)line
          if (ierr /= 0) exit main
          if (index(line,"--") > 0) then
            exit main
          else
            read(line,*)iatm, i, nw, pos(1:3,iatm)
            lab(iatm) = atom(i) % lab
          end if

        end do
      
      end if

    end do main

    if (ierr /= 0) go = .false.

  end subroutine readCoordinatesGaussian

  subroutine getNumberOfAtomsXYZ(iounit,natoms,hmat)
    use moduleVariables, only: real64
    implicit none
    integer, intent(in) :: iounit
    integer, intent(out) :: natoms
    real(real64), dimension(3,3), intent(inout) :: hmat
    character(len=2000) :: line
    character(len=500) :: cellString
    integer :: ierr, i,j
    
    natoms = 0

    read(iounit,*,iostat=ierr)natoms
    if (ierr/=0) return

    read(iounit,'(a500)',iostat=ierr)line
    if (ierr/=0 .or. len_trim(line)==0) return
    
    i = index(line,"Lattice=")
    if (i > 0) then
      read(line(i+8:),*)cellString
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1:3,1), hmat(1:3,2), hmat(1:3,3)
    end if

    i = index(line,"Cell=")
    if (i > 0) then
      read(line(i+5:),*)cellString
      hmat = 0.0_real64
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1,1), hmat(2,2), hmat(3,3)
    end if

  end subroutine getNumberOfAtomsXYZ

  subroutine readCoordinatesXYZ(iounit,n,pos,lab,chg,hmat,go)
    use moduleVariables, only: real64
    use moduleStrings
    implicit none
    integer, intent(in) :: iounit
    integer, intent(in) :: n
    real(real64), dimension(3,n), intent(out) :: pos
    real(real64), dimension(n), intent(out) :: chg
    character(len=4), dimension(n), intent(out) :: lab
    real(real64), dimension(3,3), intent(inout) :: hmat
    logical, intent(out) :: go

    integer :: i, j, nn, ierr
    character(len=500) :: line

    character(len=2000) :: cellString
    character(len=500) :: fieldString
    integer :: nfields
    character(STRLEN), dimension(100) :: fields

    integer :: nf, ilab, ipos, ichg

    integer :: typeXYZ

    go = .true. 

    ilab=1
    ipos=2
    ichg=0

    read(iounit,*,iostat=ierr)nn
    if (ierr/=0 .or. nn/=n) then 
      go = .false.
      return
    end if

    read(iounit,'(a500)',iostat=ierr)line
    if (ierr/=0) then 
      go = .false.
      return
    end if
    
    i = index(line,"Lattice=")
    if (i > 0) then
      read(line(i+8:),*)cellString
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1:3,1), hmat(1:3,2), hmat(1:3,3)
    end if

    i = index(line,"Cell=")
    if (i > 0) then
      read(line(i+5:),*)cellString
      hmat = 0.0_real64
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1,1), hmat(2,2), hmat(3,3)
    end if

    i = index(line,"Properties=")
    if (i > 0) then
      read(line(i+10:),*)fieldString
      
      j = len_trim(fieldString)
      call parse(fieldString(2:j), ":", fields, nfields)
      nf = 0
      do i=1,nfields,3
        select case (fields(i))
        case ("species")
          ilab = nf + 1
        case ("pos")
          ipos = nf + 1
        case ("final_charges" , "charges")
          ichg = nf + 1
        end select
        read(fields(i+2),*) j 
        nf = nf + j
      end do
    end if

    ! looking for special cases to speed up the reading
    if (ilab==1 .and. ipos==2 .and. ichg==0)then
      typeXYZ = 1
    else if (ilab==1 .and. ipos==2 .and. ichg==5)then
      typeXYZ = 2
    else
      typeXYZ = 0
    endif

    select case (typeXYZ)
  
      ! most general (slow extended xyz)
      case default
        do i=1,nn
          read(iounit,'(a200)',iostat=ierr)line
          if (ierr/=0) then 
            go = .false.
            return
          end if
          read(line,*,iostat=ierr)fields(1:nf)
          if (ierr/=0) then
            go = .false.
            return
          end if

          read(fields(ilab),       *)lab(i)
          read(fields(ipos:ipos+2),*)pos(1:3,i)
          read(fields(ichg),       *)chg(i)
        end do
    
      ! standard XYZ no charges
      case (1)
        do i=1,n
          read(iounit,'(a200)',iostat=ierr)line
          if (ierr/=0) then
            go = .false.
            return
          end if
          read(line,*,iostat=ierr)lab(i), pos(1:3,i)
          if (ierr/=0) then
            go = .false.
            return
          end if
        enddo

      ! standard XYZ with charges
      case (2)
        do i=1,n
          read(iounit,'(a200)',iostat=ierr)line
          if (ierr/=0) then
            go = .false.
            return
          end if
          read(line,*,iostat=ierr)lab(i), pos(1:3,i), chg(i)
          if (ierr/=0) then
            go = .false.
            return
          end if
        enddo

    end select

  end subroutine readCoordinatesXYZ

  subroutine readCoordinatesXYZ_basic(iounit,n,pos,lab,chg,hmat,go)
    use moduleVariables, only: real64
    use moduleStrings
    implicit none
    integer, intent(in) :: iounit
    integer, intent(in) :: n
    real(real64), dimension(3,n), intent(out) :: pos
    real(real64), dimension(n), intent(out) :: chg
    character(len=4), dimension(n), intent(out) :: lab
    real(real64), dimension(3,3), intent(inout) :: hmat
    logical, intent(out) :: go

    integer :: i, j, nn, ierr
    character(len=500) :: line

    character(len=2000) :: cellString
    character(len=500) :: fieldString
    integer :: nfields
    character(STRLEN), dimension(100) :: fields

    integer :: nf
    integer :: typeXYZ

    go = .true. 

    read(iounit,*,iostat=ierr)nn
    if (ierr/=0 .or. nn/=n) then 
      go = .false.
      return
    end if

    read(iounit,'(a500)',iostat=ierr)line
    if (ierr/=0) then 
      go = .false.
      return
    end if
    
    i = index(line,"Lattice=")
    if (i > 0) then
      read(line(i+8:),*)cellString
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1:3,1), hmat(1:3,2), hmat(1:3,3)
    end if

    i = index(line,"Cell=")
    if (i > 0) then
      read(line(i+5:),*)cellString
      hmat = 0.0_real64
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1,1), hmat(2,2), hmat(3,3)
    end if

    do i=1,n
      read(iounit,*,iostat=ierr)lab(i), pos(1:3,i)
      if (ierr/=0) then
        go = .false.
        return
      end if
    enddo

  end subroutine readCoordinatesXYZ_basic

  subroutine getNumberOfAtomsARC(iounit,natoms,hmat)
    use moduleVariables, only: real64
    use moduleVariables, only : identityMatrix, real64
    implicit none
    integer, intent(in) :: iounit
    integer, intent(out) :: natoms
    real(real64), dimension(3,3), intent(inout) :: hmat
    character(len=500) :: line
    integer :: ierr
    
    natoms = 0

    read(iounit,'(a500)',iostat=ierr)line
    if (ierr/=0) return

    read(line,*) natoms

    read(line,*,iostat=ierr)natoms,hmat(1:3,1), hmat(1:3,2), hmat(1:3,3)
    if (ierr/=0)then
      hmat = identityMatrix
    endif

  end subroutine getNumberOfAtomsARC

  subroutine readCoordinatesARC(iounit,n,pos,lab,chg,hmat,go)
    use moduleVariables, only: real64
    use moduleStrings
    use moduleVariables, only : identityMatrix
    implicit none
    integer, intent(in) :: iounit
    integer, intent(in) :: n
    real(real64), dimension(3,n), intent(out) :: pos
    real(real64), dimension(n), intent(out) :: chg
    character(len=4), dimension(n), intent(out) :: lab
    real(real64), dimension(3,3), intent(inout) :: hmat
    logical, intent(out) :: go

    integer :: i, j, nn, ierr
    character(len=500) :: line

    go = .true. 

    read(iounit,'(a500)',iostat=ierr)line
    if (ierr/=0 .or. len_trim(line)==0) then 
      go = .false.
      return
    end if
    read(line,*,iostat=ierr) nn
    if (ierr/=0 .or. nn/=n) then 
      go = .false.
      return
    end if

    read(line,*,iostat=ierr)nn,hmat(1:3,1), hmat(1:3,2), hmat(1:3,3)
    if (ierr/=0)then
      hmat = identityMatrix
    endif

    do i=1,n
      read(iounit,'(a200)',iostat=ierr)line
      if (ierr/=0) then
        go = .false.
        return
      end if
      read(line,*,iostat=ierr)j,lab(j), pos(1:3,j)
      if (ierr/=0) then
        go = .false.
        return
      end if
    enddo
    
    chg = 0.0_real64
    
  end subroutine readCoordinatesARC

