module moduleGULP
  use moduleSymmetry
  use moduleStrings
  use moduleResizeArrays
  implicit none

  private

  public :: getNumberOfAtomsGULP, readCoordinatesGULP, writeGulpCoordinates, writeGulpCoordinatesFractional

  CHARACTER*20  :: SPG
  INTEGER*4     :: LAUE,SGINV,SGLATT,SGUNIQ,IERR_LOCAL
  REAL*4        :: RT(5,4,25)
  INTEGER*4     :: JRT(3,5,24)
  INTEGER*4     :: SGPOL

  integer, parameter :: cp = 4
  integer, parameter :: lineLength=200

  integer :: numberOfAtomsGULP

  integer :: maxNumberOfAtoms = 1000
  real(real64), dimension(3,3) :: localMetricMatrix
  logical :: fractionalCoordinates
  character(len=cp), allocatable, dimension(:) :: localLabels
  real(real64), allocatable, dimension(:,:) :: localPositions
  real(real64), allocatable, dimension(:) :: localCharges
  logical :: lerr

contains

  subroutine getNumberOfAtomsGULP(uinp,natoms,hmat)
    use moduleMessages
    use moduleVariables, only : pi
    implicit none

    integer, intent(in) :: uinp
    integer, intent(out) :: natoms
    real(real64), dimension(3,3), intent(out) :: hmat

    integer :: ios, i
    character(len=lineLength) :: line
    character(len=lineLength) :: line2
    character(len=20)  :: str
    character(len=1), external :: uppercase

    character(len=lineLength) :: key

    character(cp) :: ctmp
    real(real64), dimension(3) :: rtmp
    real(real64) :: rqq, cell(6)
    
    allocate(localPositions(3,maxNumberOfAtoms))
    allocate(localLabels     (maxNumberOfAtoms))
    allocate(localCharges    (maxNumberOfAtoms))

    numberOfAtomsGULP = 0
    natoms = numberOfAtomsGULP

    spg = "P 1"

    ! Read the first line
    call readline(uinp,line,ios)
    if (ios/=0) return

    main : do 
      call lowercase(line,line2)
      read(line2,*)key

      select case (key(1:4))
        ! read cell line
        case ("cell")
          read(line,*,iostat=ios)str,cell  
          if (ios/=0) then
            call readline(uinp,line,ios)
            if (ios/=0) call message(-1,"Error reading GULP cell 1")
          end if

          if (ios/=0) call message(-1,"Error reading GULP cell 2")
          cell(4) = cell(4) * pi / 180.0_real64
          cell(5) = cell(5) * pi / 180.0_real64
          cell(6) = cell(6) * pi / 180.0_real64
          call cell2hmat(cell,hmat)

        case ("vect") 
          do i=1,3
            call readline(uinp,line,ios)
            if (ios/=0) call message(-1,"Error reading GULP vectors")
            read(line,*,iostat=ios)hmat(:,i)
            if (ios/=0) call message(-1,"Error reading GULP vectors")
          end do

        case ("svec") 
          do i=1,2
            call readline(uinp,line,ios)
            if (ios/=0) call message(-1,"Error reading GULP vectors")
            read(line,*,iostat=ios)hmat(:,i)
            if (ios/=0) call message(-1,"Error reading GULP vectors")
          end do
          hmat(1:2,3) = 0.0_real64
          hmat(3,3) = 1000.0_real64

        ! read the coordinates
        case ("frac" , "sfra" , "cart")
          do
            call readline(uinp,line,ios)
            if (ios==-1) then
              cycle main
            else if (ios/=0) then
              call message(-1,"Error reading GULP coordinates")
            end if
            call parseGulpCoordinates(line,ctmp,rtmp,lerr,rqq)
            if (lerr) cycle main
            
            numberOfAtomsGULP = numberOfAtomsGULP + 1

            if (numberOfAtomsGULP > maxNumberOfAtoms) then
              maxNumberOfAtoms = maxNumberOfAtoms * 2
              call resizeArray(localLabels,maxNumberOfAtoms)
              call resizeArray(localPositions,maxNumberOfAtoms)
              call resizeArray(localCharges,maxNumberOfAtoms)
            end if

            localLabels(numberOfAtomsGULP) = ctmp
            localPositions(:,numberOfAtomsGULP) = rtmp
            localCharges(numberOfAtomsGULP) = rqq

            if (key(1:4) == 'frac') then
              fractionalCoordinates = .true.
            else
              fractionalCoordinates = .false.
            end if

          enddo

        case ("spac")
          do
            call readline(uinp,line,ios)
            if (ios /= 0) call message(-1,"unknown space group "//trim(spg))
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
          
      end select

      ! read next line
      call readline(uinp,line,ios)
      if (ios/=0) exit

    end do main

    localMetricMatrix = hmat

    sgnops_tot = 0
    call sgroupnp(spg,laue,sguniq,sginv,sglatt,sgnops,sgpol,jrt,cen,sgncen,rt,IERR_LOCAL)

    block
      integer :: j, k, l
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
      
    call message(0,"Parsing GULP coordinates file")
    call message(0,"...Space group name",str=spg)
    call message(0,"......Number of representative symmetry operations",i=sgnops)
    call message(0,"......Space group has inversion centre (0=no; 1=yes)",i=sginv)
    call message(1,"......Total number of symmetry operations",i=sgnops_tot)
      
    natoms = numberOfAtomsGULP * sgnops_tot * sgncen
    
  end subroutine getNumberOfAtomsGULP

            ! call readCoordinatesGULP(numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)

  subroutine readCoordinatesGULP(natoms,pos,label,chg,hmat,go)
    integer, intent(inout) :: natoms
    real(real64), dimension(3,natoms), intent(out) :: pos
    character(cp), dimension(natoms), intent(out) :: label
    real(real64), dimension(natoms), intent(out) :: chg
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

    hmat = localMetricMatrix
    call getInverseCellMatrix(hmat,hinv,volume)

    label=''
    allocate(dij(3,natoms))
  
    iatm = 0
    do idx=1,numberOfAtomsGULP
      
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
          
          if (fractionalCoordinates) then
            ptmp = ptmp - floor(ptmp)
            pos(1,iatm) = hmat(1,1)*ptmp(1) + hmat(1,2)*ptmp(2) + hmat(1,3)*ptmp(3)
            pos(2,iatm) = hmat(2,1)*ptmp(1) + hmat(2,2)*ptmp(2) + hmat(2,3)*ptmp(3)
            pos(3,iatm) = hmat(3,1)*ptmp(1) + hmat(3,2)*ptmp(2) + hmat(3,3)*ptmp(3)
          else
            pos(1,iatm) = ptmp(1)
            pos(2,iatm) = ptmp(2)
            pos(3,iatm) = ptmp(3)
          end if
          label(iatm) = localLabels(idx)
          chg(iatm) = localCharges(idx)
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

  end subroutine readCoordinatesGULP

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
  
  subroutine parseGulpCoordinates(line,lab,rtmp,lerr1,chg)
    use moduleVariables, only : cp
    implicit none

    character(len=100), intent(inout) :: line
    character(cp) :: lab
    real(real64), dimension(3), intent(out) :: rtmp
    real(real64), intent(out) :: chg
    logical, intent(out) :: lerr1

    integer :: i, iw, nw, ia, ib, ic, ix, iy, ipos
    logical :: noword
    integer, dimension(2,20) :: iword

    chg=0.0_real64

    lerr1=.false.

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
        lerr1=.true.
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
    write(uout,'(3f22.8)') frame % hmat(:,1)
    write(uout,'(3f22.8)') frame % hmat(:,2)
    write(uout,'(3f22.8)') frame % hmat(:,3)
  
    write(uout,'("cartesian")')
  
    if (sum(abs(frame % chg)) > 0) then
      do iatom=1,frame % natoms
        write(uout,'(a8,4f22.8)') adjustl(trim(frame % lab(iatom))),frame % pos(1:3,iatom), frame % chg(iatom)
      enddo
    else
      do iatom=1,frame % natoms
        write(uout,'(a8,3f22.8)') adjustl(trim(frame % lab(iatom))),frame % pos(1:3,iatom)
      enddo
    end if
  
    if (numberOfMolecules < frame % natoms) then
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
    write(uout,'(3f22.8)') frame % hmat(:,1)
    write(uout,'(3f22.8)') frame % hmat(:,2)
    write(uout,'(3f22.8)') frame % hmat(:,3)
    
    write(uout,'("fractional")')
    
    do iatom=1,frame % natoms
      write(uout,'(a8,3f22.8)') adjustl(trim(frame % lab(iatom))),frame % frac(1:3,iatom)
    enddo
  
    if (numberOfMolecules < frame % natoms) then
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

end module moduleGULP
