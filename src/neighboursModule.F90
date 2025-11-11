!disclaimer
module moduleNeighbours 
  use moduleVariables
  use moduleSystem 
  use moduleMessages 
  use moduleActions
  use moduleDistances
  implicit none

  ! Linked cells
  integer :: iextra = 1
  integer :: nnear = 14 ! 63 ! ((iextra*2 + 1)**2 * (iextra+1) - 4*(2*(iextra-1)+1))
  integer :: ncell(3)
  real(real64) :: rcell(3,3)

  ! Timing
  real(real64) :: neigh_start, neigh_end
  real(real64) :: neighboursListTime

  ! type(allActionsProcedure) :: allActions(100)

  interface
    subroutine neighboursList() 
      implicit none
    end subroutine neighboursList
  end interface
  procedure (neighboursList), pointer :: computeNeighbours => null()

contains

  subroutine initialiseNeighboursList()
    implicit none
    integer :: iact

#ifdef DEBUG
    write(0,*) " --> Entering initialiseNeighboursList <-- "
#endif

    do iact=1,numberOfActions
      if (action(iact) % requiresNeighboursList) then
        computeNeighboursList = (computeNeighboursList .or. action(iact) % requiresNeighboursList)
        if (action(iact) % requiresNeighboursListUpdates) computeNeighboursListOnce = .false.

        computeNeighboursListDouble = (computeNeighboursListDouble .or. action(iact) % requiresNeighboursListDouble)
        cutoffNeighboursListGlobal = max(action(iact) % cutoffNeighboursList, cutoffNeighboursListGlobal)
      end if
    enddo

    ! estimated max number of neighbours as about 4 times the number density of the water
    neighmax = int(0.6_real64 * cutoffNeighboursListGlobal**3)

    neighboursListUpdates = 0
    neighboursListTime = 0.0_real64

  end subroutine initialiseNeighboursList

  subroutine setUpNeigboursList()
    implicit none
    integer :: ntmp, idx
    character(len=20) :: algorithm

#ifdef DEBUG
    write(0,*) " --> Entering setUpNeigboursList <-- "
#endif

    ! Check if I can use the linked cell
    computeNeighboursListVerlet = .true.
    if (.not. forceVerletList) then
      if (pbc_type == "tri" .or. pbc_type == "ortho") then
        do idx=1,3
          ncell(idx) = iextra*floor(frame % hmat(idx,idx)/cutoffNeighboursListGlobal)
          if (ncell(idx)>=3) computeNeighboursListVerlet = .false.
          rcell(:,idx) = frame % hmat(:,idx) / ncell(idx)
        enddo

        ! the linked cells method is faster with only one extra cell
        ! hoe however, with two extra cells is still faster than the verlet
        if (computeNeighboursListVerlet) then
          iextra = 2 
          nnear = 63
          do idx=1,3
            ncell(idx) = iextra*floor(frame % hmat(idx,idx)/cutoffNeighboursListGlobal)
            if (ncell(idx)>=3) computeNeighboursListVerlet = .false.
            rcell(:,idx) = frame % hmat(:,idx) / ncell(idx)
          enddo
        end if

      end if
    end if  

    if (computeNeighboursListVerlet) then
      algorithm = "Verlet List"
      computeNeighbours => computeVerletList
    else
      if (cutoffNeighboursListGlobal > 5.5_real64) then
        algorithm = "Linked Cell Sorting"
        computeNeighbours => computeLikedCellNeighboursSort
      else
        algorithm = "Linked Cell"
        computeNeighbours => computeLikedCellNeighbours
      end if
    end if

    ntmp = frame % natoms
! not sure about this    if (frame % natoms == 0) then
! not sure about this      ntmp = numberOfAtoms
! not sure about this    else
! not sure about this      ntmp = frame % natoms
! not sure about this    end if

    if (allocated(nneigh)) then 
      if ( size(lneigh)/ntmp /= neighmax) then
        deallocate(nneigh)
        deallocate(lneigh)
        deallocate(rneigh)
      else
        return
      end if  
    end if

    allocate(nneigh(          ntmp))
    allocate(lneigh(neighmax, ntmp))
    allocate(rneigh(neighmax, ntmp))

    call message(0,"Neighbours' list setup")
    call message(0,"...Algorithm",str=algorithm)
    if (.not. computeNeighboursListVerlet) then
      call message(0,"......Number of cells",iv=ncell)
    end if
    call message(0,"...Double list",l=computeNeighboursListDouble)
    call message(0,"...Cutoff radius",r=cutoffNeighboursListGlobal)
    call message(1,"...Maximum number of neighbours",i=neighmax)

  end subroutine setUpNeigboursList

  subroutine updateNeighboursList(lflag)
    use moduleResizeArrays 
    implicit none
    logical, intent(in) :: lflag
    integer :: iatm, jatm, ineigh
    integer, allocatable, dimension(:) :: neighTmp

#ifdef DEBUG
    write(0,*) " --> Entering updateNeighboursList <-- "
#endif

    if (.not. lflag) then
      if (updatedListFrameNumber == frame % nframe) return
    end if

    nneigh = 0
    ! Resize arrays if some atoms have been removed
    if (frame % natoms /= size(nneigh)) then
      call resizeArray(nneigh,frame % natoms)
      call resizeArray(lneigh,frame % natoms)
      call resizeArray(rneigh,frame % natoms)
    end if

    neigh_start = timing()
    call computeNeighbours()

    if (computeNeighboursListDouble) then
      ineigh = size(nneigh)
      allocate(neighTmp(ineigh) , source=nneigh)
      do iatm=1,frame % natoms
        do ineigh=1,nneigh(iatm)
          jatm = lneigh(ineigh,iatm)
          neighTmp(jatm) = neighTmp(jatm) + 1
          if (neighTmp(jatm) > neighmax) call extend_nlist(20)
          lneigh(neighTmp(jatm),jatm) = iatm
          rneigh(neighTmp(jatm),jatm) = rneigh(ineigh,iatm)
        enddo
      enddo
      call move_alloc(neighTmp, nneigh)
    end if

    updatedListFrameNumber = frame % nframe

    neigh_end = timing()
    neighboursListTime = neighboursListTime + (neigh_end-neigh_start)
    neighboursListUpdates = neighboursListUpdates + 1

  end subroutine updateNeighboursList

  subroutine computeVerletList()
    implicit none

    integer :: natm
    integer :: iatm, jatm
    real(real64) :: rcut2
    real(real64) :: dist, dij(3)
    integer :: ierr = 0

    natm = frame % natoms
    rcut2 = cutoffNeighboursListGlobal**2

    ierr = 1
    do while(ierr /= 0)
      ierr = 0
      nneigh = 0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (iatm,jatm,dij,dist)
!$OMP DO SCHEDULE(static) REDUCTION(+:ierr)
      do iatm=1,natm
        do jatm=iatm+1,natm
          dij(1:3) = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
          dist = computeDistanceSquaredPBC(dij)
          if (dist < rcut2) then
            if (nneigh(iatm) == neighmax) then
              ierr = ierr + 1
              exit
            end if
          
            nneigh(iatm) = nneigh(iatm) + 1
            lneigh(nneigh(iatm),iatm) = jatm
            rneigh(nneigh(iatm),iatm) = sqrt(dist)
          end if
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
      if (ierr>0) call extend_nlist(20)
    enddo

  end subroutine computeVerletList

  subroutine computeLikedCellNeighbours()
    integer :: ierr

    integer :: iatm, jatm, natoms
    real(real64), allocatable, dimension(:,:) :: cartesianCoord

    real(real64) :: rcut2, rij
    integer :: i, j, k, l 
    integer :: ntot, ntot2, nn
    integer, dimension(3) :: iij, iv1, iv2, ncell2
    integer, allocatable, dimension(:) :: nlink
    integer, allocatable, dimension(:,:) :: ilink
    real(real64), allocatable, dimension(:,:,:) :: rlink

    integer :: i1, i2, i3, idx, ii, jj
    real(real64), dimension(3,3) :: hmat
    real(real64), dimension(3) :: p0, dij

    natoms = frame % natoms
    
    allocate(cartesianCoord(3,natoms))
    call fractionalToCartesian(natoms, frame % frac, cartesianCoord)

    rcut2  = cutoffNeighboursListGlobal**2
    ncell2 = ncell+iextra
    ntot   = product(ncell)
    ntot2  = product(ncell2)
    nn     = int(natoms/ntot)+10

    allocate(nlink(ntot2))
    allocate(ilink(nn,ntot2))
    allocate(rlink(3,nn,ntot2))

    nlink = 0
    do i=1,natoms
      iij(1:3) = int(frame % frac(1:3,i)*ncell)+1
      j = vec2int2(iij,ncell2)
      nlink(j) = nlink(j) + 1
      
      ! Reallocate array if necessary
      if (nlink(j) > nn) then
        call extend_ivec(ilink,nn,ntot2,nn+20,ntot2)
        call extend_rmat(rlink,3,nn,ntot2,nn+20,ntot2)
        nn=nn+20
      end if

      ilink(nlink(j),j) = i
      rlink(1:3,nlink(j),j) = cartesianCoord(1:3,i)
    enddo

    hmat = frame % hmat
    do i1=1,ncell(1)
      do i2=1,ncell(2)
        do i3=1,ncell(3)

          if (i1<=iextra .or. i2<=iextra .or. i3<=iextra) then
            j = vec2int2((/i1,i2,i3/),ncell2)

            if (i1<=iextra .and. i2<=iextra .and. i3<=iextra) then
              k = vec2int2((/i1+ncell(1),i2+ncell(2),i3+ncell(3)/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,1)+hmat(:,2)+hmat(:,3)
              enddo
            end if

            if (i1<=iextra .and. i2<=iextra) then
              k = vec2int2((/i1+ncell(1),i2+ncell(2),i3/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,1)+hmat(:,2)
              enddo
            end if

            if (i1<=iextra .and. i3<=iextra) then
              k = vec2int2((/i1+ncell(1),i2,i3+ncell(3)/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,1)+hmat(:,3)
              enddo
            end if

            if (i2<=iextra .and. i3<=iextra) then
              k = vec2int2((/i1,i2+ncell(2),i3+ncell(3)/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,2)+hmat(:,3)
              enddo
            end if

            if (i1<=iextra) then
              k = vec2int2((/i1+ncell(1),i2,i3/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,1)
              enddo
            end if

            if (i2<=iextra) then
              k = vec2int2((/i1,i2+ncell(2),i3/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,2)
              enddo
            end if

            if (i3<=iextra) then
              k = vec2int2((/i1,i2,i3+ncell(3)/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,3)
              enddo
            end if

          end if

        enddo
      enddo
    enddo
    
    ierr = 1
    do while (ierr/=0)
      nneigh = 0
      ierr = 0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,iv1,ii,iatm,p0,i1,i2,i3,iij,j,jj,jatm,dij,rij,idx,iv2)
!$OMP DO SCHEDULE(static,nnear) REDUCTION(+:ierr)
      main: do idx=1,ntot2
        iv1 = int2vec2(idx,ncell2)
        iv2 = iv1
        if (iv1(1) > ncell(1)) cycle
        if (iv1(2) > ncell(2)) cycle
        if (iv1(3) > ncell(3)) cycle

        do i1=0,iextra
          iij(1) = iv1(1)+i1
          if (iij(1) <= 0) then
            iij(1) = iij(1)+ncell(1)
            iv2(1) = iv1(1)+ncell(1)
          end if
        
          do i2=-iextra,iextra
            if (i1==0 .and. i2<0) cycle
            iij(2) = iv1(2)+i2
            if (iij(2) <= 0) then
              iij(2) = iij(2)+ncell(2)
              iv2(2) = iv1(2)+ncell(2)
            end if
          
            do i3=-iextra,iextra
              if (i1==0 .and. i2==0 .and. i3<0) cycle
              iij(3) = iv1(3)+i3
              if (iij(3) <= 0) then
                iij(3) = iij(3)+ncell(3)
                iv2(3) = iv1(3)+ncell(3)
              end if
            
              i = vec2int2(iv2,ncell2)
              j = vec2int2(iij,ncell2)
            
              if (i == j) then
                do ii=1,nlink(i)
                  iatm=ilink(ii,i)
                  p0(1:3) = rlink(1:3,ii,i)
                
                  do jj=ii+1,nlink(i)
                    jatm=ilink(jj,j)
                  
                    dij = rlink(1:3,jj,j) - p0(1:3)
                    rij = sum(dij*dij)
                  
                    if (rij > rcut2) cycle

                    nneigh(iatm) = nneigh(iatm) + 1
                    if (nneigh(iatm) > neighmax) then
                      ierr = ierr + 1 
                      cycle
                    end if
                    lneigh(nneigh(iatm),iatm) = jatm
                    rneigh(nneigh(iatm),iatm) = sqrt(rij)

                  enddo
                enddo

              else
                do ii=1,nlink(i)
                  iatm=ilink(ii,i)
                  p0(1:3) = rlink(1:3,ii,i)
                
                  do jj=1,nlink(j)
                    jatm=ilink(jj,j)
                  
                    dij = rlink(1:3,jj,j) - p0(1:3)
                    rij = sum(dij*dij)
                  
                    if (rij > rcut2) cycle

                    nneigh(iatm) = nneigh(iatm) + 1
                    if (nneigh(iatm) > neighmax) then
                      ierr = ierr + 1 
                      cycle
                    end if
                    lneigh(nneigh(iatm),iatm) = jatm
                    rneigh(nneigh(iatm),iatm) = sqrt(rij)

                  enddo
                enddo

              end if

              iv2(3) = iv1(3)
            enddo
            iv2(2) = iv1(2)
          enddo
          iv2(1) = iv1(1)
        enddo

      enddo main
!$OMP END DO
!$OMP END PARALLEL
      if (ierr/=0) call extend_nlist(20)
    enddo

  end subroutine computeLikedCellNeighbours

  subroutine computeLikedCellNeighboursSort()
    use m_mrgrnk
    integer :: ierr

    integer :: iatm, jatm, natoms
    real(real64), allocatable, dimension(:,:) :: cartesianCoord

    real(real64) :: rcut2, rij
    integer :: i, j, k, l 
    integer :: ntot, ntot2, nn
    integer, dimension(3) :: iij, iv1, iv2, ncell2
    integer, allocatable, dimension(:) :: nlink
    integer, allocatable, dimension(:,:) :: ilink
    real(real64), allocatable, dimension(:,:,:) :: rlink

    integer :: i1, i2, i3, idx, ii, jj, kk
    real(real64), dimension(3,3) :: hmat
    real(real64), dimension(3) :: p0, dij, dik

    real(real64) :: rcut, rtmp
    integer :: ix2
    integer :: nmax
    real(real64), allocatable, dimension(:) :: allProjections
    integer, allocatable, dimension(:) :: sortedIndices  

    natoms = frame % natoms
    
    allocate(cartesianCoord(3,natoms))
    call fractionalToCartesian(natoms, frame % frac, cartesianCoord)

    rcut2  = cutoffNeighboursListGlobal**2
    ncell2 = ncell+iextra
    ntot   = product(ncell)
    ntot2  = product(ncell2)
    nn     = int(natoms/ntot)+10

    allocate(nlink(ntot2))
    allocate(ilink(nn,ntot2))
    allocate(rlink(3,nn,ntot2))

    nlink = 0
    do i=1,natoms
      iij(1:3) = int(frame % frac(1:3,i)*ncell)+1
      j = vec2int2(iij,ncell2)
      nlink(j) = nlink(j) + 1
      
      ! Reallocate array if necessary
      if (nlink(j) > nn) then
        call extend_ivec(ilink,nn,ntot2,nn+20,ntot2)
        call extend_rmat(rlink,3,nn,ntot2,nn+20,ntot2)
        nn=nn+20
      end if

      ilink(nlink(j),j) = i
      rlink(1:3,nlink(j),j) = cartesianCoord(1:3,i)
    enddo

    hmat = frame % hmat
    do i1=1,ncell(1)
      do i2=1,ncell(2)
        do i3=1,ncell(3)

          if (i1<=iextra .or. i2<=iextra .or. i3<=iextra) then
            j = vec2int2((/i1,i2,i3/),ncell2)

            if (i1<=iextra .and. i2<=iextra .and. i3<=iextra) then
              k = vec2int2((/i1+ncell(1),i2+ncell(2),i3+ncell(3)/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,1)+hmat(:,2)+hmat(:,3)
              enddo
            end if

            if (i1<=iextra .and. i2<=iextra) then
              k = vec2int2((/i1+ncell(1),i2+ncell(2),i3/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,1)+hmat(:,2)
              enddo
            end if

            if (i1<=iextra .and. i3<=iextra) then
              k = vec2int2((/i1+ncell(1),i2,i3+ncell(3)/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,1)+hmat(:,3)
              enddo
            end if

            if (i2<=iextra .and. i3<=iextra) then
              k = vec2int2((/i1,i2+ncell(2),i3+ncell(3)/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,2)+hmat(:,3)
              enddo
            end if

            if (i1<=iextra) then
              k = vec2int2((/i1+ncell(1),i2,i3/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,1)
              enddo
            end if

            if (i2<=iextra) then
              k = vec2int2((/i1,i2+ncell(2),i3/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,2)
              enddo
            end if

            if (i3<=iextra) then
              k = vec2int2((/i1,i2,i3+ncell(3)/),ncell2)
              nlink(k) = nlink(j)
              ilink(:,k) = ilink(:,j)
              do l=1,nlink(j)
                rlink(:,l,k) = rlink(:,l,j)+hmat(:,3)
              enddo
            end if

          end if

        enddo
      enddo
    enddo

    nmax = maxval(nlink)
    allocate(allProjections(nmax))
    allocate(sortedIndices(nmax))

    rcut = cutoffNeighboursListGlobal
    ierr = 1
    do while (ierr/=0)
      nneigh = 0
      ierr = 0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,iv1,ii,iatm,p0,i1,i2,i3,iij,j,jj,jatm,dij,rij,idx,iv2,dik) &
!$OMP PRIVATE(allProjections, sortedIndices, ix2, kk, rtmp) 
!$OMP DO SCHEDULE(static,nnear) REDUCTION(+:ierr)
      main: do idx=1,ntot2
        iv1 = int2vec2(idx,ncell2)
        iv2 = iv1
        if (iv1(1) > ncell(1)) cycle
        if (iv1(2) > ncell(2)) cycle
        if (iv1(3) > ncell(3)) cycle
  
        do i1=0,iextra
          iij(1) = iv1(1)+i1
          if (iij(1) <= 0) then
            iij(1) = iij(1)+ncell(1)
            iv2(1) = iv1(1)+ncell(1)
          end if
        
          do i2=-iextra,iextra
            if (i1==0 .and. i2<0) cycle
            iij(2) = iv1(2)+i2
            if (iij(2) <= 0) then
              iij(2) = iij(2)+ncell(2)
              iv2(2) = iv1(2)+ncell(2)
            end if
          
            do i3=-iextra,iextra
              if (i1==0 .and. i2==0 .and. i3<0) cycle
              iij(3) = iv1(3)+i3
              if (iij(3) <= 0) then
                iij(3) = iij(3)+ncell(3)
                iv2(3) = iv1(3)+ncell(3)
              end if
            
              i = vec2int2(iv2,ncell2)
              j = vec2int2(iij,ncell2)
              
              if (i==j) then
                do ii=1,nlink(i)
                  iatm = ilink(ii,i)
                  do jj=ii+1,nlink(i)
                    jatm = ilink(jj,i)
  
                    dij = rlink(1:3,jj,i) - rlink(1:3,ii,i)
                    rij = sum(dij*dij)
  
                    if (rij > rcut2) cycle

                    nneigh(iatm) = nneigh(iatm) + 1
                    if (nneigh(iatm) > neighmax) then
                      ierr = ierr + 1 
                      cycle
                    end if

                    lneigh(nneigh(iatm),iatm) = jatm
                    rneigh(nneigh(iatm),iatm) = sqrt(rij)
  
                  end do
                end do
  
              else
                dik(:) = (iij(1)-iv2(1)) * rcell(:,1) &
                       + (iij(2)-iv2(2)) * rcell(:,2) &
                       + (iij(3)-iv2(3)) * rcell(:,3)
                dik = dik / sqrt(sum(dik*dik))
  
                kk = nlink(j)
                do jj=1,kk
                  allProjections(jj) = sum(rlink(1:3,jj,j)*dik(1:3))
                end do
                call mrgrnk(allProjections(1:kk), sortedIndices(1:kk))
  
                do ii=1,nlink(i)
                  iatm = ilink(ii,i)
                  rtmp = sum(rlink(1:3,ii,i)*dik(1:3))
  
                  do jj=1,kk
                    ix2 = sortedIndices(jj)
                    if ( (allProjections(ix2) - rtmp) > rcut ) exit
                    
                    dij = rlink(1:3,ix2,j) - rlink(1:3,ii,i)
                    rij = sum(dij*dij)
  
                    if (rij > rcut2) cycle

                    jatm = ilink(ix2,j)

                    nneigh(iatm) = nneigh(iatm) + 1
                    if (nneigh(iatm) > neighmax) then
                      ierr = ierr + 1 
                      cycle
                    end if
                    lneigh(nneigh(iatm),iatm) = jatm
                    rneigh(nneigh(iatm),iatm) = sqrt(rij)
  
                  end do
                end do
              end if 
    
              iv2(3) = iv1(3)
            enddo
            iv2(2) = iv1(2)
          enddo
          iv2(1) = iv1(1)
        enddo
  
      enddo main
!$OMP END DO
!$OMP END PARALLEL
      if (ierr/=0) call extend_nlist(20)
    enddo

  end subroutine computeLikedCellNeighboursSort

  pure function vec2int2(iij,ncell)
    implicit none
    integer, dimension(3), intent(in) :: iij, ncell
    integer :: vec2int2
    vec2int2 = iij(1) + (iij(2)-1)*ncell(1) + (iij(3)-1)*ncell(1)*ncell(2)
  end function vec2int2

  pure function int2vec2(j,ncell)
    implicit none
    integer, intent(in) :: j
    integer, dimension(3), intent(in) ::  ncell(3)
    integer, dimension(3) :: int2vec2
    integer :: ii, jj, kk
    kk = int( (j-1) / ncell(1) / ncell(2) )+1
    jj = int( ((j-1) - (kk-1)*ncell(1)*ncell(2)) / ncell(1) ) +1
    ii = (j-1) - (kk-1)*ncell(1)*ncell(2) - (jj-1)*ncell(1)+1
    int2vec2(1) = ii
    int2vec2(2) = jj
    int2vec2(3) = kk
  end function int2vec2

  subroutine extend_ivec(a,i,j,n,m)
    implicit none

    integer, allocatable, dimension(:,:) :: a
    integer, intent(in) :: i, j

    integer, allocatable, dimension(:,:) :: b
    integer, intent(in) :: n, m

    allocate(b(n,m))
    b(1:i,1:j) = a(1:i,1:j)
    deallocate(a)
    allocate(a(n,m))
    a = b
    deallocate(b)

  end subroutine extend_ivec

  subroutine extend_rmat(a,i,j,k,n,m)
    implicit none

    real(real64), allocatable, dimension(:,:,:) :: a
    integer, intent(in) :: i, j, k

    real(real64), allocatable, dimension(:,:,:) :: b
    integer, intent(in) :: n, m

    allocate(b(i,n,m))
    b(1:i,1:j,1:k) = a(1:i,1:j,1:k)
    deallocate(a)
    allocate(a(i,n,m))
    a = b
    deallocate(b)

  end subroutine extend_rmat

  subroutine extend_nlist(nn)
    implicit none
    integer, intent(in) :: nn
    integer, allocatable, dimension(:,:) :: lneigh_tmp
    real(real64)   , allocatable, dimension(:,:) :: rneigh_tmp
    integer :: itmp

    allocate(lneigh_tmp(neighmax,frame % natoms))
    allocate(rneigh_tmp(neighmax,frame % natoms))

    lneigh_tmp(1:neighmax,1:frame % natoms) = lneigh(1:neighmax,1:frame % natoms)
    rneigh_tmp(1:neighmax,1:frame % natoms) = rneigh(1:neighmax,1:frame % natoms)

    deallocate(lneigh)
    deallocate(rneigh)

    itmp=neighmax
    neighmax=neighmax+nn

    allocate(lneigh(neighmax,frame % natoms))
    allocate(rneigh(neighmax,frame % natoms))

    lneigh(1:itmp,1:frame % natoms) = lneigh_tmp(1:itmp,1:frame % natoms)
    rneigh(1:itmp,1:frame % natoms) = rneigh_tmp(1:itmp,1:frame % natoms)
    lneigh(itmp+1:neighmax,1:frame % natoms) = 0
    rneigh(itmp+1:neighmax,1:frame % natoms) = 0.0_real64

    deallocate(lneigh_tmp)
    deallocate(rneigh_tmp)

  end subroutine extend_nlist

end module moduleNeighbours 
