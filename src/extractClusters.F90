!disclaimer
module moduleExtractClustersAction
  use moduleVariables
  use moduleFiles
  use moduleMessages
  use moduleSystem
  use moduleStrings
  use moduleDistances
  
  implicit none

  public :: extractClusters, extractClustersHelp
  private

  character(:), pointer :: actionCommand
  type(fileTypeDef), pointer :: outputFile
  logical, pointer :: firstAction

  integer, pointer :: ID
  real(real64), pointer :: clusterRadiusSQ
  real(real64), pointer, dimension(:) :: clusterCentre
  integer, pointer :: cluster_type

  real(real64), allocatable, dimension(:,:) :: hkl
  real(real64), allocatable, dimension(:) :: asymm
  real(real64), allocatable, dimension(:) :: thickness
  real(real64) :: rotation(4)
contains

  subroutine extractClustersHelp()
    implicit none
    call message(0,"This action extracts a spherical cluster of radius R around the selected atom.")
    call message(0,"Examples:")
    call message(0,"  gpta --i coord.pdb --cluster +r 10 +id 10 +rot 1,1,1,38")
    call message(0,"  gpta --i coord.pdb --cluster +r 15 +centre 10,10,10 +rot 1,1,1,38")
    call message(0,"  gpta --i coord.pdb --cluster +hkl 1,0,0 1,-1,0 0,1,0 0,0,1 +thick 0.3,0.3,0.3,0.2 +asym 0.2,0.2,0.2,0 +rot 1,1,1,32")
  end subroutine extractClustersHelp

  subroutine extractClusters(a)
    implicit none
    type(actionTypeDef), target :: a

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then

      if (firstAction) then
        call dumpScreenInfo()
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if
      
      call computeAction()
    end if

  end subroutine extractClusters

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a
    character(STRLEN) :: flagString = ""
    integer :: i, ntmp, nvec
    character(len=STRLEN), dimension(100) :: listOfStrings
    character(len=STRLEN), dimension(3) :: string_l

    a % actionInitialisation = .false.

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    outputFile           => a % outputFile

    ID                   => a % integerVariables(1)
    cluster_type         => a % integerVariables(2)

    clusterRadiusSQ      => a % doubleVariables(1)
    clusterCentre        => a % doubleVariables(2:4)


    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .true.
    a % requiresNeighboursListDouble = .true.
    a % cutoffNeighboursList = 3.2_real64

    call assignFlagValue(actionCommand,"+out ",outputFile % fname,'clusters.xyz')

    ! Spherical cluster parameters
    ! Parse centre (or atoms ID) and radius
    call assignFlagValue(actionCommand,"+id ",ID,0)
    call assignFlagValue(actionCommand,"+centre ",clusterCentre,[0.0_real64,0.0_real64,0.0_real64])
    call assignFlagValue(actionCommand,"+r ",clusterRadiusSQ,0.0_real64)

    call extractFlag(actionCommand,"+hkl ",flagString)
    if (len_trim(flagString) > 0) then
      call parse(flagString," ",listOfStrings,nvec)
      flagString = ''
      if (nvec==0) then
        call message(-1,"--new_cell : no vectors specified")
      else 
        allocate(hkl(3,nvec), source=0.0_real64)
      end if

      do i=1,nvec
        call parse(listOfStrings(i),",",string_l,ntmp)
        if (ntmp /= 3) then
          call message(-1,"--new_cell : vector must have 3 components")
        end if
        read(string_l(1),*) hkl(1,i)
        read(string_l(2),*) hkl(2,i)
        read(string_l(3),*) hkl(3,i)
      end do

      allocate(asymm(nvec), source=0.0_real64)
      allocate(thickness(nvec), source=0.4_real64)
      call extractFlag(actionCommand,"+asym",flagString)
      if (len_trim(flagString) > 0) then
        call parse(flagString,",",listOfStrings,ntmp)
        if (ntmp /= nvec) then
          call message(-1,"--new_cell : number of asymmetry parameters does not match number of vectors")
        end if
        do i=1,nvec
          read(listOfStrings(i),*) asymm(i)
        end do
      end if

      call extractFlag(actionCommand,"+thick ",flagString)
      if (len_trim(flagString) > 0) then
        call parse(flagString,",",listOfStrings,ntmp)
        if (ntmp /= nvec) then
          call message(-1,"--new_cell : number of thickness parameters does not match number of vectors")
        end if

        do i=1,nvec
          read(listOfStrings(i),*) thickness(i)
        end do
      end if
    end if

    ! Spherical cluster
    if (clusterRadiusSQ > 0.0_real64) then
      clusterRadiusSQ = clusterRadiusSQ**2
      if (ID == 0) then
        cluster_type = 1
      else  
        cluster_type = 2
      end if

    ! Wulff construction
    else
        cluster_type = 3
    end if

    call assignFlagValue(actionCommand,"+rot",rotation,[0.0_real64,0.0_real64,0.0_real64,0.0_real64])

    return

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    implicit none
    integer :: i

    if (cluster_type == 1 .or. cluster_type == 2) then
      call message(0,"Extract Spherical Cluster Action")
      call message(0,"...Radius",r=sqrt(clusterRadiusSQ))

      if (cluster_type == 1) then
        call message(0,"...Centre",rv=clusterCentre(1:3))
      else if (cluster_type == 2) then
        call message(0,"...Central atom ID",i=ID)
      end if

    else if (cluster_type == 3) then
      call message(0,"Extract Wulff Cluster Action")
      do i=1,size(hkl)/3
        call message(0,"...hkl vector",rv=hkl(:,i))
      end do
      call message(0,"...thickness parameters",rv=thickness)
      call message(0,"...asymmetry parameters",rv=asymm)
    else
      call message(-1,"Unknown cluster type")
    end if
    if ( abs(rotation(4)) > 0.0_real64 ) then
      call message(0,"...Rotation axis",rv=rotation(1:3))
      call message(0,"...Rotation angle",r=rotation(4))
    end if
    call message(2)

  end subroutine dumpScreenInfo

  subroutine computeAction()
    implicit none
    real(real64) :: rtmp(3)

    if (cluster_type == 1) then
      rtmp = clusterCentre(1:3)
      call cutSphericalCluster(rtmp)
    else if (cluster_type == 2) then
      rtmp(1:3) = frame % pos(1:3,ID)
      call cutSphericalCluster(rtmp)
    else if (cluster_type == 3) then
      if (allocated(hkl)) then
        call cutWulffCluster()
      else
        call message(-1,"Wulff cluster requires --hkl flag")
      end if
    else
      call message(-1,"Unknown cluster type")
    end if

  end subroutine computeAction

  subroutine cutSphericalCluster(centre)
    implicit none
    real(real64), intent(in) :: centre(3)

    integer :: idx
    integer :: iatm, jatm
    real(real64) :: dij(3), dist
    logical, allocatable, dimension(:) :: lmol
    integer :: imol, nsel
    character(cp), allocatable, dimension(:) :: llab
    real(real64), allocatable, dimension(:,:) :: lpos
    character(len=6) :: str
    type(fileTypeDef) :: localfile
    real(real64), dimension(3) :: shift

    nsel = 0

    if (numberOfMolecules > 0) then
      allocate(lmol(numberOfMolecules), source=.false.)
      do iatm=1,frame % natoms
        dij(1:3) = frame % pos(1:3,iatm) - centre(1:3)
        dist = computeDistanceSquaredPBC(dij)
        if (dist < clusterRadiusSQ) then
          imol = atomToMoleculeIndex(iatm)
          if (lmol(imol)) cycle
          lmol(imol) = .true.
          nsel = nsel + listOfMolecules(imol) % numberOfAtoms
        end if
      end do

      allocate(llab(nsel))
      allocate(lpos(3,nsel))

      idx = 0
      do imol=1,numberOfMolecules
        if (lmol(imol)) then
          do iatm=1,listOfMolecules(imol) % numberOfAtoms
            idx = idx + 1
            jatm = listOfMolecules(imol) % listOfAtoms(iatm)
            lpos(:,idx) =  frame % pos(:,jatm)
            llab(idx) = frame % lab(jatm)
          end do
        end if
      end do

    else
      allocate(lmol(frame % natoms), source=.false.)
      do iatm=1,frame % natoms
        dij(1:3) = frame % pos(1:3,iatm) - centre(1:3)
        dist = computeDistanceSquaredPBC(dij)
        if (dist < clusterRadiusSQ) then
          lmol(iatm) = .true.
          nsel = nsel + 1
        end if
      end do
      allocate(llab(nsel))
      allocate(lpos(3,nsel))

      idx = 0
      do iatm=1,frame % natoms
        if (lmol(iatm)) then
          idx = idx + 1
          lpos(:,idx) =  frame % pos(:,iatm)
          llab(idx) = frame % lab(iatm)
        end if
      end do

    end if

    ! translate the cluster to the centre of the cell
    shift = - centre(1:3)
    shift = shift + 0.5_real64 * (frame % hmat(:,1) + frame % hmat(:,2) + frame % hmat(:,3))

    do iatm=1,nsel
      lpos(1:3,iatm) = lpos(1:3,iatm) + shift(1:3)
    end do
    call wrapCoordinates(nsel,lpos)
    do iatm=1,nsel
      lpos(1:3,iatm) = lpos(1:3,iatm) - shift(1:3)
    end do

    if (abs(rotation(4)) > 0.0_real64) then
      call rotateMolecule(rotation(1:3), rotation(4), nsel, lpos, 0)
    end if

    block
      integer :: idx, ilen
      character(len=4) :: ftype
      ilen = len_trim(outputFile % fname)
      idx = index(outputFile % fname,".",back=.true.)
      if (ilen-idx>0) then
        ftype = outputFile % fname(idx+1:ilen)
      else
        ftype='NULL'
      end if
      
      if (ftype == "xyz") then
        localFile % fname = outputFile % fname(1:idx-1)
      else 
        localFile % fname = trim(outputFile % fname)
      end if
    end block
    
    write(str,'(i0)') frame % nframe
    localFile % fname = trim(localFile % fname)//"."//trim(str)//".xyz"
    call initialiseFile(localFile, localFile % fname)
    ! call dumpXYZ(localFile % funit, nsel, lpos, llab)

    call dumpCoordinates("xyz", localFile % funit, nsel, lpos, llab)
 
    close(localFile % funit)

    deallocate(llab)
    deallocate(lpos)
    
  end subroutine cutSphericalCluster

  subroutine cutWulffCluster()
    implicit none

    integer :: iatm, jatm, imol, idx
    integer :: nsel, nvec, ivec, nmols

    real(real64), dimension(3,3) :: hinv
    real(real64), dimension(3) :: xcom, dij, lvector

    logical, allocatable, dimension(:) :: lmol
    character(cp), allocatable, dimension(:) :: llab
    real(real64), allocatable, dimension(:,:) :: lpos

    real(real64), allocatable, dimension(:,:) :: ftmp, ptmp
    real(real64) :: sij

    hinv(:,1) = cross_product(frame % hmat(:,2), frame % hmat(:,3))
    hinv(:,2) = cross_product(frame % hmat(:,3), frame % hmat(:,1))
    hinv(:,3) = cross_product(frame % hmat(:,1), frame % hmat(:,2))
    hinv(:,1) = hinv(:,1) / sqrt(dot_product(hinv(:,1),hinv(:,1)))
    hinv(:,2) = hinv(:,2) / sqrt(dot_product(hinv(:,2),hinv(:,2)))
    hinv(:,3) = hinv(:,3) / sqrt(dot_product(hinv(:,3),hinv(:,3)))

    xcom(1:3) = sum(frame % pos,dim=2) / real(frame % natoms,real64)

    nvec = size(hkl) / 3

    allocate(llab(frame % natoms))
    allocate(lpos(3,frame % natoms))
    nsel = 0
    main2: do iatm=1,frame % natoms
      dij(1:3) = frame % pos(1:3,iatm) - xcom(1:3)
      do ivec=1,nvec
        lvector = hkl(1,ivec) * hinv(:,1) + hkl(2,ivec) * hinv(:,2) + hkl(3,ivec) * hinv(:,3)
        sij = dot_product(dij,lvector)
        if (sij >  thickness(ivec)+asymm(ivec)) cycle main2
        if (sij < -thickness(ivec)+asymm(ivec)) cycle main2
      end do
      nsel = nsel + 1
      llab(nsel) = frame % lab(iatm)
      lpos(1:3,nsel) = frame % pos(1:3,iatm)
    end do main2

    if (abs(rotation(4)) > 0.0_real64) then
      call rotateMolecule(rotation(1:3), rotation(4), nsel, lpos, 0)
    end if

    call initialiseFile(outputFile, outputFile % fname)
    call dumpCoordinates(outputFile % ftype, outputFile % funit, nsel, lpos, llab)
    close(outputFile % funit)

    return
    ! nsel = 0

    ! nvec = size(hkl) / 3
    ! allocate(ptmp(3,frame % natoms))
    ! call cartesianToFractionalNINT(frame % natoms, frame % pos , ftmp)
    ! call fractionalToCartesian(frame % natoms, ftmp, ptmp)
    ! frame % pos = ptmp

    ! if (numberOfMolecules > 0) then
    !   nmols = numberOfMolecules
    !   allocate(lmol(numberOfMolecules), source=.false.)
    !   main: do iatm=1,frame % natoms
    !     do ivec=1,nvec
    !       sij = dot_product(ftmp(:,iatm),hkl(:,ivec))
    !       if (sij >  thickness(ivec)+asymm(ivec)) cycle main
    !       if (sij < -thickness(ivec)+asymm(ivec)) cycle main
    !     end do

    !     imol = atomToMoleculeIndex(iatm)
    !     if (lmol(imol)) cycle main
    !     lmol(imol) = .true.
    !     nsel = nsel + listOfMolecules(imol) % numberOfAtoms
    !   end do main


    !   idx = 0
    !   do imol=1,nmols
    !     if (lmol(imol)) then
    !       do iatm=1,listOfMolecules(imol) % numberOfAtoms
    !         idx = idx + 1
    !         jatm = listOfMolecules(imol) % listOfAtoms(iatm)
    !         lpos(:,idx) =  frame % pos(:,jatm)
    !         llab(idx) = frame % lab(jatm)
    !       end do
    !     end if
    !   end do

    ! else

    !   allocate(lmol(frame % natoms), source=.false.)
    !   nmols = frame % natoms 

    !   main1: do iatm=1,frame % natoms
    !     do ivec=1,nvec
    !       sij = dot_product(ftmp(:,iatm),hkl(:,ivec))
    !       if (sij >  thickness(ivec)+asymm(ivec)) cycle main1
    !       if (sij < -thickness(ivec)+asymm(ivec)) cycle main1
    !     end do

    !     lmol(iatm) = .true.
    !     nsel = nsel + 1
        
    !   end do main1

    !   allocate(llab(nsel))
    !   allocate(lpos(3,nsel))

    !   idx = 0
    !   do iatm=1,frame % natoms
    !     if (lmol(iatm)) then
    !       idx = idx + 1
    !       lpos(:,idx) =  frame % pos(:,iatm)
    !       llab(idx) = frame % lab(iatm)
    !     end if
    !   end do

    ! end if

    ! if (abs(rotation(4)) > 0.0_real64) then
    !   call rotateMolecule(rotation(1:3), rotation(4), nsel, lpos, 0)
    ! end if

    ! call initialiseFile(outputFile, outputFile % fname)
    ! call dumpCoordinates(outputFile % ftype, outputFile % funit, nsel, lpos, llab)
    ! close(outputFile % funit)

    deallocate(llab)
    deallocate(lpos)

  end subroutine cutWulffCluster

end module moduleExtractClustersAction
