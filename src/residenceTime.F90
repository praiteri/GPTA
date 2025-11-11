!disclaimer
module moduleResidenceTime

  use moduleVariables
  use moduleFiles
  use moduleSystem
  use moduleMessages

  implicit none

  public :: computeResidenceTime, computeResidenceTimeHelp
  private
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  
  integer, pointer :: numberOfBins
  integer, pointer :: tallyExecutions
  real(real64), pointer :: cutoffRadius
  integer, pointer, dimension(:,:,:) :: coordinationList

  integer, pointer :: numberOfCentres
  integer, pointer :: MAX_NEIGH
  integer, pointer :: threshold

  logical, pointer :: lsolute

contains

  subroutine computeResidenceTimeHelp()
    implicit none
    call message(0,"This action computes the solvent residence time around selected solute atoms.")
    call message(0,"For efficiency reasons it requires to provide the number of frames in the trajectory.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --restime +s Ca OW +rcut 3.2 +thres 1 +ntraj 1000")
  end subroutine computeResidenceTimeHelp

  subroutine computeResidenceTime(a)
    use moduleSystem 
    implicit none
    type(actionTypeDef), target :: a

    call associatePointers(a)

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      if (firstAction) then
        call dumpScreenInfo()

        if (resetFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if

        ! create a list of the atoms' indices for each group
        call selectAtoms(2,actionCommand,a)
        call createSelectionList(a,2)

        numberOfCentres = count(a % isSelected(:,1))
        if (numberOfCentres == 0) call message(-1,"--restime - no solute atoms selected")

        if (count(a % isSelected(:,2)) == 0) call message(-1,"--restime - no solvent atoms selected")

        allocate(a % Iarray3D(MAX_NEIGH,numberOfCentres,numberOfBins) , source=0)
        call associatePointers(a)

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else

        if (a % updateAtomsSelection) call selectAtoms(1,actionCommand,a)

      end if

      call computeAction(a)

    end if

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) call finaliseAction()

  end subroutine computeResidenceTime

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile
    numberOfBins         => a % numberOfBins
    cutoffRadius         => a % doubleVariables(1)
    numberOfCentres      => a % integerVariables(1)
    MAX_NEIGH            => a % integerVariables(2)
    threshold            => a % integerVariables(3)
    coordinationList     => a % Iarray3D
    lsolute              => a % logicalVariables(1)

  end subroutine associatePointers

  subroutine initialiseAction(a)
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a 

    a % actionInitialisation = .false.

    call assignFlagValue(actionCommand,"+rcut",cutoffRadius,4.0_real64)
    call assignFlagValue(actionCommand,"+ntraj",numberOfBins,1000)
    call assignFlagValue(actionCommand,"+nmax",MAX_NEIGH,10)
    call assignFlagValue(actionCommand,"+thres",threshold,1)
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'restime.out')

    call assignFlagValue(actionCommand,"+solute",lsolute,.false.)
    if (lsolute) then
      a % requiresNeighboursList = .false.
    else
      a % requiresNeighboursList = .true.
      a % requiresNeighboursListUpdates = .true.
      a % requiresNeighboursListDouble = .true.
    end if

    a % cutoffNeighboursList = cutoffRadius

    tallyExecutions = 0 

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    use moduleMessages 
    implicit none
    call message(0,"Residence time calculation setup")
    call message(0,"...Maximum number of neighbours",i=MAX_NEIGH)
    call message(0,"...Maximum trajectory length",i=numberOfBins)
    call message(0,"...Frames threshold",i=threshold)
    call message(0,"...Cutoff distance",r=cutoffRadius)
    call message(1,"...Output file",str=outputFile % fname)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleVariables
    use moduleSystem 
    use moduleDistances
    use moduleMessages
    implicit none
    type(actionTypeDef), target :: a

    integer :: iatm, jatm, ineigh
    integer, allocatable, dimension(:,:) :: localList
    integer :: icentre, nshell, nsel, nsolvent
    real(real64) :: dij(3), rdist, r2
    
    allocate(localList(MAX_NEIGH,numberOfCentres), source=0)
    
    icentre = 0
    nsel = numberOfCentres

    if (lsolute) then
      nsolvent = count(a % isSelected(:,2))
      r2 = cutoffRadius**2
      do icentre = 1 , nsel
        iatm = a % idxSelection(icentre,1)
        nshell = 0
        do ineigh = 1 , nsolvent
          jatm = a % idxSelection(ineigh,2)
          dij = frame % pos(:,icentre) - frame % pos(:,jatm)
          rdist = computeDistanceSquaredPBC(dij)
          if (rdist > r2) cycle
          nshell = nshell + 1
          localList(nshell,icentre) = jatm
        end do
        if (nshell > MAX_NEIGH) call message(-1,"--restime | too many neighbours (+nmax)",i=nshell)
      end do

    else
      do icentre=1,nsel
        iatm = a % idxSelection(icentre,1)
        nshell = 0
        do ineigh = 1 , nneigh(iatm)
          if (rneigh(ineigh,iatm)>cutoffRadius) cycle
          jatm = lneigh(ineigh,iatm)
          if (a % isSelected(jatm,2)) then
            nshell = nshell + 1
            localList(nshell,icentre) = jatm
          end if
        end do
        if (nshell > MAX_NEIGH) call message(-1,"--restime | too many neighbours (+nmax)",i=nshell)
      end do
    end if
    if (frame % nframe > numberOfBins) call message(-1,"--restime | trajctory is too long (+ntraj)")
    coordinationList(:,:,frame % nframe) = localList

  end subroutine computeAction

  subroutine finaliseAction()
    implicit none
    integer :: idx, nf
    integer, allocatable, dimension(:,:) :: localList
    real(real64), allocatable, dimension(:) :: residenceTime
    real(real64), allocatable, dimension(:) :: exchangeProbability
#ifdef GPTA_MPI
    integer :: nn
#endif
    
    nf = frame % nframe

#ifdef GPTA_MPI
    call MPI_allreduce(MPI_IN_PLACE, nf, 1, MPI_INT, MPI_MAX, MPI_Working_Comm, ierr_mpi)
#endif

    allocate(residenceTime(0:nf), source=0.0_real64)
    allocate(exchangeProbability(0:nf), source=0.0_real64)
    allocate(localList(MAX_NEIGH,nf), source=0)

    do idx=1,numberOfCentres
      ! Coordination shell for one atoms in the list
      localList(1:MAX_NEIGH, 1:nf) = coordinationList(1:MAX_NEIGH, idx, 1:nf)

#ifdef GPTA_MPI
      nn = MAX_NEIGH * nf
      if (me == 0) then
        call MPI_reduce(MPI_IN_PLACE, localList, nn, MPI_INT, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)
      else 
        call MPI_reduce(localList, localList, nn, MPI_INT, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)
      end if
#endif

      call rest_original(MAX_NEIGH, nf, localList, threshold, residenceTime, exchangeProbability)
    end do

    if (me /= 0) return
    residenceTime = residenceTime / numberOfCentres
    residenceTime = residenceTime / residenceTime(0)

!    exchangeProbability = exchangeProbability / numberOfCentres
!    exchangeProbability = exchangeProbability / sum(exchangeProbability)

    call initialiseFile(outputFile,outputFile % fname)

    write(outputFile % funit,'("# Frames   | Survival func | Exchange probability")')
    do idx=0,nf
      write(outputFile % funit,'(i6,2(1x,f16.5))')idx, residenceTime(idx), exchangeProbability(idx)
    end do
  end subroutine finaliseAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rest_original(MAX_NEIGH,frames_read,local_list,thres,survivalFunction,distribution)
    implicit none
    integer, intent(in) :: frames_read, MAX_NEIGH

    integer, intent(in) :: thres
    integer, dimension(MAX_NEIGH,frames_read) :: local_list
    real(real64), dimension(0:frames_read) :: survivalFunction
    real(real64), dimension(0:frames_read) :: distribution

    integer, allocatable, dimension(:,:) :: coord

    integer :: i, j, n, n0

    ! coord(1,:) Current coordiantion shell
    ! coord(2,:) Time since the molecule has entered the coordination shell
    ! coord(3,:) Time since the molecule has left the coordination shell
    allocate(coord(3,MAX_NEIGH*2))
    coord=0

    n0 = count(local_list(:,1) > 0)

    ! Initial coordination shell
    coord(1,1:n0) = local_list(1:n0,1)

    do i=2,frames_read

      ! Number of atoms that were in the coordination shell in the previous step
      n0 = count(coord(1,:) > 0)
      ! Check if they are still alive
      do j=1,n0

        ! Yes
        if (any(local_list(:,i)==coord(1,j))) then
          if (coord(3,j)>0) then
            coord(2,j) = coord(2,j) + coord(3,j) ! Just passed out for a little while
            coord(3,j) = 0
          endif
          coord(2,j) = coord(2,j) + 1

          ! Not there anymore
        else
          coord(3,j) = coord(3,j) + 1
        endif

        ! It's been dead for too long - bury it
        ! -> add to survival function
        if ( coord(3,j) > thres) then
          if (coord(2,j) > thres) then
            survivalFunction(0:coord(2,j)) = survivalFunction(0:coord(2,j)) + 1.0_real64
            distribution(coord(2,j)) = distribution(coord(2,j)) + 1
          end if
          coord(:,j) = 0
        endif

      enddo

      ! Compact the list
      n=0
      do j=1,n0
        if (coord(1,j)>0) then
          n = n + 1
          coord(:,n) = coord(:,j)
        endif
      enddo

      ! Add new neighbours at the end of the list
      do j=1,count(local_list(:,i)>0)
        if (all(local_list(j,i)/=coord(1,:))) then
          n = n + 1
          coord(1,n) = local_list(j,i)
          coord(2:3,n) = 0
        endif
      enddo

      ! Clear the end of the array
      do j=n+1,MAX_NEIGH
        coord(:,j) = 0
      enddo

    enddo

    ! Add the molecules still in the list to the survival function
    do j=1,n
      survivalFunction(0:coord(2,j)) = survivalFunction(0:coord(2,j)) + 1.0_real64
    enddo

  end subroutine rest_original

end module moduleResidenceTime
