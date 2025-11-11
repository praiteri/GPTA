!disclaimer
program gpta
  use moduleVariables
  use moduleRead 
  use moduleActions
  use moduleRandomNumbers
  use moduleSystem 
  use moduleElements 
  use moduleMessages 
  use moduleNeighbours , only : initialiseNeighboursList, updateNeighboursList, neighboursListTime
  use moduleProperties
  use moduleDistances, only : cartesianToFractional
  
  implicit none

! Timing
  real(real64) :: startTimer, endTimer
  real(real64) :: read_start, read_end, readingTime, initialisationTime
  real(real64) :: actionTime, commTime
  character(len=20) :: str
  
  integer :: iact
  
#ifdef GPTA_MPI  
  integer :: iproc
  integer :: islave, jslave
  integer, allocatable, dimension(:) :: workginParty
  real(real64) :: comm_start, comm_end
#endif

#ifdef GPTA_MPI
  ! Initialise MPI
  call MPI_Init (ierr_mpi)
  call MPI_Comm_rank (MPI_COMM_WORLD, me, ierr_mpi)
  call MPI_Comm_size (MPI_COMM_WORLD, numberOfMpiProcesses, ierr_mpi)

  ! readingCPU is the rank of the MPI process reading the trajectory
  ! readingCPU writes all the output from the message routine
  readingCPU = numberOfMpiProcesses - 1
  numberOfActiveProcesses = numberOfMpiProcesses
  islave = 0
  jslave = readingCPU-1
  if (numberOfMpiProcesses > 4) then
    masterReadingOnly = .true.
    numberOfActiveProcesses = numberOfActiveProcesses - 1
  end if
  allocate(workginParty(0:numberOfActiveProcesses-1))
  do iact=0,numberOfActiveProcesses-1
    workginParty(iact) = iact
  end do
  
  ! all properties are written from MPI rank 0, which is always computing the properties
  ! need a communicator for all the processes that compute properties
  call MPI_Comm_group(MPI_COMM_WORLD, MPI_world, ierr_mpi)
  call MPI_Group_incl(MPI_world, numberOfActiveProcesses, workginParty, MPI_Working_Group, ierr_mpi)
  call MPI_Comm_create_group(MPI_COMM_WORLD, MPI_Working_Group, 0, MPI_Working_Comm, ierr_mpi)

#else
  me = 0
  readingCPU = 0
  numberOfMpiProcesses = 1
  numberOfActiveProcesses = 1
#endif

  nProgress = 100 * numberOfActiveProcesses

! Start timing
  startTimer = timing()
  readingTime = 0.0_real64
  commTime = 0.0_real64

  ! Initialise atoms' masses, covalent radii...
  call initialisePeriodicTable()

  call parseActionLine()

  call associateProperties()
  
  ! Process and remove one off commands (like frames)
  call executeOneOffActions()

  ! Initialise random number generator
  if (randomNumberSeed <= 0) then
    ! If no seed is given, use the current time as a seed
    call get_time_seed(randomNumberSeed)
    randomNumberSeed = randomNumberSeed + me
  end if
  call init_genrand(randomNumberSeed)

#ifdef GPTA_OMP
  ! Initialise opemMP
  if (ompNumThreads==0) then
    if (numberOfMpiProcesses>1) then
      ompNumThreads = 1
    else
      ompNumThreads = min(8,omp_get_max_threads())
    end if
  end if
  call omp_set_num_threads(ompNumThreads)
#endif
  
  ! Write time stamp out after parsing the command line for --log
  if (me == 0) call timestamp(0)

  ! Open input coordinates files
  call openCoordinatesInputFiles()
  
#ifdef GPTA_MPI
  call MPI_Bcast(numberInputFiles, 1, MPI_INT, readingCPU, MPI_COMM_WORLD, ierr_mpi)
#endif 

  ! If there are no input coordinates files to read spin off here
  if (numberInputFiles == 0) then
  
    if (numberOfMpiProcesses > 1) call message(-1,"Cannot run MPI version without an input trajectory")

    frame % nframe = -1
    call initialiseActions()
    frame % nframe = 0
    frameReadSuccessfully = .true.
    call runAllActions()
    endOfCoordinatesFiles = .true.

  else

    ! Get number of atoms from first file
    if (me == readingCPU) then
      call getNumberOfAtoms(inputFileNames(currentInputFile), inputNumberOfAtoms, frame % hmat)
    end if
    
#ifdef GPTA_MPI
    call MPI_Bcast(inputNumberOfAtoms, 1, MPI_INT, readingCPU, MPI_COMM_WORLD, ierr_mpi)
    call MPI_Bcast(frame % hmat, 9, MPI_DOUBLE, readingCPU, MPI_COMM_WORLD, ierr_mpi)
#endif 
    if (me == 0) call dumpCellInfo(frame % hmat)

    ! Initialise actions
    call initialiseActions()
    
    ! Initialise neighbours list
    call initialiseNeighboursList()
    
    ! Allocate system arrays
    frame % nframe = 0
    call createSystemArrays(frame, inputNumberOfAtoms)
  end if

#ifdef GPTA_MPI
  if (numberOfMpiProcesses > 1) then
    if (me == readingCPU) then
      allocate(allFrames(islave:jslave))
      allocate(allFramesReadSuccessfully(islave:jslave), source=.true.)
      do iproc=islave,jslave
        call createSystemArrays(allFrames(iproc), inputNumberOfAtoms)
      end do
    end if
  end if
#endif

  initialisationTime = timing() - startTimer

  ! Start processing coordinates files
  do while (.not. endOfCoordinatesFiles)

#ifdef GPTA_MPI
    ! CPU 0 reads all the frames
    if (me == readingCPU) then

      ! Read coordinates
      read_start = timing()
      do iproc=islave,jslave
        call readCoordinates(allFramesReadSuccessfully(iproc), allFrames(iproc))
        if (.not. allFramesReadSuccessfully(iproc)) endOfCoordinatesFiles = .true.
      end do
      if (.not. masterReadingOnly) then
        call readCoordinates(frameReadSuccessfully, frame)
        endOfCoordinatesFiles = (endOfCoordinatesFiles .or. (.not. frameReadSuccessfully))
      end if
      read_end = timing()
      readingTime = readingTime + (read_end-read_start)

      ! Send coordinates to slave nodes
      do iproc=islave,jslave
        call MPI_Send(allFramesReadSuccessfully(iproc), 1, MPI_LOGICAL, iproc, 0, MPI_COMM_WORLD, ierr_mpi)
        if (allFramesReadSuccessfully(iproc)) call sendFrameToProcessor(allFrames(iproc), iproc)
      end do

    else

      ! Slaves receive coordinates
      comm_start = timing()
      call MPI_Recv(frameReadSuccessfully, 1, MPI_LOGICAL, readingCPU, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)
      if (frameReadSuccessfully) call receiveFrameFromMaster(frame)
      comm_end = timing()
      commTime = commTime + (comm_end-comm_start)

    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

    ! Tell all processors if we're at the end of the file
    call MPI_Bcast(endOfCoordinatesFiles, 1, MPI_LOGICAL, readingCPU, MPI_COMM_WORLD, ierr_mpi)
    if (masterReadingOnly .and. me == readingCPU) cycle

#else

    read_start = timing()
    call readCoordinates(frameReadSuccessfully,frame)
    if (.not. frameReadSuccessfully) endOfCoordinatesFiles = .true.
    read_end = timing()
    readingTime = readingTime + (read_end-read_start)
    
#endif

    if (io == 6) then
      if (mod(frame % nframe - first_frame , nProgress)==0) then
        if (frame % nframe > nProgress) call message(-2) !??
        call message(0,"---Frames read",frame % nframe)
!        write(io,'("---Frames read",i0)')frame % nframe
      end if
    end if

    call cellProcessing() 
    call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
    
    if (frameReadSuccessfully) then
      if (lastFrameOnly) cycle
    else
      if (lastFrameOnly) frameReadSuccessfully = .true.
    end if
      
    if (frameReadSuccessfully) numberOfFramesProcessed = numberOfFramesProcessed + 1
    
    if (computeNeighboursListOnce .and. computeNeighboursList) then
      call updateNeighboursList(.false.)
      computeNeighboursList = .false.
    end if

    call runAllActions()

  end do ! End processing coordinates files

  endTimer = timing()
    
#ifdef GPTA_MPI
  call MPI_allreduce(MPI_IN_PLACE, numberOfFramesProcessed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr_mpi)
  call MPI_allreduce(MPI_IN_PLACE, neighboursListUpdates, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr_mpi)
  call MPI_allreduce(MPI_IN_PLACE, neighboursListTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr_mpi)
  neighboursListTime = neighboursListTime / numberOfActiveProcesses

  do iact=1,numberOfActions
    call MPI_allreduce(MPI_IN_PLACE, action(iact) % timeTally, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr_mpi)
  end do
  action % timeTally = action % timeTally / numberOfActiveProcesses

  call MPI_Bcast(readingTime, 1, MPI_DOUBLE, readingCPU, MPI_COMM_WORLD, ierr_mpi)
  call MPI_allreduce(MPI_IN_PLACE, commTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr_mpi)
  commTime = commTime / numberOfActiveProcesses
#endif

  if (me == 0) then
    actionTime = 0.0_real64
    do iact=1,numberOfActions
      actionTime = actionTime + action(iact) % timeTally
    end do
  
    call message(0,"Frames processing summary")
    call message(0,"...Number of frames read",i=numberOfFramesRead)
    call message(1,"...Number of frames processed",i=numberOfFramesProcessed)
    call message(0,"Timing summary (seconds)")
    call message(0,"...Initialisation time",r=initialisationTime)
    call message(0,"...Total reading time",r=readingTime)
    call message(0,"...Communication time",r=commTime)
  
    call message(0,"...Total time in actions",r=actionTime)
    do iact=1,numberOfActions
      actionTime = actionTime + action(iact) % timeTally
      call message(0,"......Total time in "//trim(actionType(iact)),r=action(iact) % timeTally)
    end do
  
    write(str,'("[n=",i0,"]")') neighboursListUpdates
    call message(0,"...Neighbours list time "//trim(str),r=neighboursListTime)
    call message(1,"Total elapsed time",r=(endTimer-startTimer))
    call timestamp(1)
    call debugTiming(-1)
    call debugTiming_old(-1)
  end if

#ifdef GPTA_MPI
  call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)
  call MPI_Finalize (ierr_mpi)
#endif

  call execute_command_line("rm -f _tmp-gpta_*")
  call exit(0)

contains

  subroutine timestamp(idx)
    use moduleMessages 
#ifdef GPTA_MPI
    use moduleVariables, only : numberOfMpiProcesses
#endif
    implicit none
    integer, intent(in) :: idx

    character(len=STRLEN) :: str

    integer ( kind = 8 ) :: d
    integer ( kind = 8 ) :: h
    integer ( kind = 8 ) :: m
    integer ( kind = 8 ) :: mm(8)
    character ( len = 9 ), parameter, dimension(12) :: month = (/ &
      'January  ', 'February ', 'March    ', 'April    ', &
      'May      ', 'June     ', 'July     ', 'August   ', &
      'September', 'October  ', 'November ', 'December ' /)
    integer ( kind = 8 ) :: n
    integer ( kind = 8 ) :: s
    integer ( kind = 8 ) :: values(8)
    integer ( kind = 8 ) :: y

    character(len=128) :: pwd
    call get_environment_variable('PWD',pwd)

    call date_and_time ( values = values )

    y  = values(1) ! year
    m  = values(2) ! month [1-12]
    d  = values(3) ! day
    h  = values(5) ! hours
    n  = values(6) ! minutes
    s  = values(7) ! seconds
    mm = values(8) ! milliseconds

    if (idx == 0) then
      call message(-10)
      write(str, '("GPTA started on ",i2.2,1x,a,1x,i4," at ",i2,":",i2.2,":",i2.2)' ) d, trim ( month(m) ), y, h, n, s
      call message(0,str)
      call message(0,"Working directory: "//trim(pwd))
      call message(2)

      call message(0,"Number of CPUs used",i=numberOfMpiProcesses*ompNumThreads)
#ifdef GPTA_OMP
      call message(0,"Number of openMP threads used",i=ompNumThreads)
#endif
#ifdef GPTA_MPI
      call message(0,"Number of MPI processes",i=numberOfMpiProcesses)
      if (masterReadingOnly) then
        call message(0,"Master CPU reading frames only")
      else
        call message(0,"Master CPU reading and processing frames")
      end if
#endif
      call message(2)

      call dumpVariablesInfo()

      if (lastFrameOnly) then
        call message(1,"Processing only the last frame")
      else
        call message(0,"First frame to process",i=first_frame)
        if (last_frame /= huge(1)) then
          call message(0,"Last frame to process",i=last_frame)
        end if
        call message(1,"Stride for processing frames",i=stride_frame)
      end if
      
    else
      write(str, '("GPTA finished on ",i2.2,1x,a,1x,i4," at ",i2,":",i2.2,":",i2.2)' ) d, trim ( month(m) ), y, h, n, s
      call message(0,str)

    end if
    call flush(io)

  end subroutine timestamp

end program
