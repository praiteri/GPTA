!disclaimer
module moduleMeanSquareDisplacement
  use moduleVariables
  use moduleFiles
  use moduleProperties
  use moduleMessages
  use moduleSystem
  use moduleStrings

  implicit none

  public :: computeMSD, computeMSDHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions

  integer, pointer :: numberOfBins
  integer, pointer :: ID
  real(real64), pointer :: timeInterval

contains

  subroutine computeMSDHelp()
    implicit none
    call message(0,"This action computes the Mean Square Displacement (MSD) of the selected species.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --msd +s Ca +out msd.out")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --msd +s O1,O2 +nt 1000")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --msd +mol M2 +dt 0.5")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --msd +i 1:300:3")
  end subroutine computeMSDHelp

  subroutine computeMSD(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: nsel, i

    call associatePointers(a)
    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if
    
    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1
      
      if (firstAction) then
        call dumpScreenInfo()

        if (resetFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if

        call selectAtoms(1,actionCommand,a)
        nsel = count(a % isSelected(:,1))
        allocate(a % localPositions(3,nsel))
        call createSelectionList(a,1)
        
        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else

        if (a % updateAtomsSelection) then
          call selectAtoms(1,actionCommand,a)
          call createSelectionList(a,1)
        end if

      end if

      if (tallyExecutions > numberOfBins+1) tallyExecutions = 1
      if (tallyExecutions == 1) then
        nsel = count(a % isSelected(:,1))
        call createReferenceCoordinatesFromList(nsel, a % idxSelection(:,1), a % localPositions)
      else
        if (a % updateAtomsSelection) then
          call computeActionReax(a)
        else
          call computeAction(a)
        end if
      end if

    end if

    if (endOfCoordinatesFiles) then
      call finaliseAction()
    end if 
  end subroutine computeMSD


  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile

    numberOfBins         => a % numberOfBins    
    ID                   => a % integerVariables(1)
    timeInterval         => a % doubleVariables(1)
    
  end subroutine associatePointers

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    a % actionInitialisation = .false.

    ! Local pointers

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.

    call assignFlagValue(actionCommand,"+out",outputFile % fname,'msd.out')
    
    call assignFlagValue(actionCommand,"+nt",numberOfBins,1000)
    call assignFlagValue(actionCommand,"+dt",timeInterval,1.0_real64)
    
    call workData % initialise(ID, "avgDist ", numberOfBins=[numberOfBins], lowerLimits=[0.0_real64], upperLimits=[dble(numberOfBins)])

    tallyExecutions = 0
  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Computing mean square displacement")
    call message(0,"...Maximum number of frames",i=numberOfBins)
    call message(0,"...Output file",str=outputFile % fname)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a

    integer :: i, nsel, iatm
    real(real64) :: msd, dij(3)
    real(real64) :: time

    nsel = count(a % isSelected(:,1))
    msd = 0.0_real64
    do i=1,nsel
      iatm = a % idxSelection(i,1)
      dij = frame % pos(1:3,iatm) - a % localPositions(1:3,i)
      msd = msd + sum(dij*dij)
    end do
    msd = msd / dble(nsel)

    time = tallyExecutions - 1.01
    call workData % compute(ID, numberOfValues=1, xValues=[msd], yValues=[time])

  end subroutine computeAction

  subroutine computeActionReax(a)
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    logical, save :: firstTimeIn = .true.

    integer :: i, nsel, iatm
    real(real64) :: msd, dij(3), plocal(3)
    real(real64), save :: psave(3)
    real(real64) :: time, rtmp

    nsel = count(a % isSelected(:,1))

    if (firstTimeIn) then
      firstTimeIn = .false.
      iatm = a % idxSelection(i,1)
      rtmp = computeDistanceSquaredPBC(dij)
      psave = frame % pos(1:3,iatm)
    end if    
    
    msd = 0.0_real64
    do i=1,nsel
      iatm = a % idxSelection(i,1)

      dij = frame % pos(1:3,iatm) - psave
      plocal = psave + dij

      dij = plocal - a % localPositions(1:3,i)
      msd = msd + sum(dij*dij)
      psave = plocal
    end do
    msd = msd / dble(nsel)

    time = tallyExecutions - 1.01
    call workData % compute(ID, numberOfValues=1, xValues=[msd], yValues=[time])

  end subroutine computeActionReax

  subroutine finaliseAction()
    implicit none
    real(real64), allocatable, dimension(:) :: msd
    integer, allocatable, dimension(:) :: tally
    real(real64), allocatable, dimension(:) :: time
    integer :: i, idx, nbin, n
    real(real64) :: alpha, beta, rfact

    ! call workData % extract(ID, msd)

    call workData % extract(ID, values=msd, tally=tally)

    if (me/=0) return

    nbin = size(tally)
    do i=1,nbin
      if (tally(i) == 0) then
        nbin = i-1
        exit
      end if
      msd(i) = msd(i) / tally(i)
    end do

    allocate(time(nbin))
    do i=1,nbin
      time(i) = i * timeInterval
    end do

    time = time * numberOfActiveProcesses

    idx = int(0.2*nbin)
    n = nbin - idx + 1

    ! y = alpha + beta*x
    call linearRegression(n,time(idx:nbin),msd(idx:nbin),alpha,beta,rfact)

    beta = beta / 6.0_real64
    beta = beta * 10.0_real64 ! conversion from angstrom^2/ps to 10^-5 cm^2/s

    call message(2)
    call message(0,"Self diffusion coefficient calculation")
    call message(0,"  Average D0 [10^-5 cm^2/s]",r=beta)
    call message(1,"  Average correlation factor",r=rfact)
    call initialiseFile(outputFile,outputFile % fname)
    write(outputFile % funit,'("# Time [ps] | MSD [A^3]")')
    do i=1,nbin
      write(outputFile % funit,'(f7.3,1x,f12.5)') time(i), msd(i)
    enddo
    close(outputFile % funit)

  end subroutine finaliseAction

  subroutine createReferenceCoordinatesFromList(natoms,isel,pos)
    implicit none
    integer, intent(in) :: natoms
    integer, intent(in), dimension(natoms) :: isel
    real(real64), intent(out), dimension(:,:) :: pos

    integer :: iatm

    do iatm=1,natoms
      pos(1:3,iatm) = frame % pos(1:3,isel(iatm))
    end do

  end subroutine createReferenceCoordinatesFromList

end module moduleMeanSquareDisplacement
