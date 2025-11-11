!disclaimer
module moduleTemplateAction
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleFiles
  use moduleProperties
  use moduleMessages
  use moduleDistances
  
  implicit none

  public :: testAction,testActionHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions
  
  integer, pointer :: ID
  integer, pointer :: numberOfBins
  real(8), pointer :: cutoffRadius
  integer, pointer :: maxNeigh

contains

subroutine testActionHelp()
  use moduleMessages
  implicit none
  call message(0,"This is a template action for developing new features.")
  call message(0,"It action computes the Radial Pair Distribution function")
  call message(0,"Examples:")
  call message(0,"  gpta.x --i coord.pdb traj.dcd --gofr +s Ca O1,O2 +out gofr.out +nbin 100")
end subroutine testActionHelp

subroutine testAction(a)
  implicit none
  type(actionTypeDef), target :: a

  if (a % actionInitialisation) then
    call initialiseAction(a)
    return
  end if

  if (frameReadSuccessfully) then
    tallyExecutions = tallyExecutions + 1

    ! Atoms' selection
    if (firstAction) then
      ! dump info about the action on the screen
      call dumpScreenInfo()

      ! select two groups of atoms
      call selectAtoms(2,actionCommand,a)

      ! update selection if "reactive" trajectory
      ! i.e. if the atom may change name 
      if (resetFrameLabels) then
        a % updateAtomsSelection = .false.
      else
        a % updateAtomsSelection = .true.
      end if

      ! create a list of the atoms' indices for each group
      call createSelectionList(a,2)

      ! Throw a warning for unused flags
      call checkUsedFlags(actionCommand)
      firstAction = .false.

    else
      
      ! Repeat selection for reactive trajectories
      if (a % updateAtomsSelection) then
        ! select two groups of atoms
        call selectAtoms(2,actionCommand,a)
        ! create a list of the atoms' indices for each group
        call createSelectionList(a,2)
      end if 
        
    end if

    call computeAction(a)
  end if

  if (endOfCoordinatesFiles) then
    call finaliseAction()
  end if 

end subroutine testAction

  subroutine initialiseAction(a)

    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile

    ID                   => a % integerVariables(1)
    numberOfBins         => a % numberOfBins    
    cutoffRadius         => a % doubleVariables(1)

    maxNeigh             => a % integerVariables(2)

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .true.
    a % requiresNeighboursListDouble = .true.
    a % cutoffNeighboursList = 3.2d0

    ! get output file name from the command line, if present
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'test.out')
    
    ! get cutoff radius from the command line, if present
    call assignFlagValue(actionCommand,"+rcut",cutoffRadius,6.d0)
    
    ! get number of bins for the distribution from the command line, if present
    call assignFlagValue(actionCommand,"+nbin",numberOfBins,100)
    
    ! Initiliase container for the 1D distribution
    call workData % initialise(ID, "probability", numberOfBins=[numberOfBins], lowerLimits=[0.d0], upperLimits=[cutoffRadius])
    
    maxNeigh = 0
    tallyExecutions = 0

  end subroutine initialiseAction

  subroutine finaliseAction()
    implicit none
!    real(8), allocatable, dimension(:) :: localDistribution
!    integer :: ntmp

    call initialiseFile(outputFile,outputFile % fname)
    
    ! thi routine will write the average number of neighbours vs distance
    call workData % dump(ID, outputFile % funit, normalisation='none')

    ! this routine extracts the distribution, which could the be used to normalise the g(r) properly
    ! call workData % extract(ID, localDistribution, numberOfCounts=ntmp)

    close(outputFile % funit)
  end subroutine finaliseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Test interface")
    call message(0,"...Output file",str=outputFile % fname)
    call message(0,"...Cutoff radius",r=cutoffRadius)
    call message(0,"...Number of bins",i=numberOfBins)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: idx, jdx
    integer :: iatm, jatm
    integer :: nsel1, nsel2, ntmp
    real(8) :: dij(3), dist2, rcut2

    real(8), allocatable, dimension(:) :: allDistances

    rcut2 = cutoffRadius**2

    nsel1 = count(a % isSelected(:,1))
    nsel2 = count(a % isSelected(:,2))

    ! this array is probably too large and it makes the calculation inefficient
    if (maxNeigh == 0) then
      maxNeigh = 100*(nsel1+nsel2)
    else 
      maxNeigh = int(maxNeigh * 1.2)
    end if
    allocate(allDistances(maxNeigh))

    ntmp = 0
    ! loop over over the first group of atoms
    do idx=1,nsel1
      ! loop over over the second group of atoms
      iatm = a % idxSelection(idx,1)
      do jdx=1,nsel2
        jatm = a % idxSelection(jdx,2)

        if (iatm == jatm) cycle

        ! calculation of the distance with PBC
        dij = frame % pos(:,iatm) - frame % pos(:,jatm)
        dist2 = computeDistanceSquaredPBC(dij)

        ! store the distance if shorter thant the cutoff radius
        if (dist2 >= rcut2) cycle
        ntmp = ntmp + 1
        allDistances(ntmp) = sqrt(dist2)

      end do
    end do
    maxNeigh = ntmp

    ! add all distances to the calculation of the distribution
    call workData % compute(ID, numberOfValues=maxNeigh, xValues=allDistances)

  end subroutine computeAction

end module moduleTemplateAction
