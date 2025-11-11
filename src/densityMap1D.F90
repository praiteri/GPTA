!disclaimer
module moduleDensityProfile
  use moduleVariables
  use moduleFiles
  use moduleMessages
  use moduleProperties

  implicit none

  public :: computeDensityProfile, computeDensityProfileHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions

  integer, pointer :: numberOfBins
  integer, pointer :: iAxis
  integer, pointer, dimension(:) :: cAxes
  integer, pointer :: ID

  real(real64), pointer :: averageSize, averageVolume
  real(real64), pointer, dimension(:) :: centre
  real(real64), pointer :: radiusSquared

contains

  subroutine computeDensityProfileHelp()
    implicit none
    call message(0,"This action computes the 1D density profile along one of the crystallographic directions.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --dmap1D +x +s O")
    call message(0,"  gpta.x --i coord.pdb --dmap1D +z +s O +out dmap.out")
    call message(0,"  gpta.x --i coord.pdb --dmap1D +z +s Ca +nbin 200 +out dmap.out")
    call message(0,"  gpta.x --i coord.pdb --dmap1D +z +i 1:300:3 +out dmap.out")
    call message(0,"  gpta.x --i coord.pdb --dmap1D +r 10,10,10 +s O")
    call message(0,"  gpta --i coord.pdb trajectory.0.dcd --dmap1D +cz 75,75,50 +s Kr")
  end subroutine computeDensityProfileHelp

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile

    numberOfBins         => a % numberOfBins
    ID                   => a % integerVariables(1)
    iAxis                => a % integerVariables(2)
    cAxes                => a % integerVariables(3:4)

    averageVolume        => a % doubleVariables(1)
    averageSize          => a % doubleVariables(2)
    centre(1:3)          => a % doubleVariables(3:5)
    radiusSquared        => a % doubleVariables(6)

  end subroutine 

  subroutine initialiseAction(a)
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a

    integer :: i
    logical :: lflag(7)
    real(real64) :: rtmp(4)

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0_real64

    call assignFlagValue(actionCommand,"+out",outputFile % fname,'dmap1D.out')

    call assignFlagValue(actionCommand,"+nbin",numberOfBins,100)
    allocate(a % array1D(0:numberOfBins+1) , source=0.0_real64)

    call assignFlagValue(actionCommand,"+x",lflag(1),.false.)
    call assignFlagValue(actionCommand,"+y",lflag(2),.false.)
    call assignFlagValue(actionCommand,"+z",lflag(3),.false.)
    call assignFlagValue(actionCommand,"+r",lflag(4),.false.)
    call assignFlagValue(actionCommand,"+cx",lflag(5),.false.)
    call assignFlagValue(actionCommand,"+cy",lflag(6),.false.)
    call assignFlagValue(actionCommand,"+cz",lflag(7),.false.)

    if (count(lflag)==0) call message(-1,"--dmap1D no direction specified +x/+y/+z/+r/+cx/+cy/+cz")
    if (count(lflag)>1) call message(-1,"--dmap1D more than one direction specified +x/+y/+z/+r/+cx/+cy/+cz")

    do i=1,size(lflag)
      if (lflag(i)) iAxis = i
    enddo

    if (lflag(4)) then
      call assignFlagValue(actionCommand,"+r",rtmp,[0.0_real64,0.0_real64,0.0_real64,1.0_real64])
      centre = rtmp(1:3)
      radiusSquared = rtmp(4)**2
    end if

    if (lflag(5)) then
      call assignFlagValue(actionCommand,"+cx",rtmp(1:3),[0.0_real64,0.0_real64,0.0_real64])
      cAxes = [2,3]
      centre(1:2) = rtmp(1:2)
      radiusSquared = rtmp(3)**2
    end if

    if (lflag(6)) then
      call assignFlagValue(actionCommand,"+cy",rtmp(1:3),[0.0_real64,0.0_real64,0.0_real64])
      cAxes = [1,3]
      centre(1:2) = rtmp(1:2)
      radiusSquared = rtmp(3)**2
    end if

    if (lflag(7)) then
      call assignFlagValue(actionCommand,"+cz",rtmp(1:3),[0.0_real64,0.0_real64,0.0_real64])
      cAxes = [1,2]
      centre(1:2) = rtmp(1:2)
      radiusSquared = rtmp(3)**2
    end if

    tallyExecutions = 0 
    averageSize = 0.0_real64
    averageVolume = 0.0_real64

    call workData % initialise(ID, "histogram", numberOfBins=[numberOfBins], lowerLimits=[0.0_real64], upperLimits=[1.0_real64])

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    use moduleMessages 
    implicit none
    call message(0,"Computing 1D density profile")
    call message(0,"...Output file",str=outputFile % fname)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleVariables
    use moduleSystem 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a

    integer :: iatm, nlocal
    real(real64), allocatable, dimension(:), save :: localCoord

    real(real64) :: dij(3), dist2, height
    integer :: axes(2)
    integer :: i, n
    
    if (iAxis == 4) then
      ! Local array with the coordinates of only the selected atoms
      nlocal = count(a % isSelected(:,1))
      allocate(localCoord(nlocal))
      nlocal = 0

      do iatm=1,frame % natoms
        if (.not. a % isSelected(iatm,1)) cycle
        dij = frame % pos(:,iatm) - centre
        dist2 = computeDistanceSquaredPBC(dij) 
        if (dist2 < radiusSquared) then
          nlocal = nlocal + 1
          localCoord(nlocal) = sqrt(dist2 / radiusSquared)
        end if
      enddo
      call workData % compute(ID, numberOfValues=nlocal, xValues=localCoord)
      
      deallocate(localCoord)
  
    else if (iAxis == 5 .or. iAxis == 6 .or. iAxis == 7) then
      averageSize = averageSize + frame % hmat(iAxis-4,iAxis-4)
      averageVolume = averageVolume + frame % volume

      ! Local array with the coordinates of only the selected atoms
      nlocal = count(a % isSelected(:,1))
      allocate(localCoord(nlocal))
      nlocal = 0

      dij = 0.d0
      do iatm=1,frame % natoms
        if (.not. a % isSelected(iatm,1)) cycle
        do i=1,2
          dij(i) = frame % pos(cAxes(i),iatm) - centre(i)
        end do
        dist2 = computeDistanceSquaredPBC(dij) 
        if (dist2 < radiusSquared) then
          nlocal = nlocal + 1
          localCoord(nlocal) = sqrt(dist2 / radiusSquared)
        end if
      enddo
      call workData % compute(ID, numberOfValues=nlocal, xValues=localCoord)
      
      deallocate(localCoord)

    else
      averageSize = averageSize + frame % hmat(iAxis,iAxis)
      averageVolume = averageVolume + frame % volume

      ! Local array with the coordinates of only the selected atoms
      nlocal = count(a % isSelected(:,1))
      allocate(localCoord(nlocal))
      nlocal = 0
      do iatm=1,frame % natoms
        if (.not. a % isSelected(iatm,1)) cycle
        nlocal = nlocal + 1
        localCoord(nlocal) = frame % frac(iAxis,iatm)
      enddo
      call workData % compute(ID, numberOfValues=nlocal, xValues=localCoord)
      deallocate(localCoord)
    end if

  end subroutine computeAction

  subroutine finaliseAction()
    implicit none
    real(real64) :: dVolume, upperLimit
    character(len=30) :: str

    call initialiseFile(outputFile,outputFile % fname)
    write(outputFile % funit,'("# Distance  | Density [atom/angstrom^3]")')

    if (iAxis == 4) then
      write(str,'("sphere=",e20.15)') dble(tallyExecutions)
      upperLimit = sqrt(radiusSquared)

    else if (iAxis == 5 .or. iAxis == 6 .or. iAxis == 7) then
      averageSize = averageSize / tallyExecutions
      write(str,'("cylinder=",e20.15)') dble(tallyExecutions * averageSize)
      upperLimit = sqrt(radiusSquared)

    else
      averageVolume = averageVolume / tallyExecutions
      dVolume =  averageVolume / dble(numberOfBins)
      write(str,'("custom=",e20.15)') dVolume * tallyExecutions 
      upperLimit = averageSize / tallyExecutions

    end if

    call workData % dump(ID, outputFile % funit, upperLimits=[upperLimit], normalisation=str)
    close(outputFile % funit)

  end subroutine finaliseAction

  subroutine computeDensityProfile(a)
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

        call selectAtoms(1,actionCommand,a)
      
        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else

        if (a % updateAtomsSelection) call selectAtoms(1,actionCommand,a)

      end if
      
      call computeAction(a)

    end if

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) then
      call finaliseAction()
    end if

  end subroutine computeDensityProfile

end module moduleDensityProfile
