!disclaimer
module moduleExtractSystemProperties
  use moduleVariables
  use moduleFiles
  use moduleProperties
  use moduleMessages 
  use moduleStrings
  use moduleSystem 

  implicit none

  public :: extractSystemProperties, extractSystemPropertiesHelp
  private

  character(:), pointer :: actionCommand
  type(fileTypeDef), pointer :: outputFile
  logical, pointer :: firstAction

  logical, pointer :: dumpCoordinatesLocal
  type(fileTypeDef), pointer :: coordinatesFile

  integer, pointer :: ID

  integer, pointer :: numberOfBins
  integer, pointer :: numberOfSelectedAtoms
  integer, dimension(:), pointer :: numberOfBins2D
  real(real64), dimension(:), pointer :: distributionLimits

  logical, pointer :: atomSelection        ! Property requires atoms' selection
  logical, pointer :: dumpProperty         ! Output dump
  logical, pointer :: averageMultiProperty ! Average
  logical, pointer :: distProperty         ! Distribution
  logical, pointer :: averagePositions
  integer, pointer :: distributionType

  integer, pointer :: numberOfProperties
  real(real64), pointer, dimension(:) :: localProperty

  procedure(stuff) :: extractcell, extractHMatrix, extractVolume, extractDensity
  procedure(stuff) :: computeInertiaTensor
  procedure(stuff) :: dielectricConstant
  procedure(stuff) :: bondLength
  procedure(stuff) :: position, averageSystem, averageSystemCell
  procedure(stuff) :: test

contains


  subroutine extractSystemPropertiesHelp()
    implicit none
    call message(0,"This action extracts a property of the system, which can then be written to a file (+dump),")
    call message(0,"averaged (+avg), the histogram (+histo) or used to compute the normalised distribution (+prob).")
    call message(0,"Several properties are available:")
    call message(0,"  +coord    : Position of a few atoms")
    call message(0,"  +volume   : cell volume")
    call message(0,"  +density  : system density")
    call message(0,"  +cell     : cell lenghts and angles")
    call message(0,"  +hmat     : metric matrix")
    call message(0,"  +system   : average cell and positions (if required)")
    call message(0,"  +inertia  : inertia tensor")
    call message(0,"  +diel     : dielectric constant")
    call message(0,"  +distance : distance between two atoms")
    call message(0,"")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +coord +i 10,11")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +volume +dist +out volume.dist")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +density +avg")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +cell +out cell.out")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +hmat +avg")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +system")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +system +s Na,Cl")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +system +out averageSystem.pdb")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +inertia +dump +out inertia.dat")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +diel +out diel.dat")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --extract +distance +i 1,3")
    end subroutine extractSystemPropertiesHelp

  subroutine extractSystemProperties(a)
    implicit none
    type(actionTypeDef), target :: a

    call associatePointers(a)
    
    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      
      if (firstAction) then       

        call message(0,"System property")
        call message(0,"...Extracting "//trim(a % sysprop % name))
        if (dumpCoordinatesLocal) then
          call message(0,"......Output file for coordinates",str=coordinatesFile % fname)
        else if (outputFile % fname /= "NULL") then
          call message(0,"......Output file",str=outputFile % fname)
        end if
        call message(2)
  
        if (atomSelection) then
          if (resetFrameLabels) then
            a % updateAtomsSelection = .false.
          else
            a % updateAtomsSelection = .true.
          end if
         
          call selectAtoms(1,actionCommand,a)
          call createSelectionList(a,1)
          numberOfSelectedAtoms = count(a % isSelected(:,1))
        end if
        
        ! get the number of properties computed and initialise the working arrays
        call initialiseCompute(a)
 
        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else

        if (atomSelection) then
          if (a % updateAtomsSelection) then
            call selectAtoms(1,actionCommand,a)
            call createSelectionList(a,1)
            numberOfSelectedAtoms = count(a % isSelected(:,1))
          end if
        end if

      end if
      localProperty => a % array1D       

      call computeAction(numberOfProperties, a)
    end if

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) then

      if (.not. dumpCoordinatesLocal) then
        if (averageMultiProperty) call message(0,"......Average(s) and standard deviation(s)")
        if (distProperty        ) then
          if (distributionType == 1) then 
            call message(0,"......Computing histogram")
          else if (distributionType == 2) then 
            call message(0,"......Computing probability")
          else if (distributionType == 3) then 
            call message(0,"......Computing scaled distributino")
          else if (distributionType == 4) then 
            call message(0,"......Kernel based histogram")
          end if
        end if
      end if
      call message(2)

      call finaliseAction(a)
    end if

  end subroutine extractSystemProperties

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand           => a % actionDetails
    firstAction             => a % firstAction
    outputFile              => a % outputFile
    
    numberOfBins            => a % numberOfBins
    numberOfBins2D(1:2)     => a % numberOfBins2D
    
    ID                      => a % integerVariables(1)
    numberOfProperties      => a % integerVariables(2)
    numberOfSelectedAtoms   => a % integerVariables(3)
    distributionType        => a % integerVariables(4)

    distributionLimits(1:2) => a % doubleVariables(1:2)
    ! localProperty(1:)       => a % doubleVariables(3:)

    atomSelection           => a % logicalVariables(1) ! Atoms selection for property
    dumpProperty            => a % logicalVariables(2) ! Output dump
    averageMultiProperty    => a % logicalVariables(3) ! Average more than one property
    distProperty            => a % logicalVariables(4) ! Distribution

    dumpCoordinatesLocal    => a % logicalVariables(5) 
    averagePositions        => a % logicalVariables(6) 
    coordinatesFile         => a % auxiliaryFile 

  end subroutine associatePointers

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: i
    logical, dimension(4) :: lflag = .false.

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.

    numberOfSelectedAtoms = 0
    atomSelection = .false. 
    distProperty = .false.
    
    call getPropertyType(a)

    ! Extract to do with the properties
    if (.not. averageMultiProperty) then
      call assignFlagValue(actionCommand,"+avg " ,averageMultiProperty, .false.) ! Average
      call assignFlagValue(actionCommand,"+dump ",dumpProperty,         .false.) ! Dump on same line

      call assignFlagValue(actionCommand,"+histo "  ,lflag(1), .false.) ! Histogram - non-normalised
      call assignFlagValue(actionCommand,"+prob "  ,lflag(2), .false.) ! Probability distribution - normalised to 1
      call assignFlagValue(actionCommand,"+kernel",lflag(4), .false.) ! Probability with Gaussian Kernel
      if (count(lflag) == 1) then
        distProperty = .true.
      else if (count(lflag) > 1) then
        call message(-1,"Only one distribution can be computed") 
      end if
      distributionType = 0
      do i=1,4
        if (lflag(i)) distributionType = i
      end do
    end if

    if ( count([distProperty,averageMultiProperty]) == 0 ) dumpProperty = .true.

    call assignFlagValue(actionCommand,"+out",outputFile % fname,"properties.out")
    if (.not. dumpCoordinatesLocal) then
      call initialiseFile(outputFile, outputFile % fname)
      write(outputFile % funit,"(a)")"#System property"
      write(outputFile % funit,'(a)')"#...Extracting "//trim(a % sysprop % name)
      if (averageMultiProperty) write(outputFile % funit,'(a)')"#......Average and standard deviation"
      if (distProperty        ) write(outputFile % funit,'(a)')"#......Distribution"
    end if

  end subroutine initialiseAction

  subroutine finaliseAction(a)
    use moduleDistances
    use moduleFiles
    use moduleMessages
    
    implicit none
    type(actionTypeDef), target :: a
    real(real64), allocatable, dimension(:) :: avgValues
    real(real64) :: rtmp
    character(len=30) :: str

    if (dumpCoordinatesLocal) then
      call workData % extract(ID, avgValues)
      frame % hmat = reshape(avgValues(1:9),[3,3])
      call hmat2cell (frame % hmat, frame % cell, "DEG")

      if (averagePositions) then
        block 
          integer :: iatm, idx
          idx = 9
          do iatm=1,frame % natoms
            if (a % isSelected(iatm,1)) then
              frame % frac(1:3,iatm) = avgValues(idx+1:idx+3)
              idx = idx + 3
            end if
          end do
        end block
      end if

      ! scale the coordinated to the average cell
      call fractionalToCartesian(frame % natoms, frame % frac, frame % pos)
      if (numberOfMolecules > 0) call reassembleAllMolecules()

      call runInternalAction("dump",trim(coordinatesFile % fname)//" +internal")

      close (coordinatesFile % funit)
      return
    end if

    if (distProperty) then
      select case (distributionType)
        case (1)
          str = 'none'
        case (2)
          str = 'probability'
        case (4)
          call assignFlagValue(actionCommand,"+sigma ",rtmp,-1.0_real64) ! Probability with Gaussian Kernel
          if (rtmp < 0) call message(-1,"Kernel sigma must be greater than zero")
          write(str,'("kernel=",e20.15)') rtmp
      end select
      call workData % dump(ID, outputFile % funit, normalisation=str)

    else
      call workData % dump(ID, outputFile % funit)
    end if
    call message(2)

  end subroutine finaliseAction

  subroutine computeAction(n, a)
    use moduleSystem
    implicit none
    integer, intent(in) :: n
    type(actionTypeDef), target :: a

    integer :: ntmp

    ntmp = n
    if (atomSelection) then
      call a % sysprop % property(numberOfProperties, localProperty, numberOfSelectedAtoms, a % idxSelection(:,1))
    else
      call a % sysprop % property(numberOfProperties, localProperty)
    end if

    if (ntmp>0) then
      call workData % compute(ID, numberOfValues=numberOfProperties, xValues=localProperty)
    else
      allocate(a % array1D(numberOfProperties), source=0.0_real64)
    end if

  end subroutine computeAction

  subroutine initialiseCompute(a)
    implicit none
    type(actionTypeDef), target :: a
    character(len=30) :: str
    numberOfProperties = -1
    call computeAction(numberOfProperties, a)

    if (dumpProperty) then 
      call workData % initialise(ID, "dump", iounit=outputFile % funit)
      
    else if (averageMultiProperty) then
      call workData % initialise(ID, "multiAverage", numberOfValues=numberOfProperties)

    else if (distProperty) then
      call assignFlagValue(actionCommand,"+nbin ",numberOfBins,100)

      if (distributionType == 1) then
        str="hist "
      else if (distributionType == 2) then
        str="prob "
      else if (distributionType == 4) then
        str="kernel "
      end if

      call assignFlagValue(actionCommand,"+"//trim(str), distributionLimits(1:2),[1.0_real64,2.0_real64])
      call workData % initialise(ID, "histogram", numberOfBins=[numberOfBins], lowerLimits=[distributionLimits(1)], upperLimits=[distributionLimits(2)])
    end if

  end subroutine initialiseCompute

  subroutine getPropertyType(a)
    implicit none
    type(actionTypeDef), target :: a
    integer :: n
    character(len=2) :: str
         
    n = 0
    if (flagExists(actionCommand,"+test")) then
      n = n + 1
      atomSelection = .true.
      a % sysprop % property => test
      a % sysprop % name = "Test Extract Routine"
    end if
         
    if (flagExists(actionCommand,"+coord")) then
      n = n + 1
      atomSelection = .true.
      a % sysprop % property => position
      a % sysprop % name = "Atoms Positions"
    end if

    if (flagExists(actionCommand,"+cell")) then
      n = n + 1
      a % sysprop % property => extractCell
      a % sysprop % name = "Cell Parameters"
    end if
    
    if (flagExists(actionCommand,"+hmat")) then
      n = n + 1
      a % sysprop % property => extractHMatrix
      a % sysprop % name = "Cell Matrix"
    end if
    
    if (flagExists(actionCommand,"+system")) then
      n = n + 1
      averageMultiProperty = .true.
      dumpCoordinatesLocal = .true.

      call assignFlagValue(actionCommand,"+out",coordinatesFile % fname,"averageSystem.pdb")
      call assignFlagValue(actionCommand,"+s ",averagePositions,.false.)
      if (averagePositions) then
        atomSelection = .true.
        a % sysprop % property => averageSystem
        a % sysprop % name = "Average cell and average coordinates"
      else
        atomSelection = .false.
        a % sysprop % property => averageSystemCell
        a % sysprop % name = "Average cell and scaled last coordinates"
      end if

    end if
    
    if (flagExists(actionCommand,"+volume")) then
      n = n + 1
      a % sysprop % property => extractVolume
      a % sysprop % name = "Cell Volume"
    end if
    
    if (flagExists(actionCommand,"+inertia")) then
      atomSelection = .true.
      n = n + 1
      a % sysprop % property => computeInertiaTensor
      a % sysprop % name = "Intertia tensor (eigenvalues)"
    end if

    if (flagExists(actionCommand,"+diel")) then
      n = n + 1
      a % sysprop % property => dielectricConstant
      a % sysprop % name = "Dielectric constant"
    end if
    
    if (flagExists(actionCommand,"+density")) then
      n = n + 1
      a % sysprop % property => extractDensity
      a % sysprop % name = "Density"
    end if

    if (flagExists(actionCommand,"+distance")) then
      n = n + 1
      atomSelection = .true.
      a % sysprop % property => bondLength
      a % sysprop % name = "Distance between atoms"
    end if

    if (n==0)  call message(-1,"--extract | no property selected")
    write(str,'(i0)') n
    if (n/=1)  call message(-1,"--extract | only one property per action can be specified "//str)
    
  end subroutine getPropertyType

end module moduleExtractSystemProperties

subroutine extractCell(n, cell, numberOfSelectedAtoms, selectionList)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: cell
  integer, intent(in), optional :: numberOfSelectedAtoms
  integer, dimension(:), intent(in), optional :: selectionList

  ! set the number of properties computed
  if (n<0) then
    n = 6
    return
  end if
  cell(1:6) = frame % cell

end subroutine extractCell

subroutine extractHMatrix(n, hmat, numberOfSelectedAtoms, selectionList)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: hmat
  integer, intent(in), optional :: numberOfSelectedAtoms
  integer, dimension(:), intent(in), optional :: selectionList

  ! set the number of properties computed
  if (n<0) then
    n = 9
    return
  end if
  hmat(1:9) = reshape(frame % hmat, [9])

end subroutine extractHMatrix

subroutine extractVolume(n, volume, numberOfSelectedAtoms, selectionList)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: volume
  integer, intent(in), optional :: numberOfSelectedAtoms
  integer, dimension(:), intent(in), optional :: selectionList

  ! set the number of properties computed
  if (n<0) then
    n = 1
    return
  end if
  volume(1) = frame % volume

end subroutine extractVolume

subroutine extractDensity(n, density, numberOfSelectedAtoms, selectionList)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: density
  integer, intent(in), optional :: numberOfSelectedAtoms
  integer, dimension(:), intent(in), optional :: selectionList

  ! set the number of properties computed
  if (n<0) then
    n = 1
    return
  end if
  density(1) = 1.6605388_real64 * totalMass / frame % volume

end subroutine extractDensity

subroutine dielectricConstant(n, eps, numberOfSelectedAtoms, selectionList)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: eps
  integer, intent(in), optional :: numberOfSelectedAtoms
  integer, dimension(:), intent(in), optional :: selectionList

  real(real64), parameter :: e0 = 0.00552642629948221_real64 ! e^2/eV/angstrom
  real(real64), parameter :: kbt = 0.025851_real64 ! in eV
  
  integer, save :: np = 0
  real(real64), dimension(3), save :: m_avg = 0.0_real64
  real(real64), save :: m2_avg = 0.0_real64
  real(real64), dimension(3) :: cellDipole
  real(real64) :: norm, a3(3), a, b
  
  integer :: iatm
  
  ! set the number of properties computed
  if (n<0) then
    n = 1
    return
  end if
  
  cellDipole = 0.0_real64
  do iatm=1,frame % natoms
    cellDipole = cellDipole + frame % chg(iatm) * frame % pos(:,iatm) 
  end do
  ! if wanted in Debye ! cellDipole = cellDipole * 4.8032513_real64    
  
  m_avg = m_avg + cellDipole
  m2_avg = m2_avg + sum(cellDipole**2)
  
  np = np + 1
  
  norm = 3.0_real64* e0 * kbt * frame % volume
  a3 = m_avg / np
  a = sum(a3*a3)
  b = m2_avg / np
  eps(1) = 1.0_real64 + (b-a) / norm
  eps(2) = 1.0_real64 + (m2_avg) / np / norm
  
end subroutine dielectricConstant

subroutine bondLength(n, val, nat, l)
  use moduleSystem
  use moduleDistances
  use moduleMessages
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: val
  integer, intent(in), optional :: nat
  integer, dimension(:), intent(in), optional :: l
  
  real(real64), dimension(3) :: dij
  
  ! set the number of properties computed
  if (n<0) then
    n = 1
    return
  end if
  
  if (nat/=2) call message(0,"bondLength - wrong number of atoms selected")
  dij = frame % pos(:,l(1)) - frame % pos(:,l(2))
  val(1) = sqrt(computeDistanceSquaredPBC(dij))
  
end subroutine bondLength
  
subroutine position(n, val, nat, l)
  use moduleSystem
  use moduleDistances
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: val
  integer, intent(in), optional :: nat
  integer, dimension(:), intent(in), optional :: l
  
  integer :: i, ntmp

  ! set the number of properties computed
  if (n < 0) then
    n = nat *3
    return
  end if
  
  ntmp = 0
  do i=1,nat
    val(ntmp+1:ntmp+3) = frame % pos(:,l(i))
    ntmp = ntmp + 3
  end do
  
end subroutine position

subroutine averageSystem(n, val, nat, l)
  use moduleSystem
  use moduleDistances
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: val
  integer, intent(in), optional :: nat
  integer, dimension(:), intent(in), optional :: l
  
  integer :: i, ntmp
  real(real64), allocatable, dimension(:,:) :: localPositions

  ! set the number of properties computed
  if (n < 0) then
    n = 9 + nat *3
    return
  end if
  ntmp = 0

  ! call fixCellJumps()

  do i=1,3
    val(ntmp+1:ntmp+3) = frame % hmat(:,i)
    ntmp = ntmp + 3
  end do

!  if (n>9) then
    allocate(localPositions(3,frame % natoms))
    call cartesianToFractionalNoWrap(frame % natoms, frame % pos, localPositions)
    do i=1,nat
      val(ntmp+1:ntmp+3) = localPositions(:,l(i))
      ntmp = ntmp + 3
    end do
!  end if
  
end subroutine averageSystem

subroutine averageSystemCell(n, val, numberOfSelectedAtoms, selectionList)
  use moduleSystem
  use moduleDistances
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: val
  integer, intent(in), optional :: numberOfSelectedAtoms
  integer, dimension(:), intent(in), optional :: selectionList

  integer :: i, ntmp

  ! set the number of properties computed
  if (n < 0) then
    n = 9
    return
  end if
  ntmp = 0

  ! call fixCellJumps()
  do i=1,3
    val(ntmp+1:ntmp+3) = frame % hmat(:,i)
    ntmp = ntmp + 3
  end do
  
end subroutine averageSystemCell

subroutine test(n, val, nat, l)
  use moduleSystem
  use moduleDistances
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: val
  integer, intent(in), optional :: nat
  integer, dimension(:), intent(in), optional :: l
  
  ! set the number of properties computed
  if (n < 0) then
    n = 3
    return
  end if
  
  val(1) = 1.0_real64
  val(2) = 2.0_real64
  val(3) = 3.0_real64
  
end subroutine
