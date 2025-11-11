module addElementsModule

  use moduleVariables
  use moduleSystem
  use moduleMessages
  use moduleDistances

  public :: addElements

  private
  
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  integer, allocatable :: newSpecies
  integer, allocatable, dimension(:) :: numberOfNewElements
  character(cp), allocatable, dimension(:) :: newElementsLabels

  logical :: fillSphere, fillBox
  real(real64), dimension(3) :: origin
  real(real64), dimension(9) :: boxSize
  real(real64) :: sphereRadius
  real(real64) :: minimumDistance
  real(real64) :: maximumOverlap

contains

  subroutine addElements(a)
    use moduleFiles
    use moduleRead
    implicit none
    type(actionTypeDef), target :: a

    integer :: n

#ifdef DEBUG
    write(0,*) " --> Entering addElements <-- "
#endif

    call associatePointers(a)

    if (a % actionInitialisation) call initialiseAction(a)

    if (frameReadSuccessfully) then
      if (firstAction) then
        ! extend coordinates arrays
        n = sum(numberOfNewElements)
        call extendFrame(n) 

        !extend molecules arrays
        call extendMolecules(n)

        call dumpScreenInfo()

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      end if
    
      call computeAction()
    end if

  end subroutine addElements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand            => a % actionDetails
    firstAction              => a % firstAction

  end subroutine associatePointers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialiseAction(a)
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a

    character(len=STRLEN), allocatable, dimension(:) :: localString

    character(len=100) :: stringCell
    real(real64) :: hmat(3,3)

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.

    call assignFlagValue(actionCommand,"+rmin",minimumDistance,1.5_real64)

    call assignFlagValue(actionCommand,"+atom",newElementsLabels)
    
    ! call assignFlagValue(actionCommand,"+rclash",maximumOverlap)

    call assignFlagValue(actionCommand,"+n",numberOfNewElements)
    if (.not. allocated(numberOfNewElements) ) \
      call message(-1,"--add | numner of new species must be specified using the +n flag")
    if (size(newElementsLabels) /= size(numberOfNewElements)) \
      call message(-1,"--add different nr of atoms' labels and nr of atoms tto add",iv=[size(localString) , numberOfNewElements])

    newSpecies = size(newElementsLabels)

    call assignFlagValue(actionCommand,"+box ",fillBox,.false.)
    call assignFlagValue(actionCommand,"+sphere ",fillSphere,.false.)

    call assignFlagValue(actionCommand,"+origin ", origin, [0.0_real64,0.0_real64,0.0_real64])

    if (fillBox) then
      call assignFlagValue(actionCommand,"+box ", stringCell, "NONE")
      call readCellFreeFormat(stringCell, hmat)
      boxSize = reshape(hmat,[9])

    else if (fillSphere) then
      call assignFlagValue(actionCommand,"+sphere ", sphereRadius, 0.0_real64)

    else
      fillBox = .true.
      boxSize = reshape(frame % hmat,[9])
    end if

  end subroutine initialiseAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine computeAction()
    use moduleNeighbours
    implicit none

    if (fillSphere) then
      call createRandomPositionsSphere()
    else
      call createRandomPositionsBox()
    end if

    ! call setUpNeigboursList()
    ! call updateNeighboursList(.true.)

  end subroutine computeAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dumpScreenInfo()
    implicit none
    integer :: i
 
    call message(0,"Add new atomic species")
    do i=1,newSpecies
      call message(0,"...Element",str=newElementsLabels(i))
      call message(0,"...Number of new atoms",i=numberOfNewElements(i))
    end do

    call message(2)

  end subroutine dumpScreenInfo

  subroutine createRandomPositionsBox()
    use moduleRandomNumbers
    use moduleDistances
    implicit none
    integer :: iElement, imol
    integer :: addedMolecules
    real(real64), dimension(3) :: centre
    real(real64), dimension(3,3) :: hmat
    real(real64), allocatable, dimension(:) :: localSize
    integer :: initialAtoms, currentAtoms
    integer :: initialMolecules, currentMolecules
    integer :: nx, ny, nz
    integer :: ix, iy, iz

    character(len=2) :: str

    hmat = reshape(boxSize,[3,3])

    nx = int( sqrt( sum(hmat(1,1:3)**2) ) / minimumDistance )
    ny = int( sqrt( sum(hmat(2,1:3)**2) ) / minimumDistance )
    nz = int( sqrt( sum(hmat(3,1:3)**2) ) / minimumDistance )
    if ( sum(numberOfNewElements) > nx*ny*nz ) then
      call message(0,"Max number of atoms that can fit in a,b,c", iv=[nx,ny,nz])
      call message(0,"Max number of atoms that can fit in the box", i=nx*ny*nz)
      call message(-1,"--add too many atoms to fit in the box",iv=[sum(numberOfNewElements),nx*ny*nz])
    end if

    allocate(localSize(sum(numberOfNewElements)))

    addedMolecules = 0
    initialAtoms = frame % natoms
    currentAtoms = frame % natoms
    initialMolecules = numberOfMolecules
    currentMolecules = numberOfMolecules

    do iElement=1,newSpecies
      imol=0
      add : do while (imol<numberOfNewElements(iElement))
        ! generate position inside the box
        ix = int(nx * grnd()) + 1.0
        iy = int(ny * grnd()) + 1.0
        iz = int(nz * grnd()) + 1.0
        ! map to box    
        centre = \
          ix * hmat(1,1:3) / real(nx,real64) + \
          iy * hmat(2,1:3) / real(ny,real64) + \
          iz * hmat(3,1:3) / real(nz,real64)

        ! add some noise
        centre = centre + \
          [ 0.25 * minimumDistance * ( grnd() - 0.5_real64 ) , \
            0.25 * minimumDistance * ( grnd() - 0.5_real64 ) , \
            0.25 * minimumDistance * ( grnd() - 0.5_real64 ) ]
        centre = centre + origin
      
        imol = imol + 1

        addedMolecules = addedMolecules + 1

        ! add atom to molecule's list
        currentMolecules = currentMolecules + 1
        listOfMolecules(currentMolecules) % ID = numberOfUniqueMoleculeTypes + iElement
        listOfMolecules(currentMolecules) % numberOfAtoms = 1

        write(str,'(i0)')iElement
        listOfMolecules(currentMolecules) % resname = "N"//trim(str)
        allocate(listOfMolecules(currentMolecules) % listOfAtoms(1))
        listOfMolecules(currentMolecules) % listOfAtoms(1) = currentAtoms + 1
        allocate(listOfMolecules(currentMolecules) % listOfLabels(1))
        listOfMolecules(currentMolecules) % listOfLabels(1) = newElementsLabels(iElement)
        listOfMolecules(currentMolecules) % centreOfMass = centre

        ! add atom to frame
        currentAtoms = currentAtoms + 1
        frame % pos(1:3,currentAtoms) = centre
        frame % lab(currentAtoms) = newElementsLabels(iElement)
        atomToMoleculeIndex(currentAtoms) = currentMolecules

      end do add
    end do

    frame % natoms = currentAtoms
    numberOfMolecules = currentMolecules

  end subroutine createRandomPositionsBox

  subroutine createRandomPositionsSphere()
    use moduleRandomNumbers
    implicit none
    integer :: iElement, imol
    integer :: addedMolecules
    real(real64), dimension(3) :: centre
    real(real64), dimension(3,3) :: hmat
    real(real64), allocatable, dimension(:) :: localSize

    real(real64) :: theta, phi
    integer :: initialAtoms, currentAtoms
    integer :: initialMolecules, currentMolecules

    character(len=2) :: str

    allocate(localSize(sum(numberOfNewElements)))

    addedMolecules = 0
    initialAtoms = frame % natoms
    currentAtoms = frame % natoms
    initialMolecules = numberOfMolecules
    currentMolecules = numberOfMolecules

    do iElement=1,newSpecies
      imol=0
      add : do while (imol<numberOfNewElements(iElement))
        ! generate position inside a sphere

        centre = 2.0_real64 * [grnd() , grnd() , grnd()] - 1
        if ( sqrt(sum(centre**2)) > 1.0_real64 ) cycle add
        centre = centre * sphereRadius + origin
        ! theta  = grnd() * pi
        ! phi    = grnd() * twopi
        ! centre = [sin(theta) * cos(phi) , sin(theta) * sin(phi) , cos(theta)] * grnd() * sphereRadius
        ! centre = origin + centre
    
        block
          integer :: i
          real(real64) :: dist, dij(3)
          do i=initialAtoms+1,currentAtoms
            dij = frame % pos(1:3,i) - centre
            dist = sqrt(computeDistanceSquaredPBC(dij))
            if (dist < minimumDistance) cycle add
          end do
        end block
      
        imol = imol + 1

        addedMolecules = addedMolecules + 1

        ! add atom to molecule's list
        currentMolecules = currentMolecules + 1
        listOfMolecules(currentMolecules) % ID = numberOfUniqueMoleculeTypes + iElement
        listOfMolecules(currentMolecules) % numberOfAtoms = 1

        write(str,'(i0)')iElement
        listOfMolecules(currentMolecules) % resname = "N"//trim(str)
        allocate(listOfMolecules(currentMolecules) % listOfAtoms(1))
        listOfMolecules(currentMolecules) % listOfAtoms(1) = currentAtoms + 1
        allocate(listOfMolecules(currentMolecules) % listOfLabels(1))
        listOfMolecules(currentMolecules) % listOfLabels(1) = newElementsLabels(iElement)
        listOfMolecules(currentMolecules) % centreOfMass = centre

        ! add atom to frame
        currentAtoms = currentAtoms + 1
        frame % pos(1:3,currentAtoms) = centre
        frame % lab(currentAtoms) = newElementsLabels(iElement)
        atomToMoleculeIndex(currentAtoms) = currentMolecules

      end do add
    end do

    frame % natoms = currentAtoms
    numberOfMolecules = currentMolecules

  end subroutine createRandomPositionsSphere

end module addElementsModule
