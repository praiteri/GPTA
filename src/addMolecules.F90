module addMoleculesModule

  use moduleVariables
  use moduleSystem
  use moduleMessages

  public :: addMoleculesFromFile

  private
  
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  integer, pointer :: numberOfFiles
  character(len=STRLEN), pointer, dimension(:) :: fileNames
  integer, pointer, dimension(:) :: numberOfNewMolecules

  type newMolecule
    integer :: natoms
    character(cp), allocatable, dimension(:) :: lab
    real(real64), allocatable, dimension(:,:) :: pos
    real(real64), allocatable, dimension(:) :: chg
    real(real64) :: centre(3)
  end type 
  type(newMolecule), allocatable, dimension(:) :: localMolecules

  logical, pointer :: fillSphere, fillBox
  real(real64), pointer, dimension(:) :: origin
  real(real64), pointer, dimension(:) :: boxSize
  real(real64), pointer :: sphereRadius
  real(real64), pointer :: minimumDistance

contains

  subroutine addMoleculesFromFile(a)
    use moduleFiles
    use moduleRead
    implicit none
    type(actionTypeDef), target :: a

    type(fileTypeDef), allocatable, dimension(:) :: inputFile
    logical :: frameRead
    integer :: n, iFile
    real(real64) :: hmat(3,3)

    call associatePointers(a)

    if (a % actionInitialisation) call initialiseAction(a)

    if (frameReadSuccessfully) then
      ! if (firstAction) then

        allocate(localMolecules(numberOfFiles))
        allocate(inputFile(numberOfFiles))
        
        do iFile=1,numberOfFiles
          call initialiseFile(inputFile(iFile), fileNames(iFile))
        
          ! -> get number of atoms
          call getNumberOfAtoms(inputFile(iFile), n, hmat)
          rewind(inputFile(iFile) % funit)
        
          call createSystemArrays(a % localFrame, n)
          call readCoordinates(frameRead, a % localFrame, inputFile(iFile))

          ! translate input molecule's centre to the origin
          block
            integer :: i
            real(real64), dimension(3) :: dij
            real(real64) :: rgyr=0.0_real64
            dij = a % localFrame % pos(:,1)
            do i=2,a % localFrame % natoms
              dij = dij + a % localFrame % pos(:,i)
            end do
            dij = dij / dble(a % localFrame % natoms)
            do i=1,a % localFrame % natoms
              a % localFrame % pos(1:3,i) = a % localFrame % pos(1:3,i) - dij(1:3)
              rgyr = rgyr + sum(a % localFrame % pos(1:3,i)**2)
            end do
            rgyr = sqrt(rgyr / dble(a % localFrame % natoms))
          end block

          localMolecules(iFile) % natoms = a % localFrame % natoms

          allocate(localMolecules(iFile) % lab(  localMolecules(iFile) % natoms))
          allocate(localMolecules(iFile) % pos(3,localMolecules(iFile) % natoms))
          allocate(localMolecules(iFile) % chg(  localMolecules(iFile) % natoms))
          localMolecules(iFile) % lab = a % localFrame % lab
          localMolecules(iFile) % pos = a % localFrame % pos
          localMolecules(iFile) % chg = a % localFrame % chg

          call deleteSystemArrays(a % localFrame)
        end do
        
        ! extend coordinates arrays
        n = sum(localMolecules(1:numberOfFiles) % natoms * numberOfNewMolecules(1:numberOfFiles))
        call extendFrame(n) 

        !extend molecules arrays
        n = sum(numberOfNewMolecules(1:numberOfFiles))
        call extendMolecules(n)

        call dumpScreenInfo()

        call checkUsedFlags(actionCommand)
        firstAction = .false.

        ! end if
      call computeAction()
        
      deallocate(localMolecules)
      deallocate(inputFile)

    end if

  end subroutine addMoleculesFromFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand            => a % actionDetails
    firstAction              => a % firstAction

    numberOfFiles            => a % integerVariables(1)
    numberOfNewMolecules(1:) => a % integerVariables(2:)
    fileNames(1:)            => a % stringVariables(1:)

    fillBox                  => a % logicalVariables(1)
    fillSphere               => a % logicalVariables(2)

    origin(1:3)              => a % doubleVariables(1:3)
    boxSize(1:9)             => a % doubleVariables(4:12)
    sphereRadius             => a % doubleVariables(13)
    minimumDistance          => a % doubleVariables(14)

  end subroutine associatePointers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialiseAction(a)
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a

    character(len=STRLEN), allocatable, dimension(:) :: localString
    character(len=STRLEN), allocatable, dimension(:) :: localString2

    character(len=100) :: stringCell
    real(real64) :: hmat(3,3)

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.
    
    call assignFlagValue(actionCommand,"+rmin",minimumDistance,1.5_real64)

    call assignFlagValue(actionCommand,"+f",localString)
    numberOfFiles = size(localString)

    fileNames(1:numberOfFiles) = localString(1:numberOfFiles)
    deallocate(localString)

    call assignFlagValue(actionCommand,"+n",localString2)
    if (size(localString2) /= numberOfFiles) call message(-1,"--add different nr of files and nr of each molecule type",iv=[size(localString2) , numberOfFiles])

    block
      integer :: i
      do i=1,numberOfFiles
        read(localString2(i),*) numberOfNewMolecules(i)
      end do
    end block

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
    implicit none

    if (fillSphere) then
      call createRandomPositionsSphere()
    else
      call createRandomPositionsBox()
    end if

  end subroutine computeAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dumpScreenInfo()
    implicit none
    integer :: i, j
    character(len=STRLEN) :: str
 
    call message(0,"Add new molecules from file(s)")
    do i=1,numberOfFiles
      call message(0,"...File name",str=fileNames(i))
      call message(0,"......Number of molecules",i=numberOfNewMolecules(i))
      call message(0,"......Number of atoms in the molecule",i=localMolecules(i) % natoms)
      str = trim(localMolecules(i) % lab(1))
      do j=2,localMolecules(i) % natoms
        str = trim(str)//","//trim(localMolecules(i) % lab(j))
      end do
      call message(0,"......Atoms' labels",str=str)
    end do

    call message(2)

  end subroutine dumpScreenInfo

  subroutine createRandomPositionsBox
    use moduleRandomNumbers
    use moduleDistances
    implicit none
    integer :: iFile, imol
    integer :: addedMolecules
    real(real64), dimension(3) :: centre
    real(real64), dimension(3,3) :: hmat
    real(real64) :: theta, phi, vec(3)
    integer :: initialAtoms, currentAtoms
    integer :: initialMolecules, currentMolecules
    integer :: nx, ny, nz
    integer :: ix, iy, iz

    logical :: overlapFlag

    hmat = reshape(boxSize,[3,3])
    
    nx = int( sqrt( sum(hmat(1,1:3)**2) ) / minimumDistance )
    ny = int( sqrt( sum(hmat(2,1:3)**2) ) / minimumDistance )
    nz = int( sqrt( sum(hmat(3,1:3)**2) ) / minimumDistance )
    if ( sum(numberOfNewMolecules) > nx*ny*nz ) then
      call message(0,"Max number of atoms that can fit in a,b,c", iv=[nx,ny,nz])
      call message(0,"Max number of atoms that can fit in the box", i=nx*ny*nz)
      call message(-1,"--add too many atoms to fit in the box",iv=[sum(numberOfNewMolecules),nx*ny*nz])
    end if

    addedMolecules = 0
    initialAtoms = frame % natoms
    currentAtoms = frame % natoms
    initialMolecules = numberOfMolecules
    currentMolecules = numberOfMolecules

    do iFile=1,numberOfFiles
      imol=0
      add : do while (imol<numberOfNewMolecules(iFile))
        ! centre = [grnd() , grnd() , grnd()]

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
        ! centre = centre + origin
        ! localMolecules(iFile) % centre = matmul(hmat,centre) + origin
        localMolecules(iFile) % centre = centre + \
          [ 0.25 * minimumDistance * ( grnd() - 0.5_real64 ) , \
            0.25 * minimumDistance * ( grnd() - 0.5_real64 ) , \
            0.25 * minimumDistance * ( grnd() - 0.5_real64 ) ]
        
        ! random rotation around COM
        theta  = grnd() * pi
        phi    = grnd() * twopi
        vec(1) = sin(theta) * cos(phi)
        vec(2) = sin(theta) * sin(phi)
        vec(3) = cos(theta)
        theta  = grnd() * twopi  

        call rotateMolecule(vec(1:3), theta, localMolecules(iFile) % natoms, localMolecules(iFile) % pos, 0)

        ! call checkMoleculesOverlap(1+initialMolecules, currentMolecules, localMolecules(iFile), overlapFlag)
        ! if (overlapFlag) cycle add

        imol = imol + 1

        addedMolecules = addedMolecules + 1

        block ! add molecule to global array
          integer :: i
          character(len=2) :: str
          currentMolecules = currentMolecules + 1
          listOfMolecules(currentMolecules) % ID = numberOfUniqueMoleculeTypes + iFile
          listOfMolecules(currentMolecules) % numberOfAtoms = localMolecules(iFile) % natoms
          write(str,'(i0)')iFile
          listOfMolecules(currentMolecules) % resname = "N"//trim(str)

          allocate(listOfMolecules(currentMolecules) % listOfAtoms(localMolecules(iFile) % natoms))
!          allocate(listOfMolecules(currentMolecules) % listOfLabels , source=localMolecules(iFile) % lab)
          allocate(listOfMolecules(currentMolecules) % listOfLabels(localMolecules(iFile) % natoms))
          listOfMolecules(currentMolecules) % listOfLabels = localMolecules(iFile) % lab
          do i=1,localMolecules(iFile) % natoms
            listOfMolecules(currentMolecules) % listOfAtoms(i) = currentAtoms + i
          end do
          listOfMolecules(currentMolecules) % centreOfMass = localMolecules(iFile) % centre
        end block

        block ! add atoms to frame
          integer :: i, j
          j = currentAtoms
          do i=1,localMolecules(iFile) % natoms
            j = j + 1
            frame % pos(1:3,j) = localMolecules(iFile) % pos(1:3,i) + localMolecules(iFile) % centre(1:3)
          end do
          j = currentAtoms
          do i=1,localMolecules(iFile) % natoms
            j = j + 1
            frame % lab(j) = localMolecules(iFile) % lab(i)
          end do
          do i=currentAtoms+1,currentAtoms+localMolecules(iFile) % natoms
            atomToMoleculeIndex(i) = currentMolecules
          end do
          currentAtoms = currentAtoms + localMolecules(iFile) % natoms
        end block

      end do add
    end do

    frame % natoms = currentAtoms
    numberOfMolecules = currentMolecules
    
  end subroutine createRandomPositionsBox

  subroutine createRandomPositionsSphere()
    use moduleRandomNumbers
    implicit none
    integer :: iFile, imol
    integer :: addedMolecules
    real(real64), dimension(3) :: centre

    real(real64) :: theta, phi, vec(3)
    integer :: initialAtoms, currentAtoms
    integer :: initialMolecules, currentMolecules

    logical :: overlapFlag

    addedMolecules = 0
    initialAtoms = frame % natoms
    currentAtoms = frame % natoms
    initialMolecules = numberOfMolecules
    currentMolecules = numberOfMolecules

    do iFile=1,numberOfFiles
      imol=0
      add : do while (imol<numberOfNewMolecules(iFile))

        ! generate position inside a sphere
        theta  = grnd() * pi
        phi    = grnd() * twopi
        centre = [sin(theta) * cos(phi) , sin(theta) * sin(phi) , cos(theta)]
        localMolecules(iFile) % centre = origin + centre * grnd() * sphereRadius

        ! random rotation around COM
        theta  = grnd() * pi
        phi    = grnd() * twopi
        vec(1) = sin(theta) * cos(phi)
        vec(2) = sin(theta) * sin(phi)
        vec(3) = cos(theta)
        theta  = grnd() * twopi  
        call rotateMolecule(vec(1:3), theta, localMolecules(iFile) % natoms, localMolecules(iFile) % pos,-1)
      
        call checkMoleculesOverlap(1+initialMolecules, currentMolecules, localMolecules(iFile), overlapFlag)
        if (overlapFlag) cycle add

        imol = imol + 1

        addedMolecules = addedMolecules + 1

        block ! add molecule to global array
          integer :: i
          character(len=2) :: str
          currentMolecules = currentMolecules + 1
          listOfMolecules(currentMolecules) % ID = numberOfUniqueMoleculeTypes + iFile
          listOfMolecules(currentMolecules) % numberOfAtoms = localMolecules(iFile) % natoms
          write(str,'(i0)')iFile
          listOfMolecules(currentMolecules) % resname = "N"//trim(str)

          allocate(listOfMolecules(currentMolecules) % listOfAtoms(localMolecules(iFile) % natoms))
!          allocate(listOfMolecules(currentMolecules) % listOfLabels , source=localMolecules(iFile) % lab)
          allocate(listOfMolecules(currentMolecules) % listOfLabels(localMolecules(iFile) % natoms))
          listOfMolecules(currentMolecules) % listOfLabels = localMolecules(iFile) % lab
          do i=1,localMolecules(iFile) % natoms
            listOfMolecules(currentMolecules) % listOfAtoms(i) = currentAtoms + i
          end do
        end block

        block ! add atoms to frame
          integer :: i, j
          j = currentAtoms
          do i=1,localMolecules(iFile) % natoms
          j = j + 1
            frame % pos(1:3,j) = localMolecules(iFile) % pos(1:3,i) + localMolecules(iFile) % centre(1:3)
          end do
          j = currentAtoms
          do i=1,localMolecules(iFile) % natoms
            j = j + 1
            frame % lab(j) = localMolecules(iFile) % lab(i)
          end do
          do i=currentAtoms+1,currentAtoms+localMolecules(iFile) % natoms
            atomToMoleculeIndex(i) = currentMolecules
          end do
          currentAtoms = currentAtoms + localMolecules(iFile) % natoms
        end block
      
      end do add
    end do

    frame % natoms = currentAtoms
    numberOfMolecules = currentMolecules

  end subroutine createRandomPositionsSphere

  subroutine checkMoleculesOverlap(n0,n1,mol,overlap)
    use moduleDistances
    implicit none
    integer, intent(in) :: n0, n1
    type(newMolecule), intent(in) :: mol
    logical, intent(out) :: overlap
    integer :: imol, iatm, idx, jatm
    real(real64), allocatable, dimension(:,:) :: pos
    real(real64), dimension(3) :: dij
    real(real64) :: dist

    allocate(pos(3,mol % natoms))
    pos = mol % pos
    do iatm=1,mol % natoms
      pos(1:3,iatm) = pos(1:3,iatm) + mol % centre(1:3)
    end do
    
    overlap = .false.
    mainLoop : do imol=n0,n1
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        iatm = listOfMolecules(imol) % listOfAtoms(idx)
        do jatm=1,mol % natoms
          dij(1:3) = frame % pos(1:3,iatm) - pos(1:3,jatm)
          dist = sqrt(computeDistanceSquaredPBC(dij))
          if (dist < minimumDistance) then
            overlap = .true.
            exit mainLoop
          end if
        end do
      end do
    end do mainLoop

  end subroutine checkMoleculesOverlap

end module addMoleculesModule
