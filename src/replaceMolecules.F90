!disclaimer
module moduleReplaceMolecule
  use moduleVariables
  use moduleStrings
  use moduleFiles
  use moduleMessages 
  use moduleNeighbours
  use moduleRead
  use moduleSuperimposeMolecules
  use moduleDistances
  use moduleProperties

  implicit none

  public :: replaceMolecules, replaceMoleculesHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  character(len=STRLEN), pointer :: targetMolecule
  integer, pointer :: numberOfLocalMolecules
  integer, pointer :: mapSize
  integer, pointer, dimension(:) :: moleculesIndices

  type(fileTypeDef), pointer :: inputFile
  logical, pointer :: frameRead, mapUsed
  integer, pointer, dimension(:) :: currentIndices
  integer, pointer, dimension(:) :: referenceIndices

  logical, pointer :: computeOnlyRMSD
  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: ID

contains

  subroutine replaceMoleculesHelp()
    implicit none
    call message(0,"This action can be used to replace selected molecules with a different molecule read from a file.")
    call message(0,"Alternatively this action could be used to simply compute the RMSD of the selected molecules with")
    call message(0,"respect to a reference structure.")
    call message(0,"At the moment the RMSD is simply written to a file, but it would be easy to compute its average or distribution.")
    call message(0,"Currently only molecules smaller than 50 atoms can be replaced.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i co3.pdb --top --subs +id M1 +map 1,2,3,4 1,4,3,2 +f hco3.pdb --o new.pdb")
    call message(0,"  gpta.x --i slab.pdb --top --subs +id M1 +f co3.pdb +rmsd rmsd.out")
  end subroutine replaceMoleculesHelp

  subroutine initialiseAction(a)
    use moduleStrings
    use moduleSystem
    implicit none
    type(actionTypeDef), target :: a
    integer :: i
    real(real64), dimension(3,3) :: hmat

    character(STRLEN) :: flagString

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0_real64

    ! molecule to replace
    call assignFlagValue(actionCommand,"+id",targetMolecule,"NULL")

    ! file with the reference molecule
    call assignFlagValue(actionCommand,"+f",inputFile % fname,'NULL')
    if (inputFile % fname /= 'NULL') then
      call initialiseFile(inputFile, inputFile % fname)

      ! -> get number of atoms
      call getNumberOfAtoms(inputFile, mapSize, hmat)
      rewind(inputFile % funit)

      call createSystemArrays(a % localFrame, mapSize)

      call readCoordinates(frameRead, a % localFrame, inputFile)

      if (.not. frameRead) call message(-1,"--replace | cannot read coordinates file")

    else
      call message(-1,"--replace - missing reference structure file")

    end if

    call assignFlagValue(actionCommand,"+map",mapUsed,.false.) 
    if (mapUsed) then
      call extractFlag(actionCommand,"+map",flagString)

      block
        integer :: nFields, n
        character(len=STRLEN), dimension(2) :: fields
        character(len=200), dimension(100) :: tokens

        call parse(flagString," ",fields,nFields)
        if (nFields /= 2) call message(-1,"--replace - +map requires two group of indices")

        call parse(fields(1),",",tokens,mapSize)
        currentIndices(1:mapSize) => a % integerVariables(101:100+mapSize)
        do i=1,mapSize
          read(tokens(i),*) currentIndices(i)
        end do

        call parse(fields(2),",",tokens,n)
        if (n /= mapSize) call message(-1,"--replace - different number of atoms for in mapping")
        referenceIndices(1:mapSize) => a % integerVariables(151:150+mapSize)
        do i=1,mapSize
          read(tokens(i),*) referenceIndices(i)
        end do

      end block

    else
      mapSize = a % localFrame % natoms
      currentIndices(1:mapSize) => a % integerVariables(101:100+mapSize)
      referenceIndices(1:mapSize) => a % integerVariables(151:150+mapSize)
      do i=1,mapSize
        currentIndices(i) = i
        referenceIndices(i) = i
      enddo

    end if
    
    call assignFlagValue(actionCommand,"+rmsd",computeOnlyRMSD,.false.)
    if (computeOnlyRMSD) then
      call assignFlagValue(actionCommand,"+rmsd",outputFile % fname,'rmsd.out')
      call initialiseFile(outputFile, outputFile % fname)
      call workData % initialise(ID, "dumpSeq", iounit=outputFile % funit)
    endif

    if ( size(referenceIndices) > 50 .or. size(currentIndices) > 50 ) call message(-1,"--replace - too many atoms in the molecule")

  end subroutine initialiseAction

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand           => a % actionDetails
    firstAction             => a % firstAction

    inputFile               => a % inputFile
    outputFile              => a % outputFile

    targetMolecule          => a % stringVariables(1)

    ID                      => a % integerVariables(1)
    numberOfLocalMolecules  => a % integerVariables(2)
    mapSize                 => a % integerVariables(3)

    moleculesIndices        => a % localIndices

    frameRead               => a % logicalVariables(1)
    computeOnlyRMSD         => a % logicalVariables(2)
    mapUsed                 => a % logicalVariables(3)

  end subroutine associatePointers

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Replace Molecules Action")
    call message(0,"...Target molecule",str=targetMolecule)
    call message(0,"...New molecule read from",str=inputFile % fname)
    if (mapUsed) then
      call message(0,"...Atoms used for the alignment",mapSize)
      call message(0,"......Current molecule",iv=currentIndices(1:mapSize))
      call message(0,"......New molecule",iv=referenceIndices(1:mapSize))
    else
      call message(0,"...All atoms used for the alignment")
    end if
    call message(2)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleSystem
    implicit none
    type(actionTypeDef), target :: a

    integer :: nref, ncur, ilocal, imol
    real(real64) :: rmsd, rvec(9), rotmat(3,3), dij(3)
    real(real64), allocatable, dimension(:,:) :: referencePositions
    real(real64), allocatable, dimension(:,:) :: currentPositions
    real(real64), allocatable, dimension(:) :: weights

    integer :: idx, jdx, iatm

    logical, allocatable, dimension(:) :: selectedMolecules
    
    integer :: newNumberOfAtoms
    real(real64), allocatable, dimension(:,:) :: newPositions
    character(len=cp), allocatable, dimension(:) :: newLabels

    real(real64), allocatable, dimension(:) :: localProperty

    ! Number of atoms in the new molecule
    nref = a % localFrame % natoms

    ! number of atoms in the current molecule
    imol = moleculesIndices(1)
    ncur = listOfMolecules(imol) % numberOfAtoms

    ! allocating arrays for the alignment
    allocate(referencePositions(3,mapSize))
    allocate(currentPositions(3,mapSize))

    ! extracting the atoms for the reference structure to match
    ! custom weights can be added here
    allocate(weights(mapSize), source=1.0_real64)
    do idx=1,mapSize
      iatm = referenceIndices(idx)
      referencePositions(1:3,idx) = a % localFrame % pos(1:3,iatm)
    end do

    ! shifting reference (new) molecule to the origin
    call CenterCoords(dij, mapSize, referencePositions, weights) 
    do idx=1, a % localFrame % natoms
      a % localFrame % pos(1:3,idx) = a % localFrame % pos(1:3,idx) - dij(1:3)
    end do

    if (computeOnlyRMSD) then
      allocate(localProperty(numberOfLocalMolecules))
      do ilocal=1,numberOfLocalMolecules
        imol = moleculesIndices(ilocal)
        do idx=1,mapSize
          jdx = currentIndices(idx)
          iatm = listOfMolecules(imol) % listOfAtoms(jdx)
          currentPositions(1:3,idx) = frame % pos(1:3,iatm)
        end do

        ! extracting the rotation matrix to rotate referencePositions onto currentPositions
        ! current positions are centered to the origin
        call Superimpose(mapSize, referencePositions, currentPositions, rmsd, rvec, weights)
        localProperty(ilocal) = rmsd
      end do
      call workData % compute(ID, numberOfValues=numberOfLocalMolecules, xValues=localProperty)

    else

      ! allocating temporary arrays for the new positions
      allocate(selectedMolecules(numberOfMolecules), source=.false.)
      do ilocal=1,numberOfLocalMolecules
        imol = moleculesIndices(ilocal)
        selectedMolecules(imol) = .true.
      end do
      newNumberOfAtoms = frame % natoms + count(selectedMolecules)*(nref-ncur)
      allocate(newPositions(3,newNumberOfAtoms))
      allocate(newLabels(newNumberOfAtoms))

      ! aligning and replacing the molecules
      newNumberOfAtoms = 0
      do imol=1,numberOfMolecules
      
        if (selectedMolecules(imol)) then
          ! extracting the atoms for the reference structure to match
          do idx=1,mapSize
            jdx = currentIndices(idx)
            iatm = listOfMolecules(imol) % listOfAtoms(jdx)
            currentPositions(1:3,idx) = frame % pos(1:3,iatm)
          end do
          ! This allows for a better centering of the molecule
          listOfMolecules(imol) % centreOfMass = computecenter(mapSize,currentPositions)

          ! extracting the rotation matrix to rotate referencePositions onto currentPositions
          ! current positions don't need to be centered to the origin
          call Superimpose(mapSize, referencePositions, currentPositions, rmsd, rvec, weights)
          rotmat = reshape(rvec,[3,3])

          ! rotating the reference molecule + shifting its original centre of mass
          do idx=1,a % localFrame % natoms
            dij = matmul(rotmat , a % localFrame % pos(1:3,idx))
            newPositions(1:3,newNumberOfAtoms+idx) = dij + listOfMolecules(imol) % centreOfMass
          end do

          do idx=1,a % localFrame % natoms
            newLabels(newNumberOfAtoms+idx) = a % localFrame % lab(idx)
          end do
          newNumberOfAtoms = newNumberOfAtoms + a % localFrame % natoms

        else

          ! ignoring this molecule, just copy it over the temporary arrays
          do idx=1,listOfMolecules(imol) % numberOfAtoms
            iatm = listOfMolecules(imol) % listOfAtoms(idx)
            newPositions(1:3,newNumberOfAtoms+idx) = frame % pos(1:3,iatm)
          end do
          do idx=1,listOfMolecules(imol) % numberOfAtoms
            iatm = listOfMolecules(imol) % listOfAtoms(idx)
            newLabels(newNumberOfAtoms+idx) = frame % lab(iatm)
          end do
          newNumberOfAtoms = newNumberOfAtoms + listOfMolecules(imol) % numberOfAtoms
        end if

      end do

      ! swapping allocations over
      frame % natoms = newNumberOfAtoms
      call move_alloc(newPositions, frame % pos)
      call move_alloc(newLabels, frame % lab)

      ! extending other frame arrays
      if (size(frame % chg) /= frame % natoms) then
        deallocate(frame % frac)
        deallocate(frame % chg)
        allocate(frame % frac(3,frame % natoms))
        allocate(frame % chg(frame % natoms), source=0.0_real64)
      end if

      ! recomputing topology
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
      call updateNeighboursList(.true.)
      call runInternalAction("topology","+update")
      
    end if ! RMSD only

  end subroutine computeAction

  subroutine replaceMolecules(a)
    use moduleVariables
    use moduleSystem 
    ! use moduleOpenMM
    use moduleDistances
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
        call dumpScreenInfo()
        
        if (numberOfMolecules == 0) call runInternalAction("topology","NULL")
        
        ! get the indices of the molecules to replace
        call getLocalIndices(targetMolecule, a, numberOfLocalMolecules)
        moleculesIndices => a % localIndices
        
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if

      call associatePointers(a)

      call computeAction(a)
    end if

  end subroutine replaceMolecules

end module moduleReplaceMolecule
