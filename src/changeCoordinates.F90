!disclaimer
module moduleModifyCoordinates

contains

  subroutine rescaleCellHelp()
    use moduleMessages
    implicit none
    call message(0,"This action changes the system cell and rescales the atomic coordinates, if required.")
    call message(0,"The new cell can be given using the +cell flag or as a scale version of the input one. (+scale)")
    call message(0,"If +scale is used the lenght of the lattice vectors are scaled, but the angles are kept constant.")
    call message(0,"The +nopos option can be used to leave the coordinates unchanged, useful to create a vaccum gap.")
    call message(0,"The cell can be expressed as:")
    call message(0,"  1 number -> cubic")
    call message(0,"  3 number -> orthorhombic")
    call message(0,"  6 number -> a, b, c, alpha, beta, gamma")
    call message(0,"  9 number -> metric matrix")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --rescale +cell 10")
    call message(0,"  gpta.x --i coord.pdb --rescale +cell 10,20,30")
    call message(0,"  gpta.x --i coord.pdb --rescale +cell 10,10,10,90,90,120")
    call message(0,"  gpta.x --i coord.pdb --rescale +cell 10,0,0,0,10,0,0,0,10")
    call message(0,"  gpta.x --i coord.pdb --rescale +scale 1,1,3 +nopos")
  end subroutine rescaleCellHelp

  subroutine shiftCoordinatesHelp()
    use moduleMessages
    implicit none
    call message(0,"This action shifts the atomic coordinates by a 3D vector.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --shift 1,2,3")
  end subroutine shiftCoordinatesHelp

  subroutine shiftCOMHelp()
    use moduleMessages
    implicit none
    call message(0,"This actions shifts the centre of mass of the selected atoms to:")
    call message(0,"  * the centre of the cell (+centre)")
    call message(0,"  * the same position as the first frame (+initial)")
    call message(0,"  * a user defined location (+loc)")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --fixcom +s Ca +centre")
    call message(0,"  gpta.x --i coord.pdb --fixcom +s Ca +initial")
    call message(0,"  gpta.x --i coord.pdb --fixcom +s Ca +loc 10,20,30")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --fixcom +s Ca +initial")
  end subroutine shiftCOMHelp

  subroutine applyPeriodicboundaryConditionsHelp()
    use moduleMessages
    implicit none
    call message(0,"This actions wraps the atoms/molecules inside the cell.")
    call message(0,"If the topology is known the molecues are not broken.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --pbc")
    call message(0,"  gpta.x --i coord.pdb --top --pbc")
    call message(0,"  gpta.x --i coord.pdb --top --pbc +nint")
  end subroutine applyPeriodicboundaryConditionsHelp

  subroutine unwrapCoordinatesHelp()
    use moduleMessages
    implicit none
    call message(0,"This actions unwraps the atoms/molecules across the PBC to generate a continuous trajectory.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb trajectory.dcd --unwrap")
    call message(0,"  gpta.x --i coord.pdb trajectory.dcd --top --unwrap")
  end subroutine unwrapCoordinatesHelp

  subroutine replicateCellHelp()
    use moduleMessages
    implicit none
    call message(0,"This action creates a supercell of the initial system.")
    call message(0,"If the topology is known the replication is done by molecule; the molecues must not be broken.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --repl 2,2,2")
    call message(0,"  gpta.x --i coord.pdb --top --repl 2,2,2")
  end subroutine replicateCellHelp

  subroutine mirrorCellHelp()
    use moduleMessages
    implicit none
    call message(0,"This action replicates the system by making its mirror image.")
    call message(0,"The mirror plane is one of the faces of the cell normal to x, y or z")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --mirror +x")
    call message(0,"  gpta.x --i coord.pdb --mirror +y")
    call message(0,"  gpta.x --i coord.pdb --mirror +z")
  end subroutine mirrorCellHelp

  subroutine removeOverlappingMoleculesHelp()
    use moduleMessages
    implicit none
    call message(0,"This action remove molecules that overlap. The atoms with the larger index are removed")
    call message(0,"The minimum distance between atoms can be specified by the +rmin flag.")
    call message(0,"The topology is calculated by default and the the entire molecules are removed.")
    call message(0,"If the selectoin flag is used the overla with only the selected atoms is considered.")
    call message(0,"It is also possible to do a dry run to see how many atoms/molecules would be removed with certain distances.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --noclash +rmin 1.1")
    call message(0,"  gpta.x --i coord.pdb --top --noclash +rmin 1.1")
    call message(0,"  gpta.x --i coord.pdb --top --noclash +rmin 1.1 +s K,Br")
    call message(0,"  gpta.x --i coord.pdb --top --noclash +rmin 1.1 +dry")
  end subroutine removeOverlappingMoleculesHelp

  subroutine fixCellHelp()
    use moduleMessages
    implicit none
    call message(0,"This action tries to straighthen the cell by changing the selected axys.")
    call message(0,"Examples:")
    call message(0,"  gpta --i calcite.pdb --repl 1,2,1 --fixcell +y --top +rebuild +reorder --pbc --o c.pdb")
    call message(0,"  gpta --i calcite.pdb --repl 2,1,1 --fixcell +x --top +rebuild +reorder --pbc --o c.pdb")
  end subroutine fixCellHelp

  subroutine defineNewCellHelp()
    use moduleMessages
    implicit none
    call message(0,"This action defines a new cell using the specified vectors.")
    call message(0,"The vectors can be given in the form of hkl or as a metric matrix.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --newcell +vec 1,0,0 0,1,0 0,0,1")
    call message(0,"  gpta.x --i coord.pdb --newcell +hmat 1,2,3 4,5,6 7,8,9")
  end subroutine defineNewCellHelp  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rescaleCell(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    integer :: numberOfWords
    character(len=STRLEN), dimension(100) :: listOfWords

    logical, pointer :: actionInitialisation
    character(len=100) :: stringCell
    real(real64), dimension(3,3), save :: hmat, hinv
    real(real64), save :: volume, cell(6)
    real(real64), dimension(3) :: cScaling

    logical :: lflag
    logical, pointer :: cellOnly

    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    cellOnly => a % logicalVariables(1)
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      call parse(actionCommand," ",listOfWords,numberOfWords)

      if (numberOfWords==0) call message(-1,"--rescale : nothing to do")

      hmat = frame % hmat
      call hmat2cell(hmat,cell,"DEG")

      call assignFlagValue(actionCommand,"+nopos ",cellOnly,.false.)

      call assignFlagValue(actionCommand,"+cell ", stringCell, "NONE")
      if (trim(stringCell) /= "NONE") then
        call readCellFreeFormat(stringCell, hmat)
        call getInverseCellMatrix(hmat,hinv,volume)
        call hmat2cell (hmat, cell, "DEG")
      end if

      call assignFlagValue(actionCommand,"+scale ", stringCell, "NONE")
      if (trim(stringCell) /= "NONE") then
        read(stringCell,*) cScaling(1:3)
        cell(1) = cell(1) * cScaling(1)
        cell(2) = cell(2) * cScaling(2)
        cell(3) = cell(3) * cScaling(3)
        call cell2hmat(cell,hmat)
        call getInverseCellMatrix(hmat,hinv,volume)
      end if

      call message(0,"Rescaled cell")
      call message(0,"...Cell vector A",rv=hmat(1:3,1))
      call message(0,"...Cell vector B",rv=hmat(1:3,2))
      call message(0,"...Cell vector C",rv=hmat(1:3,3))
      call message(0,"...Cell lengths",rv=cell(1:3))
      call message(0,"...Cell angles",rv=cell(4:6))
      call message(1,"...Volume",r=volume)

      call checkUsedFlags(actionCommand)
      return
    end if
    
    if (frameReadSuccessfully) then 
      frame % hmat = hmat
      frame % hinv = hinv
      frame % volume = volume
      frame % cell = cell
    
      if (.not. cellOnly) then
        call fractionalToCartesian(frame % natoms, frame % frac, frame % pos)
      end if
    end if

    if (endOfCoordinatesFiles) return

  end subroutine rescaleCell

  subroutine defineNewCell(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    logical, pointer :: actionInitialisation

    integer :: i, j

    character(STRLEN) :: flagString = ""
    integer :: numberOfWords, ntmp
    character(len=STRLEN), dimension(100) :: listOfWords
    character(len=STRLEN), dimension(3) :: listOfVectors
    character(len=STRLEN), dimension(3) :: string_l

    real(real64), save,  dimension(3,3) :: hkl, hmat_new(3,3)
    real(real64) :: hinv_new(3,3), cell_new(6), volume_new

    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      call parse(actionCommand," ",listOfWords,numberOfWords)

      if (numberOfWords==0) call message(-1,"--newcell : nothing to do")

      hkl = 0.0_real64
      hmat_new = 0.0_real64

      call extractFlag(actionCommand,"+hkl ",flagString)
      if (len_trim(flagString) > 0) then
        call parse(flagString," ",listOfVectors,ntmp)
        flagString = ''
        if (ntmp==0) then
          call message(-1,"--newcell : no vectors specified")
        else if (ntmp>3) then
          call message(-1,"--newcell : too many vectors specified")
        end if

        do i=1,3
          call parse(listOfVectors(i),",",string_l,ntmp)
          if (ntmp /= 3) then
            call message(-1,"--newcell : vector must have 3 components")
          end if
          read(string_l(1),*) hkl(1,i)
          read(string_l(2),*) hkl(2,i)
          read(string_l(3),*) hkl(3,i)
        end do
        do i=1,3
          do j=1,3
            hmat_new(:,i) = hmat_new(:,i) + frame % hmat(:,j) * hkl(j,i)
          end do
        end do
        where (abs(hmat_new) < 1.0e-6_real64) hmat_new = 0.0_real64
      end if

      call extractFlag(actionCommand,"+hmat",flagString)
      
      if (len_trim(flagString) > 0) then
        call parse(flagString," ",listOfVectors,ntmp)
        flagString = ''
        if (ntmp==0) then
          call message(-1,"--newcell : no vectors specified")
        else if (ntmp>3) then
          call message(-1,"--newcell : too many vectors specified")
        end if

        do i=1,3
          call parse(listOfVectors(i),",",string_l,ntmp)
          if (ntmp /= 3) then
            call message(-1,"--newcell : vector must have 3 components")
          end if
          read(string_l(1),*) hmat_new(1,i)
          read(string_l(2),*) hmat_new(2,i)
          read(string_l(3),*) hmat_new(3,i)
        end do
      end if

      call checkUsedFlags(actionCommand)
      return
    end if
    
    if (frameReadSuccessfully) then 
      call create_new_cell(hmat_new)

      call getInverseCellMatrix(hmat_new,hinv_new,volume_new)
      call message(0,"New cell defined as")
      call message(0,"...Cell vector A",rv=hmat_new(1:3,1))
      call message(0,"...Cell vector B",rv=hmat_new(1:3,2))
      call message(0,"...Cell vector C",rv=hmat_new(1:3,3))
      call message(1,"...Volume",r=volume_new)

      call message(0,"Final cell")
      call message(0,"...Cell vector A",rv=frame % hmat(1:3,1))
      call message(0,"...Cell vector B",rv=frame % hmat(1:3,2))
      call message(0,"...Cell vector C",rv=frame % hmat(1:3,3))
      call message(0,"...Cell lengths",rv=frame % cell(1:3))
      call message(0,"...Cell angles",rv=frame % cell(4:6))
      call message(1,"...Volume",r=frame % volume)
    end if

    if (endOfCoordinatesFiles) return


  end subroutine defineNewCell

  subroutine shiftCoordinates_old(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    integer :: numberOfWords
    character(len=STRLEN), dimension(100) :: listOfWords

    logical, pointer :: actionInitialisation

    real(real64), pointer, dimension(:) :: dpos
    integer :: i

    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    dpos(1:3)            => a % doubleVariables(1:3)
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      dpos = 0.0_real64
      call parse(actionCommand," ",listOfWords,numberOfWords)

      if (numberOfWords==0) call message(-1,"--shift : nothing to do")

      if (numberOfWords==3) then
        read(actionCommand,*)dpos(1:3)
      else
        call parse(actionCommand,",",listOfWords,numberOfWords)
        if (numberOfWords==3) read(actionCommand,*)dpos(1:3)
      end if

      if (sum(abs(dpos(1:3)))<1.0e-6_real64) call message(-1,"--shift : zero displacement")
      
      call message(1,"Translating atoms' coordinates by",rv=dpos)

      call checkUsedFlags(actionCommand)
      return
    end if
    
    if (frameReadSuccessfully) then 
      do i=1,frame % natoms
        frame % pos(1:3,i) = frame % pos(1:3,i) + dpos
      enddo
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
    end if

!    if (numberOfMolecules>0) call computeMoleculesCOM(0)
    if (numberOfMolecules>0) then
      block
        integer :: imol
        do imol=1,numberOfMolecules
          listOfMolecules(imol) % centreOfMass = listOfMolecules(imol) % centreOfMass + dpos
        end do
      end block
    end if
    if (endOfCoordinatesFiles) return

  end subroutine shiftCoordinates_old

  subroutine shiftCoordinates(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    integer :: numberOfFlags, numberOfWords, numberOfStrings
    character(len=STRLEN), dimension(100) :: listOfWords
    character(len=STRLEN), dimension(100) :: listOfFlags

    logical, pointer :: actionInitialisation

    real(real64), pointer, dimension(:) :: dpos
    integer :: i, iatm, nsel1

    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    dpos(1:3)            => a % doubleVariables(1:3)
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.

      dpos = 0.0_real64

      ! Check for flags
      call parse2(actionCommand,"+",listOfFlags,numberOfFlags)

      ! If no falgs default to old behaviour
      if (numberOfFlags == 0) then
        call parse(actionCommand," ",listOfWords,numberOfWords)
        if (numberOfWords == 0) call message(-1,"--shift : nothing to do")

        if (numberOfWords==3) then
          read(actionCommand,*)dpos(1:3)
        else
          call parse(actionCommand,",",listOfWords,numberOfWords)
          if (numberOfWords==3) read(actionCommand,*)dpos(1:3)
        end if

        actionCommand = "+s all +vec "//actionCommand
      end if

      call assignFlagValue(actionCommand,"+vec",dpos,[0.d0,0.d0,0.d0])

      if (sum(abs(dpos(1:3)))<1.0e-6_real64) call message(-1,"--shift : zero displacement")
      
      call message(1,"Translating atoms' coordinates by",rv=dpos)
      return
    end if
    
    if (frameReadSuccessfully) then 
      call selectAtoms(1,actionCommand,a)
      call createSelectionList(a,1)

      ! do i=1,frame % natoms
      nsel1 = count(a % isSelected(:,1))
      do i=1,nsel1
        iatm = a % idxSelection(i,1)
        frame % pos(1:3,iatm) = frame % pos(1:3,iatm) + dpos
      enddo
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

      call checkUsedFlags(actionCommand)
    end if

    if (numberOfMolecules>0) call computeMoleculesCOM(0)
    ! if (numberOfMolecules>0) then
    !   block
    !     integer :: imol
    !     do imol=1,numberOfMolecules
    !       listOfMolecules(imol) % centreOfMass = listOfMolecules(imol) % centreOfMass + dpos
    !     end do
    !   end block
    ! end if
    if (endOfCoordinatesFiles) return

  end subroutine shiftCoordinates

  subroutine shiftCOM(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand
    logical, pointer :: firstAction
    logical, pointer :: actionInitialisation
    logical, pointer :: lcentre
    logical, pointer :: linitial
    real(real64), pointer, dimension(:) :: dpos
    
    real(real64), dimension(3) :: xcom, xcentre, delta

    integer :: i, idx
    integer, save :: nsel

    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    dpos(1:3)            => a % doubleVariables(1:3)
    firstAction          => a % firstAction
    lcentre              => a % logicalVariables(1)
    linitial             => a % logicalVariables(2)

    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      ! fix the centre of mass to the cell centre      
      call assignFlagValue(actionCommand,"+centre",lcentre,.false.)
      
      ! fix the centre of mass intial position
      call assignFlagValue(actionCommand,"+initial",linitial,.false.)

      if (.not. lcentre) call assignFlagValue(actionCommand,"+loc",dpos,[0.0_real64,0.0_real64,0.0_real64])

      ! call checkUsedFlags(actionCommand)

      return
    end if
    
    if (frameReadSuccessfully) then 

      if (firstAction) then
        if (lcentre) then
          call message(0,"Translating centre of mass to the centre of the cell")

        else if (linitial) then
          call message(0,"Translating centre of mass to its initial position")
        
        else
          call message(0,"Translating centre of mass to",rv=dpos)
        end if

        call selectAtoms(1,actionCommand,a)
        call createSelectionList(a,1)
        if (linitial) then
          block
            nsel = count(a % isSelected(:,1))

            dpos = 0.0_real64
            do i=1,nsel
              idx = a % idxSelection(i,1)
              dpos(1:3) = dpos(1:3) + frame % pos(1:3,idx)
            end do
            dpos(1:3) = dpos(1:3) / dble(nsel)
          end block
        end if

        if (resetFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if

        call checkUsedFlags(actionCommand)

        firstAction = .false.

        if (linitial) return

      else

        if (a % updateAtomsSelection) then
          call selectAtoms(1,actionCommand,a)
          call createSelectionList(a,1)
        end if

      end if
      nsel = count(a % isSelected(:,1))

      xcom = 0.0_real64
      do i=1,nsel
        idx = a % idxSelection(i,1)
        xcom(1:3) = xcom(1:3) + frame % pos(1:3,idx)
      end do
      xcom(1:3) = xcom(1:3) / dble(nsel)
      
      if (lcentre) then
        xcentre(1:3) = frame % hmat(1:3,1) + frame % hmat(1:3,2) + frame % hmat(1:3,3)
        xcentre = xcentre / 2.0_real64
        delta(1:3) = xcentre(1:3) - xcom(1:3)
      else
        delta(1:3) = dpos(1:3) - xcom(1:3)
      end if

      do i=1,frame % natoms
        frame % pos(1:3,i) = frame % pos(1:3,i) + delta
      enddo
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
    end if

!    if (numberOfMolecules>0) call computeMoleculesCOM(0)
    if (numberOfMolecules>0) then
      block
        integer :: imol
        do imol=1,numberOfMolecules
          listOfMolecules(imol) % centreOfMass = listOfMolecules(imol) % centreOfMass + delta
        end do
      end block
    end if

    if (endOfCoordinatesFiles) return

  end subroutine shiftCOM

  subroutine applyPeriodicboundaryConditions(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    integer :: i
    logical, pointer :: lnint
    real(real64), allocatable, dimension(:,:) :: tmpArray

    lnint => a % logicalVariables(1)

    if (a % actionInitialisation) then
      a % actionInitialisation = .false.
      a % requiresNeighboursList = .false.

      ! fix the centre of mass to the cell centre
      call assignFlagValue(a % actionDetails,"+nint",lnint,.false.)
      if (lnint) then
        call message(1,"Applying periodic boundary conditions",str='[-0.5:+0.5]')
      else
        call message(1,"Applying periodic boundary conditions",str='[0:1]')
      endif
      return
    end if
    
    if (frameReadSuccessfully) then 
      if (lnint) then
        allocate(tmpArray(3,frame % natoms))
        do i=1,frame % natoms
          tmpArray(1:3,i) = frame % frac(1:3,i) - nint(frame % frac(1:3,i))
        end do
        call fractionalToCartesian(frame % natoms, tmpArray, frame % pos)
      else
        call fractionalToCartesian(frame % natoms,frame % frac,frame % pos)
      endif
      if (numberOfMolecules > 0) call reassembleAllMolecules()
    end if

    if (endOfCoordinatesFiles) return

  end subroutine applyPeriodicboundaryConditions

  subroutine unwrapCoordinates(a)
    use moduleVariables
    use moduleSystem 
    use moduleDistances
    use moduleMessages
    implicit none
    type(actionTypeDef), target :: a

    integer :: iatm
    real(real64) :: dij(3), dist

    ! Associate variables
    
    if (a % actionInitialisation) then
      a % actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      return
    end if
    
    if (frameReadSuccessfully) then 

      if (a % firstAction) then
        call message(1,"Unwrapping trajectory")
        allocate(a % localPositions(3,frame % natoms), source=frame % pos)
        call checkUsedFlags(a % actionDetails)
        a % firstAction = .false.
        return
      endif

      do iatm=1,frame % natoms
        dij = frame % pos(:,iatm) - a % localPositions(:,iatm)
        dist = computeDistanceSquaredPBC(dij)
        frame % pos(:,iatm) = a % localPositions(:,iatm) + dij
      end do
      a % localPositions = frame % pos

      if (numberOfMolecules > 0) call reassembleAllMolecules()
    end if

    if (endOfCoordinatesFiles) return

  end subroutine unwrapCoordinates

  subroutine replicateCell(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    use moduleNeighbours
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand
    integer, pointer :: tallyExecutions
    logical, pointer :: actionInitialisation

    real(real64), allocatable, dimension(:,:) :: cartesianCoord
    real(real64), allocatable, dimension(:) :: saveCharges
    character(cp), allocatable, dimension(:) :: saveLabels
    real(real64), allocatable, dimension(:) :: saveMasses
    
    integer, pointer, dimension(:) :: idx, jdx, kdx
    integer :: nargs
    character(len=10) :: args(10)
    integer :: iatm, i, j, k, nrepl, nn, natoms_old
    
    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    tallyExecutions      => a % tallyExecutions

    idx(1:2)             => a % integerVariables(1:2)
    jdx(1:2)             => a % integerVariables(3:4)
    kdx(1:2)             => a % integerVariables(5:6)
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      call message(1,"Replicating cell")
      
      call parse(actionCommand," ",args,nargs)
      if (nargs == 1) then
        call parse(args(1),",",args,nargs)
      endif
      if (nargs/=3) call message(-1,"--repl - wrong number of arguments")
      idx = 0
      jdx = 0
      kdx = 0
      read(args(1),*)idx(2)
      read(args(2),*)jdx(2)
      read(args(3),*)kdx(2)
      idx(2) = idx(2) - 1
      jdx(2) = jdx(2) - 1
      kdx(2) = kdx(2) - 1

      call checkUsedFlags(actionCommand)

      tallyExecutions = 0 

      return
    end if
    
    if (frameReadSuccessfully) then 
      if (tallyExecutions > 0) then
        call message(-1,"Replicating cell can be done only on one frame")
      else
        tallyExecutions = tallyExecutions + 1
      end if

      if (numberOfMolecules > 0) call runInternalAction("pbc","NULL")

      nrepl = (idx(2)-idx(1)+1) * (jdx(2)-jdx(1)+1) * (kdx(2)-kdx(1)+1)

      allocate(cartesianCoord(3, frame % natoms))
      allocate(    saveLabels(   frame % natoms))
      allocate(   saveCharges(   frame % natoms))
      allocate(    saveMasses(   frame % natoms))
      cartesianCoord = frame % pos
          saveLabels = frame % lab
         saveCharges = frame % chg
          saveMasses = frame % mass
  
      natoms_old = frame % natoms
      call deleteSystemArrays(frame)
      call createSystemArrays(frame, natoms_old * nrepl ) ! also sets frame % natoms

      nn = 0 
      do k=kdx(1),kdx(2)
        do j=jdx(1),jdx(2)
          do i=idx(1),idx(2)
            do iatm=1,natoms_old
              nn = nn + 1
              frame % pos(1:3,nn) = cartesianCoord(1:3,iatm) & 
                                  + i*frame % hmat(1:3,1) &
                                  + j*frame % hmat(1:3,2) &
                                  + k*frame % hmat(1:3,3)
              frame % lab(nn) = saveLabels(iatm)
              frame % chg(nn) = saveCharges(iatm)
              frame % mass(nn)= saveMasses(iatm)
            end do
          end do
        end do
      end do
      
      frame % hmat(1:3,1) = (idx(2)-idx(1)+1)*frame % hmat(1:3,1)
      frame % hmat(1:3,2) = (jdx(2)-jdx(1)+1)*frame % hmat(1:3,2)
      frame % hmat(1:3,3) = (kdx(2)-kdx(1)+1)*frame % hmat(1:3,3)

      call getInverseCellMatrix(frame % hmat, frame % hinv, frame % volume)
      call hmat2cell (frame % hmat, frame % cell, "DEG")

      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

      call setUpNeigboursList()
      call updateNeighboursList(.true.)
      
      call message(0,"New cell")
      call message(0,"...Cell vector A",rv=frame % hmat(1:3,1))
      call message(0,"...Cell vector B",rv=frame % hmat(1:3,2))
      call message(0,"...Cell vector C",rv=frame % hmat(1:3,3))
      call message(0,"...Cell lengths",rv=frame % cell(1:3))
      call message(0,"...Cell angles",rv=frame % cell(4:6))
      call message(1,"...Volume",r=frame % volume)

      if (numberOfMolecules > 0) call runInternalAction("topology","+update +reorder +rebuild")

    end if

    if (endOfCoordinatesFiles) return

  end subroutine replicateCell

  subroutine mirrorCell(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    use moduleNeighbours
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand  
    logical, pointer :: firstAction

    logical, pointer :: actionInitialisation
    
    integer, pointer :: idx
    integer :: nargs
    character(len=10) :: args(10)
    real(real64), dimension(3) :: rtmp, shift
    integer :: iatm, nrepl, nn, i
    logical :: lflag(3)
    
    ! Associate variables
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    actionInitialisation => a % actionInitialisation

    idx                  => a % integerVariables(1)
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      call message(1,"Mirroring the cell")
      
      call parse(actionCommand," ",args,nargs)
      call assignFlagValue(actionCommand,"+x",lflag(1),.false.)
      call assignFlagValue(actionCommand,"+y",lflag(2),.false.)
      call assignFlagValue(actionCommand,"+z",lflag(3),.false.)
      if (count(lflag) == 0) then
        call message(-1,"--mirror | no direction specified")
      else if (count(lflag) > 1) then
        call message(-1,"--mirror | only one direction can be specified at one time")
      else
        do i=1,3
          if (lflag(i)) idx = i
        end do
      end if

      call checkUsedFlags(actionCommand)

      return
    end if
    
    if (frameReadSuccessfully) then 

      if (numberOfMolecules > 0) call runInternalAction("pbc","NULL")

      nrepl = 2
     
      if (a % firstAction) then
        call extendFrame(frame % natoms * nrepl)
        call checkUsedFlags(actionCommand)
        a % firstAction = .false.
      end if

      rtmp = 2.0_real64 * frame % hmat(:,idx)
      ! nn = frame % natoms 
      ! do iatm=1,frame % natoms
      !   frame % pos(1:3,nn+iatm) = frame % pos(1:3,iatm) 
      ! end do

      nn = frame % natoms 

      frame % natoms = frame % natoms * nrepl
      frame % hmat(1:3,idx) = nrepl * frame % hmat(1:3,idx)
      frame % lab(nn+1:nn+nn) = frame % lab(1:nn)
      frame % chg(nn+1:nn+nn) = frame % chg(1:nn)

      rtmp = [1.d0 , 1.d0 , 1.d0] ; rtmp(idx) = -1.d0
      do iatm=1,nn
        frame % pos(1:3,nn+iatm) = rtmp(1:3) * frame % pos(1:3,iatm)
      end do
      
      ! Shift the atoms back in the cell
      shift = frame % hmat(:,idx) / 2.0_real64
      do iatm=1,frame % natoms
        frame % pos(1:3,iatm) = frame % pos(1:3,iatm) + shift(1:3)
      end do

      call getInverseCellMatrix(frame % hmat, frame % hinv, frame % volume)
      call hmat2cell (frame % hmat, frame % cell, "DEG")

      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

      call setUpNeigboursList()
      call updateNeighboursList(.true.)
      
      if (numberOfMolecules > 0) call runInternalAction("topology","+update +rebuild")

    end if

    if (endOfCoordinatesFiles) return

  end subroutine mirrorCell

  subroutine removeOverlappingMolecules(a)
    use moduleDistances
    use moduleSystem
    use moduleVariables
    use moduleMessages 
    use moduleNeighbours
    use moduleStrings

    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand  
    logical, pointer :: firstAction
    logical, pointer :: dryRun, lselect

    real(real64), pointer :: overlapDistance
    integer ::natm0, natm1
    integer :: idx

    ! Associate variables
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction

    overlapDistance      => a % doubleVariables(1)
    dryRun               => a % logicalVariables(1)
    lselect              => a % logicalVariables(2)

    if (a % actionInitialisation) then
      a % actionInitialisation = .false.
      a % requiresNeighboursList = .true.
      a % requiresNeighboursListUpdates = .true.
      a % requiresNeighboursListDouble = .true.
      
      call assignFlagValue(actionCommand,"+rmin",overlapDistance,1.0_real64)
      a % cutoffNeighboursList = overlapDistance

      call assignFlagValue(actionCommand,"+dry",dryRun,.false.)

      call assignFlagValue(actionCommand,"+s",lselect,.false.)

      call checkUsedFlags(actionCommand)

      return
    end if

    if (frameReadSuccessfully) then 

      if (firstAction) then
        call message(0,"Removing overlaping molecules")
        call checkUsedFlags(actionCommand)
        firstAction = .false.

        if (lselect) then
          call selectAtoms(1,actionCommand,a)
          call createSelectionList(a,1)
        endif
  
      end if
      
      block
        integer :: iatm, ineigh
        integer :: ndist0 = 0
        integer :: ndist1 = 0
        integer :: ndist2 = 0
        integer :: ndist3 = 0
        real(real64) :: rdist1 = 0.75_real64
        real(real64) :: rdist2 = 0.50_real64
        real(real64) :: rdist3 = 0.10_real64
        character(len=5) :: str
        do iatm=1,frame % natoms
          do ineigh=1,nneigh(iatm)
            if (rneigh(ineigh,iatm) < overlapDistance) ndist0 = ndist0 + 1
            if (rneigh(ineigh,iatm) < rdist1) ndist1 = ndist1 + 1
            if (rneigh(ineigh,iatm) < rdist2) ndist2 = ndist2 + 1
            if (rneigh(ineigh,iatm) < rdist3) ndist3 = ndist3 + 1
          end do
        end do

        write(str,'(f5.3)')overlapDistance
        call message(0,"...Number of contacts below "//trim(str)//" A",i=ndist0)
        write(str,'(f5.2)')rdist1
        call message(0,"...Number of contacts below "//trim(str)//" A",i=ndist1)
        write(str,'(f5.2)')rdist2
        call message(0,"...Number of contacts below "//trim(str)//" A",i=ndist2)
        write(str,'(f5.2)')rdist3
        call message(0,"...Number of contacts below "//trim(str)//" A",i=ndist3)
      end block
        
      if (dryRun) then
        call message(1,"...Dry run - no molecules removed !!")
        return
      else
        call message(0,"...Minimum overlap distance",r=overlapDistance)
      end if

      natm0 = frame % natoms
      ! if (numberOfMolecules == 0) call runInternalAction("topology","NULL")

      if (lselect) then
        call deleteOverlapAllMolecules(overlapDistance, &
                                        size(a % idxSelection(:,1)), &
                                        a % idxSelection(:,1))
      else
        call deleteOverlapAllMolecules(overlapDistance)
      end if
                                                                      
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

      ! update neighbours' list
      call updateNeighboursList(.true.)
      
      natm1 = frame % natoms
      call message(0,"...Number of atoms removed",i=natm0-natm1)
      call message(2)
      !if (numberOfMolecules > 0) then
      !  call runInternalAction("topology","+update")
      !end if

    end if

    if (endOfCoordinatesFiles) return
    
  end subroutine removeOverlappingMolecules

  subroutine fixCell(a)
    use moduleSystem
    use moduleStrings
    use moduleVariables
    use moduleDistances
    use moduleMessages
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand  
    logical, pointer :: firstAction
    logical, pointer :: actionInitialisation

    character(len=1) :: dirs(3)
    integer :: nargs, i, j
    integer, save :: idx
    character(len=10) :: args(10)
    logical :: lflag(3)
    real(real64), save :: box(3,3), tmp(3,3)
  
    ! Associate variables
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    actionInitialisation => a % actionInitialisation

    dirs = ["x","y","z"]

    if (a % actionInitialisation) then
      a % actionInitialisation = .false.
      a % requiresNeighboursList = .false.

      call parse(actionCommand," ",args,nargs)
      call assignFlagValue(actionCommand,"+x",lflag(1),.false.)
      call assignFlagValue(actionCommand,"+y",lflag(2),.false.)
      call assignFlagValue(actionCommand,"+z",lflag(3),.false.)


      if (count(lflag) == 0) then
        call message(-1,"--fixcell | no direction specified")
      else if (count(lflag) > 1) then
        call message(-1,"--fixcell | only one direction can be specified at one time")
      else
        do i=1,3
          if (lflag(i)) idx = i
        end do
      end if

      call checkUsedFlags(actionCommand)

      return
    end if

    if (frameReadSuccessfully) then
      call message(0,"Cleaning up the cell...")
      call straightenCell(frame % hmat,dirs(idx),box)

      do i=1,3
        do j=1,3
          tmp(i,j) = abs(frame % hmat(i,j) - box(i,j))
        end do
      end do

      if (sum(tmp) < 1.0e-6_real64) then
        call message(1,"...Cell already clean")
        return
      end if

      if (abs(box(2,1)) + abs(box(3,1)) + abs(box(3,2)) .gt. 1.0e-6_real64) then
        call makeUpperTriangularCell(box, frame % pos, frame % natoms)
      end if

      do i=1,3
        do j=1,3
          if (abs(box(i,j)) < 1.0e-6_real64) then
            box(i,j) = 0.0_real64
          end if
        end do
      end do

      frame % hmat = box
      call getInverseCellMatrix(frame % hmat, frame % hinv, frame % volume)
      call hmat2cell (frame % hmat, frame % cell, "DEG")

      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

      call message(0,"...New cell")
      call message(0,"...Cell vector A",rv=frame % hmat(1:3,1))
      call message(0,"...Cell vector B",rv=frame % hmat(1:3,2))
      call message(0,"...Cell vector C",rv=frame % hmat(1:3,3))
      call message(0,"...Cell lengths",rv=frame % cell(1:3))
      call message(0,"...Cell angles",rv=frame % cell(4:6))
      call message(1,"...Volume",r=frame % volume)

    end if
    ! rtmp = (box - frame % hmat)

    ! iVal = 0
    ! iVal(1) = nint(sqrt(sum(rtmp(1:1,2)**2)) / frame % cell(1))
    ! ! iVal(2) = nint(sqrt(sum(rtmp(1,3)**2)) / frame % cell(1))
    ! ! iVal(3) = nint(sqrt(sum(rtmp(1:2,3)**2)) / frame % cell(2))

    ! if (iVal(1) > 0) then
    !   rVec = sign(dble(iVal(1)), box(1,2))
    !   frame % hmat(:,2) = frame % hmat(:,2) + rVec * frame % hmat(:,1)
    ! end if

    ! ! if (iVal(2) > 0) then
    ! !   rVec = sign(dble(iVal(2)), box(1,3))
    ! !   frame % hmat(:,3) = frame % hmat(:,3) + rVec * frame % hmat(:,1)
    ! ! end if

    ! ! if (iVal(3) > 0) then
    ! !   rVec = sign(dble(iVal(2)), box(2,3))
    ! !   frame % hmat(:,3) = frame % hmat(:,3) + rVec * frame % hmat(:,2)
    ! ! end if

    ! if (any(iVal>0)) then
    !   call hmat2cell (frame % hmat, frame % cell, "DEG")
    !   call getInverseCellMatrix(frame % hmat,frame % hinv,frame % volume)
    !   call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
    ! end if
    
    ! box = frame % hmat

  end subroutine fixCell

end module moduleModifyCoordinates

subroutine checkForBrokenMolecules(brokenMolecule)
  use moduleSystem 
  use moduleDistances
  use moduleMessages
  implicit none

  integer :: ibond, iatm, jatm, imol
  real(real64) :: dij(3), distance2
  logical, intent(out) :: brokenMolecule

  brokenMolecule = .false.
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (ibond, iatm, jatm, dij, distance2, imol)
!$OMP DO
  do ibond=1,numberOfUniqueBonds 
    iatm = listOfUniqueBonds(1,ibond)
    jatm = listOfUniqueBonds(2,ibond)
    dij = frame % pos(1:3,iatm) - frame % pos(1:3,jatm)
    distance2 = sum(dij*dij)
    if (distance2 > 9.0_real64) then
      brokenMolecule = .true.
      imol = atomToMoleculeIndex(iatm)
      listOfMolecules(imol) % brokenBonds = .true.
    end if
  enddo
!$OMP END DO
!$OMP END PARALLEL
  
end subroutine checkForBrokenMolecules

subroutine reassembleAllMolecules()
  use moduleVariables
  use moduleSystem 
  use moduleDistances
  implicit none

  integer :: imol, idx, ib, iatm, jatm
  integer :: n, m, itmp(1)
  integer, allocatable, dimension(:) :: atomList
  real(real64) :: dij(3), rtmp

  do imol=1,numberOfMolecules
    n = listOfMolecules(imol) % numberOfAtoms
    allocate(atomList(n),source=0)
    atomList(1) = 1
    do while (count(atomList==2) < n)
      idx = getatom(n,atomList)
      iatm = listOfMolecules(imol) % listOfAtoms(idx)
      do ib=1,numberOfCovalentBondsPerAtom(iatm)
        jatm = listOfCovalentBondsPerAtom(ib,iatm)
        itmp = minloc( (listOfMolecules(imol) % listOfAtoms(:)-jatm)**2 )
        if (atomList(itmp(1)) ==2) cycle

        dij(1:3) = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
        rtmp = computeDistanceSquaredPBC(dij)
        frame % pos(1:3,jatm) = frame % pos(1:3,iatm) + dij(1:3)
  
        atomList(itmp(1)) = 1
      end do
      atomList(idx) = 2
    end do
    deallocate(atomList)
  end do

  contains
  integer function getatom(n,l) result(i)
    implicit none
    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: l
    do i=1,n
      if (l(i)==1) exit
    end do
  end function getatom
end subroutine reassembleAllMolecules

subroutine reassembleAllMolecules2(n,pos)
  use moduleVariables
  use moduleSystem 
  use moduleDistances
  implicit none

  integer :: n
  real(real64), dimension(3,n) :: pos

  integer :: imol, idx, iatm, jatm
  real(real64) :: dij(3), rtmp

  do imol=1,numberOfMolecules
    iatm = listOfMolecules(imol) % listOfAtoms(1)
    do idx=2,listOfMolecules(imol) % numberOfAtoms
      jatm = listOfMolecules(imol) % listOfAtoms(idx)
      dij(1:3) = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
      rtmp = computeDistanceSquaredPBC(dij)
      pos(1:3,jatm) = pos(1:3,iatm) + dij(1:3)
    enddo
  enddo

end subroutine reassembleAllMolecules2

subroutine reassembleBrokenMolecules()
  use moduleVariables
  use moduleSystem 
  use moduleDistances
  implicit none

  integer :: imol, idx, iatm, jatm
  real(real64) :: dij(3), rtmp

  do imol=1,numberOfMolecules
    if (.not. listOfMolecules(imol) % brokenBonds) cycle

    iatm = listOfMolecules(imol) % listOfAtoms(1)
    do idx=2,listOfMolecules(imol) % numberOfAtoms
      jatm = listOfMolecules(imol) % listOfAtoms(idx)
      dij(1:3) = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
      rtmp = computeDistanceSquaredPBC(dij)
      frame % pos(1:3,jatm) = frame % pos(1:3,iatm) + dij(1:3)
    enddo
    listOfMolecules(imol) % brokenBonds = .false.
  enddo

end subroutine reassembleBrokenMolecules

subroutine computeMoleculesCOM(mol0)
  use moduleSystem 
  use moduleElements 
  use moduleSuperimposeMolecules, only : computeCenter

  implicit none

  integer, intent(in) :: mol0
  integer :: imol, idx, iatm, n

  real(real64), allocatable, dimension(:,:) :: localPositions

  idx = maxval(listOfMolecules(mol0+1:numberOfMolecules) % numberOfAtoms)
  allocate(localPositions(3,idx))
  do imol=mol0+1,numberOfMolecules
    n = listOfMolecules(imol) % numberOfAtoms
    do idx=1,n
      iatm = listOfMolecules(imol) % listOfAtoms(idx)
      localPositions(1:3,idx) = frame % pos(1:3,iatm)
    end do
    listOfMolecules(imol) % centreOfMass(1:3) = computeCenter(n, localPositions(:,1:n))
  end do

end subroutine computeMoleculesCOM

subroutine deleteOverlapAllMolecules(rcut,n,alist)
  use moduleSystem
  use moduleVariables
  use moduleDistances
  implicit none
  real(real64), intent(in) :: rcut
  integer, intent(in), optional :: n
  integer, dimension(*), intent(in), optional :: alist

  logical, allocatable, dimension(:) :: lremove
  integer :: imol, jmol
  integer :: idx
  integer :: iatm, jatm, ineigh
  real(real64) :: dij(3), rtmp, rmax
  integer :: n0
  integer :: kdx

  integer, allocatable, dimension(:) :: tmpIDX

  if (.not. present(n)) then
    allocate(tmpIDX(frame % natoms))
    do idx=1,frame % natoms
      tmpIDX(idx) = idx
    end do
  else
    allocate(tmpIDX(n),source=alist(1:n))
  end if

  allocate(lremove(frame % natoms), source=.false.)
  do idx=1,size(tmpIDX)
    iatm = tmpIDX(idx)
    if (lremove(iatm)) cycle
    do ineigh=1,nneigh(iatm)
      jatm = lneigh(ineigh,iatm)
      if (lremove(jatm)) cycle
      if (rneigh(ineigh,iatm) < rcut) lremove(jatm) = .true.
    end do
  end do
  
  if (numberOfMolecules > 0) then
    do jatm=1,frame % natoms
      if (lremove(jatm)) then
        jmol = atomToMoleculeIndex(jatm)
        do idx=1,listOfMolecules(jmol) % numberOfAtoms
          iatm = listOfMolecules(jmol) % listOfAtoms(idx)
          imol = atomToMoleculeIndex(iatm)
          lremove(iatm) = .true.
        end do
      end if
    end do
  end if

  n0 = 0
  do iatm=1,frame % natoms
    if (lremove(iatm)) cycle
    n0 = n0 + 1
    frame % lab(n0) = frame % lab(iatm)
    frame % chg(n0) = frame % chg(iatm)
    frame % pos(:,n0) = frame % pos(:,iatm)
  end do
  frame % natoms = n0

end subroutine deleteOverlapAllMolecules

! subroutine deleteOverlapAllAtoms(rcut,removeSecond)
!   use moduleSystem
!   use moduleVariables
!   use moduleDistances
!   implicit none
!   real(real64), intent(in) :: rcut
!   logical, intent(in) :: removeSecond

!   logical, allocatable, dimension(:) :: lremove
!   integer, target :: iatm, jatm, ineigh
!   real(real64) :: dij(3), rtmp, rmax
!   integer :: n0
!   integer, pointer :: idx

!   rmax = rcut**2

!   n0 = frame % natoms
!   allocate(lremove(n0), source=.false.)

!   if (removeSecond) then
!     idx => jatm
!   else
!     idx => iatm
!   end if

!   m1 : do iatm=1,n0
!     if (lremove(iatm)) cycle

!     m2 : do ineigh=1,nneigh(iatm)
!       jatm=lneigh(ineigh,iatm)
!       if (lremove(jatm)) cycle

!       dij = frame % pos(1:3,iatm) - frame % pos(1:3,jatm)
!       rtmp = computeDistanceSquaredPBC(dij) 
!       if (rtmp < rmax) then
!         lremove(idx) = .true.
!         if (idx == iatm) then
!           cycle m1
!         else
!           cycle m2
!         end if
!       end if

!     end do m2
!   end do m1

!   n0 = 0
!   do iatm=1,frame % natoms
!     if (lremove(iatm)) cycle
!     n0 = n0 + 1
!     frame % lab(n0) = frame % lab(iatm)
!     frame % chg(n0) = frame % chg(iatm)
!     frame % pos(:,n0) = frame % pos(:,iatm)
!   end do
!   frame % natoms = n0

! end subroutine deleteOverlapAllAtoms


! subroutine deleteOverlapAllMolecules(rcut,removeSecond)
!   use moduleSystem
!   use moduleVariables
!   use moduleDistances
!   implicit none
!   real(real64), intent(in) :: rcut
!   logical, intent(in) :: removeSecond

!   logical, allocatable, dimension(:) :: lremove
!   integer :: imol, jmol
!   integer :: idx
!   integer :: iatm, jatm, ineigh
!   real(real64) :: dij(3), rtmp, rmax
!   integer :: n0
!   integer :: kdx

!   allocate(lremove(numberOfMolecules), source=.false.)

!   rmax = rcut**2
!   m1 : do imol=1,numberOfMolecules
!     if (lremove(imol)) cycle
!     do idx=1,listOfMolecules(imol) % numberOfAtoms
!       iatm = listOfMolecules(imol) % listOfAtoms(idx)
!       m2 : do ineigh=1,nneigh(iatm)
!         jatm = lneigh(ineigh,iatm)
!         jmol = atomtomoleculeindex(jatm)
!         if (imol == jmol) cycle
!         if (lremove(jmol)) cycle
!         dij = frame % pos(1:3,iatm) - frame % pos(1:3,jatm)
!         rtmp = computeDistanceSquaredPBC(dij)
!         if (rtmp < rmax) then
!           kdx = imol
!           if (removeSecond) then
!             if (jmol > imol) kdx = jmol
!           else
!             if (jmol < imol) kdx = jmol
!           end if
!           lremove(kdx) = .true.
!           if (kdx == imol) then
!             cycle m1
!           else
!             cycle m2
!           end if
!         end if
!       end do m2
!     end do
!   end do m1

!   n0 = 0
!   do imol=1,numberOfMolecules
!     if (lremove(imol)) cycle
!     do idx=1,listOfMolecules(imol) % numberOfAtoms
!       iatm = listOfMolecules(imol) % listOfAtoms(idx)

!       n0 = n0 + 1
!       frame % lab(n0) = frame % lab(iatm)
!       frame % chg(n0) = frame % chg(iatm)
!       frame % pos(:,n0) = frame % pos(:,iatm)
!     end do
!   end do
!   frame % natoms = n0

! end subroutine deleteOverlapAllMolecules
