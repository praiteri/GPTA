!disclaimer
module moduleCreateSurface
  use moduleVariables
  use moduleDistances
  use moduleFiles
  use moduleMessages
  use moduleSystem
  use moduleStrings
  use m_mrgrnk
  
  implicit none

  public :: createSurface, createSurfaceHelp
  private

  character(:), pointer :: actionCommand
  type(fileTypeDef), pointer :: outputFile
  logical, pointer :: firstAction
  logical, pointer :: computeDipole

  integer, pointer, dimension(:) :: imiller

contains

  subroutine createSurfaceHelp()
    implicit none
    call message(0,"This action creates a new system unit cell oriented with the chosen Miller index parallel to the z axis.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --surface +hkl 1,0,4 +out surface.pdb")
  end subroutine createSurfaceHelp

  subroutine initialiseAction(a)

    implicit none
    type(actionTypeDef), target :: a

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0_real64

    ! Local pointers
    actionCommand          => a % actionDetails
    firstAction            => a % firstAction

    outputFile             => a % outputFile
    imiller(1:3)           => a % integerVariables(1:3)
    computeDipole          => a % logicalVariables(1)

    call assignFlagValue(actionCommand,"+out",outputFile % fname,'surface.pdb')
    call initialiseFile(outputFile, outputFile % fname)

    call assignFlagValue(actionCommand,"+hkl",imiller,[0,0,0])

    call assignFlagValue(actionCommand,"+dipole",computeDipole,.false.)
  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0,"Create Surface")
    call message(0,"...Miller indices",iv=imiller)

  end subroutine dumpScreenInfo

  subroutine computeAction()
    implicit none
    integer :: ii, idx
    integer :: iatm
    character(50) :: str
    real(real64), dimension(3,3) :: hmat_new
    real(real64), dimension(3,3) :: hmat_fin
    real(real64), pointer, dimension(:,:) :: pos_fin, frac_fin
    real(real64), dimension(3) :: shift
    real(real64), dimension(3) :: dipole
    real(real64) :: rtmp

    ! shift the atoms to avoid numeriacal error for the atoms on the box sides
    shift = [0.0_real64,0.0_real64,0.0_real64]
    do iatm=1,frame%natoms
      do idx=1,3
        if (abs(frame%pos(idx,iatm)) < 1e-5_real64) shift(idx) = 1e-4_real64
      end do
    end do
    do iatm=1,frame%natoms
      frame%pos(:,iatm) = frame%pos(:,iatm) + shift
    end do

    ! Determine the vector normal to the surface from the reciprocal space vectors
    call defineSurfaceVectors(imiller, frame % hmat, hmat_new)

    call create_new_cell(hmat_new)

    call straightenCell(frame % hmat,"z",hmat_new)
    frame % hmat = hmat_new
    call getInverseCellMatrix(frame % hmat, frame % hinv, frame % volume)
    call hmat2cell (frame % hmat, frame % cell, "DEG")

    allocate(pos_fin(3,frame % natoms))
    allocate(frac_fin(3,frame % natoms))

    pos_fin = frame % pos
    ! call dumpCoordinates(outputFile % ftype, outputFile % funit, &
    !   frame % natoms, pos_fin, frame % lab, frame % hmat)
    ! stop

    do ii=1,numberOfMolecules
      shift(1:3) = listOfMolecules(ii) % centreOfMass(1:3) 

      do iatm=1,frame % natoms
        pos_fin(1:3,iatm) = frame % pos(1:3,iatm) - shift(1:3) -0.01 + frame % hmat(:,3)
      enddo
      call cartesianToFractional(frame % natoms, pos_fin, frac_fin)
      call fractionalToCartesian(frame % natoms, frac_fin, pos_fin)
      call reassembleAllMolecules2(frame % natoms, pos_fin)

      ! Calculate the cell dipole moment
      if (computeDipole) then
        dipole=0.0_real64
        do iatm=1,frame % natoms
          dipole = dipole + frame % chg(iatm)*(pos_fin(1:3,iatm))
        enddo
        where(abs(dipole)<1e-5) dipole=0.0_real64

        ! Convert to Debye when writing  
        write(str,'("  Cell dipole (Debye) for configuration ",i0)')ii
        call message(0,str,rv=dipole)

        if (abs(dipole(3))>1e-2) cycle
      end if

      call dumpCoordinates(outputFile % ftype, outputFile % funit, &
        frame % natoms, pos_fin, frame % lab, frame % hmat)
    enddo
    call message(2)
  end subroutine computeAction

  subroutine createSurface(a)
    implicit none
    type(actionTypeDef), target :: a

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then

      if (firstAction) then
        if (numberOfMolecules == 0) call runInternalAction("topology","NULL")

        call dumpScreenInfo()
        call computeAction()
  
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if

    end if

  end subroutine createSurface

end module moduleCreateSurface
