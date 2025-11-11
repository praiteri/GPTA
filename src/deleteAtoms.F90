!disclaimer
module moduleDeleteAtoms
  use moduleMessages 

contains


  subroutine deleteAtomsHelp()
    implicit none
    call message(0,"This action deletes the selected atoms from the input configuration.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --delete +s O")
    call message(0,"  gpta.x --i coord.pdb --delete +s O,H")
    call message(0,"  gpta.x --i coord.pdb --delete +i 1:30")
    call message(0,"  gpta.x --i coord.pdb --delete +i 1:30:3")
    call message(0,"  gpta.x --i coord.pdb --delete +i 1,34,6789")
    call message(0,"  gpta.x --i coord.pdb --top --delete +mol M2")
  end subroutine deleteAtomsHelp

  subroutine deleteAtoms(a)
    use moduleVariables
    use moduleStrings
    use moduleSystem 
    use moduleNeighbours

    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    logical, pointer :: actionInitialisation
    logical, pointer :: firstAction

    logical, pointer :: labelSelection
    logical, pointer :: indexSelection
    logical, pointer :: moleculeSelection
    logical, pointer :: positionSelection

    integer :: iatm

    actionCommand        => a % actionDetails

    ! Associate variables
    actionInitialisation => a % actionInitialisation
    firstAction          => a % firstAction

    labelSelection       => a % logicalVariables(1)
    indexSelection       => a % logicalVariables(2)
    moleculeSelection    => a % logicalVariables(3)
    positionSelection    => a % logicalVariables(4)

    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.

      call assignFlagValue(actionCommand,"+s",labelSelection,.false.)
      call assignFlagValue(actionCommand,"+i",indexSelection,.false.)
      call assignFlagValue(actionCommand,"+mol",moleculeSelection,.false.)
      call assignFlagValue(actionCommand,"+pos",positionSelection,.false.)
      if ( count([labelSelection,indexSelection,moleculeSelection,positionSelection]) /= 1) then
          write(0,*)labelSelection,indexSelection,moleculeSelection,positionSelection
          call message(-1,"--set | mixed selection not available")
      endif

      call checkUsedFlags(actionCommand)
      return
    end if

    if (frameReadSuccessfully) then

      ! Atoms' selection
      if (firstAction) then
        call message(0,"Deleting Atoms")

        if (resetFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if
        call selectAtoms(1,actionCommand,a)
        ! call createInvertedSelectionList(a,1)

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      ! Repeat selection for reactive trajectories
      else
        if (a % updateAtomsSelection) call selectAtoms(1,actionCommand,a)
      end if

      block
        integer ::  icentre
        icentre = 0
        do iatm=1,frame % natoms
          if (a % isSelected(iatm, 1)) cycle
          icentre = icentre + 1
          frame % pos(1:3, icentre) = frame % pos(1:3, iatm)
          frame % frac(1:3,icentre) = frame % frac(1:3,iatm)
          frame % lab(     icentre) = frame % lab(     iatm)
          frame % chg(     icentre) = frame % chg(     iatm)
          frame % element( icentre) = frame % element( iatm)
        end do
        frame % natoms = icentre
      end block
      
      if (computeNeighboursList) call updateNeighboursList(.true.)

    end if

    if (endOfCoordinatesFiles) return

  end subroutine deleteAtoms

end module moduleDeleteAtoms
