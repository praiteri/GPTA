!disclaimer
module moduleAtomsAttributes
  use moduleMessages 
  
contains

  subroutine setAtomAttributesHelp()
    implicit none
    call message(0,"This action overrides some of the basic properties of the atoms in the system that can be used by other actions.")
    call message(0,"Examples:")
    call message(0,"gpta.x --i coord.pdb --top --set +s O2,H2 +l OW,HW")
    call message(0,"gpta.x --i coord.pdb --top --set +i 1:300:3 +l OW")
    call message(0,"gpta.x --i coord.pdb --top --set +mol M1 +l OW,HW,HW")
    call message(0,"gpta.x --i coord.pdb --top --set +s O2,H2 +tq 1,2")
  end subroutine setAtomAttributesHelp

  subroutine setAtomAttributes(a)
    use moduleVariables
    use moduleStrings
    use moduleSystem 
    use moduleElements

    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    logical, pointer :: actionInitialisation
    logical, pointer :: firstAction

    logical, pointer :: uniqueLabels
    logical, pointer :: elementLabels
    logical, pointer :: labelSelection
    logical, pointer :: indexSelection
    logical, pointer :: moleculeSelection
    logical, pointer :: renameAtoms
    logical, pointer :: rechargeAtoms
    logical, pointer :: tinkerAtoms

    integer :: iatm, natoms, idx, imol, itmp, jdx
    integer :: indices(3)
    integer, allocatable, dimension(:) :: atomFlag
    character(cp), allocatable, dimension(:) :: oldLabels, moleculeLabels
    character(len=40), allocatable, dimension(:) :: newLabels
    character(len=STRLEN) :: ctmp 
    character(len=STRLEN) :: str

    actionCommand        => a % actionDetails

    ! Associate variables
    actionInitialisation => a % actionInitialisation
    firstAction          => a % firstAction

    uniqueLabels         => a % logicalVariables(1) 
    elementLabels        => a % logicalVariables(2) 
    labelSelection       => a % logicalVariables(3)
    indexSelection       => a % logicalVariables(4)
    moleculeSelection    => a % logicalVariables(5)
    renameAtoms          => a % logicalVariables(6)
    rechargeAtoms        => a % logicalVariables(7)
    tinkerAtoms          => a % logicalVariables(8)

    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.

      call assignFlagValue(actionCommand,"+l",renameAtoms,.false.)
      call assignFlagValue(actionCommand,"+q",rechargeAtoms,.false.)
      call assignFlagValue(actionCommand,"+tq",tinkerAtoms,.false.)

      call assignFlagValue(actionCommand,"+s",labelSelection,.false.)
      call assignFlagValue(actionCommand,"+i",indexSelection,.false.)
      call assignFlagValue(actionCommand,"+mol",moleculeSelection,.false.)
      
      call assignFlagValue(actionCommand,"+ul",uniqueLabels,.false.)

      call assignFlagValue(actionCommand,"+el",elementLabels,.false.)

      if ( count([labelSelection,indexSelection,moleculeSelection]) /= 1) call message(-1,"--set | mixed selection not available")

      call checkUsedFlags(actionCommand)
      return
    end if

    if (frameReadSuccessfully) then
      if (firstAction) then
        allocate(atomFlag(frame % natoms), source=0)

        if (labelSelection) then
          call assignFlagValue(actionCommand,"+s",oldLabels)
          if (oldLabels(1) == "all") then
            deallocate(oldLabels)
            allocate(oldLabels(numberOfUniqueLabels), source=listOfUniqueLabels)
          end if
          do iatm=1,frame % natoms
            do idx=1,size(oldLabels)
              if (frame % lab(iatm) == oldLabels(idx)) atomFlag(iatm) = idx
            enddo
          enddo

        else if (moleculeSelection) then
          call assignFlagValue(actionCommand,"+mol",moleculeLabels)
          idx = 0
          do itmp=1,size(moleculeLabels)
            do imol=1,numberOfMolecules
              if (listOfMolecules(imol) % resname == moleculeLabels(itmp) .or. moleculeLabels(itmp) == "all") then
                natoms = listOfMolecules(imol) % numberOfAtoms
                do iatm=1,natoms
                  jdx = listOfMolecules(imol) % listOfAtoms(iatm)
                  atomFlag(jdx) = idx + iatm
                end do
              end if
            end do
            idx = idx + natoms
          end do

        else if (indexSelection) then
          call assignFlagValue(actionCommand,"+i",newLabels)
          do idx=1,size(newLabels)
            indices = extractSelectionIndices(newLabels(idx))
            do jdx=indices(1),indices(2),indices(3)
              atomFlag(jdx) = 1
            end do
          end do

        else 
          call message(-1,"--set | Unknown selection")
        end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! call message(0,"Modify atoms' Attributes")
        if (labelSelection) then
          str = trim(oldLabels(1))
          do idx=2,size(oldLabels)
            str = trim(str)//" "//trim(oldLabels(idx))
          end do
          call message(0,"...Atoms selected by label",str=str)

        else if (indexSelection) then
          str = trim(newLabels(1))
          do idx=2,size(newLabels)
            str = trim(str)//" "//trim(newLabels(idx))
          end do
          call message(0,"...Atoms selected by index",str=str)

        else if (moleculeSelection) then
          str = trim(moleculeLabels(1))
          do idx=2,size(moleculeLabels)
            str = trim(str)//" "//trim(moleculeLabels(idx))
          end do
          call message(0,"...Atoms selected by molecule",str=str)

        end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (renameAtoms) then
          allocate(a % localLabels(frame % natoms))
          call assignFlagValue(actionCommand,"+l",newLabels)
          ! if (moleculeSelection) then
          !   if (size(newLabels) /= maxval(atomFlag)) &
          !     call message(-1,"--set | +mol and +l have different number of elements",iv=[size(oldLabels),size(newLabels)])
          ! else if (indexSelection) then
          !     if (size(newLabels) /= 1) call message(-1,"--set | only one new label is allowed with +is")
          ! else
          !   if (size(newLabels) /= size(oldLabels)) then
          !     if (size(newLabels) /= 1) call message(-1,"--set | +s and +l have different number of elements",iv=[size(oldLabels),size(newLabels)])
          !     ctmp = newLabels(1)
          !     deallocate(newLabels)
          !     allocate(newLabels(size(oldLabels)), source=ctmp)
          !     newLabels = ctmp
          !   end if
          ! end if

          if (newLabels(1) == "unique") then
            do iatm=1,frame % natoms
              if (atomFlag(iatm)==0) then
                a % localLabels(iatm) = frame % lab(iatm)
              else
                write(str,"(a,i0)")trim(getElement(frame % lab(iatm))),atomFlag(iatm)
                
                a % localLabels(iatm) = trim(str)
              end if
            enddo
            call message(0,"...Giving unique names to residue atoms")

          else if (newLabels(1) == "element") then
              do iatm=1,frame % natoms
                if (atomFlag(iatm)==0) then
                  a % localLabels(iatm) = frame % lab(iatm)
                else
                  write(str,"(a)")trim(getElement(frame % lab(iatm)))
                  a % localLabels(iatm) = trim(str)
                end if
              enddo
              call message(0,"...Replacing atoms names with elements")
          else
            do iatm=1,frame % natoms
              if (atomFlag(iatm)==0) then
                a % localLabels(iatm) = frame % lab(iatm)
              else
                a % localLabels(iatm) = newLabels(atomFlag(iatm))(1:4)
              end if
            enddo
            str = trim(newLabels(1))
            do idx=2,size(newLabels)
              str = trim(str)//" "//trim(newLabels(idx))
            enddo
            call message(0,"...Renaming selected atoms as",str=str)
          end if
        end if

        if (rechargeAtoms) then
          allocate(a % localCharges(frame % natoms))
          call assignFlagValue(actionCommand,"+q",newLabels)
          if (size(newLabels) /= maxval(atomFlag)) then
            if (size(newLabels) /= 1) call message(-1,"select +s and +q have different number of elements")
            ctmp = newLabels(1)
            deallocate(newLabels)
            allocate(newLabels(size(oldLabels)))
            newLabels = ctmp
          end if

          do iatm=1,frame % natoms
            if (atomFlag(iatm)==0) then
              a % localCharges(iatm) = frame % chg(iatm)
            else
              read(newLabels(atomFlag(iatm)),*) a % localCharges(iatm) 
            end if
          enddo

          str = trim(newLabels(1))
          do idx=2,size(newLabels)
            str = trim(str)//" "//trim(newLabels(idx))
          enddo
          call message(0,"...Setting atomic charges to",str=str)
        end if

        call message(2)

        call checkUsedFlags(actionCommand)
        firstAction = .false.
      
      end if

      if (renameAtoms) frame % lab = a % localLabels

      if (rechargeAtoms) frame % chg = a % localCharges

    end if

    if (endOfCoordinatesFiles) return

  end subroutine setAtomAttributes

end module moduleAtomsAttributes
