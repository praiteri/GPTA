!disclaimer
module XMLtopology
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleFiles
  use moduleProperties
  use moduleMessages
  use moduleDistances
  use moduleElements  

  implicit none

  public :: writeXML
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions
  
contains

  subroutine initialiseAction(a)

    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.2_real64

    ! get output file name from the command line, if present
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'test.out')
    
  end subroutine initialiseAction

  subroutine finaliseAction()
    implicit none
    close(outputFile % funit)
  end subroutine finaliseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0, "Topology XML template ")
  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: funit = 0
    integer :: i, j, k, kdx, n, m, n0, iatm
    integer :: idx(4)
    character(len=cp) :: l1,l2,l3,l4,e
    real(real64) :: mass

    integer :: nlist
    character(cp), allocatable, dimension(:,:) :: list
    character(cp), dimension(4) :: labels
    character(len=200) :: str, ldx
    integer, allocatable, dimension(:,:) :: blist

    call computeConnectivity()

    write(funit,'(a/)') "<ForceField>"
    write(funit,'(a)') " <AtomTypes>"
    do i=1,numberOfUniqueLabels
      mass = getElementMass(listOfUniqueLabels(i))
      l1 = listOfUniqueLabels(i)
      e = getElement(l1)
      write(str,'("  <Type name=""",a,""" class=""",a,""" element=""",a,""" mass=""",f8.4"""/>")') trim(l1),trim(l1),trim(e),mass
      write(funit,'(a)') trim(str)
    enddo
    write(funit,'(a/)')  " </AtomTypes>"

    write(funit,'(a)') " <Residues>"
    do k=1,numberOfUniqueMoleculeTypes
      write(str,'("  <Residue name=""",a,"""/>")') trim(listOfMolecules(k) % resname)
      write(funit,'(a)') trim(str)
      kdx = uniqueMoleculesID(k)
      do i=1,listOfMolecules(kdx) % numberOfAtoms
        l1 = listOfMolecules(kdx) % listOfLabels(i)
        write(ldx,'(i0)')i-1
        write(str,'("   <Atom name=""",a,""" type=""",a,""" charge=""000000""/>""")')trim(l1)//trim(ldx),trim(l1)
        write(funit,'(a)') trim(str)
      end do

      allocate(blist(2,sum(numberOfCovalentBondsPerAtom)))
      n0 = listOfMolecules(kdx) % listOfAtoms(1)
      m = 0
      do iatm=1,listOfMolecules(kdx) % numberOfAtoms
        i = listOfMolecules(kdx) % listOfAtoms(iatm)
        n = numberOfCovalentBondsPerAtom(i)
        do j=1,numberOfCovalentBondsPerAtom(i)
          m = m + 1
          blist(:,m) = [i-n0,listOfCovalentBondsPerAtom(j,i)-n0]
        end do
      end do

      n = 1
      do i=2,m
        do j=1,i-1
          if ( (blist(1,j) == blist(1,i)) .and. (blist(2,j) == blist(2,i)) ) exit
          if ( (blist(2,j) == blist(1,i)) .and. (blist(1,j) == blist(2,i)) ) exit
        end do
        if (j==i) then
          n = n + 1
          blist(:,n) = blist(:,i)
        end if
      end do
      
      do i=1,n 
        write(str,'("   <Bond from=""",i0,""" to=""",i0,"""/>""")')blist(:,i)
        write(funit,'(a)') trim(str)
      end do
      deallocate(blist)
      write(funit,'(a)') "  </Residue>"

    end do

    write(funit,'(a/)') " </Residues>"

    ! Bonds
    if (numberOfUniqueBonds > 0) then

      allocate(list(2,numberOfUniqueBonds))
      nlist = 1
      labels(1:2) = [frame%lab(listOfUniqueBonds(1,nlist)), \
                    frame%lab(listOfUniqueBonds(2,nlist))]
      list(1:2,nlist) = labels(1:2)
      ax : do i=2,numberOfUniqueBonds
        labels(1:2) = [frame%lab(listOfUniqueBonds(1,i)), \
                      frame%lab(listOfUniqueBonds(2,i))]

        do j=1,nlist
          if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) ) cycle ax
          if ( (list(1,j) == labels(2)) .and. (list(2,j) == labels(1)) ) cycle ax
        end do
        nlist = nlist + 1 
        labels(1:2) = [frame%lab(listOfUniqueBonds(1,i)), \
                      frame%lab(listOfUniqueBonds(2,i))]
        list(1:2,nlist) = labels(1:2)
      end do ax

      write(funit,'(a)') " <HarmonicBondForce>"
      do i=1,nlist
        write(str,'("  <Bond type1=""",a,"""  type2=""",a,"""  k=""000000"" length=""000000"" />""")') (trim(list(j,i)),j=1,2)
        write(funit,'(a)') trim(str)
      enddo
      deallocate(list)
      write(funit,'(a/)') " </HarmonicBondForce>"
    end if

    ! Angles
    if (numberOfUniqueAngles > 0) then
      allocate(list(3,numberOfUniqueAngles))

      nlist = 1
      labels(1:3) = [frame%lab(listOfUniqueAngles(1,nlist)), \
                    frame%lab(listOfUniqueAngles(2,nlist)), \
                    frame%lab(listOfUniqueAngles(3,nlist))]
      list(1:3,nlist) = labels(1:3)
      b : do i=2,numberOfUniqueAngles
        labels(1:3) = [frame%lab(listOfUniqueAngles(1,i)), \
                      frame%lab(listOfUniqueAngles(2,i)), \
                      frame%lab(listOfUniqueAngles(3,i))]

        do j=1,nlist
          if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(3)) ) cycle b
          if ( (list(1,j) == labels(3)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(1)) ) cycle b
        end do
        nlist = nlist + 1 
        labels(1:3) = [frame%lab(listOfUniqueAngles(1,i)), \
                      frame%lab(listOfUniqueAngles(2,i)), \
                      frame%lab(listOfUniqueAngles(3,i))]
        list(1:3,nlist) = labels(1:3)
      end do b

      write(funit,'(a)') " <HarmonicAngleForce>"
      do i=1,nlist
        write(str,'("  <Angle type1=""",a,""" type2=""",a,""" type3=""",a,""" k=""000000"" angle=""000000"" />""")') (trim(list(j,i)),j=1,3)
        write(funit,'(a)') trim(str)
      enddo
      deallocate(list) 
      write(funit,'(a/)') " </HarmonicAngleForce>"
    end if

    ! Torsions
    if (numberOfUniqueTorsions > 0) then
      allocate(list(4,numberOfUniqueTorsions))

      nlist = 1
      labels(1:4) = [frame%lab(listOfUniqueTorsions(1,nlist)), \
                    frame%lab(listOfUniqueTorsions(2,nlist)), \
                    frame%lab(listOfUniqueTorsions(3,nlist)), \
                    frame%lab(listOfUniqueTorsions(4,nlist))]
      list(1:4,nlist) = labels(1:4)
      c : do i=2,numberOfUniqueTorsions
        labels(1:4) = [frame%lab(listOfUniqueTorsions(1,i)), \
                      frame%lab(listOfUniqueTorsions(2,i)), \
                      frame%lab(listOfUniqueTorsions(3,i)), \
                      frame%lab(listOfUniqueTorsions(4,i))]

        do j=1,nlist
          if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(3)) .and. (list(4,j) == labels(4)) ) cycle c
          if ( (list(1,j) == labels(4)) .and. (list(2,j) == labels(3)) .and. (list(3,j) == labels(2)) .and. (list(4,j) == labels(1)) ) cycle c
        end do
        nlist = nlist + 1 
        labels(1:4) = [frame%lab(listOfUniqueTorsions(1,i)), \
                      frame%lab(listOfUniqueTorsions(2,i)), \
                      frame%lab(listOfUniqueTorsions(3,i)), \
                      frame%lab(listOfUniqueTorsions(4,i))]
        list(1:4,nlist) = labels(1:4)
      end do c

      write(funit,'(a)') " <RBTorsionForce>"
      do i=1,nlist
        write(str,'(" <Proper class1=""",a,""" class2=""",a,""" class3=""",a,""" class4=""",a,""" c0=0.0 c1=0.0 c2=0.0 c3=0.0 c4=0.0 c5=0.0 />""")') (trim(list(j,i)),j=1,4)
        write(funit,'(a)') trim(str)
      enddo
      deallocate(list) 
      write(funit,'(a/)') " </RBTorsionForce>"
    end if

    ! Impropers
    if (numberOfUniqueOutOfPlane > 0) then
      allocate(list(4,numberOfUniqueOutOfPlane))

      nlist = 1
      labels(1:4) = [frame%lab(listOfUniqueOutOfPlane(1,nlist)), \
                    frame%lab(listOfUniqueOutOfPlane(2,nlist)), \
                    frame%lab(listOfUniqueOutOfPlane(3,nlist)), \
                    frame%lab(listOfUniqueOutOfPlane(4,nlist))]
      list(1:4,nlist) = labels(1:4)
      d : do i=2,numberOfUniqueOutOfPlane
        labels(1:4) = [frame%lab(listOfUniqueOutOfPlane(1,i)), \
                      frame%lab(listOfUniqueOutOfPlane(2,i)), \
                      frame%lab(listOfUniqueOutOfPlane(3,i)), \
                      frame%lab(listOfUniqueOutOfPlane(4,i))]

        do j=1,nlist
          if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(3)) .and. (list(4,j) == labels(4)) ) cycle d
          if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(2)) .and. (list(3,j) == labels(4)) .and. (list(4,j) == labels(3)) ) cycle d
          if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(3)) .and. (list(3,j) == labels(2)) .and. (list(4,j) == labels(4)) ) cycle d
          if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(3)) .and. (list(3,j) == labels(4)) .and. (list(4,j) == labels(2)) ) cycle d
          if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(4)) .and. (list(3,j) == labels(2)) .and. (list(4,j) == labels(3)) ) cycle d
          if ( (list(1,j) == labels(1)) .and. (list(2,j) == labels(4)) .and. (list(3,j) == labels(3)) .and. (list(4,j) == labels(2)) ) cycle d
        end do
        nlist = nlist + 1 
        labels(1:4) = [frame%lab(listOfUniqueOutOfPlane(1,i)), \
                      frame%lab(listOfUniqueOutOfPlane(2,i)), \
                      frame%lab(listOfUniqueOutOfPlane(3,i)), \
                      frame%lab(listOfUniqueOutOfPlane(4,i))]
        list(1:4,nlist) = labels(1:4)
      end do d

      write(funit,'(a)') " <PeriodicTorsionForce>"
      do i=1,nlist
        write(str,'("  <Improper class1=""",a,""" class2=""",a,""" class3=""",a,""" class4=""",a,""" periodicity1=0 phase1=""000000"" k1=""000000""/>""")') (trim(list(j,i)),j=1,4)
        write(funit,'(a)') trim(str)
      enddo
      deallocate(list) 
      write(funit,'(a/)') " </PeriodicTorsionForce>"
    end if

    write(funit,'(a)') " <NonbondedForce coulomb14scale=1.0 lj14scale=0.5>"
    write(funit,'(a)') "  <UseAttributeFromResidue name='charge'/>"""
    do i=1,numberOfUniqueLabels
      mass = getElementMass(listOfUniqueLabels(i))
      l1 = listOfUniqueLabels(i)
      e = getElement(l1)
      write(str,'("  <Atom type=""",a,""" sigma=""1.0"" epsilon=""0.0""/>""")') trim(l1)
      write(funit,'(a)') trim(str)
    enddo

    write(funit,'(a/)') " </NonbondedForce>"

    write(funit,'(a)') " <Script>"

    write(funit,'(a/)') " </Script>"

    write(funit,'(a)') "</ForceField>"

  end subroutine computeAction

  subroutine writeXML(a)
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
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if

      call computeAction(a)
    end if

    if (endOfCoordinatesFiles) then
      call finaliseAction()
    end if 

  end subroutine writeXML

end module XMLtopology
