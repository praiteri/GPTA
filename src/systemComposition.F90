!disclaimer
subroutine systemComposition(localFrame)
  use moduleVariables, only : fileTypeDef
  use moduleSystem , only : cp, frameTypeDef, numberOfUniqueLabels, listOfUniqueLabels, listOfElements
  use moduleMessages 
  use moduleElements, only : getElement
  implicit none

  type(frameTypeDef), intent(in) :: localFrame

  integer :: i, iatm
  integer :: nspecies
  integer, parameter :: maxSpecies = 1000
  integer, dimension(maxSpecies) ::  ncount
  character(cp), dimension(maxSpecies) :: species
  character(2), dimension(maxSpecies) :: element
  character(5)  :: str


  call message(0,"System composition")
  call message(0,'...Number of atoms found',i=localFrame % natoms)

  ncount = 0
  nspecies = 1
  species(1) = localFrame % lab(1)
  element(1) = localFrame % element(1)

  ncount(1) = 1
  do iatm=2,localFrame % natoms
    if (all(species(1:nspecies) /= localFrame % lab(iatm))) then
      nspecies = nspecies + 1
      species(nspecies) = localFrame % lab(iatm)
      element(nspecies) = localFrame % element(iatm)

    end if
    where (species(1:nspecies) == localFrame % lab(iatm)) ncount(1:nspecies) = ncount(1:nspecies) + 1
  enddo
  do i=1,nspecies
    element(i) = getElement(species(i))
    write(str,'("[",a2,"] ")') adjustl(element(i))
    call message(0,"......Number of unique species "//str//trim(species(i)),i=ncount(i))
  enddo

  call message(3)

  numberOfUniqueLabels = nspecies
  allocate(listOfUniqueLabels(numberOfUniqueLabels))
  allocate(listOfElements(numberOfUniqueLabels))
  listOfUniqueLabels = species(1:numberOfUniqueLabels)
  listOfElements = element(1:numberOfUniqueLabels)

end subroutine systemComposition
