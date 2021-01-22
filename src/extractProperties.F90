! Copyright (c) 2021, Paolo Raiteri, Curtin University.
! All rights reserved.
! 
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License as published by the 
! Free Software Foundation; either version 3 of the License, or 
! (at your option) any later version.
!  
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! * Redistributions of source code must retain the above copyright notice, 
!   this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright notice, 
!   this list of conditions and the following disclaimer in the documentation 
!   and/or other materials provided with the distribution.
! * Neither the name of the <ORGANIZATION> nor the names of its contributors 
!   may be used to endorse or promote products derived from this software 
!   without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
module moduleExtractSystemProperties
  use moduleVariables
  use moduleFiles
  use moduleProperties
  use moduleMessages 
  use moduleStrings
  use moduleSystem 

  implicit none

  public :: extractSystemProperties
  private

  character(:), pointer :: actionCommand
  type(fileTypeDef), pointer :: outputFile
  logical, pointer :: firstAction

  logical, pointer :: dumpCoordinates
  type(fileTypeDef), pointer :: coordinatesFile

  integer, pointer :: ID

  integer, pointer :: numberOfBins
  integer, pointer :: numberOfSelectedAtoms
  integer, dimension(:), pointer :: numberOfBins2D
  real(8), dimension(:), pointer :: distributionLimits

  logical, pointer :: atomSelection        ! Property requires atoms' selection
  logical, pointer :: dumpProperty         ! Output dump
  logical, pointer :: averageMultiProperty ! Average
  logical, pointer :: distProperty         ! Distribution

  integer, pointer :: numberOfProperties
  real(8), pointer, dimension(:) :: localProperty

  interface
    subroutine stuff(np, property, numberOfSelectedAtoms, selectionList)
      integer, intent(inout) :: np
      real(8), dimension(*), intent(out) :: property
      integer, intent(in), optional :: numberOfSelectedAtoms
      integer, dimension(:), intent(in), optional :: selectionList
    end subroutine stuff
  end interface
  procedure(stuff) :: extractcell, extractHMatrix, extractVolume, extractDensity
  procedure(stuff) :: computeInertiaTensor
  procedure(stuff) :: dielectricConstant
  procedure(stuff) :: bondLength
  procedure(stuff) :: position
  procedure(stuff) :: test

contains

  subroutine extractSystemProperties(a)
    implicit none
    type(actionTypeDef), target :: a

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      
      if (firstAction) then       

        call message(0,"System property")
        call message(0,"...Extracting "//trim(a % sysprop % name))
        
        if (atomSelection) then
          if (keepFrameLabels) then
            a % updateAtomsSelection = .false.
          else
            a % updateAtomsSelection = .true.
          end if
         
          call selectAtoms(1,actionCommand,a)
          call createSelectionList(a,1)
          numberOfSelectedAtoms = count(a % isSelected(:,1))
        end if
        
        call checkUsedFlags(actionCommand)
        firstAction = .false.

        ! get the number of properties computed and initialise the working arrays
        call initialiseCompute(a)

      else

        if (atomSelection) then
          if (a % updateAtomsSelection) then
            call selectAtoms(1,actionCommand,a)
            call createSelectionList(a,1)
            numberOfSelectedAtoms = count(a % isSelected(:,1))
          end if
        end if

      end if

      call computeAction(numberOfProperties, a)
    end if

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) then
      if (outputFile % fname /= "NULL") call message(0,"......Output file",str=outputFile % fname)
      if (dumpCoordinates     ) call message(0,"......Output file for scaled coordinates",str=coordinatesFile % fname)
      if (averageMultiProperty) call message(0,"......Average(s) and standard deviation(s)")
      if (distProperty        ) call message(0,"......Distribution")
      call finaliseAction()
      call message(2)
    end if

  end subroutine extractSystemProperties

  subroutine finaliseAction()
    use moduleDistances
    use moduleFiles
    use moduleMessages

    implicit none
    real(8), allocatable, dimension(:) :: cell

    call workData % dump(ID, outputFile % funit)
    call message(4)

    if (dumpCoordinates) then
      call workData % extract(ID, cell)
      frame % hmat = reshape(cell,[3,3])
      call fractionalToCartesian(frame % natoms, frame % frac, frame % pos)
      if (numberOfMolecules > 0) call reassembleAllMolecules()

      call runInternalAction("dump",trim(coordinatesFile % fname)//" +internal")

      close (coordinatesFile % funit)
    end if

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

    if (ntmp>0) call workData % compute(ID, numberOfValues=numberOfProperties, xValues=localProperty)

  end subroutine computeAction

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .false.

    ! Local pointers
    actionCommand           => a % actionDetails
    firstAction             => a % firstAction
    outputFile              => a % outputFile
    
    numberOfBins            => a % numberOfBins
    numberOfBins2D(1:2)     => a % numberOfBins2D
    
    ID                      => a % integerVariables(1)
    numberOfProperties      => a % integerVariables(2)
    numberOfSelectedAtoms   => a % integerVariables(3)

    distributionLimits(1:2) => a % doubleVariables(1:2)
    localProperty(1:)       => a % doubleVariables(3:)

    atomSelection           => a % logicalVariables(1) ! Atoms selection for property
    dumpProperty            => a % logicalVariables(2) ! Output dump
    averageMultiProperty    => a % logicalVariables(3) ! Average more than one property
    distProperty            => a % logicalVariables(4) ! Distribution

    dumpCoordinates         => a % logicalVariables(5) 
    coordinatesFile         => a % auxiliaryFile 

    numberOfSelectedAtoms = 0
    atomSelection = .false. 
    
    call getPropertyType(a)

    ! Extract to do with the properties
    if (.not. averageMultiProperty) then
      call assignFlagValue(actionCommand,"+avg " ,averageMultiProperty, .false.) ! Average
      call assignFlagValue(actionCommand,"+dist ",distProperty,         .false.) ! Distribution
      call assignFlagValue(actionCommand,"+dump ",dumpProperty,         .false.) ! Dump on same line
    end if

    if ( count([distProperty,averageMultiProperty]) == 0 ) dumpProperty = .true.

    call assignFlagValue(actionCommand,"+out",outputFile % fname,"NULL")
    if (outputFile % fname /= "NULL") then
      call initialiseFile(outputFile, outputFile % fname)
      write(outputFile % funit,"(a)")"#System property"
      write(outputFile % funit,'(a)')"#...Extracting "//trim(a % sysprop % name)
      if (averageMultiProperty) write(outputFile % funit,'(a)')"#......Average and standard deviation"
      if (distProperty        ) write(outputFile % funit,'(a)')"#......Distribution"
    else
      outputFile % funit = 6
    end if
    
  end subroutine initialiseAction

  subroutine initialiseCompute(a)
    implicit none
    type(actionTypeDef), target :: a

    numberOfProperties = -1
    call computeAction(numberOfProperties, a)

    if (dumpProperty) then 
      call workData % initialise(ID, "dump", iounit=outputFile % funit)
      
    else if (averageMultiProperty) then
      call workData % initialise(ID, "multiAverage", numberOfValues=numberOfProperties)

    else if (distProperty) then
      call assignFlagValue(actionCommand,"+nbin",numberOfBins,100)
      call assignFlagValue(actionCommand,"+dist ", distributionLimits(1:2),[1.d0,2.d0])
      call workData % initialise(ID, "prob1D", numberOfBins=[numberOfBins], lowerLimits=[distributionLimits(1)], upperLimits=[distributionLimits(2)])
    end if

  end subroutine initialiseCompute

  subroutine getPropertyType(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: n
         
    if (flagExists(actionCommand,"+test")) then
      n = n +1
      atomSelection = .true.
      a % sysprop % property => test
      a % sysprop % name = "Test Extract Routine"
    end if
         
    if (flagExists(actionCommand,"+pos")) then
      n = n +1
      atomSelection = .true.
      a % sysprop % property => position
      a % sysprop % name = "Atoms Positions"
    end if

    if (flagExists(actionCommand,"+cell")) then
      n = n +1
      a % sysprop % property => extractCell
      a % sysprop % name = "Cell Parameters"
    end if
    
    if (flagExists(actionCommand,"+hmat")) then
      n = n +1
      a % sysprop % property => extractHMatrix
      a % sysprop % name = "Cell Matrix"
    end if
    
    if (flagExists(actionCommand,"+system")) then
      n = n +1
      averageMultiProperty = .true.
      dumpCoordinates = .true.
      a % sysprop % property => extractHMatrix
      a % sysprop % name = "Average cell and scaled coordinates"
      if (countFlagArguments(actionCommand,"+system") == 0) then
        coordinatesFile % fname = "averageSystem.pdb"
      else 
        call assignFlagValue(actionCommand,"+system",coordinatesFile % fname,"NULL")
      end if
    end if
    
    if (flagExists(actionCommand,"+volume")) then
      n = n +1
      a % sysprop % property => extractVolume
      a % sysprop % name = "Cell Volume"
    end if
    
    if (flagExists(actionCommand,"+inertia")) then
      n = n +1
      a % sysprop % property => computeInertiaTensor
      a % sysprop % name = "Intertia tensor (eigenvalues)"
    end if

    if (flagExists(actionCommand,"+diel")) then
      n = n +1
      a % sysprop % property => dielectricConstant
      a % sysprop % name = "Dielectric constant"
    end if
    
    if (flagExists(actionCommand,"+density")) then
      n = n +1
      a % sysprop % property => extractDensity
      a % sysprop % name = "Density"
    end if

    if (flagExists(actionCommand,"+distance")) then
      n = n +1
      atomSelection = .true.
      a % sysprop % property => bondLength
      a % sysprop % name = "Distance between atoms"
    end if

    if (n==0)  call message(-1,"--extract | no property selected")
    if (n/=1)  call message(-1,"--extract | only one property per action can be specified")
    
  end subroutine getPropertyType

end module moduleExtractSystemProperties

subroutine extractCell(n, cell)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(8), dimension(*), intent(out) :: cell
  
  ! set the number of properties computed
  if (n<0) then
    n = 6
    return
  end if
  cell(1:6) = frame % cell

end subroutine extractCell

subroutine extractHMatrix(n, hmat)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(8), dimension(*), intent(out) :: hmat
  
  ! set the number of properties computed
  if (n<0) then
    n = 9
    return
  end if
  hmat(1:9) = reshape(frame % hmat, [9])

end subroutine extractHMatrix

subroutine extractVolume(n, volume)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(8), dimension(*), intent(out) :: volume
  
  ! set the number of properties computed
  if (n<0) then
    n = 1
    return
  end if
  volume(1) = frame % volume

end subroutine extractVolume

subroutine extractDensity(n, density)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(8), dimension(*), intent(out) :: density
  
  ! set the number of properties computed
  if (n<0) then
    n = 1
    return
  end if
  density(1) = 1.6605388d0 * totalMass / frame % volume

end subroutine extractDensity

subroutine dielectricConstant(n, eps)
  use moduleSystem
  implicit none
  integer, intent(inout) :: n
  real(8), dimension(*), intent(out) :: eps
  
  real(8), parameter :: e0 = 0.00552642629948221d0 ! e^2/eV/angstrom
  real(8), parameter :: kbt = 0.025851d0 ! in eV
  
  integer, save :: np = 0
  real(8), dimension(3), save :: m_avg = 0.d0
  real(8), save :: m2_avg = 0.d0
  real(8), dimension(3) :: cellDipole
  real(8) :: norm, a3(3), a, b
  
  integer :: iatm
  
  ! set the number of properties computed
  if (n<0) then
    n = 1
    return
  end if
  
  cellDipole = 0.d0
  do iatm=1,frame % natoms
    cellDipole = cellDipole + frame % chg(iatm) * frame % pos(:,iatm) 
  end do
  ! if wanted in Debye ! cellDipole = cellDipole * 4.8032513d0    
  
  m_avg = m_avg + cellDipole
  m2_avg = m2_avg + sum(cellDipole**2)
  
  np = np + 1
  
  norm = 3.d0* e0 * kbt * frame % volume
  a3 = m_avg / np
  a = sum(a3*a3)
  b = m2_avg / np
  eps(1) = 1.d0 + (b-a) / norm
  eps(2) = 1.d0 + (m2_avg) / np / norm
  
end subroutine dielectricConstant

subroutine bondLength(n, val, nat, l)
  use moduleSystem
  use moduleDistances
  use moduleMessages
  implicit none
  integer, intent(inout) :: n
  real(8), dimension(*), intent(out) :: val
  integer, intent(in), optional :: nat
  integer, dimension(:), intent(in) :: l
  
  real(8), dimension(3) :: dij
  
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
  real(8), dimension(*), intent(out) :: val
  integer, intent(in), optional :: nat
  integer, dimension(:), intent(in) :: l
  
  integer :: i, ntmp

  ! set the number of properties computed
  if (n < 0) then
    n = nat*3
    return
  end if
  
  ntmp = 0
  do i=1,nat
    val(ntmp+1:ntmp+3) = frame % pos(:,l(i))
    ntmp = ntmp + 3
  end do
  
end subroutine position
  
subroutine test(n, val, nat, l)
  use moduleSystem
  use moduleDistances
  implicit none
  integer, intent(inout) :: n
  real(8), dimension(*), intent(out) :: val
  integer, intent(in), optional :: nat
  integer, dimension(:), intent(in) :: l
  
  ! set the number of properties computed
  if (n < 0) then
    n = 3
    return
  end if
  
  val(1) = 1.d0
  val(2) = 2.d0
  val(3) = 3.d0
  
end subroutine
