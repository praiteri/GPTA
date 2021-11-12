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
module moduleExtractFramesByProperty
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleFiles
  use moduleMessages
  use moduleProperties
  use moduleResizeArrays
  
  implicit none

  public :: extractFrames
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: inputfile
  type(fileTypeDef), pointer :: outputfile
  integer, pointer :: tallyExecutions
  real(8), pointer, dimension(:) :: limits

  integer, pointer :: nProperties
  real(8), pointer, dimension(:,:) :: properties
  logical, pointer :: header 

contains

  subroutine extractFrames(a)
    implicit none
    type(actionTypeDef), target :: a

    call associatePointers(a)

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      ! Atoms' selection
      if (firstAction) then
        ! write info 
        call dumpScreenInfo()
 
        ! Throw a warning for unused flags
        call checkUsedFlags(actionCommand)
        firstAction = .false.          
      end if

      call computeAction(a)
    end if

    if (endOfCoordinatesFiles) then
      ! call finaliseAction(a)
    end if 

  end subroutine extractFrames

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    inputfile            => a % inputfile
    outputfile           => a % outputfile
    limits(1:)           => a % doubleVariables(1:)
    nProperties          => a % integerVariables(1)
    properties           => a % array2D
    header               => a % logicalVariables(1)

  end subroutine associatePointers

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a
    character(STRLEN) :: flagString
    character(STRLEN), dimension(100) :: str
    integer :: i, nval
    
    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .true.
    a % cutoffNeighboursList = 3.0d0

    call assignFlagValue(actionCommand,"+f",inputfile % fname,'properties')
    call assignFlagValue(actionCommand,"+out",outputfile % fname,'out.pdb')
    call initialiseFile(outputFile, outputFile % fname)

    block
      integer :: io, ierr, nl, i
      character(len=200) :: line, words(10)
      open(newunit=io, file= a % inputfile % fname, status='old')

      ! count the number of lines
      nl = 0
      do
        read(io,'(a200)',iostat=ierr)line
        if (index(line,"#") > 0) cycle ! skip comments
        if (len_trim(line) == 0) cycle ! skip empty lines
        if (ierr/=0) exit ! exit if EOF or EOL
        nl = nl + 1
      end do
      rewind(io)

      ! count the number of properties in the line
      do nProperties=1,20
        read(line,*,iostat=ierr)words(1:nProperties)
        if (ierr/=0) exit
      end do
      nProperties = nProperties - 1
      allocate(a % array2D(nProperties,nl))

      read(io,'(a200)',iostat=ierr)line
      do i=1,nl
        read(io,*)a % array2D(:,i)
      end do

    end block

    call extractFlag(actionCommand,"+min",flagString)
    call parse(flagString, ",", str, nval)
    if (nval /= nProperties) call message(-1,"Number of values /= from number of columns",iv=[nval,nProperties])
    do i=1,nProperties
      if (str(i) == "*" .or. str(i) == "INF") then
        limits(i) = minval(a % array2D(i,:)) - 1.d0
      else
        read(str(i),*)limits(i)
      end if
    end do

    call extractFlag(actionCommand,"+max",flagString)
    call parse(flagString, ",", str, nval)
    if (nval /= nProperties) call message(-1,"Number of values /= from number of columns",iv=[nval,nProperties])
    do i=1,nProperties
      if (str(i) == "*" .or. str(i) == "INF") then
        limits(i+nProperties) = maxval(a % array2D(i,:)) + 1.d0
      else
        read(str(i),*)limits(i+nProperties)
      end if
    end do
    
    header = .true.
    tallyExecutions = 0 

  end subroutine initialiseAction

  subroutine finaliseAction(a)
    implicit none
    type(actionTypeDef), target :: a

  end subroutine finaliseAction

  subroutine dumpScreenInfo()
    implicit none
    integer :: i
    character(len=5) :: str
    call message(0,"Extract frames by properties")
    call message(0,"...Properties file",str=inputfile % fname)
    call message(0,"...Number of properties in file",i=nProperties)
    do i=1,nProperties
      write(str,'(i0)')i
      call message(0,"......Limits for property #"//trim(str),rv=[limits(i),limits(i+nProperties)])
    end do
    call message(0,"...Output file",str=outputfile % fname)
    call message(2)
  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleSystem

    implicit none
    type(actionTypeDef), target :: a

    integer :: i
    logical, allocatable, dimension(:) :: dumpFrame

    allocate(dumpFrame(nProperties), source=.false.)
    do i=1,nProperties
      if ( properties(i,tallyExecutions) > limits(i) .and. &
           properties(i,tallyExecutions) < limits(i+nProperties) ) dumpFrame(i) = .true.
    end do
    
    if (count(dumpFrame) == nProperties) then
      call dumpCoordinates(outputFile % ftype, &
                           outputFile % funit, &
                           frame % natoms, &
                           frame % pos, &
                           frame % lab, &
                           hmat = frame % hmat, &
                           header = header)
      header = .false.
    end if
 
  end subroutine computeAction

end module moduleExtractFramesByProperty
