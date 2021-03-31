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
module moduleDumpCoordinates

#ifdef GPTA_XDR
  use moduleXDR, only: xtcfile, trrfile
#endif

contains

  subroutine writeCoordinates(a)
    use moduleVariables
    use moduleStrings
    use moduleMessages 
    use moduleSystem 
    use moduleNeighbours 
    use moduleFiles

    use moduleLammps

    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    logical, pointer :: actionInitialisation
    logical, pointer :: firstAction
    type(fileTypeDef), pointer :: outputFile
    logical :: lappend

    character(len=STRLEN), pointer :: forcefieldFile
    logical, pointer :: header

    integer :: i, nw
    character(len=STRLEN), dimension(MAXWORDS) :: words  
    character(len=STRLEN) :: str

    logical, pointer :: internalOutput

#ifdef GPTA_XDR
    type(xtcfile), save :: xtc_out
#endif

    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    firstAction          => a % firstAction
    outputFile           => a % outputFile
    forcefieldFile       => a % stringVariables(1)
    header               => a % logicalVariables(1)
    internalOutput       => a % logicalVariables(2)

    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      a % requiresNeighboursListDouble = .false.
      a % cutoffNeighboursList = 3.0d0

      ! output file name must immediately follow the --o command
      outputFile % fname = "NULL"
      call parse(actionCommand," ",words,nw)
      if (nw == 0) call message(-1,"--o : missing output filename")
      outputFile % fname = trim(words(1))

      ! remove the file name from the list of flags
      actionCommand = ''
      do i=2,nw
        actionCommand = trim(actionCommand)//" "//trim(words(i))
      end do

      if (outputFile % fname(1:1) == "+") call message(-1,"--o : output filename cannot start with '+'",str=outputFile % fname)

      call assignFlagValue(actionCommand,"+append ",lappend,.false.)
      if (lappend) then
        call initialiseFile(outputFile, outputFile % fname, fstatus='unknown', fposition='append')
        header = .false. ! Flag for formats that require a header
      else
        call initialiseFile(outputFile, outputFile % fname, fstatus='unknown', fposition='rewind')
        header = .true. ! Flag for formats that require a header
      endif
      if (len_trim(a % name) > 3) outputFile % ftype = a % name(4:3+fp)

      ! variable to force output from other actions
      call assignFlagValue(actionCommand,"+internal ",internalOutput,.false.)

#ifdef GPTA_XDR
      if (outputFile % ftype == 'xtc') then
        close(outputFile % funit)
        call xtc_out % init(trim(outputFile % fname),'w')
      endif
#endif
      
      if (outputFile % ftype == "lmp") then
      
        ! LAMMPS forcefield file
        call assignFlagValue(actionCommand,"+f ",forcefieldFile,'forcefield.lmp')

        ! flags for missing forcefield terms
        call assignFlagValue(actionCommand,"+ignore ",str,'check')
        call parse(str,",",words,nw)
        do i=1,nw
          select case(words(i))
            case ('check')
              checkMissingTerms = .true.
            case ('all')
              ignoreMissingBonds = .true.
              ignoreMissingAngles = .true.
              ignoreMissingTorsions = .true.
              ignoreMissingImpropers = .true.
            case ('bond', 'bonds')
              ignoreMissingBonds = .true.
            case ('angle', 'angles')
              ignoreMissingAngles = .true.
            case ('torsion', 'torsions')
              ignoreMissingTorsions = .true.
            case ('improper', 'impropers')
              ignoreMissingImpropers = .true.
            case default
              call message(-1,"LAMMPS forcefield - Unknown +ignore option",str=str)
          end select
        end do
      end if

      call assignFlagValue(actionCommand,"+bohr ",outputCoordInBohr,.false.)

      ! The neigbours' list is require for the topology
      if (any(outputFile % ftype == ["lmp ","pdb2","psf ","arc "])) then
        a % requiresNeighboursList = .true.
      end if

      call checkUsedFlags(actionCommand)

      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully .or. internalOutput) then

      if (firstAction) then
        call message(1,"Writing coordinates in "//trim(outputFile % ftype)//" format",str=trim(outputFile % fname))
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if

      ! Compute topology and connectivity before convestion to Bohr, if required
      if (any(outputFile % ftype == ["lmp ","pdb2","psf ","arc "])) then
        if (numberOfMolecules == 0) then
          call runInternalAction("topology","NULL")
        end if
        call computeConnectivity() 
      end if

      !! --> define local frame here <--- !!

      ! Convert coordinates to Bohr
      if (outputCoordInBohr) then
        frame % hmat = frame % hmat / rbohr
        frame % pos = frame % pos / rbohr
      end if

      if (outputFile % ftype == "xyz") then
        call writeCoordinatesXYZ(outputFile % funit)

      else if (outputFile % ftype == "pdb") then
        call writeCoordinatesPDB(outputFile % funit,.false.)

      else if (outputFile % ftype == "pdb2") then
        call writeCoordinatesPDB(outputFile % funit,.true.)

      else if (outputFile % ftype == "lmp") then
        call readLammpsForcefieldFile(forcefieldFile)
        call writeLammpsCoordinates(outputFile % funit)

      else if (outputFile % ftype == "gin") then
        call writeGulpCoordinates(outputFile % funit)

      else if (outputFile % ftype == "gin2") then
        call writeGulpCoordinatesFractional(outputFile % funit)

       else if (outputFile % ftype == "dcd") then
        call writeCoordinatesDCD(outputFile % funit, header)

      else if (outputFile % ftype == "psf") then
        call writeCoordinatesPSF(outputFile % funit)

      else if (outputFile % ftype == "arc") then
        call writeCoordinatesARC(outputFile % funit)

#ifdef GPTA_XDR
      else if (outputFile % ftype == "xtc") then
          ! this breaks the read xtc
          ! xtc_out % natoms = frame % natoms
          ! xtc_out % step   = frame % nframe
          ! xtc_out % time   = real(frame % nframe)
          ! xtc_out % box    = real(frame % hmat / 10.d0)
          ! xtc_out % pos    = real(frame % pos / 10.d0)
          ! xtc_out % prec   = 1000.
          call xtc_out % write(frame % natoms, &
                               frame % nframe,   &
                               real(frame % nframe),   &
                               real(frame % hmat / 10.d0),    &
                               real(frame % pos / 10.d0),    &
                               1000.e0)
#endif
      else
        call message(-1,"--o : unknown file format "//trim(outputFile % ftype)//" for file "//trim(outputFile % fname))

      end if

      ! Convert coordinates back from Bohr
      if (outputCoordInBohr) then
        frame % hmat = frame % hmat * rbohr
        frame % pos = frame % pos * rbohr
      end if

    end if

    if (endOfCoordinatesFiles) return

  end subroutine writeCoordinates

end module moduleDumpCoordinates

subroutine writeCoordinatesARC(io)
  use moduleOpenMM, only : qType
  use moduleSystem, only : numberOfAtoms, frame, numberOfCovalentBondsPerAtom, listOfCovalentBondsPerAtom
  implicit none
  integer, intent(in) :: io
  integer :: i, itmp
  write(io,*)numberOfAtoms
  do i=1,numberOfAtoms
    itmp = numberOfCovalentBondsPerAtom(i)
    write(io,'(i5,2x,a4,2x,3(f8.3,2x),6(i7,2x))') i \
      , frame % lab(i) \
      , frame % pos(1:3,i) \
      , qType(i) \
      , listOfCovalentBondsPerAtom(1:itmp,i)
  enddo

end subroutine writeCoordinatesARC

subroutine writeCoordinatesXYZ(io)
  use moduleSystem, only : frame
  use moduleStrings
  implicit none
  integer, intent(in) :: io

  character(len=200) :: latticeString
  character(len=200) :: propertiesString
  integer :: i

  write(latticeString,"(8(f15.5,1x),f12.5)") frame % hmat
  call compact(latticeString)
  latticeString = 'Lattice="'//trim(latticeString)//'"'

  propertiesString = 'Properties=species:S:1:pos:R:3:charges:R:1 pbc="T T T"'

  write(io,'(i0)') frame % natoms
  write(io,'(a,1x,a)') trim(latticeString) , trim(propertiesString)

  do i=1,frame % natoms
    write(io,'(a4,1x,4(f12.5,1x))')frame % lab(i), frame % pos(1:3,i), frame % chg(i)
  end do

end subroutine writeCoordinatesXYZ
