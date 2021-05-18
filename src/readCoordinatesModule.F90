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
module moduleRead 

#ifdef GPTA_XDR
  use moduleXDR, only: xtcfile, trrfile
#endif
  contains

  subroutine openCoordinatesInputFiles()
    use moduleVariables
    use moduleStrings
    use moduleFiles
    use moduleSystem 
    use moduleActions
    use moduleMessages 
    use moduleStrings
    implicit none

    integer :: iact, jcmd
    integer :: i, iw, nw
    character(len=STRLEN), dimension(MAXWORDS) :: words  
    character(len=STRLEN) :: stringCell

    integer :: idx, ilen
    character(len=STRLEN), dimension(MAXWORDS) :: fname
    character(len=fp), dimension(MAXWORDS) :: ftype

    jcmd = 0
    if (me /= readingCPU) then
      do iact=1,numberOfActions
        if (actionType(iact)(1:3) == "--i") then
          call parse(actionDetails(iact)," ",words,nw)
          do iw=1,nw
            call compact(words(iw))
            if (words(iw)(1:1)=="+") exit
            numberInputFiles = numberInputFiles + 1
            fname(numberInputFiles) = trim(words(iw))
            if (len_trim(actionType(iact)) > 3) then
              ftype(numberInputFiles) = actionType(iact)(4:3+fp)
            else 
              ilen = len_trim(fname(numberInputFiles))
              idx = index(fname(numberInputFiles),".",back=.true.)
              if (ilen-idx>0) then
                ftype(numberInputFiles) = fname(numberInputFiles)(idx+1:ilen)
              else
                ftype(numberInputFiles) = 'NULL'
              end if
            end if
          end do
          cycle
        end if
        jcmd = jcmd + 1
        actionType(jcmd) = actionType(iact)
        actionDetails(jcmd) = actionDetails(iact)
      enddo
      numberOfActions = jcmd
      
      if (numberInputFiles == 0) return
      ! call message(0,"Opening input coordinates files")
      ! do i=1,numberInputFiles
      !   call message(0,"..."//trim(ftype(i))//" file",str=fname(i))
      ! enddo
      ! call message(2)

      return
      
    end if

    currentInputFile = 1
    numberInputFiles = 0
    numberOfFramesRead = 0
    numberOfFramesProcessed = 0

    do iact=1,numberOfActions
      if (actionType(iact)(1:3) == "--i") then

        call parse(actionDetails(iact)," ",words,nw)
        do iw=1,nw
          call compact(words(iw))
          if (words(iw)(1:1)=="+") exit
          numberInputFiles = numberInputFiles + 1
          call initialiseFile(inputFileNames(numberInputFiles), words(iw), fstatus='old')
          ! file type from extension
          if (len_trim(actionType(iact)) > 3) inputFileNames(numberInputFiles) % ftype = actionType(iact)(4:3+fp)

! #ifdef GPTA_XDR
!           if ( inputFileNames(numberInputFiles) % ftype == "xtc") then
!             close(inputFileNames(numberInputFiles) % funit)
!           else if ( inputFileNames(numberInputFiles) % ftype == "trr") then
!             close(inputFileNames(numberInputFiles) % funit)
!           end if
! #endif
        enddo
 
        call assignFlagValue(actionDetails(iact),"+nm ",inputCoordInNM,.false.)
        call assignFlagValue(actionDetails(iact),"+bohr ",inputCoordInBohr,.false.)
        call assignFlagValue(actionDetails(iact),"+cell ", stringCell, "NONE")
        if (stringCell == "NONE") then
          userDefinedCell = .false.
        else
          userDefinedCell = .true.
          call readCellFreeFormat(stringCell, userDefinedHMatrix)
        end if

      ! Remove --i commands from list
      else
        jcmd = jcmd + 1
        actionType(jcmd) = actionType(iact)
        actionDetails(jcmd) = actionDetails(iact)

      end if

    enddo
    call message(0,"Opening input coordinates files")
    do i=1,numberInputFiles
      call message(0,"..."//trim(inputFileNames(i) % ftype)//" file",str=inputFileNames(i) % fname)
    enddo
    call message(2)

    numberOfActions = jcmd

  end subroutine openCoordinatesInputFiles

  subroutine getNumberOfAtoms(f, n, hmat)
    use moduleVariables, only : fileTypeDef
    use moduleMessages 
    use moduleSystem, only : userDefinedCell, userDefinedHMatrix
    implicit none
    type(fileTypeDef), intent(in) :: f
    integer, intent(out) :: n
    real(8), dimension(3,3),  intent(out) :: hmat

    hmat = reshape([1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0],[3,3])
    select case(f % ftype)
      case default
        call message(-1,"Unknown input file type (getNumberOfAtoms)",str=f % ftype)

      case ("xyz")  
        call getNumberOfAtomsXYZ(f % funit, n, hmat)

      case ("pdb" , "PDB")
        call getNumberOfAtomsPDB(f % funit, n, hmat)

      case ("dcd", "xtc", "trr")
        call message(-1,"First file type should be ASCII")

      case ("gin")  
        call getNumberOfAtomsGULP(f % funit, n, hmat)

      case ("gau")  
        call getNumberOfAtomsGaussian(f % funit, n)

      case ("gro")  
        call getNumberOfAtomsGromacs(f % funit, n, hmat)

      case ("lmp")  
        call getNumberOfAtomsLammps(f % funit, n, hmat)

      case ("lmptrj")  
        call getNumberOfAtomsLammpsTrajectory(f % funit, n, hmat)

      end select
      rewind(f % funit)
      
      if (userDefinedCell) then
        hmat = userDefinedHMatrix
      end if

  end subroutine getNumberOfAtoms

  subroutine readCoordinates(lerr,localFrame,infile)
    use moduleVariables
    use moduleSystem 
    use moduleMessages 
    implicit none
    logical, intent(out) :: lerr
    integer :: numberOfAtomsLocal
    type(frameTypeDef), intent(inout) :: localFrame
    type(fileTypeDef), optional, intent(inout) :: infile
    
    logical, save :: firstTimeIn = .true.
    character(cp), allocatable, dimension(:), save :: savedLabels
    real(8), allocatable, dimension(:), save :: savedCharges

#ifdef GPTA_XDR
    type(xtcfile), save :: xtcf
    type(trrfile), save :: trrf
#endif

    logical :: go
    integer :: iounit
    character(fp) :: ftype
    character(len=STRLEN) :: fname
    logical :: firstAccess
    integer :: now, next_frame

    ! Default to no error
    lerr = .true.

    ! Read coordinates for an action
    if (present(infile)) then
      iounit = infile % funit
      ftype = infile % ftype
      fname = infile % fname
      firstAccess = infile % first_access
      now = 0
      next_frame = 1

      numberOfAtomsLocal = size(localFrame % pos)/3

    ! Read frames for general processing
    else
      iounit = inputFileNames(currentInputFile) % funit
      ftype = inputFileNames(currentInputFile) % ftype
      fname = inputFileNames(currentInputFile) % fname
      firstAccess = inputFileNames(currentInputFile) % first_access

      now = numberOfFramesRead
      next_frame = getNextFrameNumber(numberOfFramesRead)

      numberOfAtomsLocal = numberOfAtoms

      if (next_frame < 0) then
        lerr = .false.
        call message(-2)
        return
      end if

    end if

    do while (now < next_frame)
      select case (ftype)
        case default 
          call message(-1,"Unknown input FILE",str=fname)
        
        case ("xyz")
          call readCoordinatesXYZ(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)

        case ("pdb")
          call readCoordinatesPDB(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)
        
        case ("dcd")
          call readCoordinatesDCD(iounit,firstAccess,numberOfAtomsLocal,localFrame % pos,localFrame % hmat,go)
          if (present(infile)) then
            infile % first_access = firstAccess
          else
            inputFileNames(currentInputFile) % first_access = firstAccess
          end if

        case("gin")
          call readCoordinatesGULP(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)
        
        case("gau")
          call readCoordinatesGaussian(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, go)
   
        case("gro")
          call readCoordinatesGromacs(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)
   
        case("lmp")
          call readCoordinatesLammps(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)

        case("lmptrj")
          call readCoordinatesLammpsTrajectory(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)

#ifdef GPTA_XDR
        case("xtc")
          if (inputFileNames(currentInputFile) % first_access) then
            inputFileNames(currentInputFile) % first_access = .false.
            call xtcf % init(trim(inputFileNames(currentInputFile) % fname))
          endif

          call xtcf % read
          if (xtcf%NATOMS .ne. size(localFrame % pos)/3) &
            call message(-1,"Wrong number of atoms in XTC file",xtcf % NATOMS)

          if (xtcf % STAT ==0) then
            go = .true.
            localFrame % hmat = (xtcf % box) * 10.d0
            localFrame % pos(1:3,1:localFrame % natoms) = xtcf % pos(1:3,1:localFrame % natoms) * 10.d0
          else
            go = .false.
          endif

        case("trr")
          if (inputFileNames(numberInputFiles) % first_access) then
            inputFileNames(numberInputFiles) % first_access = .false.
            call trrf % init(trim(inputFileNames(numberInputFiles) % fname))
          endif

          call trrf % read
          if (trrf%NATOMS .ne. localFrame % natoms) &
            call message(-1,"Wrong number of atoms in TRR file",trrf % NATOMS)

          if (xtcf % STAT ==0) then
            go = .true.
            localFrame % hmat = (trrf % box)* 10.d0
            localFrame % pos(1:3,1:localFrame % natoms) = trrf % pos(1:3,1:localFrame % natoms) * 10.d0
          else
            go = .false.
          endif
#endif
      end select

      if (go) then
        if (inputCoordInBohr) then
          localFrame % hmat = localFrame % hmat * rbohr
          localFrame % pos = localFrame % pos * rbohr
        else if (inputCoordInNM) then
          localFrame % hmat = localFrame % hmat * 10.d0
          localFrame % pos = localFrame % pos * 10.d0
        end if
      end if 

      if (present(infile)) then
        if (go) then
          return
        else
          call message(-1,"Error while reading auxiliary coordinated file")
        end if
      end if

      ! this part is for trajectory reading only
      if (go) then
        now = now + 1
      else
        if (currentInputFile < numberInputFiles .and. .not. present(infile)) then
          currentInputFile = currentInputFile + 1
          iounit = inputFileNames(currentInputFile) % funit
          ftype = inputFileNames(currentInputFile) % ftype
          fname = inputFileNames(currentInputFile) % fname
          firstAccess = inputFileNames(currentInputFile) % first_access
          cycle
        end if
        lerr = .false.
        call message(-2)
        return
      end if

    enddo
    
    numberOfFramesRead = now
    localFrame % nframe = numberOfFramesRead
    localFrame % natoms = numberOfAtomsLocal
    
    if (firstTimeIn) then
      firstTimeIn = .false.
      
      allocate(savedLabels(numberOfAtoms))
      allocate(savedCharges(numberOfAtoms))
      
      savedLabels = localFrame % lab
      savedCharges = localFrame % chg
    else
      if (keepFrameLabels) localFrame% lab = savedLabels
      if (keepFrameCharges) localFrame% chg = savedCharges
    end if
    
    ! if (mod(numberOfFramesRead,nProgress)==0) then
    !   if (numberOfFramesRead > nProgress) call message(-2)
    !   call message(0,"---Frames read",numberOfFramesRead)
    ! end if
    
  end subroutine readCoordinates
  
  function getNextFrameNumber(current) result(next)
    use moduleSystem , only : first_frame, last_frame, stride_frame, nlist_frames, list_frames
    implicit none
    integer, intent(in) :: current
    integer :: next
    integer, save :: iframe = 0

    if (nlist_frames == 0) then
      if (current < first_frame) then
        next = first_frame
      else
        next = current + stride_frame
      end if  
      if (next > last_frame) next = -1

    else
      iframe = iframe + 1
      if (iframe > nlist_frames) then
        next = -1
      else
        next = list_frames(iframe)
      end if

    end if

  end function getNextFrameNumber

  subroutine cellProcessing()
    use moduleVariables
    use moduleMessages 
    use moduleSystem
    use moduleNeighbours
    use moduleDistances, only : initialisePBC
    use moduleElements, only : getElementMass


    implicit none
    integer :: iatm
    real(8) :: density

    logical, save :: firstTimeIn = .true.

    if (userDefinedCell) then
      frame % hmat = userDefinedHMatrix
    end if

    call getInverseCellMatrix (frame % hmat, frame % hinv, frame % volume)
    
    if (frame % volume < 1.1d0) then
      pbc_type = "none"
      frame % cell = 0.d0
      frame % hmat = identityMatrix
      frame % hinv = identityMatrix

    else

      ! makeUpperTriangularCell(hmat,pos,nn)
      if (abs(frame % hmat(2,1)) + abs(frame % hmat(3,1)) + abs(frame % hmat(3,2)) .gt. 1.0d-6) then
        call makeUpperTriangularCell(frame % hmat, frame % pos, frame % natoms)
        call getInverseCellMatrix (frame % hmat, frame % hinv, frame % volume)
      endif
      call hmat2cell (frame % hmat, frame % cell, "DEG")

      if (abs(frame % hmat(1,2)) + abs(frame % hmat(1,3)) + abs(frame % hmat(2,3)) .lt. 1.0d-6) then
        pbc_type = "ortho"
      else
        pbc_type = "tri"
      end if
    end if
    
    if (firstTimeIn) then
      firstTimeIn = .false.
      
          call initialisePBC(pbc_type)
          if (computeNeighboursList) call setUpNeigboursList() 

      call systemComposition(frame)

      totalMass = 0.d0
      do iatm=1,frame % natoms
        totalMass = totalMass + getElementMass(frame % lab(iatm))
      enddo

      totalCharge = 0.d0
      do iatm=1,frame % natoms
        totalCharge = totalCharge + frame % chg(iatm)
      enddo

      call message(0,"...Total mass (g/mole)",r=totalMass)
      if (pbc_type /= "none") then
        density = 1.6605388d0 * totalMass / frame % volume
        call message(0,"...Density (g/cm^3)",r=density)
      end if
      call message(1,"...Total charge",r=totalCharge)

    end if

  end subroutine cellProcessing

  subroutine dumpCellInfo(hmat)
    use moduleMessages
    implicit none
    real(8), intent(in) :: hmat(3,3)
    real(8) :: cell(6), hinv(3,3), volume
    call hmat2cell (hmat, cell, "DEG")
    call getInverseCellMatrix (hmat, hinv, volume)

    if (volume > 1.01d0) then
      call message(0,"Initial cell")
      call message(0,"...Cell vector A",rv=hmat(1:3,1))
      call message(0,"...Cell vector B",rv=hmat(1:3,2))
      call message(0,"...Cell vector C",rv=hmat(1:3,3))
      call message(0,"...Cell lengths",rv=cell(1:3))
      call message(0,"...Cell angles",rv=cell(4:6))
      call message(1,"...Volume",r=volume)
    else
      call message(1,"Non periodic system")
    end if

  end subroutine dumpCellInfo

end module moduleRead 
