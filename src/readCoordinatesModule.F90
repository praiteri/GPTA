!disclaimer
module moduleRead 
  use moduleCIF
  use moduleGULP
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
        enddo
 
        call assignFlagValue(actionDetails(iact),"+nm ",inputCoordInNM,.false.)
        call assignFlagValue(actionDetails(iact),"+bohr ",inputCoordInBohr,.false.)
        call assignFlagValue(actionDetails(iact),"+frac ",inputCoordInFractional,.false.)
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

    if (numberInputFiles > 0) then
      call message(0,"Opening input coordinates files")
      do i=1,numberInputFiles
        call message(0,"..."//trim(inputFileNames(i) % ftype)//" file",str=inputFileNames(i) % fname)
      enddo
      call message(2)
    end if

    numberOfActions = jcmd

  end subroutine openCoordinatesInputFiles

  subroutine getNumberOfAtoms(f, n, hmat)
    use moduleVariables, only : fileTypeDef, real64
    use moduleMessages 
    use moduleSystem, only : userDefinedCell, userDefinedHMatrix

    implicit none
    type(fileTypeDef), intent(in) :: f
    integer, intent(out) :: n
    real(real64), dimension(3,3),  intent(out) :: hmat

    hmat = reshape([1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 1.0_real64],[3,3])
    select case(f % ftype)
      case default
        call message(-1,"Unknown input file type (getNumberOfAtoms)",str=f % ftype)

      case ("xyz" , "xyz0")  
        call getNumberOfAtomsXYZ(f % funit, n, hmat)

      case ("arc" , "txyz")  
        call getNumberOfAtomsARC(f % funit, n, hmat)

      case ("pdb" , "PDB", "pdb2")
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

      case ("cif") 
        call getNumberOfAtomsCIF(f % funit, n, hmat)

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
    use moduleElements
    use moduleDistances

    implicit none
    logical, intent(out) :: lerr
    integer :: numberOfAtomsLocal
    type(frameTypeDef), intent(inout) :: localFrame
    type(fileTypeDef), optional, intent(inout) :: infile
    
    logical, save :: firstTimeIn = .true.
    character(cp), allocatable, dimension(:), save :: savedLabels
    real(real64), allocatable, dimension(:), save :: savedCharges
    character(2), allocatable, dimension(:), save :: savedElements

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
    integer :: iatm
    integer :: n0
    real(real64) :: dij(3)

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

      numberOfAtomsLocal = inputNumberOfAtoms

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
        
        case ("xyz0")
          call readCoordinatesXYZ_basic(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)

        case ("xyz")
          call readCoordinatesXYZ(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)

        case ("arc" , "txyz")
          call readCoordinatesARC(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)

        case ("pdb", "PDB", "pdb2")
          call readCoordinatesPDB(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, localFrame % element, go)
        
        case ("dcd2")
          block 
            use dcdfort_common
            use dcdfort_reader, only: dcdfile
            type(dcdfile) :: dcd
            integer(kind=int32) :: nframes, istart, nevery, iend, natoms
            real(kind=real32) :: timestep
            real(kind=real32), allocatable :: xyz(:,:)
            real(kind=real64) :: box(6)
            
            if (inputFileNames(currentInputFile) % first_access) then
              inputFileNames(currentInputFile) % first_access = .false.
              close(iounit)
              call dcd % open(trim(inputFileNames(currentInputFile) % fname))
              call dcd % read_header(nframes, istart, nevery, iend, timestep, natoms)

              if (numberOfAtomsLocal /= natoms) call message(-1,"Wrong number of atoms in DCD file")
            end if
            
            allocate(xyz(3,numberOfAtomsLocal))
            call dcd%read_next(xyz, box,go)

            if (go) then

              if ( all(box(4:6) >= -1.0_real64) .and.  all(box(4:6) <= 1.0_real64) ) then
                box(4) = dacos(box(4)) 
                box(5) = dacos(box(5)) 
                box(6) = dacos(box(6)) 
              end if
              call cell2hmat(real(box,8),localFrame % hmat)          

              localFrame % pos = real(xyz,8)
            end if
          end block

        case ("dcd")
          call readCoordinatesDCD(iounit,firstAccess,numberOfAtomsLocal,localFrame % pos,localFrame % hmat,go)
          if (present(infile)) then
            infile % first_access = firstAccess
          else
            inputFileNames(currentInputFile) % first_access = firstAccess
          end if

        case("gin")
          call readCoordinatesGULP(numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)
       
        case("gau")
          call readCoordinatesGaussian(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, go)
   
        case("gro")
          call readCoordinatesGromacs(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)
   
        case("lmp")
          call readCoordinatesLammps(iounit, numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)

        case("lmptrj")
          call readCoordinatesLammpsTrajectory(iounit, numberOfAtomsLocal, &
            localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go)
!            localFrame % pos, localFrame % lab, localFrame % chg, localFrame % hmat, go, localFrame % forces)

        case("cif")
          call readCoordinatesCIF(numberOfAtomsLocal, localFrame % pos, localFrame % lab, localFrame % hmat, go)

#ifdef GPTA_XDR
        case("xtc")
          if (inputFileNames(currentInputFile) % first_access) then
            inputFileNames(currentInputFile) % first_access = .false.
            call xtcf % init(trim(inputFileNames(currentInputFile) % fname))
          end if

          call xtcf % read
          if (xtcf%NATOMS .ne. size(localFrame % pos)/3) &
            call message(-1,"Wrong number of atoms in XTC file",xtcf % NATOMS)

          if (xtcf % STAT ==0) then
            go = .true.
            localFrame % hmat = (xtcf % box) * 10.0_real64
            localFrame % pos(1:3,1:localFrame % natoms) = xtcf % pos(1:3,1:localFrame % natoms) * 10.0_real64
          else
            go = .false.
          end if

        case("trr")
          if (inputFileNames(numberInputFiles) % first_access) then
            inputFileNames(numberInputFiles) % first_access = .false.
            call trrf % init(trim(inputFileNames(numberInputFiles) % fname))
          end if

          call trrf % read
          if (trrf%NATOMS .ne. localFrame % natoms) &
            call message(-1,"Wrong number of atoms in TRR file",trrf % NATOMS)

          if (xtcf % STAT ==0) then
            go = .true.
            localFrame % hmat = (trrf % box)* 10.0_real64
            localFrame % pos(1:3,1:localFrame % natoms) = trrf % pos(1:3,1:localFrame % natoms) * 10.0_real64
          else
            go = .false.
          end if
#endif
      end select

      if (go) then
        if (inputCoordInBohr) then
          localFrame % hmat = localFrame % hmat * rbohr
          localFrame % pos = localFrame % pos * rbohr
        else if (inputCoordInNM) then
          localFrame % hmat = localFrame % hmat * 10.0_real64
          localFrame % pos = localFrame % pos * 10.0_real64
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

      if (ftype /= "pdb") then        
        do iatm=1,frame % natoms
          localFrame % element(iatm) = getElement(localFrame % lab(iatm))
        end do
      end if
      
      allocate(savedLabels(frame % natoms))
      allocate(savedCharges(frame % natoms))
      allocate(savedElements(frame % natoms))
      
      savedLabels = localFrame % lab
      savedCharges = localFrame % chg
      savedElements = localFrame % element

    else

      n0 = size(savedLabels)
      if (resetFrameLabels) localFrame % lab(1:n0) = savedLabels(1:n0)
      if (resetFrameCharges) localFrame % chg(1:n0) = savedCharges(1:n0)
      if (resetFrameElements) localFrame % element(1:n0) = savedElements(1:n0)

    end if
    
    if (inputCoordInFractional) then
      do iatm=1,localFrame % natoms
        dij(1:3) = localFrame % pos(1:3,iatm)
        localFrame % pos(1,iatm) = localFrame % hmat(1,1)*dij(1) + localFrame % hmat(1,2)*dij(2) + localFrame % hmat(1,3)*dij(3)
        localFrame % pos(2,iatm) = localFrame % hmat(2,1)*dij(1) + localFrame % hmat(2,2)*dij(2) + localFrame % hmat(2,3)*dij(3)
        localFrame % pos(3,iatm) = localFrame % hmat(3,1)*dij(1) + localFrame % hmat(3,2)*dij(2) + localFrame % hmat(3,3)*dij(3)
      end do
    end if

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
    real(real64) :: density

    logical, save :: firstTimeIn = .true.

    if (userDefinedCell) then
      frame % hmat = userDefinedHMatrix
    end if

    call getInverseCellMatrix (frame % hmat, frame % hinv, frame % volume)
    
    if (frame % volume < 1.1_real64) then
      pbc_type = "none"
      frame % cell = 0.0_real64
      frame % hmat = identityMatrix
      frame % hinv = identityMatrix

    else

      ! makeUpperTriangularCell(hmat,pos,nn)
      if (abs(frame % hmat(2,1)) + abs(frame % hmat(3,1)) + abs(frame % hmat(3,2)) .gt. 1.0e-6_real64) then
        call makeUpperTriangularCell(frame % hmat, frame % pos, frame % natoms)
        call getInverseCellMatrix (frame % hmat, frame % hinv, frame % volume)
      end if
      call hmat2cell (frame % hmat, frame % cell, "DEG")

      if (abs(frame % hmat(1,2)) + abs(frame % hmat(1,3)) + abs(frame % hmat(2,3)) .lt. 1.0e-6_real64) then
        pbc_type = "ortho"
      else
        pbc_type = "tri"
      end if
    end if
    
    if (firstTimeIn) then
      firstTimeIn = .false.
      
      call initialisePBC(pbc_type)
      if (computeNeighboursList) call setUpNeigboursList() 

      if (frame % natoms > 0) call systemComposition(frame)

      totalMass = 0.0_real64
      block
        real(real64) :: rmass
        do iatm=1,frame % natoms
          rmass = getElementMass(frame % lab(iatm))
          frame % mass(iatm) = rmass
          totalMass = totalMass + rmass
        enddo
      end block

      totalCharge = 0.0_real64
      do iatm=1,frame % natoms
        totalCharge = totalCharge + frame % chg(iatm)
      enddo

      call message(0,"...Total mass (g/mole)",r=totalMass)
      if (pbc_type /= "none") then
        density = frame % natoms / frame % volume
        call message(0,"...Initial density (atoms/A^3)",r=density)
        density = 1.6605388_real64 * totalMass / frame % volume
        call message(0,"...Initial density (g/cm^3)",r=density)
      end if
      call message(1,"...Total charge",r=totalCharge)

    end if

  end subroutine cellProcessing

  subroutine dumpCellInfo(hmat)
    use moduleVariables, only: real64
    use moduleMessages
    implicit none
    real(real64), intent(in) :: hmat(3,3)
    real(real64) :: cell(6), hinv(3,3), volume
    call hmat2cell (hmat, cell, "DEG")
    call getInverseCellMatrix (hmat, hinv, volume)

    if (volume > 1.1_real64) then
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
