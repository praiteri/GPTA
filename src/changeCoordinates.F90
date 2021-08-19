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
module moduleModifyCoordinates

contains

  subroutine rescaleCellHelp()
    use moduleMessages
    implicit none
    call message(0,"This action changes the system cell and rescales the atomic coordinates")
    call message(0,"The cell can be expressed as:")
    call message(0,"  1 number -> cubic")
    call message(0,"  3 number -> orthorhombic")
    call message(0,"  6 number -> a, b, c, alpha, beta, gamma")
    call message(0,"  9 number -> metric matrix")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --rescale +cell 10")
    call message(0,"  gpta.x --i coord.pdb --rescale +cell 10,20,30")
    call message(0,"  gpta.x --i coord.pdb --rescale +cell 10,10,10,90,90,120")
    call message(0,"  gpta.x --i coord.pdb --rescale +cell 10,0,0,0,10,0,0,0,10")
  end subroutine rescaleCellHelp

  subroutine shiftCoordinatesHelp()
    use moduleMessages
    implicit none
    call message(0,"This action shifts the atomic coordinates by a 3D vector")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --shift 1,2,3")
  end subroutine shiftCoordinatesHelp

  subroutine shiftCOMHelp()
    use moduleMessages
    implicit none
    call message(0,"This actions shifts the centre of mass of the selected atoms to:")
    call message(0,"  * the centre of the cell (+centre)")
    call message(0,"  * the same position as the first frame (+initial)")
    call message(0,"  * a user defined point (+pos)")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --fixcom +c Ca +centre")
    call message(0,"  gpta.x --i coord.pdb --fixcom +c Ca +initial")
    call message(0,"  gpta.x --i coord.pdb --fixcom +c Ca +pos 10,20,30")
  end subroutine shiftCOMHelp

  subroutine applyPeriodicboundaryConditionsHelp()
    use moduleMessages
    implicit none
    call message(0,"This actions wraps the atoms/molecules inside the cell.")
    call message(0,"If the topology is known the molecues are not broken.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --pbc")
    call message(0,"  gpta.x --i coord.pdb --top --pbc")
  end subroutine applyPeriodicboundaryConditionsHelp

  subroutine unwrapCoordinatesHelp()
    use moduleMessages
    implicit none
    call message(0,"This actions unwraps the atoms/molecules across the PBC to generate a continuous trajectory.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb trajectory.dcd --unwrap")
    call message(0,"  gpta.x --i coord.pdb trajectory.dcd --top --unwrap")
  end subroutine unwrapCoordinatesHelp

  subroutine replicateCellHelp()
    use moduleMessages
    implicit none
    call message(0,"This action creates a supercell of the initial system.")
    call message(0,"If the topology is known the replication is done by molecule; the molecues must not be broken.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --repl 2,2,2")
    call message(0,"  gpta.x --i coord.pdb --top --repl 2,2,2")
  end subroutine replicateCellHelp

  subroutine mirrorCellHelp()
    use moduleMessages
    implicit none
    call message(0,"This action replicates the system by making its mirror image.")
    call message(0,"The mirror plane is one of the faces of the cell normal to x, y or z")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --mirror +x")
    call message(0,"  gpta.x --i coord.pdb --mirror +y")
    call message(0,"  gpta.x --i coord.pdb --mirror +z")
  end subroutine mirrorCellHelp

  subroutine removeOverlappingMoleculesHelp()
    use moduleMessages
    implicit none
    call message(0,"This action remove molecules that overlap.")
    call message(0,"The minimum distance between atoms can be specified by the +rmin flag.")
    call message(0,"If the topology is known the the entire molecule is removed.")
    call message(0,"It is also possoble to do a dry run to see how many atoms/molecules will be removed with certain distances.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --noclash +rmin 1.1")
    call message(0,"  gpta.x --i coord.pdb --top --noclash +rmin 1.1")
    call message(0,"  gpta.x --i coord.pdb --top --noclash +rmin 1.1 +dry")
  end subroutine removeOverlappingMoleculesHelp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rescaleCell(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    integer :: numberOfWords
    character(len=STRLEN), dimension(100) :: listOfWords

    logical, pointer :: actionInitialisation
    character(len=100) :: stringCell
    real(8), dimension(3,3), save :: hmat, hinv
    real(8), save :: volume, cell(6)

    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      call parse(actionCommand," ",listOfWords,numberOfWords)

      if (numberOfWords==0) call message(-1,"--rescale : nothing to do")

      call assignFlagValue(actionCommand,"+cell ", stringCell, "NONE")
      call readCellFreeFormat(stringCell, hmat)
      call getInverseCellMatrix(hmat,hinv,volume)
      call hmat2cell (hmat, cell, "DEG")

      call message(0,"Rescaled cell")
      call message(0,"...Cell vector A",rv=hmat(1:3,1))
      call message(0,"...Cell vector B",rv=hmat(1:3,2))
      call message(0,"...Cell vector C",rv=hmat(1:3,3))
      call message(0,"...Cell lengths",rv=cell(1:3))
      call message(0,"...Cell angles",rv=cell(4:6))
      call message(1,"...Volume",r=volume)

      call checkUsedFlags(actionCommand)
      return
    end if
    
    if (frameReadSuccessfully) then 
      frame % hmat = hmat
      frame % hinv = hinv
      frame % volume = volume
      frame % cell = cell
      call fractionalToCartesian(frame % natoms, frame % frac, frame % pos)
    end if

    if (endOfCoordinatesFiles) return

  end subroutine rescaleCell

  subroutine shiftCoordinates(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    integer :: numberOfWords
    character(len=STRLEN), dimension(100) :: listOfWords

    logical, pointer :: actionInitialisation

    real(8), pointer, dimension(:) :: dpos
    integer :: i

    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    dpos(1:3)            => a % doubleVariables(1:3)
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      dpos = 0.d0
      call parse(actionCommand," ",listOfWords,numberOfWords)

      if (numberOfWords==0) call message(-1,"--shift : nothing to do")

      if (numberOfWords==3) then
        read(actionCommand,*)dpos(1:3)
      else
        call parse(actionCommand,",",listOfWords,numberOfWords)
        if (numberOfWords==3) read(actionCommand,*)dpos(1:3)
      end if

      if (sum(abs(dpos(1:3)))<1d-6) call message(-1,"--shift : zero displacement")
      
      call message(1,"Translating atoms' coordinates by",rv=dpos)

      call checkUsedFlags(actionCommand)
      return
    end if
    
    if (frameReadSuccessfully) then 
      do i=1,frame % natoms
        frame % pos(1:3,i) = frame % pos(1:3,i) + dpos
      enddo
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
    end if

!    if (numberOfMolecules>0) call computeMoleculesCOM(0)
    if (numberOfMolecules>0) then
      block
        integer :: imol
        do imol=1,numberOfMolecules
          listOfMolecules(imol) % centreOfMass = listOfMolecules(imol) % centreOfMass + dpos
        end do
      end block
    end if
    if (endOfCoordinatesFiles) return

  end subroutine shiftCoordinates

  subroutine shiftCOM(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand
    logical, pointer :: firstAction
    logical, pointer :: actionInitialisation
    logical, pointer :: lcentre
    logical, pointer :: linitial
    real(8), pointer, dimension(:) :: dpos
    
    real(8), dimension(3) :: xcom, xcentre, delta

    integer :: i, idx
    integer, save :: nsel

    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation
    dpos(1:3)            => a % doubleVariables(1:3)
    firstAction          => a % firstAction
    lcentre              => a % logicalVariables(1)
    linitial             => a % logicalVariables(2)

    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      ! fix the centre of mass to the cell centre      
      call assignFlagValue(actionCommand,"+centre",lcentre,.false.)
      
      ! fix the centre of mass intial position
      call assignFlagValue(actionCommand,"+initial",linitial,.false.)

      if (.not. lcentre) call assignFlagValue(actionCommand,"+pos",dpos,[0.d0,0.d0,0.d0])

      ! call checkUsedFlags(actionCommand)

      return
    end if
    
    if (frameReadSuccessfully) then 

      if (firstAction) then
        if (lcentre) then
          call message(0,"Translating centre of mass to the centre of the cell")

        else if (linitial) then
          call message(0,"Translating centre of mass to its initial position")
        
        else
          call message(0,"Translating centre of mass to",rv=dpos)
        end if

        call selectAtoms(1,actionCommand,a)
        call createSelectionList(a,1)
        if (linitial) then
          block
            nsel = count(a % isSelected(:,1))

            dpos = 0.d0
            do i=1,nsel
              idx = a % idxSelection(i,1)
              dpos(1:3) = dpos(1:3) + frame % pos(1:3,idx)
            end do
            dpos(1:3) = dpos(1:3) / dble(nsel)
          end block
        end if

        if (resetFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if

        call checkUsedFlags(actionCommand)

        firstAction = .false.

        if (linitial) return

      else

        if (a % updateAtomsSelection) then
          call selectAtoms(1,actionCommand,a)
          call createSelectionList(a,1)
        end if

      end if
      nsel = count(a % isSelected(:,1))

      xcom = 0.d0
      do i=1,nsel
        idx = a % idxSelection(i,1)
        xcom(1:3) = xcom(1:3) + frame % pos(1:3,idx)
      end do
      xcom(1:3) = xcom(1:3) / dble(nsel)
      
      if (lcentre) then
        xcentre(1:3) = frame % hmat(1:3,1) + frame % hmat(1:3,2) + frame % hmat(1:3,3)
        xcentre = xcentre / 2.d0
        delta(1:3) = xcentre(1:3) - xcom(1:3)
      else
        delta(1:3) = dpos(1:3) - xcom(1:3)
      end if

      do i=1,frame % natoms
        frame % pos(1:3,i) = frame % pos(1:3,i) + delta
      enddo
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)
    end if

!    if (numberOfMolecules>0) call computeMoleculesCOM(0)
    if (numberOfMolecules>0) then
      block
        integer :: imol
        do imol=1,numberOfMolecules
          listOfMolecules(imol) % centreOfMass = listOfMolecules(imol) % centreOfMass + delta
        end do
      end block
    end if

    if (endOfCoordinatesFiles) return

  end subroutine shiftCOM

  subroutine applyPeriodicboundaryConditions(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a

    if (a % actionInitialisation) then
      a % actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      call message(1,"Applying periodic boundary conditions",str='[0:1]')
      return
    end if
    
    if (frameReadSuccessfully) then 
      call fractionalToCartesian(frame % natoms,frame % frac,frame % pos)
      if (numberOfMolecules > 0) call reassembleAllMolecules()
    end if

    if (endOfCoordinatesFiles) return

  end subroutine applyPeriodicboundaryConditions

  subroutine unwrapCoordinates(a)
    use moduleVariables
    use moduleSystem 
    use moduleDistances
    use moduleMessages
    implicit none
    type(actionTypeDef), target :: a

    integer :: iatm
    real(8) :: dij(3), dist

    ! Associate variables
    
    if (a % actionInitialisation) then
      a % actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      return
    end if
    
    if (frameReadSuccessfully) then 

      if (a % firstAction) then
        call message(1,"Unwrapping trajectory")
        allocate(a % localPositions(3,frame % natoms), source=frame % pos)
        call checkUsedFlags(a % actionDetails)
        a % firstAction = .false.
      endif

      do iatm=1,frame % natoms
        if (any(frame % frac(:,iatm) < 0.d0) .or. any(frame % frac(:,iatm) > 1.d0)) then
          dij = frame % pos(:,iatm) - a % localPositions(:,iatm)
          dist = computeDistanceSquaredPBC(dij)
          frame % pos(:,iatm) = a % localPositions(:,iatm) + dij
        end if
       ! dij = frame % pos(:,iatm) - a % localPositions(:,iatm)
       ! dist = computeDistanceSquaredPBC(dij)
       ! frame % pos(:,iatm) = a % localPositions(:,iatm) + dij
      end do

      a % localPositions = frame % pos

      if (numberOfMolecules > 0) call reassembleAllMolecules()
    end if

    if (endOfCoordinatesFiles) return

  end subroutine unwrapCoordinates

  subroutine replicateCell(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    use moduleNeighbours
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand

    logical, pointer :: actionInitialisation
    real(8), allocatable, dimension(:,:) :: cartesianCoord
    real(8), allocatable, dimension(:) :: saveCharges
    character(cp), allocatable, dimension(:) :: saveLabels
    real(8), allocatable, dimension(:) :: saveMasses
    
    integer, pointer, dimension(:) :: idx, jdx, kdx
    integer :: nargs
    character(len=10) :: args(10)
    integer :: iatm, i, j, k, nrepl, nn
    
    ! Associate variables
    actionCommand        => a % actionDetails
    actionInitialisation => a % actionInitialisation

    idx(1:2)             => a % integerVariables(1:2)
    jdx(1:2)             => a % integerVariables(3:4)
    kdx(1:2)             => a % integerVariables(5:6)
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      call message(1,"Replicating cell")
      
      call parse(actionCommand," ",args,nargs)
      if (nargs == 1) then
        call parse(args(1),",",args,nargs)
      endif
      if (nargs/=3) call message(-1,"--repl - wrong number of arguments")
      idx = 0
      jdx = 0
      kdx = 0
      read(args(1),*)idx(2)
      read(args(2),*)jdx(2)
      read(args(3),*)kdx(2)
      idx(2) = idx(2) - 1
      jdx(2) = jdx(2) - 1
      kdx(2) = kdx(2) - 1

      call checkUsedFlags(actionCommand)

      return
    end if
    
    if (frameReadSuccessfully) then 
      if (numberOfMolecules > 0) call runInternalAction("pbc","NULL")

      nrepl = (idx(2)-idx(1)+1) * (jdx(2)-jdx(1)+1) * (kdx(2)-kdx(1)+1)

      allocate(cartesianCoord(3, frame % natoms))
      allocate(    saveLabels(   frame % natoms))
      allocate(   saveCharges(   frame % natoms))
      allocate(    saveMasses(   frame % natoms))
      cartesianCoord = frame % pos
          saveLabels = frame % lab
         saveCharges = frame % chg
          saveMasses = frame % mass
  
  
      ! call move_alloc(frame % pos, cartesianCoord)
      ! call move_alloc(frame % lab, saveLabels)
      ! call move_alloc(frame % chg, saveCharges)
      ! call move_alloc(frame % mass, saveMasses)
      deallocate(frame % pos)
      deallocate(frame % lab)
      deallocate(frame % chg)
      deallocate(frame % mass)
      deallocate(frame % frac)
      
      allocate(frame % pos(3,frame % natoms * nrepl))
      allocate(frame % lab(frame % natoms * nrepl))
      allocate(frame % chg(frame % natoms * nrepl))
      allocate(frame % mass(frame % natoms * nrepl))
      allocate(frame % frac(3,frame % natoms * nrepl))
      
      nn = 0 
      do iatm=1,frame % natoms
        do i=idx(1),idx(2)
          do j=jdx(1),jdx(2)
            do k=kdx(1),kdx(2)
              nn = nn + 1
              frame % pos(1:3,nn) = cartesianCoord(1:3,iatm) & 
                                  + i*frame % hmat(1:3,1) &
                                  + j*frame % hmat(1:3,2) &
                                  + k*frame % hmat(1:3,3)
              frame % lab(nn) = saveLabels(iatm)
              frame % chg(nn) = saveCharges(iatm)
              frame % mass(nn)= saveMasses(iatm)
            end do
          end do
        end do
      end do
      
      frame % natoms = nn
      frame % hmat(1:3,1) = (idx(2)-idx(1)+1)*frame % hmat(1:3,1)
      frame % hmat(1:3,2) = (jdx(2)-jdx(1)+1)*frame % hmat(1:3,2)
      frame % hmat(1:3,3) = (kdx(2)-kdx(1)+1)*frame % hmat(1:3,3)

      call getInverseCellMatrix(frame % hmat, frame % hinv, frame % volume)
      call hmat2cell (frame % hmat, frame % cell, "DEG")

      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

      call setUpNeigboursList()
      call updateNeighboursList(.true.)
      
      call message(0,"New cell")
      call message(0,"...Cell vector A",rv=frame % hmat(1:3,1))
      call message(0,"...Cell vector B",rv=frame % hmat(1:3,2))
      call message(0,"...Cell vector C",rv=frame % hmat(1:3,3))
      call message(0,"...Cell lengths",rv=frame % cell(1:3))
      call message(0,"...Cell angles",rv=frame % cell(4:6))
      call message(1,"...Volume",r=frame % volume)

      if (numberOfMolecules > 0) call runInternalAction("topology","+update +reorder +rebuild")

    end if

    if (endOfCoordinatesFiles) return

  end subroutine replicateCell

  subroutine mirrorCell(a)
    use moduleVariables
    use moduleSystem 
    use moduleStrings
    use moduleMessages 
    use moduleDistances
    use moduleNeighbours
    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand  
    logical, pointer :: firstAction

    logical, pointer :: actionInitialisation
    
    integer, pointer :: idx
    integer :: nargs
    character(len=10) :: args(10)
    real(8), dimension(3) :: rtmp
    integer :: iatm, nrepl, nn, i
    logical :: lflag(3)
    
    ! Associate variables
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    actionInitialisation => a % actionInitialisation

    idx                  => a % integerVariables(1)
    
    if (actionInitialisation) then
      actionInitialisation = .false.
      a % requiresNeighboursList = .false.
      
      call message(1,"Mirroring the cell")
      
      call parse(actionCommand," ",args,nargs)
      call assignFlagValue(actionCommand,"+x",lflag(1),.false.)
      call assignFlagValue(actionCommand,"+y",lflag(2),.false.)
      call assignFlagValue(actionCommand,"+z",lflag(3),.false.)
      if (count(lflag) == 0) then
        call message(-1,"--mirror | no direction specified")
      else if (count(lflag) > 1) then
        call message(-1,"--mirror | only one direction can be specified at one time")
      else
        do i=1,3
          if (lflag(i)) idx = i
        end do
      end if

      call checkUsedFlags(actionCommand)

      return
    end if
    
    if (frameReadSuccessfully) then 

      if (numberOfMolecules > 0) call runInternalAction("pbc","NULL")

      nrepl = 2
     
      if (a % firstAction) then
        call extendFrame(frame % natoms * nrepl)
        call checkUsedFlags(actionCommand)
        a % firstAction = .false.
      end if

      rtmp = 2.d0 * frame % hmat(:,idx)
      nn = frame % natoms 
      do iatm=1,frame % natoms
        frame % pos(1:3,nn+iatm) = frame % pos(1:3,iatm) 
        frame % lab(    nn+iatm) = frame % lab(    iatm) 
        frame % chg(    nn+iatm) = frame % chg(    iatm) 
      end do

      do iatm=1,frame % natoms
        frame % pos(idx,nn+iatm) = rtmp(idx) - frame % pos(idx,iatm) 
      end do
      
      frame % natoms = frame % natoms * nrepl
      frame % hmat(1:3,idx) = nrepl * frame % hmat(1:3,idx)

      call getInverseCellMatrix(frame % hmat, frame % hinv, frame % volume)
      call hmat2cell (frame % hmat, frame % cell, "DEG")

      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

      call setUpNeigboursList()
      call updateNeighboursList(.true.)
      
      if (numberOfMolecules > 0) call runInternalAction("topology","+update +reorder +rebuild")

    end if

    if (endOfCoordinatesFiles) return

  end subroutine mirrorCell

  subroutine removeOverlappingMolecules(a)
    use moduleDistances
    use moduleSystem
    use moduleVariables
    use moduleMessages 
    use moduleNeighbours
    use moduleStrings

    implicit none
    type(actionTypeDef), target :: a
    character(:), pointer :: actionCommand  
    logical, pointer :: firstAction
    logical, pointer :: dryRun

    real(8), pointer :: overlapDistance
    integer :: nmol0, natm0
    integer :: nmol1, natm1
    
    ! Associate variables
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction

    overlapDistance      => a % doubleVariables(1)
    dryRun               => a % logicalVariables(1)

    if (a % actionInitialisation) then
      a % actionInitialisation = .false.
      a % requiresNeighboursList = .true.
      a % requiresNeighboursListUpdates = .true.
      a % requiresNeighboursListDouble = .false.
      
      call assignFlagValue(actionCommand,"+rmin",overlapDistance,0.5d0)
      a % cutoffNeighboursList = overlapDistance

      call assignFlagValue(actionCommand,"+dry",dryRun,.false.)

      call checkUsedFlags(actionCommand)

      return
    end if

    if (frameReadSuccessfully) then 

      if (firstAction) then
        call message(0,"Removing overlap")
        call message(0,"...Minimum overlap distance",r=overlapDistance)
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if
      
      block
        integer :: iatm, ineigh
        integer :: ndist0 = 0
        integer :: ndist1 = 0
        integer :: ndist2 = 0
        integer :: ndist3 = 0
        real(8) :: rdist1 = 0.1d0
        real(8) :: rdist2 = 0.5d0
        real(8) :: rdist3 = 1.0d0
        character(len=5) :: str
        do iatm=1,frame % natoms
          do ineigh=1,nneigh(iatm)

            if (rneigh(ineigh,iatm) < rdist1) then
              ndist1 = ndist1 + 1
            else if (rneigh(ineigh,iatm) < rdist2) then
              ndist2 = ndist2 + 1
            else if (rneigh(ineigh,iatm) < rdist3) then
              ndist3 = ndist3 + 1
            end if
          
            if (rneigh(ineigh,iatm) < overlapDistance) ndist0 = ndist0 + 1
            
          end do
        end do

        write(str,'(f5.3)')overlapDistance
        call message(0,"...Number of contacts below "//trim(str)//" A",i=ndist0)
        write(str,'(f3.1)')rdist1
        call message(0,"...Number of contacts below "//trim(str)//" A",i=ndist1)
        write(str,'(f3.1)')rdist2
        call message(0,"...Number of contacts below "//trim(str)//" A",i=ndist2)
        write(str,'(f3.1)')rdist3
        call message(1,"...Number of contacts below "//trim(str)//" A",i=ndist3)
      end block
        
      if (dryRun) return

      natm0 = frame % natoms
      nmol0 = numberOfMolecules

      if (numberOfMolecules == 0) then
        call deleteOverlapAllAtoms(overlapDistance,.true.)
      else
        call deleteOverlapAllMolecules(overlapDistance,.true.)
      end if
      call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

      ! update neighbours' list
      call updateNeighboursList(.true.)
      
      natm1 = frame % natoms
      call message(0,"...Number of atoms removed",natm0-natm1)
      if (numberOfMolecules > 0) then
        call runInternalAction("topology","+update")
        nmol1 = numberOfMolecules
        call message(0,"...Number of molecules removed",nmol0-nmol1)
      end if
      call message(2)

    end if

    if (endOfCoordinatesFiles) return
    
  end subroutine removeOverlappingMolecules

end module moduleModifyCoordinates

subroutine checkForBrokenMolecules(brokenMolecule)
  use moduleSystem 
  use moduleDistances
  use moduleMessages
  implicit none

  integer :: ibond, iatm, jatm, imol
  real(8) :: dij(3), distance2
  logical, intent(out) :: brokenMolecule

  brokenMolecule = .false.
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (ibond, iatm, jatm, dij, distance2, imol)
!$OMP DO
  do ibond=1,numberOfUniqueBonds 
    iatm = listOfUniqueBonds(1,ibond)
    jatm = listOfUniqueBonds(2,ibond)
    dij = frame % pos(1:3,iatm) - frame % pos(1:3,jatm)
    distance2 = sum(dij*dij)
    if (distance2 > 9.d0) then
      brokenMolecule = .true.
      imol = atomToMoleculeIndex(iatm)
      listOfMolecules(imol) % brokenBonds = .true.
    end if
  enddo
!$OMP END DO
!$OMP END PARALLEL
  
end subroutine checkForBrokenMolecules

subroutine reassembleAllMolecules()
  use moduleVariables
  use moduleSystem 
  use moduleDistances
  implicit none

  integer :: imol, idx, iatm, jatm
  real(8) :: dij(3), rtmp

  do imol=1,numberOfMolecules
    iatm = listOfMolecules(imol) % listOfAtoms(1)
    do idx=2,listOfMolecules(imol) % numberOfAtoms
      jatm = listOfMolecules(imol) % listOfAtoms(idx)
      dij(1:3) = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
      rtmp = computeDistanceSquaredPBC(dij)
      frame % pos(1:3,jatm) = frame % pos(1:3,iatm) + dij(1:3)
    enddo
  enddo

end subroutine reassembleAllMolecules

subroutine reassembleBrokenMolecules()
  use moduleVariables
  use moduleSystem 
  use moduleDistances
  implicit none

  integer :: imol, idx, iatm, jatm
  real(8) :: dij(3), rtmp

  do imol=1,numberOfMolecules
    if (.not. listOfMolecules(imol) % brokenBonds) cycle

    iatm = listOfMolecules(imol) % listOfAtoms(1)
    do idx=2,listOfMolecules(imol) % numberOfAtoms
      jatm = listOfMolecules(imol) % listOfAtoms(idx)
      dij(1:3) = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
      rtmp = computeDistanceSquaredPBC(dij)
      frame % pos(1:3,jatm) = frame % pos(1:3,iatm) + dij(1:3)
    enddo
    listOfMolecules(imol) % brokenBonds = .false.
  enddo

end subroutine reassembleBrokenMolecules

subroutine computeMoleculesCOM(mol0)
  use moduleSystem 
  use moduleElements 

  implicit none

  integer, intent(in) :: mol0
  integer :: imol, idx, iatm
  real(8) :: xcom(3), mass, massTotal

  do imol=mol0+1,numberOfMolecules
    xcom = 0.0d0
    massTotal = 0.d0
    do idx=1,listOfMolecules(imol) % numberOfAtoms
      iatm = listOfMolecules(imol) % listOfAtoms(idx)
      ! mass = getElementMass(listOfMolecules(imol) % listOfLabels(idx))
      mass = frame % mass(iatm)
      xcom = xcom + frame % pos(1:3,iatm) * mass
      massTotal = massTotal + mass
    enddo
    ! listOfMolecules(imol) % centreOfMass(1:3) = xcom / listOfMolecules(imol) % numberOfAtoms
    listOfMolecules(imol) % centreOfMass(1:3) = xcom / massTotal
  enddo

end subroutine computeMoleculesCOM

subroutine deleteOverlapAllMolecules(rcut,removeSecond)
  use moduleSystem
  use moduleVariables
  use moduleDistances
  implicit none
  real(8), intent(in) :: rcut
  logical, intent(in) :: removeSecond

  logical, allocatable, dimension(:) :: lremove
  integer, target :: imol, jmol
  integer :: idx
  integer :: iatm, jatm, ineigh
  real(8) :: dij(3), rtmp, rmax
  integer :: n0
  integer, pointer :: kdx

  if (removeSecond) then
    kdx => jmol
  else
    kdx => imol
  end if

  allocate(lremove(numberOfMolecules), source=.false.)

  rmax = rcut**2
  m1 : do imol=1,numberOfMolecules
    if (lremove(imol)) cycle
    do idx=1,listOfMolecules(imol) % numberOfAtoms
      iatm = listOfMolecules(imol) % listOfAtoms(idx)
      m2 : do ineigh=1,nneigh(iatm)
        jatm = lneigh(ineigh,iatm)
        jmol = atomtomoleculeindex(jatm)
        if (imol == jmol) cycle
        if (lremove(jmol)) cycle
        dij = frame % pos(1:3,iatm) - frame % pos(1:3,jatm)
        rtmp = computeDistanceSquaredPBC(dij) 
        if (rtmp < rmax) then
          lremove(kdx) = .true.
          if (kdx == imol) then
            cycle m1
          else
            cycle m2
          end if
        end if
      end do m2
    end do
  end do m1

  n0 = 0
  do imol=1,numberOfMolecules
    if (lremove(imol)) cycle
    do idx=1,listOfMolecules(imol) % numberOfAtoms
      iatm = listOfMolecules(imol) % listOfAtoms(idx)

      n0 = n0 + 1
      frame % lab(n0) = frame % lab(iatm)
      frame % chg(n0) = frame % chg(iatm)
      frame % pos(:,n0) = frame % pos(:,iatm)
    end do
  end do
  frame % natoms = n0

end subroutine deleteOverlapAllMolecules

subroutine deleteOverlapAllAtoms(rcut,removeSecond)
  use moduleSystem
  use moduleVariables
  use moduleDistances
  implicit none
  real(8), intent(in) :: rcut
  logical, intent(in) :: removeSecond

  logical, allocatable, dimension(:) :: lremove
  integer, target :: iatm, jatm, ineigh
  real(8) :: dij(3), rtmp, rmax
  integer :: n0
  integer, pointer :: idx

  rmax = rcut**2

  n0 = frame % natoms
  allocate(lremove(n0), source=.false.)

  if (removeSecond) then
    idx => jatm
  else
    idx => iatm
  end if

  m1 : do iatm=1,n0
    if (lremove(iatm)) cycle

    m2 : do ineigh=1,nneigh(iatm)
      jatm=lneigh(ineigh,iatm)
      if (lremove(jatm)) cycle

      dij = frame % pos(1:3,iatm) - frame % pos(1:3,jatm)
      rtmp = computeDistanceSquaredPBC(dij) 
      if (rtmp < rmax) then
        lremove(idx) = .true.
        if (idx == iatm) then
          cycle m1
        else
          cycle m2
        end if
      end if

    end do m2
  end do m1

  n0 = 0
  do iatm=1,frame % natoms
    if (lremove(iatm)) cycle
    n0 = n0 + 1
    frame % lab(n0) = frame % lab(iatm)
    frame % chg(n0) = frame % chg(iatm)
    frame % pos(:,n0) = frame % pos(:,iatm)
  end do
  frame % natoms = n0

end subroutine deleteOverlapAllAtoms
