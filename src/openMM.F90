! ! Copyright (c) 2021, Paolo Raiteri, Curtin University.
! ! All rights reserved.
! ! 
! ! This program is free software; you can redistribute it and/or modify it 
! ! under the terms of the GNU General Public License as published by the 
! ! Free Software Foundation; either version 3 of the License, or 
! ! (at your option) any later version.
! !  
! ! Redistribution and use in source and binary forms, with or without 
! ! modification, are permitted provided that the following conditions are met:
! ! 
! ! * Redistributions of source code must retain the above copyright notice, 
! !   this list of conditions and the following disclaimer.
! ! * Redistributions in binary form must reproduce the above copyright notice, 
! !   this list of conditions and the following disclaimer in the documentation 
! !   and/or other materials provided with the distribution.
! ! * Neither the name of the <ORGANIZATION> nor the names of its contributors 
! !   may be used to endorse or promote products derived from this software 
! !   without specific prior written permission.
! ! 
! ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! ! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! ! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! ! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! ! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! ! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! ! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! ! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! ! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! ! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ! 
module moduleOpenMM

  real(8), allocatable, dimension(:) :: x, y, z
  integer, allocatable, dimension(:) :: qType
  integer, allocatable, dimension(:) :: vdwType
  ! real(8), allocatable, dimension(:) :: moleculeDipole
  ! real(8), allocatable, dimension(:) :: electricField
  ! real(8), allocatable, dimension(:) :: moleculeDipole

  type amoebaVdw_Type
    real(8) :: sigma
    real(8) :: epsilon
    real(8) :: reduction
  end type
  type(amoebaVdw_Type), dimension(1000) :: amoebaVdwType
  real(8), allocatable, dimension(:,:) :: amoebaVdwParameters
  integer, allocatable, dimension(:) :: ivIndex

  type amoebaMultipole_Type
    integer :: axisAtom(3)
    character(len=8) :: axisType
    integer :: axisID
    real(8) :: multipoles(13)
    real(8) :: thole
    real(8) :: polarisation
    integer :: pg(8)
  end type
  type(amoebaMultipole_Type), dimension(1000) :: amoebaMultipoleType

  integer, allocatable, dimension(:,:) :: axes
  integer :: numberOfPolarisationGroups
  integer, allocatable, dimension(:) :: groupID
  integer, allocatable, dimension(:,:) :: polarisationGroups

  type amoebaBondType
    integer :: atomTypes(2)
    real(8) :: kappa
    real(8) :: distance
  end type
  integer :: numberOfAmoebaBonds
  type(amoebaBondType), allocatable, dimension(:) :: amoebaBonds
  real(8), allocatable, dimension(:,:) :: amoebaBondsParameters

  type amoebaAngleType
    integer :: atomTypes(3)
    real(8) :: kappa
    real(8) :: theta
  end type
  integer :: numberOfAmoebaAngles
  type(amoebaAngleType), allocatable, dimension(:) :: amoebaAngles
  real(8), allocatable, dimension(:,:) :: amoebaAnglesParameters

  type amoebaTorsionType
    integer :: atomTypes(4)
    real(8) :: kappa(3)
    real(8) :: phi(3)
    integer :: fold(3)
  end type
  integer :: numberOfAmoebaTorsions
  type(amoebaTorsionType), allocatable, dimension(:) :: amoebaTorsions
  real(8), allocatable, dimension(:,:) :: amoebaTorsionsParameters

  type amoebaOutOfPlaneType
    integer :: atomTypes(4)
    real(8) :: kappa
  end type
  integer :: numberOfAmoebaOutOfPlane
  type(amoebaOutOfPlaneType), allocatable, dimension(:) :: amoebaOutOfPlane
  real(8), allocatable, dimension(:) :: amoebaOutOfPlaneParameters

  type amoebaStretchBendType
    integer :: atomTypes(3)
    real(8) :: k1, k2
  end type
  integer :: numberOfAmoebaStretchBend
  type(amoebaStretchBendType), allocatable, dimension(:) :: amoebaStretchBend
  real(8), allocatable, dimension(:,:) :: amoebaStretchBendParameters

  type ureyBradleyType
    integer :: atomTypes(3)
    real(8) :: kappa
    real(8) :: distance
  end type
  integer :: numberOfUreyBradley
  type(ureyBradleyType), allocatable, dimension(:) :: ureyBradley
  real(8), allocatable, dimension(:,:) :: ureyBradleyParameters

  character(len=6) :: polarisationType = "mutual"
  real(8) :: cutoffOpenMM = 10.d0
  real(8) :: ewaldPrecision = 1d-5

  contains

#ifdef GPTA_OPENMM

  subroutine initialiseOpenMM(forcefieldFile)
    use moduleSystem
    use moduleElements
    implicit none
    character(len=STRLEN), intent(in) :: forcefieldFile
    integer :: iounit
    integer :: iatm
    integer :: i

    integer, allocatable, dimension(:) :: buckList

    logical, save :: firstTimeIn = .true.
    real(8), allocatable, dimension(:) :: atomicMass
    if (.not. firstTimeIn) return
    firstTimeIn = .false.
    
    call computeConnectivity() 

    open(newunit=iounit,file=forcefieldFile,status='old')
    call readTinkerParameters(iounit)

    allocate(x(frame % natoms), y(frame % natoms), z(frame % natoms))
    ! allocate(electricField(frame % natoms))
    
    ! Atom types assigment done via --select +s O2,H2 +tq 1,2
    ! allocate(qType(frame % natoms))
    allocate(vdwType(frame % natoms))
    vdwType = qType

    call defineAmoebaBonds()
    call defineAmoebaAngles()
    call defineAmoebaTorsions()
    call defineAmoebaOutOfPlane()
    call defineAmoebaStretchBend()
    call defineUreyBradley()

    ! Define the axes for the multipoles
    allocate(axes(3,numberOfAtoms), source=0)
    call defineMultipoleAxes()
    call definePolarisationGroups()

    ! set vdw interactions
    call defineAmoebaVdw()

    ! Initilise OpenMM
    allocate(atomicMass(frame % natoms))
    do iatm=1,frame % natoms
      atomicMass(iatm) = getElementMass(frame % lab(iatm))
    enddo
    call openmm_create_atoms(frame % natoms, atomicMass)

    call openmm_setup_multipoles(cutoffOpenMM, ewaldPrecision, polarisationType//C_NULL_CHAR)
    
    do iatm=1,numberOfAtoms
      i = groupID(iatm)
      call openmm_add_amoeba_multipoles_to_atom( &
        qType(iatm), &
        amoebaMultipoleType(qType(iatm)) % multipoles, &
        amoebaMultipoleType(qType(iatm)) % polarisation, &
        amoebaMultipoleType(qType(iatm)) % thole, &
        amoebaMultipoleType(qType(iatm)) % axisID, &
        axes(:,iatm), &
        polarisationGroups(:,i) )
    end do
  
    call openmm_add_amoeba_bonds(numberOfUniqueBonds, listOfUniqueBonds, amoebaBondsParameters)
    call openmm_add_amoeba_angles(numberOfUniqueAngles, listOfUniqueAngles, amoebaAnglesParameters)
    call openmm_add_amoeba_torsions(numberOfUniqueTorsions, listOfUniqueTorsions, amoebaTorsionsParameters)
    call openmm_add_amoeba_out_of_plane(numberOfUniqueOutOfPlane, listOfUniqueOutOfPlane, amoebaOutOfPlaneParameters)
    
    call openmm_add_amoeba_stretchbend(numberOfUniqueAngles, listOfUniqueAngles, amoebaStretchBendParameters)
    call openmm_add_amoeba_ureybradley(numberOfUniqueAngles, listOfUniqueAngles, ureyBradleyParameters)

    call openmm_add_amoeba_vdw(numberOfAtoms, ivIndex, amoebaVdwParameters, &
                               shape(listOfCovalentBondsPerAtom), &
                               numberOfCovalentBondsPerAtom,    listOfCovalentBondsPerAtom, &
                               numberOf_13_InteractionsPerAtom, listOf_13_InteractionsAtom)

    call openmm_setup_exclusions(shape(listOfCovalentBondsPerAtom), &
                                 numberOfCovalentBondsPerAtom,    listOfCovalentBondsPerAtom, &
                                 numberOf_13_InteractionsPerAtom, listOf_13_InteractionsAtom, &
                                 numberOf_14_InteractionsPerAtom, listOf_14_InteractionsAtom)

    ! allocate(buckList(numberOfAtoms), source=0)
    ! buckList(1) = 1
    ! buckList(3) = 2
    ! buckList(4) = 2
    ! buckList(5) = 2
    ! call openmm_add_custom_non_bonded(numberOfAtoms, buckList)

  end subroutine initialiseOpenMM

  subroutine readTinkerParameters(iounit)
    implicit none
    integer, intent(in) :: iounit

    ! call readTinkerAtoms(iounit)
    call readTinkerVdw(iounit)

    call readTinkerBonds(iounit)
    call readTinkerAngles(iounit)
    call readTinkerTorsions(iounit)
    call readTinkerOutOfPlane(iounit)
    call readTinkerStretchBend(iounit)
    call readTinkerUreyBradley(iounit)

    call readTinkerMultipoles(iounit)
    call readTinkerPolarisation(iounit)

    close(iounit)

  end subroutine readTinkerParameters

  ! subroutine readTinkerAtoms(iounit)
  !   implicit none
  !   integer, intent(in) :: iounit
  !   integer :: ios, iatm
  !   character(len=100) :: line
  !   character(len=20) :: word

  !   rewind(iounit)
  !   do
  !     read(iounit,'(a100)',iostat=ios) line
  !     if (ios/=0) exit
      
  !     if (line(1:5) == 'atom ' .or. line(1:5) == 'ATOM ') then
  !       read(line,*) word, iatm, amoebaVdwType(iatm) % vdwType
  !     end if
  !   end do

  ! end subroutine readTinkerAtoms
  
!-------------------------------------------------------------------------------------------------------------!
  subroutine readTinkerVdw(iounit)
    implicit none
    integer, intent(in) :: iounit
    integer :: ios, iatm
    real(8) :: epsilon, sigma, reduction
    character(len=100) :: line
    character(len=20) :: word

    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      
      if (line(1:4) == 'vdw ' .or. line(1:4) == 'VDW ') then
        read(line,*,iostat=ios) word, iatm, sigma, epsilon, reduction
        if (ios/=0) then
          read(line,*) word, iatm, sigma, epsilon
          reduction = 1.d0
        endif
        amoebaVdwType(iatm) % sigma = sigma / 10.d0       ! conversion angstrom -> nm
        amoebaVdwType(iatm) % epsilon = epsilon * 4.184d0 ! conversion kcal -> kJ
        amoebaVdwType(iatm) % reduction = reduction       ! conversion kcal -> kJ
      end if
    end do
    
  end subroutine readTinkerVdw

  subroutine defineAmoebaVdw()
    use moduleSystem
    use moduleElements
    implicit none
    integer :: i, iatm

    allocate(amoebaVdwParameters(3,numberOfAtoms), source=0.d0)
    allocate(ivIndex(numberOfAtoms))

    do i=1,numberOfAtoms
      if (getAtomicNumber(frame % lab(i)) == 1) then
        ivIndex(i) = listOfCovalentBondsPerAtom(1,i) ! index of the atom to which the H is attached
      else
        ivIndex(i) = i
      end if

      iatm = vdwType(i)
      amoebaVdwParameters(1,i) = amoebaVdwType(iatm) % sigma
      amoebaVdwParameters(2,i) = amoebaVdwType(iatm) % epsilon
      amoebaVdwParameters(3,i) = amoebaVdwType(iatm) % reduction
    end do

  end subroutine defineAmoebaVdw

!-------------------------------------------------------------------------------------------------------------!
  subroutine readTinkerBonds(iounit)
    implicit none
    integer, intent(in) :: iounit
    integer :: ios, iatm, jatm, i
    real(8) :: K, b
    character(len=100) :: line
    character(len=20) :: word

    numberOfAmoebaBonds = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:5) == 'bond ' .or. line(1:5) == 'BOND ') then
        numberOfAmoebaBonds = numberOfAmoebaBonds + 1
      end if
    end do
    allocate(amoebaBonds(numberOfAmoebaBonds))

    i = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:5) == 'bond ' .or. line(1:5) == 'BOND ') then
        i = i + 1
        read(line,*) word, iatm, jatm, k, b
        amoebaBonds(i) % atomTypes = [iatm,jatm]
        amoebaBonds(i) % kappa = k * 418.4d0  ! conversion kcal/angstrom^2 -> kJ/nm^2
        amoebaBonds(i) % distance = b / 10.d0 ! conversion angstrom -> nm
      end if
    end do

  end subroutine readTinkerBonds

  subroutine defineAmoebaBonds()
    use moduleSystem
    implicit none
    integer :: ibond, i
    integer :: iatm, jatm
    character(len=8) :: str0
    character(len=8), allocatable, dimension(:) :: str1, str2

    allocate(str1(numberOfAmoebaBonds), source="        ")
    allocate(str2(numberOfAmoebaBonds), source="        ")
    do i=1,numberOfAmoebaBonds
      write(str1(i),'(2i4)') amoebaBonds(i) % atomTypes(1), amoebaBonds(i) % atomTypes(2)
      write(str2(i),'(2i4)') amoebaBonds(i) % atomTypes(2), amoebaBonds(i) % atomTypes(1)
    end do

    allocate(amoebaBondsParameters(2,numberOfUniqueBonds), source=0.d0)
    do ibond=1,numberOfUniqueBonds
      iatm = listOfUniqueBonds(1,ibond) 
      jatm = listOfUniqueBonds(2,ibond) 
      write(str0,'(2i4)') qType(iatm), qType(jatm)
      do i=1,numberOfAmoebaBonds
        if (str0==str1(i) .or. str0==str2(i)) then
          amoebaBondsParameters(1,ibond) = amoebaBonds(i) % distance
          amoebaBondsParameters(2,ibond) = amoebaBonds(i) % kappa
          exit
        end if
      end do
    end do

  end subroutine defineAmoebaBonds

!-------------------------------------------------------------------------------------------------------------!
  subroutine readTinkerAngles(iounit)
    implicit none
    integer, intent(in) :: iounit
    integer :: ios, iatm, jatm, katm, i
    real(8) :: K, b
    character(len=100) :: line
    character(len=20) :: word

    numberOfAmoebaAngles = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:6) == 'angle' .or. line(1:6) == 'ANGLE ') then
        numberOfAmoebaAngles = numberOfAmoebaAngles + 1
      end if
    end do
    allocate(amoebaAngles(numberOfAmoebaAngles))

    i = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:6) == 'angle' .or. line(1:6) == 'ANGLE ') then
        i = i + 1
        read(line,*) word, iatm, jatm, katm, k, b
        amoebaAngles(i) % atomTypes = [iatm,jatm,katm]
        amoebaAngles(i) % kappa = k * 1.27452d-3  ! conversion kcal/rad^2 -> kj/mol/deg^2
        amoebaAngles(i) % theta = b ! no conversion degrees!
      end if
    end do

  end subroutine readTinkerAngles

  subroutine defineAmoebaAngles()
    use moduleSystem
    implicit none
    integer :: iangle, i
    integer :: iatm, jatm, katm
    character(len=12) :: str0
    character(len=12), allocatable, dimension(:) :: str1, str2

    allocate(str1(numberOfAmoebaAngles), source="            ")
    allocate(str2(numberOfAmoebaAngles), source="            ")
    do i=1,numberOfAmoebaAngles
      write(str1(i),'(3i4)') amoebaAngles(i) % atomTypes(1), amoebaAngles(i) % atomTypes(2), amoebaAngles(i) % atomTypes(3)
      write(str2(i),'(3i4)') amoebaAngles(i) % atomTypes(3), amoebaAngles(i) % atomTypes(2), amoebaAngles(i) % atomTypes(1)
    end do

    allocate(amoebaAnglesParameters(2,numberOfUniqueAngles), source=0.d0)
    do iangle=1,numberOfUniqueAngles
      iatm = listOfUniqueAngles(1,iangle) 
      jatm = listOfUniqueAngles(2,iangle) 
      katm = listOfUniqueAngles(3,iangle) 
      write(str0,'(3i4)') qType(iatm), qType(jatm), qType(katm)
      do i=1,numberOfAmoebaAngles
        if (str0==str1(i) .or. str0==str2(i)) then
          amoebaAnglesParameters(1,iangle) = amoebaAngles(i) % theta
          amoebaAnglesParameters(2,iangle) = amoebaAngles(i) % kappa
          exit
        end if
      end do
    end do

  end subroutine defineAmoebaAngles

!-------------------------------------------------------------------------------------------------------------!
  subroutine readTinkerTorsions(iounit)
    implicit none
    integer, intent(in) :: iounit
    integer :: ios, iatm, jatm, katm, latm, i
    real(8) :: k(3), b(3)
    integer :: f(3)
    character(len=100) :: line
    character(len=20) :: word

    numberOfAmoebaTorsions = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:9) == 'torsion ' .or. line(1:9) == 'TORSION ') then
        numberOfAmoebaTorsions = numberOfAmoebaTorsions + 1
      end if
    end do
    allocate(amoebaTorsions(numberOfAmoebaTorsions))

    i = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:9) == 'torsion ' .or. line(1:9) == 'TORSION ') then
        i = i + 1
        read(line,*) word, iatm, jatm, katm, latm, k(1), b(1), f(1), k(2), b(2), f(2), k(3), b(3), f(3)
        amoebaTorsions(i) % atomTypes = [iatm,jatm,katm,latm]
        amoebaTorsions(i) % kappa = k * 2.092                   ! conversion kcal/rad -> kj/mol/deg (includes the 1/2)
        amoebaTorsions(i) % phi   = b / 180. * 3.14159265358979 ! onversion degrees -> radians
      end if
    end do

  end subroutine readTinkerTorsions
  
  subroutine defineAmoebaTorsions()
    use moduleSystem
    implicit none
    integer :: iangle, i
    integer :: iatm, jatm, katm, latm
    character(len=16) :: str0
    character(len=16), allocatable, dimension(:) :: str1, str2

    allocate(str1(numberOfAmoebaTorsions), source="                ")
    allocate(str2(numberOfAmoebaTorsions), source="                ")
    do i=1,numberOfAmoebaTorsions
      write(str1(i),'(4i4)') amoebaTorsions(i) % atomTypes(1), amoebaTorsions(i) % atomTypes(2), amoebaTorsions(i) % atomTypes(3), amoebaTorsions(i) % atomTypes(4)
      write(str2(i),'(4i4)') amoebaTorsions(i) % atomTypes(4), amoebaTorsions(i) % atomTypes(3), amoebaTorsions(i) % atomTypes(2), amoebaTorsions(i) % atomTypes(1)
    end do

    allocate(amoebaTorsionsParameters(6,numberOfUniqueTorsions), source=0.d0)
    do iangle=1,numberOfUniqueTorsions
      iatm = listOfUniqueTorsions(1,iangle) 
      jatm = listOfUniqueTorsions(2,iangle) 
      katm = listOfUniqueTorsions(3,iangle) 
      latm = listOfUniqueTorsions(4,iangle) 
      write(str0,'(4i4)') qType(iatm), qType(jatm), qType(katm), qType(latm)
      do i=1,numberOfAmoebaTorsions
        if (str0==str1(i) .or. str0==str2(i)) then
          amoebaTorsionsParameters(1:3,iangle) = amoebaTorsions(i) % phi
          amoebaTorsionsParameters(4:6,iangle) = amoebaTorsions(i) % kappa
          exit
        end if
      end do
    end do

  end subroutine defineAmoebaTorsions

!-------------------------------------------------------------------------------------------------------------!
  subroutine readTinkerOutOfPlane(iounit)
    implicit none
    integer, intent(in) :: iounit
    integer :: ios, iatm, jatm, katm, latm, i
    real(8) :: k
    character(len=100) :: line
    character(len=20) :: word

    numberOfAmoebaOutOfPlane = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:9) == 'opbend ' .or. line(1:9) == 'OPBEND ') then
        numberOfAmoebaOutOfPlane = numberOfAmoebaOutOfPlane + 1
      end if
    end do
    allocate(amoebaOutOfPlane(numberOfAmoebaOutOfPlane))

    i = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:9) == 'opbend ' .or. line(1:9) == 'OPBEND ') then
        i = i + 1
        read(line,*) word, iatm, jatm, katm, latm, k
        amoebaOutOfPlane(i) % atomTypes = [iatm,jatm,katm,latm]
        amoebaOutOfPlane(i) % kappa = k * 1.27452e-3   ! conversion kcal/rad^2 -> kj/mol/deg^2 (includes the 1/2)
      end if
    end do

  end subroutine readTinkerOutOfPlane
  
  subroutine defineAmoebaOutOfPlane()
    use moduleSystem
    implicit none
    integer :: iangle, i, j
    integer :: iatm, jatm, katm, latm
    character(len=16) :: str0
    character(len=16), allocatable, dimension(:,:) :: str

    allocate(str(6,numberOfAmoebaOutOfPlane), source="                ")
    do i=1,numberOfAmoebaOutOfPlane
      write(str(1,i),'(4i4)') amoebaOutOfPlane(i) % atomTypes(1), amoebaOutOfPlane(i) % atomTypes(2), amoebaOutOfPlane(i) % atomTypes(3), amoebaOutOfPlane(i) % atomTypes(4)
      write(str(2,i),'(4i4)') amoebaOutOfPlane(i) % atomTypes(1), amoebaOutOfPlane(i) % atomTypes(2), amoebaOutOfPlane(i) % atomTypes(4), amoebaOutOfPlane(i) % atomTypes(3)
      write(str(3,i),'(4i4)') amoebaOutOfPlane(i) % atomTypes(3), amoebaOutOfPlane(i) % atomTypes(2), amoebaOutOfPlane(i) % atomTypes(1), amoebaOutOfPlane(i) % atomTypes(4)
      write(str(4,i),'(4i4)') amoebaOutOfPlane(i) % atomTypes(3), amoebaOutOfPlane(i) % atomTypes(2), amoebaOutOfPlane(i) % atomTypes(4), amoebaOutOfPlane(i) % atomTypes(1)
      write(str(5,i),'(4i4)') amoebaOutOfPlane(i) % atomTypes(4), amoebaOutOfPlane(i) % atomTypes(2), amoebaOutOfPlane(i) % atomTypes(1), amoebaOutOfPlane(i) % atomTypes(3)
      write(str(6,i),'(4i4)') amoebaOutOfPlane(i) % atomTypes(4), amoebaOutOfPlane(i) % atomTypes(2), amoebaOutOfPlane(i) % atomTypes(3), amoebaOutOfPlane(i) % atomTypes(1)
    end do

    allocate(amoebaOutOfPlaneParameters(numberOfUniqueOutOfPlane), source=0.d0)
    do iangle=1,numberOfUniqueOutOfPlane
      iatm = listOfUniqueOutOfPlane(2,iangle) 
      jatm = listOfUniqueOutOfPlane(1,iangle) 
      katm = listOfUniqueOutOfPlane(3,iangle) 
      latm = listOfUniqueOutOfPlane(4,iangle) 
      write(str0,'(4i4)') qType(iatm), qType(jatm), qType(katm), qType(latm)
      do i=1,numberOfAmoebaOutOfPlane
        do j=1,6
          if (str0==str(j,i)) then
            amoebaOutOfPlaneParameters(iangle) = amoebaOutOfPlane(i) % kappa 
            exit
          end if
        end do
      end do
    end do

  end subroutine defineAmoebaOutOfPlane

!-------------------------------------------------------------------------------------------------------------!
  subroutine readTinkerStretchBend(iounit)
    implicit none
    integer, intent(in) :: iounit
    integer :: ios, iatm, jatm, katm, i
    real(8) :: k1, k2
    character(len=100) :: line
    character(len=20) :: word

    numberOfAmoebaStretchBend = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:7) == 'strbnd' .or. line(1:7) == 'STRBND ') then
        numberOfAmoebaStretchBend = numberOfAmoebaStretchBend + 1
      end if
    end do
    allocate(amoebaStretchBend(numberOfAmoebaStretchBend))

    i = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:7) == 'strbnd' .or. line(1:7) == 'STRBND ') then
        i = i + 1
        read(line,*) word, iatm, jatm, katm, k1, k2
        amoebaStretchBend(i) % atomTypes = [iatm,jatm,katm]
        amoebaStretchBend(i) % k1 = k1 * 0.730246  ! conversion kcal/mol/angstrom/rad -> kJ/mol/nm/deg
        amoebaStretchBend(i) % k2 = k2 * 0.730246  ! conversion kcal/mol/angstrom/rad -> kJ/mol/nm/deg
      end if
    end do

  end subroutine readTinkerStretchBend

  subroutine defineAmoebaStretchBend()
    use moduleSystem
    implicit none
    integer :: iangle, i, j, k, l
    integer :: iatm, jatm, katm
    character(len=12) :: str0
    character(len=12), allocatable, dimension(:) :: str1, str2

    character(len=8), allocatable, dimension(:) :: bnd1, bnd2
    character(len=12), allocatable, dimension(:) :: ang1, ang2

    ! bonds
    allocate(bnd1(numberOfAmoebaBonds), source="        ")
    allocate(bnd2(numberOfAmoebaBonds), source="        ")
    do i=1,numberOfAmoebaBonds
      write(bnd1(i),'(2i4)') amoebaBonds(i) % atomTypes(1), amoebaBonds(i) % atomTypes(2)
      write(bnd2(i),'(2i4)') amoebaBonds(i) % atomTypes(2), amoebaBonds(i) % atomTypes(1)
    end do

    ! angles
    allocate(ang1(numberOfAmoebaAngles), source="            ")
    allocate(ang2(numberOfAmoebaAngles), source="            ")
    do i=1,numberOfAmoebaAngles
      write(ang1(i),'(3i4)') amoebaAngles(i) % atomTypes(1), amoebaAngles(i) % atomTypes(2), amoebaAngles(i) % atomTypes(3)
      write(ang2(i),'(3i4)') amoebaAngles(i) % atomTypes(3), amoebaAngles(i) % atomTypes(2), amoebaAngles(i) % atomTypes(1)
    end do

    ! stretch bend
    allocate(str1(numberOfAmoebaStretchBend), source="            ")
    allocate(str2(numberOfAmoebaStretchBend), source="            ")
    do i=1,numberOfAmoebaStretchBend
      write(str1(i),'(3i4)') amoebaStretchBend(i) % atomTypes(1), amoebaStretchBend(i) % atomTypes(2), amoebaStretchBend(i) % atomTypes(3)
      write(str2(i),'(3i4)') amoebaStretchBend(i) % atomTypes(3), amoebaStretchBend(i) % atomTypes(2), amoebaStretchBend(i) % atomTypes(1)
    end do

    allocate(amoebaStretchBendParameters(5,numberOfUniqueAngles), source=0.d0)
    do iangle=1,numberOfUniqueAngles
      iatm = listOfUniqueAngles(1,iangle) 
      jatm = listOfUniqueAngles(2,iangle) 
      katm = listOfUniqueAngles(3,iangle) 
      write(str0,'(3i4)') qType(iatm), qType(jatm), qType(katm)

      ! bond AB
      do i=1,numberOfAmoebaBonds
        if (str0(1:8)==bnd1(i) .or. str0(1:8)==bnd2(i)) exit
      end do

      ! bond BC
      do j=1,numberOfAmoebaBonds
        if (str0(5:12)==bnd1(j) .or. str0(5:12)==bnd2(j)) exit
      end do

      ! angles
      do k=1,numberOfAmoebaAngles
        if (str0==ang1(k) .or. str0==ang2(k)) exit
      end do
      ! if (i>numberOfAmoebaBonds .or. j>numberOfAmoebaBonds .or. k>numberOfAmoebaAngles) then
      !   write(0,*)"Cannot find strbnd terms : "//trim(str0),iatm,jatm,katm,qType(iatm), qType(jatm), qType(katm)
      !   ! stop
      ! endif

      ! stretch bend
      do l=1,numberOfAmoebaStretchBend
        if (str0==str1(l) .or. str0==str2(l)) then
          amoebaStretchBendParameters(1,iangle) = amoebaBonds(i) % distance
          amoebaStretchBendParameters(2,iangle) = amoebaBonds(j) % distance
          amoebaStretchBendParameters(3,iangle) = amoebaAngles(k) % theta / 180. * 3.14159265358979
          amoebaStretchBendParameters(4,iangle) = amoebaStretchBend(l) % k1
          amoebaStretchBendParameters(5,iangle) = amoebaStretchBend(l) % k2
          exit
        end if
      end do
    end do

  end subroutine defineAmoebaStretchBend

!-------------------------------------------------------------------------------------------------------------!
  subroutine readTinkerUreyBradley(iounit)
    implicit  none
    integer, intent(in) :: iounit
    integer :: ios, iatm, jatm, katm, i
    real(8) :: K, b
    character(len=100) :: line
    character(len=20) :: word

    numberOfUreyBradley = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:9) == 'ureybrad ' .or. line(1:9) == 'UREYBRAD ') then
        numberOfUreyBradley = numberOfUreyBradley + 1
      end if
    end do
    allocate(ureyBradley(numberOfUreyBradley))

    i = 0
    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      if (line(1:9) == 'ureybrad ' .or. line(1:9) == 'UREYBRAD ') then
        i = i + 1
        read(line,*) word, iatm, jatm, katm, k, b
        ureyBradley(i) % atomTypes = [iatm,jatm,katm]
        ureyBradley(i) % kappa = k * 418.4  ! conversion kcal/angstrom -> kj/nm
        ureyBradley(i) % distance = b / 10. ! conversion angstrom -> nm
      end if
    end do

  end subroutine readTinkerUreyBradley

  subroutine defineUreyBradley()
    use moduleSystem
    implicit none
    integer :: iangle, i
    integer :: iatm, jatm
    character(len=8) :: str0
    character(len=8), allocatable, dimension(:) :: str1, str2

    allocate(str1(numberOfUreyBradley), source="        ")
    allocate(str2(numberOfUreyBradley), source="        ")
    do i=1,numberOfUreyBradley
      write(str1(i),'(2i4)') ureyBradley(i) % atomTypes(1), ureyBradley(i) % atomTypes(3)
      write(str2(i),'(2i4)') ureyBradley(i) % atomTypes(3), ureyBradley(i) % atomTypes(1)
    end do

    allocate(ureyBradleyParameters(2,numberOfUniqueAngles), source=0.d0)
    do iangle=1,numberOfUniqueAngles
      iatm = listOfUniqueAngles(1,iangle) 
      jatm = listOfUniqueAngles(3,iangle) 
      write(str0,'(2i4)') qType(iatm), qType(jatm)
      do i=1,numberOfUreyBradley
        if (str0==str1(i) .or. str0==str2(i)) then
          ureyBradleyParameters(1,iangle) = ureyBradley(i) % distance
          ureyBradleyParameters(2,iangle) = ureyBradley(i) % kappa
          exit
        end if
      end do
    end do

  end subroutine defineUreyBradley

!-------------------------------------------------------------------------------------------------------------!
  subroutine readTinkerMultipoles(iounit)
    implicit none
    integer, intent(in) :: iounit
    integer :: ios
    integer :: indices(4), i
    real(8) :: charge, dipole(3), quadrupole(3,3)
    character(len=100) :: line
    character(len=20) :: word
    character(len=8) :: axis

    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit

      ! Read multipoles
      if (line(1:10) == 'multipole ' .or. line(1:10) == 'MULTIPOLE ') then
        indices = 0

        i = 5
        ios = 1
        do while(ios/=0)
          i = i - 1
          read(line,*,iostat=ios) word, indices(1:i), charge
        end do

        read(iounit,*,iostat=ios) dipole(1:3)

        read(iounit,*,iostat=ios) quadrupole(1,1)
        read(iounit,*,iostat=ios) quadrupole(1:2,2)
        read(iounit,*,iostat=ios) quadrupole(1:3,3)
        quadrupole(2,1) = quadrupole(1,2)
        quadrupole(3,1) = quadrupole(1,3)
        quadrupole(3,2) = quadrupole(2,3)

        axis = 'Z-then-X'
        if (indices(2) .eq. 0                        )  axis = 'None'
        if (indices(2) .ne. 0 .and. indices(3) .eq. 0)  axis = 'Z-Only'
        if (indices(2) .lt. 0 .or.  indices(3) .lt. 0)  axis = 'Bisector'
        if (indices(3) .lt. 0 .and. indices(4) .lt. 0)  axis = 'Z-Bisect'
        if (maxval(indices(2:4)) .lt. 0              )  axis = '3-Fold'

        indices = abs(indices)
        i = indices(1)
        amoebaMultipoleType(i) % axisAtom(1:3)       = indices(2:4)
        amoebaMultipoleType(i) % axisType         = axis  
        amoebaMultipoleType(i) % multipoles(1)    = charge
        amoebaMultipoleType(i) % multipoles(2:4)  = dipole *0.0529177208590237d0
        amoebaMultipoleType(i) % multipoles(5:13) = reshape(quadrupole, [9]) *0.000933428393638d0

        ! ZThenX = 0, Bisector = 1, ZBisect = 2, ThreeFold = 3, ZOnly = 4, NoAxisType = 5, LastAxisTypeIndex = 6
        if (axis == "Z-then-X") amoebaMultipoleType(i) % axisID = 0 ! ZThenX	
        if (axis == "Bisector") amoebaMultipoleType(i) % axisID = 1 ! Bisector	
        if (axis == "Z-Bisect") amoebaMultipoleType(i) % axisID = 2 ! ZBisect	
        if (axis == "3-Fold"  ) amoebaMultipoleType(i) % axisID = 3 ! ThreeFold	
        if (axis == "Z-Only"  ) amoebaMultipoleType(i) % axisID = 4 ! ZOnly	
        if (axis == "None"    ) amoebaMultipoleType(i) % axisID = 5 ! NoAxisType	
      end if
    end do
  end subroutine readTinkerMultipoles

  subroutine defineMultipoleAxes()
    use moduleSystem
    implicit none
    integer :: iatm, jatm, katm, latm
    integer :: j, k, l

    atm:do iatm=1,numberOfAtoms

      ! assign multipole parameters via only 1-2 connected atoms
      do j=1,numberOfCovalentBondsPerAtom(iatm)
        jatm = listOfCovalentBondsPerAtom(j,iatm)

        if (qType(jatm) == amoebaMultipoleType(qType(iatm)) % axisAtom(1)) then

          do k=1,numberOfCovalentBondsPerAtom(iatm)
            katm = listOfCovalentBondsPerAtom(k,iatm)

            if (qType(katm) == amoebaMultipoleType(qType(iatm)) % axisAtom(2) .and. katm/=jatm) then

              if (amoebaMultipoleType(qType(iatm)) % axisAtom(3) == 0) then
                ! axes are jatm and katm
                axes(1:2,iatm) = [jatm, katm]
                cycle atm
                ! Done for the atom
              else
                do l=1,numberOfCovalentBondsPerAtom(iatm)
                  latm = listOfCovalentBondsPerAtom(l,iatm)
                  if (qType(latm) == amoebaMultipoleType(qType(iatm)) % axisAtom(3) .and. latm/=jatm  .and. latm/=katm) then
                    ! axes are jatm, katm and latm
                    axes(1:3,iatm) = [jatm, katm, latm]
                    cycle atm
                    ! Done for the atom
                  end if
                end do
              end if

            end if
          end do

        end if
      end do

      ! assign multipole parameters via 1-2 and 1-3 connected atoms
      do j=1,numberOfCovalentBondsPerAtom(iatm)
        jatm = listOfCovalentBondsPerAtom(j,iatm)
        if (qType(jatm) == amoebaMultipoleType(qType(iatm)) % axisAtom(1)) then

          do k=1,numberOfCovalentBondsPerAtom(jatm)
            katm = listOfCovalentBondsPerAtom(k,jatm)
            if (qType(katm) == amoebaMultipoleType(qType(iatm)) % axisAtom(2) .and. katm/=iatm) then

              if (amoebaMultipoleType(qType(iatm)) % axisAtom(3) == 0) then
                ! axes are jatm and katm
                axes(1:2,iatm) = [jatm, katm]
                cycle atm
                ! Done for the atom
              else
                do l=1,numberOfCovalentBondsPerAtom(jatm)
                  latm = listOfCovalentBondsPerAtom(l,jatm)
                  if (qType(latm) == amoebaMultipoleType(qType(iatm)) % axisAtom(3) .and. latm/=iatm  .and. latm/=katm) then
                    ! axes are jatm, katm and latm
                    axes(1:3,iatm) = [jatm, katm, latm]
                    cycle atm
                    ! Done for the atom
                  end if
                end do
              end if

            end if
          end do
          
        end if
      end do

      ! assign multipole parameters via only a z-defining atom
      do j=1,numberOfCovalentBondsPerAtom(iatm)
        jatm = listOfCovalentBondsPerAtom(j,iatm)
        if (qType(jatm) == amoebaMultipoleType(qType(iatm)) % axisAtom(1)) then
          if (amoebaMultipoleType(qType(iatm)) % axisAtom(2) == 0) then
            ! axes are jatm
            axes(1,iatm) = jatm
            cycle atm
            ! Done for the atom
          end if
        end if
      end do

      if (amoebaMultipoleType(qType(iatm)) % axisAtom(1) /=0) then
        write(0,*)"something has gone wrong" 
        stop
      endif
    end do atm

    ! do iatm=1,3
    !   write(0,*)iatm, frame % lab(iatm)
    !   write(0,*)qType(iatm)
    !   write(0,*)amoebaMultipoleType(qType(iatm)) % multipoles(1)
    !   write(0,*)amoebaMultipoleType(qType(iatm)) % multipoles(2:4)
    !   write(0,*)amoebaMultipoleType(qType(iatm)) % multipoles(5:7)
    !   write(0,*)amoebaMultipoleType(qType(iatm)) % multipoles(8:10)
    !   write(0,*)amoebaMultipoleType(qType(iatm)) % multipoles(11:13)
    !   write(0,*)amoebaMultipoleType(qType(iatm)) % polarisation
    !   write(0,*)amoebaMultipoleType(qType(iatm)) % thole
    !   write(0,*)amoebaMultipoleType(qType(iatm)) % axisType
    !   write(0,*)axes(:,iatm)
    ! enddo

  end subroutine defineMultipoleAxes

!-------------------------------------------------------------------------------------------------------------!
  subroutine readTinkerPolarisation(iounit)
    implicit none
    integer, intent(in) :: iounit
    integer :: ios, iatm, i
    character(len=100) :: line
    character(len=20) :: word
    real(8) :: polarisation, thole
    integer :: pg(8)

    rewind(iounit)
    do
      read(iounit,'(a100)',iostat=ios) line
      if (ios/=0) exit
      
      if (line(1:9) == 'polarize ' .or. line(1:9) == 'POLARIZE ') then
        i = 0
        polarisation = 0.d0
        thole = 0.d0
        pg = 0
        ios = 0
        do while (ios==0)
          i = i + 1
          read(line,*,iostat=ios) word, iatm, polarisation, thole, pg(1:i)
        end do

        amoebaMultipoleType(iatm) % polarisation =  polarisation / 1000.d0
        amoebaMultipoleType(iatm) % thole = thole
        amoebaMultipoleType(iatm) % pg = pg
      end if
    end do
  end subroutine readTinkerPolarisation

  subroutine definePolarisationGroups()
    use moduleSystem
    implicit none
    integer :: iatm, jatm
    integer :: i, j, ipg
    logical, allocatable, dimension(:) :: lused

    allocate(lused(numberOfatoms), source=.false.)

    numberOfPolarisationGroups = 0
    allocate(groupID(numberOfAtoms), source=0)
    allocate(polarisationGroups(0:8,numberOfAtoms), source=0)
    do iatm=1,numberOfAtoms
      if (.not. lused(iatm)) then
        lused(iatm) = .true.
        numberOfPolarisationGroups = numberOfPolarisationGroups + 1
        groupID(iatm) = numberOfPolarisationGroups
        polarisationGroups(0,numberOfPolarisationGroups) = 1 
        polarisationGroups(1,numberOfPolarisationGroups) = iatm
      endif
      ipg = groupID(iatm)

      do j=1,numberOfCovalentBondsPerAtom(iatm)
        jatm = listOfCovalentBondsPerAtom(j,iatm)
        if (any(qType(jatm)==amoebaMultipoleType(qType(iatm)) % pg)) then
          lused(jatm) = .true.
          groupID(jatm) = ipg
          i = polarisationGroups(0,ipg)
          if (all(polarisationGroups(1:i,ipg)/=jatm)) then
            i = i + 1
            polarisationGroups(i,ipg) = jatm
            polarisationGroups(0,ipg) = i
          end if
        end if
      end do
    end do  

  end subroutine definePolarisationGroups

!-------------------------------------------------------------------------------------------------------------!
  subroutine computeAmoebaMutipoles(outFile)
    use moduleSystem
    use moduleMessages
    use moduleVariables
    implicit none
    character(len=*), intent(in) :: outfile
    integer, save :: computedForFrame = 0

    ! Just to save time
    ! The dipoles are computed for all molecules and stored on the C++ side
    if (computedForFrame == frame % nframe) return

    computedForFrame = frame % nframe
    
    x = frame % pos(1,1:frame % natoms) / 10.d0
    y = frame % pos(2,1:frame % natoms) / 10.d0
    z = frame % pos(3,1:frame % natoms) / 10.d0
 
    call openmm_compute_multipoles(frame % nframe, frame % hmat/10.d0, x, y, z, qType, trim(outfile)//C_NULL_CHAR)
    
  end subroutine computeAmoebaMutipoles

  subroutine computeAmoebaElectricField(cartesian, electricField) 
    use moduleSystem 
    implicit none
    character(len=1) :: cartesian
    real(8), dimension(:), intent(out) :: electricField
    integer :: idir

    if (cartesian == "X" ) then
      idir = 1
    else if (cartesian == "Y" ) then
      idir = 2
    else if (cartesian == "Z" ) then
      idir = 3
    end if

    call openmm_compute_electric_field(qType, electricField, idir, z)
    electricField = electricField / 100.d0 * 14.39d0 ! conversion from [e*nm / nm^3] to [V/Angstrom]

  end subroutine computeAmoebaElectricField

  subroutine computeAmoebaDipoles(molecularDipoles)
    use moduleSystem 
    implicit none
    real(8), dimension(:,:), intent(out) :: molecularDipoles

    call openmm_compute_dipoles(qType, molecularDipoles, x, y, z)
    molecularDipoles = molecularDipoles * 48.032047d0 ! conversion from e*nm to Debye

  end subroutine computeAmoebaDipoles

#endif

end module moduleOpenMM
