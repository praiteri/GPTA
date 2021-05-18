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
module moduleLammps 

  use moduleVariables, only : cp
  implicit none

  logical :: lmpMap = .false.
  integer :: numberOfLammpsTypes
  character(len=100), dimension(100) :: mapLabelsTypes

  integer :: lammpsNumberOfAtoms = 0
  integer :: lammpsNumberOfBonds = 0
  integer :: lammpsNumberOfAngles = 0
  integer :: lammpsNumberOfTorsions = 0
  integer :: lammpsNumberOfOutOfPlane = 0
  
  character(cp), allocatable, dimension(:) :: lammpsAtomTypes
  character(cp), allocatable, dimension(:,:) :: lammpsBondTypes
  character(cp), allocatable, dimension(:,:) :: lammpsAngleTypes
  character(cp), allocatable, dimension(:,:) :: lammpsTorsionTypes
  character(cp), allocatable, dimension(:,:) :: lammpsImproperTypes
  character(len=20), allocatable, dimension(:) :: lammpsImproperForcefield
  character(len=8), dimension(6) :: lammpsImproperCentralFirst = ["distance", &
                                                                  "fourier ", &
                                                                  "harmonic", &
                                                                  "cvff    ", &
                                                                  "gulp    ", &
                                                                  "umbrella"]


  type :: lammpsForceFieldType
    character(cp), dimension(4) :: types
    integer :: id = 0
    integer :: number = 0
  end type

  type(lammpsForceFieldType), allocatable, dimension(:) :: lmp_atoms
  type(lammpsForceFieldType), allocatable, dimension(:) :: lmp_bonds
  type(lammpsForceFieldType), allocatable, dimension(:) :: lmp_angles
  type(lammpsForceFieldType), allocatable, dimension(:) :: lmp_torsions
  type(lammpsForceFieldType), allocatable, dimension(:) :: lmp_impropers

  logical :: ignoreMissingBonds = .false.
  logical :: ignoreMissingAngles = .false.
  logical :: ignoreMissingTorsions = .false.
  logical :: ignoreMissingImpropers = .false.
  logical :: checkMissingTerms = .false.

  integer, parameter :: maxIgnored = 100
  integer :: numberOfIgnoredTerms
  character(len=30), dimension(maxIgnored) :: missingTermsIgnored

end module moduleLammps 

subroutine readLammpsForcefieldFile(forcefieldFile)
  use moduleFiles
  use moduleLammps 
  use moduleStrings
  use moduleMessages 
  implicit none
  character(len=*), intent(in) :: forcefieldFile

  integer :: numberOfWords
  character(len=STRLEN), dimension(100) :: listOfWords

  integer :: iounit, ios
  integer, parameter :: lineLength = 100
  character(len=lineLength) :: line
  character(len=20) :: dummy
  character(len=8) :: fmt
  integer :: idx, jdx, itmp

  write(fmt,'("(a",i0,")")') lineLength
  open(newunit=iounit,file=forcefieldFile,status='old')

  do
    read(iounit,fmt,iostat=ios) line
    if (ios/=0) exit
    if (index(line,"#@") > 0) then
      if (index(line,"atom") > 0) then
        read (line(3:),*) lammpsNumberOfAtoms
        allocate(lammpsAtomTypes(lammpsNumberOfAtoms))
        allocate(lmp_atoms(lammpsNumberOfAtoms))

        idx = 0
        do while (idx < lammpsNumberOfAtoms)
          read(iounit,fmt,iostat=ios) line
          if (ios/=0) exit
          if (index(line,"variable") > 0) then
            idx = idx + 1
            read(line(9:),*) lammpsAtomTypes(idx)
          end if
        enddo
      
      else if (index(line,"bond") > 0) then
        read (line(3:),*) lammpsNumberOfBonds
        allocate(lammpsBondTypes(2,lammpsNumberOfBonds))
        allocate(lmp_bonds(lammpsNumberOfBonds))

        idx = 0
        jdx = 0
        do while (jdx < lammpsNumberOfBonds)
          read(iounit,fmt,iostat=ios) line
          if (ios/=0) exit
          if (len_trim(line) == 0 )cycle
          line = trim(line)
          if (line(1:1) == "#" .and. line(2:2) == "@") then
            idx = idx + 1
            call parse(line(3:),"-",listOfWords,numberOfWords)
            read(listOfWords(1),*) lammpsBondTypes(1,idx)
            read(listOfWords(2),*) lammpsBondTypes(2,idx)

          else if (line(1:10) == "bond_coeff") then
            jdx = jdx + 1
            lmp_bonds(jdx) % types(1:2) = lammpsBondTypes(1:2,idx)
            read(line,*)dummy , lmp_bonds(jdx) % ID
          end if

        enddo

      else if (index(line,"angle") > 0) then
        read (line(3:),*) lammpsNumberOfAngles
        allocate(lammpsAngleTypes(3,lammpsNumberOfAngles))
        allocate(lmp_angles(lammpsNumberOfAngles))
      
        idx = 0
        jdx = 0
        do while (jdx < lammpsNumberOfAngles)
          read(iounit,fmt,iostat=ios) line
          if (ios/=0) exit
          if (len_trim(line) == 0 )cycle
          line = trim(line)
          if (line(1:1) == "#" .and. line(2:2) == "@") then
            idx = idx + 1
            call parse(line(3:),"-",listOfWords,numberOfWords)
            read(listOfWords(1),*) lammpsAngleTypes(1,idx)
            read(listOfWords(2),*) lammpsAngleTypes(2,idx)
            read(listOfWords(3),*) lammpsAngleTypes(3,idx)

          else if (line(1:11) == "angle_coeff") then
            read(line,*)dummy , itmp
            if (any(lmp_angles(1:jdx) % ID == itmp)) cycle
            jdx = jdx + 1
            lmp_angles(jdx) % ID = itmp
            lmp_angles(jdx) % types(1:3) = lammpsAngleTypes(1:3,idx)

          end if
        enddo

      else if (index(line,"dihedral") > 0) then
        read (line(3:),*) lammpsNumberOfTorsions
        allocate(lammpsTorsionTypes(4,lammpsNumberOfTorsions))
        allocate(lmp_torsions(lammpsNumberOfTorsions))

        idx = 0
        jdx = 0
        do while (jdx < lammpsNumberOfTorsions)
          read(iounit,fmt,iostat=ios) line
          if (ios/=0) exit
          if (len_trim(line) == 0 )cycle
          line = trim(line)
          if (line(1:1) == "#" .and. line(2:2) == "@") then
            idx = idx + 1
            call parse(line(3:),"-",listOfWords,numberOfWords)
            read(listOfWords(1),*) lammpsTorsionTypes(1,idx)
            read(listOfWords(2),*) lammpsTorsionTypes(2,idx)
            read(listOfWords(3),*) lammpsTorsionTypes(3,idx)
            read(listOfWords(4),*) lammpsTorsionTypes(4,idx)
            
          else if (line(1:14) == "dihedral_coeff") then

            jdx = jdx + 1
            lmp_torsions(jdx) % types(1:4) = lammpsTorsionTypes(1:4,idx)
            read(line,*)dummy , lmp_torsions(jdx) % ID
            
          end if
        enddo

      else if (index(line,"improper") > 0) then
        read (line(3:),*) lammpsNumberOfOutOfPlane
        allocate(lammpsImproperTypes(4,lammpsNumberOfOutOfPlane))
        allocate(lammpsImproperForcefield(lammpsNumberOfOutOfPlane))
        allocate(lmp_impropers(lammpsNumberOfOutOfPlane))

        idx = 0
        jdx = 0
        do while (jdx < lammpsNumberOfOutOfPlane)
          read(iounit,fmt,iostat=ios) line
          if (ios/=0) exit
          if (len_trim(line) == 0 )cycle
          line = trim(line)
          
          ! Default Improper
          if (line(1:14) == "improper_style") then
            call parse(line," ",listOfWords,numberOfWords)
            if (numberOfWords == 2) read(listOfWords(2),*) lammpsImproperForcefield(1)
            lammpsImproperForcefield(2:) = lammpsImproperForcefield(1)
          end if
          
          if (line(1:1) == "#" .and. line(2:2) == "@") then
            idx = idx + 1
            call parse(line(3:),"-",listOfWords,numberOfWords)
            read(listOfWords(1),*) lammpsImproperTypes(1,idx)
            read(listOfWords(2),*) lammpsImproperTypes(2,idx)
            read(listOfWords(3),*) lammpsImproperTypes(3,idx)
            read(listOfWords(4),*) lammpsImproperTypes(4,idx)

          else if (line(1:14) == "improper_coeff") then

            jdx = jdx + 1
            lmp_impropers(jdx) % types(1:4) = lammpsImproperTypes(1:4,idx)
            read(line,*)dummy , lmp_impropers(jdx) % ID
                        
            call parse(line(3:)," ",listOfWords,numberOfWords)
            if (ichar(listOfWords(3)(1:1)).ge.ichar("a") .and. ichar(listOfWords(3)(1:1)).le.ichar("z")) then
              read(listOfWords(3),*) lammpsImproperForcefield(jdx)
            end if
          end if
        enddo
      
      end if

    end if
  enddo

  close(iounit)

end subroutine readLammpsForcefieldFile

subroutine writeLammpsCoordinates(iout)
  use moduleSystem 
  use moduleLammps 
  use moduleMessages
  implicit none
  integer, intent(in) :: iout
  integer :: i, itype, idx
  integer :: numberOfBondTerms
  integer :: numberOfAngleTerms
  integer :: numberOfTorsionTerms
  integer :: numberOfImproperTerms
  integer :: jtype
  character(len=100) :: str

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i=1,frame % natoms
    itype = findLammpsAtomType(frame % lab(i))
  enddo
  if (sum(lmp_atoms(:) % number) /= frame % natoms) call message(-1,"LAMMPS - not all atoms have a force field type")

  do i=1,numberOfUniqueBonds
    itype = findLammpsBondType(0,listOfUniqueBonds(1:2,i))
  enddo

  do i=1,numberOfUniqueAngles
    itype = findLammpsAngleType(0,listOfUniqueAngles(1:3,i))
  enddo

  numberOfTorsionTerms = 0
  do i=1,numberOfUniqueTorsions
    jtype = findLammpsTorsionType(0,listOfUniqueTorsions(1:4,i))
  end do
  
  do i=1,numberOfUniqueOutOfPlane
    itype = findLammpsImproperType(0,listOfUniqueOutOfPlane(1:4,i))
  enddo

  numberOfBondTerms = 0
  numberOfAngleTerms = 0
  numberOfTorsionTerms = 0
  numberOfImproperTerms = 0

  if (allocated(lmp_bonds)) numberOfBondTerms = sum(lmp_bonds(:) % number)
  if (allocated(lmp_angles)) numberOfAngleTerms = sum(lmp_angles(:) % number)
  if (allocated(lmp_torsions)) numberOfTorsionTerms = sum(lmp_torsions(:) % number)
  if (allocated(lmp_impropers)) numberOfImproperTerms = sum(lmp_impropers(:) % number)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call message(0,"LAMMPS forcefiel summary")
  call message(0,"...Number of atom types",i=lammpsNumberOfAtoms)
  do idx=1,lammpsNumberOfAtoms
    call message(0,"......Atom type -> "//trim(lammpsAtomTypes(idx)),lmp_atoms(idx) % number)
  enddo
  call message(0,"...Number of bond types",i=lammpsNumberOfbonds)
  do idx=1,lammpsNumberOfBonds
    call message(0,"......Bond type -> "//trim(lammpsBondTypes(1,idx))//"-"//trim(lammpsBondTypes(2,idx)),lmp_bonds(idx) % number)
  enddo
  call message(0,"...Number of angle types",i=lammpsNumberOfAngles)
   do idx=1,lammpsNumberOfAngles
    call message(0,"......Angle type -> "//trim(lammpsAngleTypes(1,idx))//"-"//trim(lammpsAngleTypes(2,idx))//"-"//trim(lammpsAngleTypes(3,idx)),lmp_angles(idx) % number)
  enddo
  call message(0,"...Number of dihedral types",i=lammpsNumberOfTorsions)
   do idx=1,lammpsNumberOfTorsions
    call message(0,"......Torsion type -> "//trim(lammpsTorsionTypes(1,idx))//"-"//trim(lammpsTorsionTypes(2,idx))//"-"//trim(lammpsTorsionTypes(3,idx))//"-"//trim(lammpsTorsionTypes(4,idx)),lmp_torsions(idx) % number)
  enddo
  call message(0,"...Number of improper types",i=lammpsNumberOfOutOfPlane)
   do idx=1,lammpsNumberOfOutOfPlane
    call message(0,"......Improper type -> "//trim(lammpsImproperTypes(1,idx))//"-"//trim(lammpsImproperTypes(2,idx))//"-"//trim(lammpsImproperTypes(3,idx))//"-"//trim(lammpsImproperTypes(4,idx)),lmp_impropers(idx) % number)
  enddo
  call message(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(iout,'("LAMMPS description"/)')

  write(iout,'(i8,a15)') frame % natoms        , adjustr("atoms    ")
  write(iout,'(i8,a15)') numberOfBondTerms     , adjustr("bonds    ")
  write(iout,'(i8,a15)') numberOfAngleTerms    , adjustr("angles   ")
  write(iout,'(i8,a15)') numberOfTorsionTerms  , adjustr("dihedrals")
  write(iout,'(i8,a15)') numberOfImproperTerms , adjustr("impropers")
  write(iout,'(i8,a15)')

  write(iout,'(i8,a15)') lammpsNumberOfAtoms      , adjustr("atom types    ")
  write(iout,'(i8,a15)') lammpsNumberOfBonds      , adjustr("bond types    ")
  write(iout,'(i8,a15)') lammpsNumberOfAngles     , adjustr("angle types   ")
  write(iout,'(i8,a15)') lammpsNumberOfTorsions   , adjustr("dihedral types")
  write(iout,'(i8,a15)') lammpsNumberOfOutOfPlane , adjustr("improper types")
  write(iout,'(i8,a15)')

  write(iout,'(/,g15.8,g15.8," xlo xhi")')0.0d0,frame % hmat(1,1)
  write(iout,'(  g15.8,g15.8," ylo yhi")')0.0d0,frame % hmat(2,2)
  write(iout,'(  g15.8,g15.8," zlo zhi")')0.0d0,frame % hmat(3,3)
  if (abs(frame % hmat(1,2))+abs(frame % hmat(1,3))+abs(frame % hmat(2,3))>1.d-8) &
    write(iout,'(3g15.8," xy xz yz")')frame % hmat(1,2),frame % hmat(1,3),frame % hmat(2,3)

  write(iout,'(/" Atoms"/)')
  do i=1,frame % natoms
    itype = findLammpsAtomType(frame % lab(i))
    if (itype==0) call message(-1,"Cannot find atom type for",i=findLammpsAtomType(frame % lab(i)))
    write(iout,'(3i7,4f15.7)') i, atomToMoleculeIndex(i), itype, frame % chg(i), frame % pos(:,i)
  enddo

  if (numberOfUniqueBonds>0) write(iout,'(/" Bonds"/)')
  do i=1,numberOfUniqueBonds
    itype = findLammpsBondType(iout,listOfUniqueBonds(1:2,i))
     if (itype==0) then
      str = trim(frame % lab(listOfUniqueBonds(1,i)))//" - "//trim(frame % lab(listOfUniqueBonds(2,i)))
      if (stopForMissingForceField(str,ignoreMissingBonds)) call message(-1,"Cannot find bond type for atom types",str=str)
    end if
  enddo

  if (numberOfUniqueAngles>0) write(iout,'(/" Angles"/)')
  do i=1,numberOfUniqueAngles
    itype = findLammpsAngleType(iout,listOfUniqueAngles(1:3,i))
    if (itype==0) then
      str = trim(frame % lab(listOfUniqueAngles(1,i)))
      do idx=2,3
        str = trim(str)//" - "//trim(frame % lab(listOfUniqueAngles(idx,i)))
      end do
      if (stopForMissingForceField(str,ignoreMissingAngles)) call message(-1,"Cannot find angle type for atom types",str=str)
    end if
  enddo

  if (numberOfUniqueTorsions>0) write(iout,'(/" Dihedrals"/)')
  numberOfTorsionTerms = 0
  do i=1,numberOfUniqueTorsions
    jtype = findLammpsTorsionType(iout,listOfUniqueTorsions(1:4,i))
    if (itype==0) then
      str = trim(frame % lab(listOfUniqueTorsions(1,i)))
      do idx=2,4
        str = trim(str)//" - "//trim(frame % lab(listOfUniqueTorsions(idx,i)))
      end do
      if (stopForMissingForceField(str,ignoreMissingTorsions)) call message(-1,"Cannot find torsion type for atom types",str=str)
    end if
  end do
  
  if (numberOfUniqueOutOfPlane>0) write(iout,'(/" Impropers"/)')
  do i=1,numberOfUniqueOutOfPlane
    itype = findLammpsImproperType(iout,listOfUniqueOutOfPlane(1:4,i))
    if (itype==0) then
      str = trim(frame % lab(listOfUniqueOutOfPlane(1,i)))
      do idx=2,4
        str = trim(str)//" - "//trim(frame % lab(listOfUniqueOutOfPlane(idx,i)))
      end do
      if (stopForMissingForceField(str,ignoreMissingImpropers)) call message(-1,"Cannot find improper type for atom types",str=str)
    end if
  enddo

contains
 
  integer function findLammpsAtomType(l) result(i)
    implicit none
    character(len=*), intent(in) :: l
    do i=lammpsNumberOfAtoms,1,-1
      if (l == lammpsAtomTypes(i)) then
        lmp_atoms(i) % number = lmp_atoms(i) % number + 1
        exit
      end if
    enddo
    return
  end function findLammpsAtomType

  integer function findLammpsBondType(iout,idx)
    use moduleVariables, only : cp
    implicit none
    integer, intent(in) :: iout
    integer, intent(in) :: idx(2)
    character(len=cp) :: l(2)
    integer :: i
    integer, save :: ibond = 0
    l = [ frame % lab(idx(1)) , frame % lab(idx(2)) ]
    findLammpsBondType = 0
    do i=lammpsNumberOfBonds,1,-1
      if ( all( [l(1),l(2)] == lmp_bonds(i) % types(1:2)) .or. &
           all( [l(2),l(1)] == lmp_bonds(i) % types(1:2)) ) then
            
        if (iout /= 0) then
          ibond = ibond + 1
          write(iout,'(4i8)')ibond, lmp_bonds(i) % ID, idx(1), idx(2)
        else
          lmp_bonds(i) % number = lmp_bonds(i) % number + 1
        endif

        findLammpsBondType = ibond
      end if
    enddo   
    return
  end function findLammpsBondType

  integer function findLammpsAngleType(iout,idx)
    use moduleVariables, only : cp
    implicit none
    integer, intent(in) :: iout
    integer, intent(in) :: idx(3)
    character(len=cp) :: l(3)
    integer :: i
    integer, save :: iangle = 0
    l = [ frame % lab(idx(1)) , frame % lab(idx(2)) , frame % lab(idx(3)) ]
    findLammpsAngleType = 0
    do i=lammpsNumberOfAngles,1,-1
      if ( all( [l(1),l(2),l(3)] == lmp_angles(i) % types(1:3)) .or. &
           all( [l(3),l(2),l(1)] == lmp_angles(i) % types(1:3)) ) then
        
        if (iout /= 0) then
          iangle = iangle + 1
          write(iout,'(5i8)')iangle, lmp_angles(i) % ID, idx(1), idx(2), idx(3)
        else
          lmp_angles(i) % number = lmp_angles(i) % number + 1
        end if

        findLammpsAngleType = iangle
      end if
    enddo
    return
  end function findLammpsAngleType

  integer function findLammpsTorsionType(iout,idx)
    use moduleVariables, only : cp
    implicit none
    integer, intent(in) :: iout
    integer, intent(in) :: idx(4)
    character(len=cp) :: l(4)
    integer :: i
    integer, save :: itors = 0

    l = [ frame % lab(idx(1)) , frame % lab(idx(2)) , frame % lab(idx(3)) , frame % lab(idx(4)) ]
    findLammpsTorsionType = 0
    do i=lammpsNumberOfTorsions,1,-1
      if ( all( [l(1),l(2),l(3),l(4)] == lmp_torsions(i) % types(1:4)) .or. &
           all( [l(4),l(3),l(2),l(1)] == lmp_torsions(i) % types(1:4)) ) then

        if (iout /= 0) then
          itors = itors + 1
          write(iout,'(6i8)')itors, lmp_torsions(i) % ID, idx(1), idx(2), idx(3), idx(4)
        else
          lmp_torsions(i) % number = lmp_torsions(i) % number + 1
        endif
  
        findLammpsTorsionType = itors
      end if
    enddo
    return
  end function findLammpsTorsionType

  integer function findLammpsImproperType(iout,idx)
    use moduleVariables, only : cp
    implicit none
    integer, intent(in) :: iout
    integer, intent(in) :: idx(4)
    character(len=cp) :: l(4)
    integer :: i
    integer, save :: iimp = 0
    l = [ frame % lab(idx(1)) , frame % lab(idx(2)) , frame % lab(idx(3)) , frame % lab(idx(4)) ]
    findLammpsImproperType = 0
    do i=lammpsNumberOfOutOfPlane,1,-1
      if (l(1) /= lmp_impropers(i) % types(1)) cycle
      if ( all( [l(2),l(3),l(4)] == lmp_impropers(i) % types(2:4)) .or. &
           all( [l(2),l(4),l(3)] == lmp_impropers(i) % types(2:4)) .or. &
           all( [l(3),l(2),l(4)] == lmp_impropers(i) % types(2:4)) .or. &
           all( [l(3),l(4),l(2)] == lmp_impropers(i) % types(2:4)) .or. &
           all( [l(4),l(2),l(3)] == lmp_impropers(i) % types(2:4)) .or. &
           all( [l(4),l(3),l(2)] == lmp_impropers(i) % types(2:4)) ) then

        if (iout /= 0) then
          iimp = iimp + 1
          if ( any(lammpsImproperForcefield(i) == lammpsImproperCentralFirst) ) then
            write(iout,'(6i8)')iimp, lmp_impropers(i) % ID, idx(1), idx(2), idx(3), idx(4)
          else
            write(iout,'(6i8)')iimp, lmp_impropers(i) % ID, idx(2), idx(1), idx(3), idx(4)
          end if
        else
          lmp_impropers(i) % number = lmp_impropers(i) % number + 1
        end if

        findLammpsImproperType = iimp
      end if
    end do
    return
  end function findLammpsImproperType

  logical function stopForMissingForceField(str,loverride)
    character(len=*), intent(in) :: str
    logical, intent(in) :: loverride
    character(len=1) :: ignore
    integer :: i

    do i=1,numberOfIgnoredTerms
      if (trim(missingTermsIgnored(i)) == trim(str)) then
        stopForMissingForceField = .false.
        return
      end if
    end do

    write(0,*) "LAMMPS WARNING : Missing Force Field Term For "//trim(str)

    if (loverride) then
      numberOfIgnoredTerms = numberOfIgnoredTerms + 1
      missingTermsIgnored(numberOfIgnoredTerms) = trim(str)
      stopForMissingForceField = .false.
      return
    end if

    write(0,*) "Ignore Warning And Continue? [Y/N]"
    do
      read(*,*) ignore
      select case (ignore)
        case ("Y" , "y")
          numberOfIgnoredTerms = numberOfIgnoredTerms + 1
          missingTermsIgnored(numberOfIgnoredTerms) = trim(str)
          stopForMissingForceField = .false.
          exit
        case ("N" , "n")
          stopForMissingForceField = .true.
          exit
        case default
          write(0,*) "Reply yes or no"
      end select
    end do
    
  end function stopForMissingForceField
  
end subroutine writeLammpsCoordinates

!------ reading lammps coordinates ------!
subroutine getNumberOfAtomsLammps(uinp, n, hmat)
  implicit none
  integer, intent(in)  :: uinp
  integer, intent(out) :: n
  real(8), dimension(3,3), intent(out) :: hmat
  
  integer :: ios
  character(len=20)  :: line
  real(8) :: rtmp(2)
  
  n = 0
  do
    read(uinp, '(a200)',iostat=ios)line
    if (ios /= 0) exit
    if (index(line, "atoms") >0 )then
      read(line,*) n
    end if
    
    if (index(line, "xlo xhi") >0 )then
      read(line,*,iostat=ios) rtmp
      if (ios == 0) then
        hmat(1,1) = rtmp(2) - rtmp(1)
      else
        read(line,*,iostat=ios) hmat(1,1)
      end if
    end if
    
    if (index(line, "ylo yhi") >0 )then
      read(line,*,iostat=ios) rtmp
      if (ios == 0) then
        hmat(2,2) = rtmp(2) - rtmp(1)
      else
        read(line,*,iostat=ios) hmat(2,2)
      end if
    end if
    
    if (index(line, "zlo zhi") >0 )then
      read(line,*,iostat=ios) rtmp
      if (ios == 0) then
        hmat(3,3) = rtmp(2) - rtmp(1)
      else
        read(line,*,iostat=ios) hmat(3,3)
      end if
    end if
    
    if (index(line, "xy xz yz") >0 )then
      read(line,*) hmat(1,2) , hmat(1,3) , hmat(2,3) 
    end if
    
    if (index(line, "Atoms") >0 ) exit
  end do
  
end subroutine getNumberOfAtomsLammps

subroutine readCoordinatesLammps(uinp,natoms,pos,label,charge,hmat,go)
  use moduleMessages
  use moduleLammps
  implicit none
  integer, intent(in) :: uinp
  integer, intent(in) :: natoms
  real(8), dimension(3,natoms), intent(inout) :: pos
  character(*), dimension(natoms), intent(inout) :: label
  real(8), dimension(natoms), intent(inout) :: charge
  real(8), dimension(3,3), intent(inout) :: hmat
  logical, intent(inout) :: go

  character(len=10) :: fmtline
  integer, parameter :: maxCharacters = 200
  character(len=maxCharacters)  :: line

  integer :: ios
  real(8) :: rtmp(2)
  integer :: idx, iatm, itmp

  write(fmtline,'("(a",i0,")")') maxCharacters

  go = .true.

  do
    read(uinp, fmtline,iostat=ios)line
    if (ios /= 0) exit

    if (index(line, "xlo xhi") >0 )then
      read(line,*,iostat=ios) rtmp
      if (ios == 0) then
        hmat(1,1) = rtmp(2) - rtmp(1)
      else
        read(line,*,iostat=ios) hmat(1,1)
      end if
    end if
    
    if (index(line, "ylo yhi") >0 )then
      read(line,*,iostat=ios) rtmp
      if (ios == 0) then
        hmat(2,2) = rtmp(2) - rtmp(1)
      else
        read(line,*,iostat=ios) hmat(2,2)
      end if
    end if
    
    if (index(line, "zlo zhi") >0 )then
      read(line,*,iostat=ios) rtmp
      if (ios == 0) then
        hmat(3,3) = rtmp(2) - rtmp(1)
      else
        read(line,*,iostat=ios) hmat(3,3)
      end if
    end if
    
    if (index(line, "xy xz yz") >0 )then
      read(line,*) hmat(1,2) , hmat(1,3) , hmat(2,3) 
    end if

    if (index(line, "Atoms") > 0) then
      read(uinp, fmtline,iostat=ios)
      if (ios /= 0) then
        go = .false.
        exit
      end if

      ! read all atoms
      do idx=1,natoms
        read(uinp, fmtline,iostat=ios)line
        if (ios /= 0) then
          go = .false.
          exit
        end if
        read(line,*) iatm, itmp, label(iatm), charge(iatm), pos(1:3,iatm)
      end do
      write(0,*) pos(1:3,1)
      return
    end if

  end do

  go = .false.

end subroutine readCoordinatesLammps

! ---------------- LAMMPS trajectory files ----------------
subroutine getNumberOfAtomsLammpsTrajectory(uinp, n, hmat)
  implicit none
  integer, intent(in)  :: uinp
  integer, intent(out) :: n
  real(8), dimension(3,3), intent(out) :: hmat
  
  integer :: ios
  character(len=10) :: fmtline
  integer, parameter :: maxCharacters = 200
  character(len=maxCharacters)  :: line
  real(8) :: rtmp(2)
  
  write(fmtline,'("(a",i0,")")') maxCharacters
  do
    read(uinp,fmtline,iostat=ios)line
    if (ios /= 0) return

    if (index(line,"ITEM: NUMBER OF ATOMS") > 0) then
      read(uinp,*)n
      return
    end if
  end do
  
end subroutine getNumberOfAtomsLammpsTrajectory

subroutine readCoordinatesLammpsTrajectory(uinp,natoms,pos,label,charge,hmat,go)
  use moduleMessages
  use moduleStrings
  use moduleVariables, only : MAXWORDS

  implicit none
  integer, intent(in) :: uinp
  integer, intent(in) :: natoms
  real(8), dimension(3,natoms), intent(inout) :: pos
  character(*), dimension(natoms), intent(inout) :: label
  real(8), dimension(natoms), intent(inout) :: charge
  real(8), dimension(3,3), intent(inout) :: hmat
  logical, intent(inout) :: go

  integer :: ios
  character(len=10) :: fmtline
  integer, parameter :: maxCharacters = 200
  character(len=maxCharacters)  :: line

  integer :: pbcType, n
  real(8) :: rtmp(2)

  write(fmtline,'("(a",i0,")")') maxCharacters
  go = .true.

  main : do
    read(uinp,fmtline, iostat=ios) line 
    if (ios /= 0) then
      go = .false.
      exit main
    end if

    if (index(line, "ITEM:") > 0) then
      if (index(line,"ITEM: NUMBER OF ATOMS") > 0) then
        read(uinp,*)n
        if ( n /= natoms ) call message(0,"Wrong number of atoms in LAMMPS trajectory file")
      end if

      ! read the cell
      if (index(line,"ITEM: BOX BOUNDS") > 0) then

        if ( index(line,"ITEM: BOX BOUNDS pp pp pp") > 0 ) then
          pbcType = 1
        else if ( index(line,"ITEM: BOX BOUNDS xy xz yz pp pp pp") > 0 ) then
          pbcType = 2
        else
          pbcType = 0
        end if

        ! read an orthorhombic cell
        if (pbcType == 1) then
          block
            integer :: idx
            do idx=1,3
              read(uinp,fmtline, iostat=ios) line 
              if (ios /= 0) then
                go = .false.
                exit main
              end if
              read(line,*) rtmp
              hmat (idx, idx) = rtmp(2) - rtmp(1)
            end do
          end block

        ! read a triclinic cell
        else if (pbcType == 2) then
          block
            integer :: idx
            real(8) :: xlo_bound, xhi_bound, xy
            real(8) :: ylo_bound, yhi_bound, xz
            real(8) :: zlo_bound, zhi_bound, yz
            real(8) :: xlo, xhi
            real(8) :: ylo, yhi
            real(8) :: zlo, zhi
            
            read(uinp,*) xlo_bound, xhi_bound, xy
            read(uinp,*) ylo_bound, yhi_bound, xz
            read(uinp,*) zlo_bound, zhi_bound, yz

            xlo = xlo_bound - MIN(0.0,xy,xz,xy+xz)
            xhi = xhi_bound - MAX(0.0,xy,xz,xy+xz)
            ylo = ylo_bound - MIN(0.0,yz)
            yhi = yhi_bound - MAX(0.0,yz)
            zlo = zlo_bound
            zhi = zhi_bound

            hmat(1,1) = xhi - xlo
            hmat(2,2) = yhi - ylo
            hmat(3,3) = zhi - zlo
            hmat(1,2) = xy
            hmat(1,3) = xz
            hmat(2,3) = yz

          end block
        end if
      end if

      ! read the atoms' coordinated
      if (index(line,"ITEM: ATOMS") > 0) then
        block
          integer :: numberOfFields
          character(len=STRLEN), dimension(100) :: fields
          integer :: idx, iatm, ilab, ix, iy, iz, i
          
          call parse(line," ",fields,numberOfFields)
          do idx=1,numberOfFields
            if (index(fields(idx), "id")   > 0) iatm = idx - 2
            if (index(fields(idx), "type") > 0) ilab = idx - 2
            if (index(fields(idx), "x")    > 0) ix   = idx - 2
            if (index(fields(idx), "y")    > 0) iy   = idx - 2
            if (index(fields(idx), "z")    > 0) iz   = idx - 2
          end do

          do idx=1,natoms
            read(uinp,fmtline, iostat=ios) line 
            if (ios /= 0) then
              go = .false.
              exit main
            end if
            call parse(line," ",fields,numberOfFields)
            read(fields(iatm),*) i
            read(fields(ilab),*) label(i)
            read(fields(ix),*) pos(1,i)
            read(fields(iy),*) pos(2,i)
            read(fields(iz),*) pos(3,i)
          end do
          
        end block
        return
      end if

    end if

  end do main

  go = .false.

end subroutine readCoordinatesLammpsTrajectory

subroutine writeCoordinatesLammpsTrajectory(io)
  use moduleVariables, only : cp
  use moduleSystem, only : frame
  use moduleLammps
  use moduleMessages
  use moduleStrings
  implicit none
  integer, intent(in) :: io

  integer :: nTypes
  character(len=cp), dimension(100), save :: uniqueTypes
  integer, allocatable, dimension(:), save :: atomTypes
  logical, save :: firstTimeIn = .true.
! ITEM: TIMESTEP
! 100000
! ITEM: NUMBER OF ATOMS
! 3360
! ITEM: BOX BOUNDS pp pp pp
! 0.0000000000000000e+00 3.4527000000000001e+01
! 0.0000000000000000e+00 3.4232999999999997e+01
! 0.0000000000000000e+00 3.4454999999999998e+01
! ITEM: ATOMS id type x y z
! 1 5 0.238477 0.129618 34.392

  if (firstTimeIn) then
    firstTimeIn = .false.

    allocate(atomTypes(frame % natoms), source=0)

    ! if (.false.) then
    if (lmpMap) then
      block
        integer :: idx, iatm
        integer :: nlab, mlab
        character(cp), dimension(50) :: li, lo

        call parse(mapLabelsTypes(1),",",li,nlab)
        call parse(mapLabelsTypes(2),",",lo,mlab)
        if (nlab /= mlab) call message(0,'Inconsistent number of labels for the lmptrj map',iv=[nlab,mlab])

        do iatm=1,frame % natoms
          do idx=1,nlab
            if (frame % lab(iatm) == li(idx)) then
              read(lo(idx),*) atomTypes(iatm) 
              exit
            end if
          end do
        end do
      end block

    else
      nTypes = 0
      block
        integer :: i, j
        do i=1,frame % natoms
          do j=1,nTypes
            if (uniqueTypes(j) == frame % lab(i)) exit
          end do
          if (j > nTypes) then
            nTypes = nTypes + 1
            uniqueTypes(nTypes) = frame % lab(i)
            j = nTypes
          end if
          atomTypes(i) = j
        end do
      end block
    end if

  end if

  write(io,'("ITEM: TIMESTEP")')
  write(io,'(i0)') frame % nframe
  
  write(io,'("ITEM: NUMBER OF ATOMS")')
  write(io,'(i0)') frame % natoms
  
  if (abs(frame % hmat(2,1)) + abs(frame % hmat(3,1)) + abs(frame % hmat(3,2)) .gt. 1.0d-6) then
    block
      real(8) :: xlo_bound, xhi_bound, xy
      real(8) :: ylo_bound, yhi_bound, xz
      real(8) :: zlo_bound, zhi_bound, yz
      real(8) :: xlo, xhi
      real(8) :: ylo, yhi
      real(8) :: zlo, zhi

      xlo = 0.d0
      ylo = 0.d0
      zlo = 0.d0
      xhi = frame % hmat(1,1)
      yhi = frame % hmat(2,2)
      zhi = frame % hmat(3,3)
      xy = frame % hmat(1,2)
      xz = frame % hmat(2,3)
      yz = frame % hmat(3,3)

      xlo_bound = xlo + MIN(0.d0,xy,xz,xy+xz)
      xhi_bound = xhi + MAX(0.d0,xy,xz,xy+xz)
      ylo_bound = ylo + MIN(0.d0,yz)
      yhi_bound = yhi + MAX(0.d0,yz)
      zlo_bound = zlo
      zhi_bound = zhi

      write(io,'(3(e20.15,1x))') xlo_bound, xhi_bound, xy
      write(io,'(3(e20.15,1x))') ylo_bound, yhi_bound, xz
      write(io,'(3(e20.15,1x))') zlo_bound, zhi_bound, yz
    end block      

  else
    write(io,'("ITEM: BOX BOUNDS pp pp pp")')
    write(io,'(2(e20.15,1x))') 0.d0, frame % hmat(1,1)
    write(io,'(2(e20.15,1x))') 0.d0, frame % hmat(2,2)
    write(io,'(2(e20.15,1x))') 0.d0, frame % hmat(3,3)
  end if


  write(io,'("ITEM: ATOMS id type x y z")')
  block
    integer :: i
    character(len=200) :: str
    do i=1,frame % natoms
      write(str,'(2(i0,1x),3(f13.6,1x))')i, atomTypes(i), frame % pos(:,i)
      call compact(str)
      write(io,'(a)')trim(str)
    end do
  end block

end subroutine writeCoordinatesLammpsTrajectory
