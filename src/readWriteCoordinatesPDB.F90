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
subroutine getNumberOfAtomsPDB(uinp,natoms,hmat)
  implicit none
  integer, intent(in)  :: uinp
  integer, intent(out) :: natoms
  real(8), dimension(3,3), intent(out) :: hmat
  real(8), parameter :: pi = 3.1415926535898d0
  character(len=80)  :: line
  character(len=10)  :: str
  real(8), dimension(6) :: cell

  natoms=0
  do
    read(uinp,'(a80)',end=546)line
    read(line(1:6),*)str
    
    if(str == "CRYST1") then
      read(line,*)str,cell
      cell(4) = cell(4) * pi / 180.0d0
      cell(5) = cell(5) * pi / 180.0d0
      cell(6) = cell(6) * pi / 180.0d0
      call cell2hmat(cell,hmat)
      
    else if(str=="ATOM" .or. str=="HETATM")then
      natoms=natoms+1

    else if(str=="END" .or. str=="ENDMDL")then
      exit
    else
      cycle
    end if
  enddo
546 continue

end subroutine getNumberOfAtomsPDB

subroutine readCoordinatesPDB(uinp,natoms,pos,label,charge,hmat,element,go)
  use moduleMessages 
  implicit none
  integer, intent(in) :: uinp
  integer, intent(in) :: natoms
  real(8), dimension(3,natoms), intent(inout) :: pos
  character(*), dimension(natoms), intent(inout) :: label
  real(8), dimension(natoms), intent(inout) :: charge
  real(8), dimension(3,3), intent(inout) :: hmat
  character(*), dimension(natoms), intent(inout) :: element
  logical, intent(inout) :: go

  real(8), parameter :: pi = 3.1415926535898d0

  integer :: iatom
  real(8) :: cell(6)
  character(len=100) :: line
  character(len=10)  :: str
  integer :: ios

  iatom = 0
  go=.true.
  do
    read(uinp,'(a100)',iostat=ios)line
    if (ios/=0) then
      if (is_iostat_end(ios) .and. iatom==natoms) return 
      go=.false.
      return
    end if

    ! if (len_trim(line)==0) cycle

    read(line(1:6),*)str
    if(str=="ATOM" .or. str=="HETATM")then
      iatom=iatom+1
      if (iatom>natoms) then
        call message(-1,"###### WARNING ###### too many atoms in the PDB")
        go=.false.
        return
      end if
! Labels and positions
      read(line(13:16),*)label(iatom) 
      read(line(31:54),'(3f8.3)')pos(:,iatom)
      read(line(77:78),*)element(iatom)

! charges - beyond character 80
      charge(iatom)=0.0d0
      if (len_trim(line)>80) then
        read(line(81:100),*,iostat=ios) charge(iatom) 
      end if
 
    else if(str == "CRYST1") then
      read(line,*)str,cell
      cell(4) = cell(4) * pi / 180.0d0
      cell(5) = cell(5) * pi / 180.0d0
      cell(6) = cell(6) * pi / 180.0d0
      call cell2hmat(cell,hmat)

    else if(str=="END" .or. str=="ENDMDL")then
      exit
    else
      cycle
    end if
  enddo

end subroutine readCoordinatesPDB

subroutine writeCoordinatesPDB(uout,lpdb2)
  use moduleSystem 
  use moduleElements 
  implicit none
  integer, intent(in) :: uout
  logical, intent(in) :: lpdb2

  character(len=100) :: line
  integer :: iatm
  integer :: id, jd, itmp
  character(len=4) :: fmtAtom
  character(len=5) :: chr

  write(uout,'(a18,i6)')"REMARK --- frame:  ",frame % nframe
  if(frame % volume>1.0d0)write(uout,'(a6,3f9.3,3f7.2)')"CRYST1",frame % cell

  if (frame % natoms<=99999) then
    fmtAtom = '(i5)'
  else
    fmtAtom = '(z5)'
  end if

  do iatm=1,frame % natoms
    id=getAtomicNumber(frame % lab(iatm))
    if (numberOfMolecules>0) then
      jd=atomToMoleculeIndex(iatm)
    else
      jd=iatm
    end if
    line=" "
 
    write(line( 1: 6),'(a6   )')"ATOM  "                           ! a literal "ATOM  " (note two trailing spaces).
    write(line( 7:11),fmtAtom) iatm                                ! atom serial number, e.g. "   86". 
    write(line(12:12),'(a1   )')" "                                ! space
    write(line(13:16),'(a4   )')adjustl(frame % lab(iatm))         ! Atom role name, e.g. " CG1;". 
    write(line(17:17),'(a1   )')" "                                ! atom variant, officially called the "alternate location indicator". This is usually " " 
 
    if (numberOfMolecules>0) then
      if (jd<1 .or. jd>numberOfMolecules) then
        write(line(18:22),'(a4,1x)')"UNK"                            ! amino acid abbreviation, e.g. "ARG". 
      else
        write(line(18:22),'(a4,1x)')listOfMolecules(jd) % resname    ! amino acid abbreviation, e.g. "ARG". 
      end if
    else
      write(line(18:22),'(a4,1x)')"UNK"                            ! amino acid abbreviation, e.g. "ARG". 
    end if
 
    write(chr,fmtAtom)jd
    write(line(23:27),'(a5   )')adjustl(chr)                       ! residue sequence number (I4) and insertion code (A1), e.g. "  11 " or " 256C". 
 
    write(line(28:30),'(a3   )')"   "                              ! NON STANDARD one spaces
    write(line(31:54),'(3f8.3)')frame % pos(:,iatm)                ! atom coordinates (X, Y, Z)
    write(line(55:60),'(f6.2 )')1.d0                               ! atom occupancy, usually "  1.00".
    write(line(61:66),'(f6.2 )')0.d0                               ! B value or temperature factor.
    write(line(67:72),'(a6   )')"      "                           ! Atom's species
    write(line(73:76),'(a4   )')"    "                             ! segment identifier, left-justified. [format version 2.0 and later.]
    write(line(77:78),'(a2   )')adjustr(trim(atom(id) % lab))      ! element symbol, right-justified. [format version 2.0 and later.]
    write(line(79:80),'(a2   )')"  "                               ! charge on the atom. [format version 2.0 and later.]
    if (.not.lpdb2) write(line(81:100),'(f12.8)')frame % chg(iatm) ! charges 
    write(uout,'(a100)')line
  enddo

  if (.not.lpdb2) then
    write(uout,'("END")')
  else
    write(uout,'("ENDMDL")')
    do iatm=1,frame % natoms
      itmp = numberOfCovalentBondsPerAtom(iatm)
      if (itmp>0) write(uout,'("CONECT",5i5)') iatm , listOfCovalentBondsPerAtom(1:itmp,iatm)
    enddo
  end if

end subroutine writeCoordinatesPDB
