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
subroutine dumpXYZ(iounit,nat,pos,lab,hmat)
  use moduleVariables, only : cp
  implicit none
  integer, intent(in) :: iounit
  integer, intent(in) :: nat
  real(8), dimension(3,nat), intent(in) :: pos
  character(cp), dimension(nat), optional, intent(in) :: lab
  real(8), dimension(3,3), optional, intent(in) :: hmat
  integer :: i
  
  character(len=200) :: line

  write(iounit,'(i0)')nat

  if (present(hmat)) then
    write(line,'(9(f17.8,1x))') hmat
  else
    line = ""
  end if
  write(iounit,'(a)') trim(line)

  if (present(lab)) then
    do i=1,nat
      write(iounit,'(a4,1x,4(f12.5,1x))')lab(i), pos(1:3,i)
    enddo
  else
    do i=1,nat
      write(iounit,'(a4,1x,4(f12.5,1x))')"Ar", pos(1:3,i)
    enddo
  end if
  call flush(iounit)

end subroutine

subroutine dumpPDB(uout,natoms,pos,label,hmat)
  use moduleVariables, only : cp
  use moduleElements, only : getAtomicNumber, atom

  implicit none
  integer, intent(in) :: uout
  integer, intent(in) :: natoms
  real(8), dimension(3,natoms), intent(in) :: pos
  character(cp), dimension(natoms), intent(in) :: label
  real(8), intent(in) :: hmat(3,3)

  integer :: i,nn
  character(len=80)  :: line
  real(8) :: cell(6)
  real(8)  :: occ, beta, tmp
  character(len=5) :: chr
  real(8), allocatable, dimension(:,:), target :: pos_tmp
  real(8) :: hnew(3,3), hmati(3,3), volume
  integer :: id

  allocate(pos_tmp(3,natoms), source=pos)

  call getInverseCellMatrix (hmat,hmati,volume)
  if (volume>1e-6) then
    tmp=abs(hmat(2,1)) + abs(hmat(3,1)) + abs(hmat(3,2))
    if (tmp.gt.1.d-9) then
      hnew = hmat
      call makeUpperTriangularCell(hnew,pos_tmp,natoms)
      call hmat2cell(hnew,cell,"DEG")
    else
      call hmat2cell(hmat,cell,"DEG")
    endif
    if(volume>1.0d0)write(uout,'(a6,3f9.3,3f7.2)')"CRYST1",cell
  endif

  occ=1.0d0
  beta=0.0d0

  do i=1,natoms
    id=getAtomicNumber(label(i))
    line=" "

    write(line( 1: 6),'(a6   )')"ATOM  "                    ! a literal "ATOM  " (note two trailing spaces).
    if (natoms<=99999) then
      write(line( 7:11),'(i5   )')i                         ! atom serial number, e.g. "   86".
    else
      write(line( 7:11),'(z5   )')i                         ! NON-STANDARD atom serial number, e.g. "   86".
    endif
    write(line(12:12),'(a1   )')" "                         ! space
    write(line(13:16),'(a4   )')adjustl(label(i))           ! Atom role name, e.g. " CG1;".
    write(line(17:17),'(a1   )')" "                         ! atom variant, officially called the "alternate location indicator". This is usually " "

    write(line(18:22),'(a4,1x)')"UNK"                       ! amino acid abbreviation, e.g. "ARG".

    if (nn<=99999) then
      write(chr,'(i5)')i
    else
      write(chr,'(z5)')i
    endif
    write(line(23:27),'(a5   )')adjustl(chr)                ! residue sequence number (I4) and insertion code (A1), e.g. "  11 " or " 256C".

    write(line(28:30),'(a3   )')"   "                       ! NON STANDARD one spaces
    write(line(31:54),'(3f8.3)')pos_tmp(:,i)                ! atom coordinates (X, Y, Z)
    write(line(55:60),'(f6.2 )')occ                         ! atom occupancy, usually "  1.00".
    write(line(61:66),'(f6.2 )')beta                        ! B value or temperature factor.
    write(line(67:72),'(a6   )')"      "                    ! Atom's species
    write(line(73:76),'(a4   )')"    "                      ! segment identifier, left-justified. [format version 2.0 and later.]
    write(line(77:78),'(a2   )')adjustr(trim(atom(id)%lab)) ! element symbol, right-justified. [format version 2.0 and later.]
    write(line(79:80),'(a2   )')"  "                        ! charge on the atom. [format version 2.0 and later.]
    write(uout,'(a80)')line
  enddo
  write(uout,'("END")')

  deallocate(pos_tmp)

end subroutine dumpPDB

subroutine dumpGULP(uout, natoms, pos, lab, hmat)
  use moduleVariables, only : cp
  implicit none
  integer, intent(in) :: uout
  integer, intent(in) :: natoms
  real(8), dimension(3,natoms), intent(in) :: pos
  character(cp), dimension(natoms) :: lab
  real(8), dimension(3,3), intent(in) :: hmat

  integer :: iatom

  write(uout,'("vectors")')
  write(uout,'(3f22.8)') hmat(:,1)
  write(uout,'(3f22.8)') hmat(:,2)
  write(uout,'(3f22.8)') hmat(:,3)

  write(uout,'("cartesian")')

  do iatom=1,natoms
    write(uout,'(a8,3f22.8)') adjustl(trim(lab(iatom))), pos(1:3,iatom)
  enddo

end subroutine dumpGULP

subroutine dumpCoordinates(ftype, funit, natoms, pos, label, hmat)
  use moduleVariables
  use moduleMessages
  implicit none 
  character(len=*) :: ftype
  integer, intent(in) :: funit
  integer, intent(in) :: natoms
  real(8), dimension(3,natoms), intent(in) :: pos
  character(cp), dimension(natoms), intent(in) :: label
  real(8), dimension(3,3), optional, intent(in) :: hmat

  select case (ftype)
    case default
      call message(-1,"Unknown coordinates type",str=ftype)
    case ("xyz" , "XYZ")

      if (present(hmat)) then
        call dumpXYZ(funit, natoms, pos, label, hmat)
      else
        call dumpXYZ(funit, natoms, pos, label)
      end if

    case ("pdb" , "PDB")
      call dumpPDB(funit, natoms, pos, label, hmat)
      
    case ("gin")
      call dumpGULP(funit, natoms, pos, label, hmat)

  end select

end subroutine dumpCoordinates
