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
  subroutine getNumberOfAtomsGaussian(iounit,natoms)
    integer, intent(in) :: iounit
    integer, intent(out) :: natoms
    integer :: ierr
    character(len=100) :: line

    main : do
      read(iounit,'(a100)',iostat=ierr)line
      if (ierr /= 0) exit main

      if (len_trim(line) == 0) cycle

      if (index(line,"Standard orientation") > 0) then
        do i=1,4
          read(iounit,'(a100)',iostat=ierr)line
          if (ierr /= 0) exit main
        end do
      
        natoms = 0
        do
          read(iounit,'(a100)',iostat=ierr)line
          if (ierr /= 0) exit main
          if (index(line,"--") > 0) exit main
          natoms = natoms + 1
        end do
      
      end if

    end do main

  end subroutine getNumberOfAtomsGaussian

  subroutine readCoordinatesGaussian(iounit,n,pos,lab,go)
    use moduleElements
    implicit none
    integer, intent(in) :: iounit
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(out) :: pos
    character(len=4), dimension(n), intent(out) :: lab
    logical, intent(out) :: go

    integer :: iatm, i, ierr, nw
    character(len=100) :: line

!                              Standard orientation:
!  ---------------------------------------------------------------------
!  Center     Atomic      Atomic             Coordinates (Angstroms)
!  Number     Number       Type             X           Y           Z
!  ---------------------------------------------------------------------
!       1          6           0        1.415358    0.307944    1.371401
!       ...
!      60          1           0       -5.696168   -0.945555   -1.360240
! ---------------------------------------------------------------------
!  Rotational constants (GHZ):           0.1348428           0.1085367           0.0968257

    go = .true.
    main : do
      read(iounit,'(a100)',iostat=ierr)line
      if (ierr /= 0) exit main

      if (len_trim(line) == 0) cycle

      if (index(line,"Standard orientation") > 0) then
        do i=1,4
          read(iounit,'(a100)',iostat=ierr)line
          if (ierr /= 0) exit main
        end do
      
        do
          read(iounit,'(a100)',iostat=ierr)line
          if (ierr /= 0) exit main
          if (index(line,"--") > 0) then
            exit main
          else
            read(line,*)iatm, i, nw, pos(1:3,iatm)
            lab(iatm) = atom(i) % lab
          end if

        end do
      
      end if

    end do main

    if (ierr /= 0) go = .false.

  end subroutine readCoordinatesGaussian

  subroutine getNumberOfAtomsXYZ(iounit,natoms,hmat)
    implicit none
    integer, intent(in) :: iounit
    integer, intent(out) :: natoms
    real(8), dimension(3,3), intent(inout) :: hmat
    character(len=500) :: line
    character(len=100) :: cellString
    integer :: ierr, i,j
    
    natoms = 0

    read(iounit,*,iostat=ierr)natoms
    if (ierr/=0) return

    read(iounit,'(a500)',iostat=ierr)line
    if (ierr/=0 .or. len_trim(line)==0) return
    
    i = index(line,"Lattice=")
    if (i > 0) then
      read(line(i+8:),*)cellString
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1:3,1), hmat(1:3,2), hmat(1:3,3)
    end if

    i = index(line,"Cell=")
    if (i > 0) then
      read(line(i+5:),*)cellString
      hmat = 0.d0
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1,1), hmat(2,2), hmat(3,3)
    end if

  end subroutine getNumberOfAtomsXYZ

  subroutine readCoordinatesXYZ(iounit,n,pos,lab,chg,hmat,go)
    use moduleStrings
    implicit none
    integer, intent(in) :: iounit
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(out) :: pos
    real(8), dimension(n), intent(out) :: chg
    character(len=4), dimension(n), intent(out) :: lab
    real(8), dimension(3,3), intent(inout) :: hmat
    logical, intent(out) :: go

    integer :: i, j, nn, ierr
    character(len=500) :: line

    character(len=100) :: cellString
    character(len=100) :: fieldString
    integer :: nfields
    character(STRLEN), dimension(100) :: fields

    integer :: nf, ilab, ipos, ichg

    logical :: oldXYZ

    go = .true. 

    read(iounit,*,iostat=ierr)nn
    if (ierr/=0 .or. nn/=n) then 
      go = .false.
      return
    end if

    read(iounit,'(a500)',iostat=ierr)line
    if (ierr/=0) then 
      go = .false.
      return
    end if
    
    i = index(line,"Lattice=")
    if (i > 0) then
      read(line(i+8:),*)cellString
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1:3,1), hmat(1:3,2), hmat(1:3,3)
    end if

    i = index(line,"Cell=")
    if (i > 0) then
      read(line(i+5:),*)cellString
      hmat = 0.d0
      j = len_trim(cellString)-1
      read(cellString,*)hmat(1,1), hmat(2,2), hmat(3,3)
    end if

    oldXYZ = .true.
    i = index(line,"Properties=")
    if (i > 0) then
      oldXYZ = .false.
      read(line(i+10:),*)fieldString
      
      j = len_trim(fieldString)
      call parse(fieldString(2:j), ":", fields, nfields)
      nf = 0
      do i=1,nfields,3
        select case (fields(i))
        case ("species")
          ilab = nf + 1
        case ("pos")
          ipos = nf + 1
        case ("final_charges" , "charges")
          ichg = nf + 1
        end select
        read(fields(i+2),*) j 
        nf = nf + j
      end do
    end if

    if (oldXYZ) then

      do i=1,n
        read(iounit,'(a200)',iostat=ierr)line
        if (ierr/=0) then
          go = .false.
          return
        end if
        read(line,*,iostat=ierr)lab(i), pos(1:3,i)
        if (ierr/=0) then
          go = .false.
          return
        end if
      enddo


    else
      do i=1,nn
        read(iounit,'(a200)',iostat=ierr)line
        if (ierr/=0) then 
          go = .false.
          return
        end if
        read(line,*,iostat=ierr)fields(1:nf)
        if (ierr/=0) then
          go = .false.
          return
        end if

        read(fields(ilab),       *)lab(i)
        read(fields(ipos:ipos+2),*)pos(1:3,i)
        read(fields(ichg),       *)chg(i)
      end do
    end if

  end subroutine readCoordinatesXYZ

  subroutine getNumberOfAtomsARC(iounit,natoms,hmat)
    use moduleVariables, only : identityMatrix
    implicit none
    integer, intent(in) :: iounit
    integer, intent(out) :: natoms
    real(8), dimension(3,3), intent(inout) :: hmat
    character(len=500) :: line
    integer :: ierr
    
    natoms = 0

    read(iounit,'(a500)',iostat=ierr)line
    if (ierr/=0) return

    read(line,*) natoms

    read(line,*,iostat=ierr)natoms,hmat(1:3,1), hmat(1:3,2), hmat(1:3,3)
    if (ierr/=0)then
      hmat = identityMatrix
    endif

  end subroutine getNumberOfAtomsARC

  subroutine readCoordinatesARC(iounit,n,pos,lab,chg,hmat,go)
    use moduleStrings
    use moduleVariables, only : identityMatrix
    implicit none
    integer, intent(in) :: iounit
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(out) :: pos
    real(8), dimension(n), intent(out) :: chg
    character(len=4), dimension(n), intent(out) :: lab
    real(8), dimension(3,3), intent(inout) :: hmat
    logical, intent(out) :: go

    integer :: i, j, nn, ierr
    character(len=500) :: line

    go = .true. 

    read(iounit,'(a500)',iostat=ierr)line
    if (ierr/=0 .or. len_trim(line)==0) then 
      go = .false.
      return
    end if
    read(line,*,iostat=ierr) nn
    if (ierr/=0 .or. nn/=n) then 
      go = .false.
      return
    end if

    read(line,*,iostat=ierr)nn,hmat(1:3,1), hmat(1:3,2), hmat(1:3,3)
    if (ierr/=0)then
      hmat = identityMatrix
    endif

    do i=1,n
      read(iounit,'(a200)',iostat=ierr)line
      if (ierr/=0) then
        go = .false.
        return
      end if
      read(line,*,iostat=ierr)j,lab(j), pos(1:3,j)
      if (ierr/=0) then
        go = .false.
        return
      end if
    enddo
    
    chg = 0.d0
    
  end subroutine readCoordinatesARC

