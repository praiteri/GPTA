!
! This file is part of libdcdfort
! https://github.com/wesbarnett/dcdfort
!
! Copyright (c) 2017,2018 by James W. Barnett
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

!> @file
!> @author James W. Barnett, Columbia University
!
!> @brief Module that contains dcdreader class

module dcdfort_reader

    use dcdfort_common
    use iso_c_binding, only: C_NULL_CHAR

    implicit none
    private

    !> @brief dcdwriter class
    type, public :: dcdfile
        integer(kind=int32), private :: u
        integer(kind=int64) :: filesize, framesize
        character(len=80), allocatable :: titles(:)
    contains
        !> Opens file to read from
        procedure :: open => dcdfile_open
        !> Reads header of open DCD file
        procedure :: read_header => dcdfile_read_header
        !> Closes DCD file
        procedure :: close => dcdfile_close
        !> Reads next frame into memory
        procedure :: read_next => dcdfile_read_next
        !> Skips reading this frame into memory
        procedure :: skip_next => dcdfile_skip_next
    end type dcdfile

contains

    !> @brief Opens file to read from
    !> @details Also detects if the DCD file is compatible with this library and swaps endinness if necessary.
    !> @param[inout] this dcdreader object
    !> @param[in] filename name of of DCD file to read from
    subroutine dcdfile_open(this, filename)

        implicit none
        character(len=*), parameter :: magic_string = "CORD"
        integer(int32), parameter :: magic_number = 84
        character(len=*), intent(in) :: filename
        class(dcdfile), intent(inout) :: this
        integer(kind=int32) :: line1, charmm_version, has_extra_block, four_dimensions
        character(len=4) :: line2
        logical :: ex

        ! Does file exist?
        inquire(file=trim(filename), exist=ex, size=this%filesize)
        if (ex .eqv. .false.) then
            call error_stop_program("The specified DCD file '"//trim(filename)//"' does not exist.")
        end if

        ! Open file in native endinness
        open(newunit=this%u, file=trim(filename), form="unformatted", access="stream", status="old")

        ! Read in magic number of magic string
        read(this%u,pos=1) line1
        read(this%u) line2

        ! Ensure the magic number and string are correct, if not we'll swap the endinness
        if (line1 .ne. magic_number .or. line2 .ne. magic_string) then

            ! Try converting to the reverse endianness
            close(this%u)
            open(newunit=this%u, file=trim(filename), form="unformatted", access="stream", status="old", convert="swap")

            read(this%u,pos=1) line2
            read(this%u) line2

            ! We tried both native and reverse endiness and didn't have magic number or string
            if (line1 .ne. magic_number .or. line2 .ne. magic_string) then
                call error_stop_program("This DCD file format is not supported, or the file header is corrupt.")
            end if

        end if

        ! Check if the file identifies as CHARMM (LAMMPS pretends to be CHARMM v. 24)
        read(this%u, pos=85) charmm_version
        if (charmm_version .eq. 0) then
            call error_stop_program("DCD file indicates it is not CHARMM. Only CHARMM-style DCD files are supported.")
        end if

        ! We only support files with the extra unitcell block
        read(this%u, pos=49) has_extra_block
        if (has_extra_block .ne. 1) then
            call error_stop_program("DCD file indicates it does not have unit cell information. Only DCD files with&
                & unit cell information are supported.")
        end if

        ! We don't support files with four dimensions
        read(this%u) four_dimensions
        if (four_dimensions .eq. 1) then
            call error_stop_program("DCD file indicates it has four dimensions. Only DCD files with three dimensions&
                & are supported.")
        end if

    end subroutine dcdfile_open

    !> @brief Reads header of open DCD file
    !> @details Should be called after open()
    !> @param[inout] this dcdreader object
    !> @param[in] nframes rnumber of frames (snapshots) in file
    !> @param[out] istart first timestep of trajectory file
    !> @param[out] nevery how often (in timesteps) file was written to
    !> @param[out] iend last timestep of trajectory file
    !> @param[out] timestep timestep of simulation
    !> @param[out] natoms number of atoms in each snapshot
    subroutine dcdfile_read_header(this, nframes, istart, nevery, iend, timestep, natoms)

        implicit none

        integer(kind=int32), intent(out) :: nframes, istart, nevery, iend, natoms
        integer(kind=int32) :: i, ntitle, n, dummy, pos
        integer(kind=int64) :: nframes2
        character(len=80) :: title_string
        real(kind=real32), intent(out) :: timestep
        class(dcdfile), intent(inout) :: this

        read(this%u, pos=9) nframes, istart, nevery, iend

        read(this%u, pos=45) timestep

        read(this%u, pos=97) ntitle
        if (ntitle > 0) then
            write(error_unit,'(a)') prompt//"The following titles were found:"
        end if
        allocate(this%titles(ntitle))
        do i = 1, ntitle
            read(this%u) title_string

            n = 1
            do while (n .le. 80) ! .and. title_string(n:n) .ne. C_NULL_CHAR)
                n = n + 1
            end do

            this%titles(i) = trim(title_string(1:n-1))
            write(error_unit,'(a)') prompt//"  "//this%titles(i)
        end do

        read(this%u) dummy, dummy
        if (dummy .ne. 4) then
            call error_stop_program("This DCD file format is not supported, or the file header is corrupt.")
        end if

        ! Number of atoms in each snapshot
        read(this%u) natoms, dummy
        if (dummy .ne. 4) then
            call error_stop_program("This DCD file format is not supported, or the file header is corrupt.")
        end if

        inquire(unit=this%u, pos=pos)
        pos = pos - 1

        ! Each frame has natoms*3 (4 bytes each) = natoms*12
        ! plus 6 box dimensions (8 bytes each) = 48
        ! Additionally there are 32 bytes of file information in each frame
        this%framesize = natoms*12 + 80
        ! Header is typically 276 bytes, but inquire gives us exact size
        nframes2 = (this%filesize-pos)/this%framesize
        if ( nframes2 .ne. nframes) then
            write(error_unit,'(a,i0,a,i0,a)') prompt//"WARNING: Header indicates ", nframes, &
                &" frames, but file size indicates ", nframes2, "." 
            nframes = int(nframes2)
        end if

    end subroutine dcdfile_read_header

    !> @brief Closes DCD file
    !> @param[inout] this dcdreader object
    subroutine dcdfile_close(this)

        implicit none
        class(dcdfile), intent(inout) :: this

        deallocate(this%titles)
        close(this%u)

    end subroutine dcdfile_close

    !> @brief Reads next frame into memory
    !> @param[inout] this dcdreader object
    !> @param[inout] xyz coordinates of all atoms in snapshot. First dimension is atom number, 
    !! second dimension is x, y, or z !coordinate.
    !> @param[inout] box dimensions of snapshot. Array containing 6 elements, ordered as A, B, C, alpha, beta, gamma.
    !! A = length of unit cell vector along x-axis;
    !! B = length of unit cell vector in xy-plane;
    !! C = length of unit cell vector in yz-plane;
    !! alpha = cosine of angle between B and C;
    !! beta = cosine of angle between A and C;
    !! gamma = cosine of angle between A and B;

    subroutine dcdfile_read_next(this, xyz, box, error)

        implicit none
        real(kind=real32), allocatable, intent(inout) :: xyz(:,:)
        real(kind=real64), intent(inout) :: box(6)
        integer(kind=int32) :: dummy(6), nbytes, ndims, i
        class(dcdfile), intent(inout) :: this
        logical :: error
        integer(kind=int32), save :: saveIO
        logical, save :: firstTimein = .true.
        integer :: ios

        if (firstTimeIn) then
            firstTimeIn = .false.
            saveIO = this%u
        else
            this%u = saveIO
        end if

        error = .true.
        nbytes = size(xyz,2)*4
        ndims = size(xyz,1)

        if (ndims /= 3) then
            error = .false.
            call error_stop_program("Number of dimensions of xyz array is incorrect.")
            return
        end if
    
        read(this%u,iostat=ios) dummy(1)
        if (ios/=0) then
            error = .false.
            return
        end if

        if (dummy(1) /= 48) then
            error = .false.
            call error_stop_program("Problem reading in DCD snapshot.")
            return
        end if

        !            A       gamma   B       beta    alpha   C
        read(this%u) box(1), box(6), box(2), box(5), box(4), box(3)
        if (box(1) < 0 .or. box(2) < 0 .or. box(3) < 0) then
            error = .false.
            call error_stop_program("Problem reading in DCD snapshot box dimensions.")
            return
        end if

        ! 48, then no. of bytes for x coordinates, x coordinates (repeat for y and z coordinates)
        read(this%u) dummy(1:2), xyz(1,:), dummy(3:4), xyz(2,:), dummy(5:6), xyz(3,:)

        if (dummy(1) /= 48) then
            error = .false.
            call error_stop_program("Problem reading in DCD snapshot coordinates.")
            return
        end if

        do i = 2, 6
            if (dummy(i) /= nbytes) then
                error = .false.
                call error_stop_program("Number of bytes in DCD snapshot is incorrect for size of xyz array passed.")
                return
            end if
        end do

        read(this%u) dummy(1)
        if (dummy(1) .ne. nbytes) then
            error = .false.
            call error_stop_program("Problem reading in DCD snapshot.")
            return
        end if

    end subroutine dcdfile_read_next

    !> @brief Skips reading this frame into memory
    !> @param[inout] this dcdreader object
    !> @param[in] n number of frames to skip
    subroutine dcdfile_skip_next(this, n)

        implicit none
        integer(kind=int32) :: dummy
        integer(kind=int64) :: pos, newpos
        integer(kind=int32), intent(in), optional :: n
        class(dcdfile), intent(inout) :: this
   
        ! Where are we?
        inquire(unit=this%u, pos=pos)

        ! We subtract 4 bytes so that the next read of the 4-byte integer will line things up properly for the next read
        if (.not. present(n)) then
            newpos = pos + this%framesize - 4
        else
            newpos = pos + this%framesize*n - 4
        end if
        
        read(this%u, pos=newpos) dummy

    end subroutine dcdfile_skip_next

end module dcdfort_reader
