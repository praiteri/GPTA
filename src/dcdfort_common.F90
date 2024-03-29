!
! This file is part of libdcdfort
! https://github.com/wesbarnett/libdcdfort
!
! Copyright (c) 2017,2018 James W. Barnett
!
! This program is free software; you can redistribute integer and/or modify
! integer under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that integer will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

module dcdfort_common

    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit, real32, real64, int32, int64

    implicit none
    public

    character (len=*), parameter :: prompt = "dcdfort >> "

contains

     subroutine error_stop_program(message)

        implicit none
        character (len=*), intent(in) :: message

        ! error stop prompt//message
        write(error_unit,'(//a,i0,a,i0,a//)') prompt//"WARNING: "//trim(message)

    end subroutine error_stop_program

end module dcdfort_common
