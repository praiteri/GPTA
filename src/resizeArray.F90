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
module moduleResizeArrays 
  
  interface resizeArray
    module procedure resizeIntegerArray1D, resizeIntegerArray2D, resizeDoubleArray1D, resizeDoubleArray2D, resizeCharacterArray1D
  end interface

contains
  
  subroutine resizeCharacterArray1D(array,n)
    implicit none
    integer, parameter :: cp=4
    integer, intent(in) :: n
    character(len=cp), allocatable, dimension(:) :: array
    character(len=cp), allocatable, dimension(:) :: buffer
    character(len=cp), allocatable, dimension(:) :: zeroes
    integer :: m, o
    if (n==0) then
      deallocate(array)
      return
    end if
    m = size(array)
    if (m==n) return
    if (n>m) then
      o=n-m
      allocate(zeroes(o) , source="    ")
      allocate(buffer(n) , source=[array(1:m),zeroes])
      deallocate(zeroes)
    else
      allocate(buffer(n) , source=array(1:n))
    end if
    call move_alloc(buffer,array)
    return
  end subroutine resizeCharacterArray1D

  subroutine resizeIntegerArray1D(array,n)
    implicit none
    integer, intent(in) :: n
    integer, allocatable, dimension(:) :: array
    integer, allocatable, dimension(:) :: buffer
    integer, allocatable, dimension(:) :: zeroes
    integer :: m, o
    if (n==0) then
      deallocate(array)
      return
    end if
    m = size(array)
    if (m==n) return
    if (n>m) then
      o=n-m 
      allocate(zeroes(o) , source=0)
      allocate(buffer(n) , source=[array(1:m),zeroes])
      deallocate(zeroes)
    else
      allocate(buffer(n) , source=array(1:n))
    end if
    call move_alloc(buffer,array)
    return
  end subroutine resizeIntegerArray1D

  subroutine resizeIntegerArray2D(array,n)
    implicit none
    integer, intent(in) :: n
    integer, allocatable, dimension(:,:) :: array
    integer, allocatable, dimension(:,:) :: buffer
    integer, allocatable, dimension(:,:) :: zeroes
    integer :: m(2), o
    if (n==0) then
      deallocate(array)
      return
    end if
    m = shape(array)
    if (m(2)==n) return
    if (n>m(2)) then
      o=n-m(2)
      allocate(zeroes(m(1),o) , source=0)
      allocate(buffer(m(1),n))
      buffer(1:m(1),1:m(2)) = array(1:m(1),1:m(2))
      buffer(1:m(1),m(2)+1:) = zeroes
      deallocate(zeroes)
    else
      allocate(buffer(m(1),n) , source=array(1:m(1),1:n))
    end if
    call move_alloc(buffer,array)
    return
  end subroutine resizeIntegerArray2D

  subroutine resizeDoubleArray1D(array,n)
    implicit none
    integer, intent(in) :: n
    real(8), allocatable, dimension(:) :: array
    real(8), allocatable, dimension(:) :: buffer
    real(8), allocatable, dimension(:) :: zeroes
    integer :: m, o
    if (n==0) then
      deallocate(array)
      return
    end if
    m = size(array)
    if (m==n) return
    if (n>m) then
      o=n-m 
      allocate(zeroes(o) , source=0.d0)
      allocate(buffer(n) , source=[array(1:m),zeroes])
      deallocate(zeroes)
    else
      allocate(buffer(n) , source=array(1:n))
    end if
    call move_alloc(buffer,array)
    return
  end subroutine resizeDoubleArray1D

  subroutine resizeDoubleArray2D(array,n)
    implicit none
    integer, intent(in) :: n
    real(8), allocatable, dimension(:,:) :: array
    real(8), allocatable, dimension(:,:) :: buffer
    real(8), allocatable, dimension(:,:) :: zeroes
    integer :: m(2), o
    if (n==0) then
      deallocate(array)
      return
    end if
    m = shape(array)
    if (m(2)==n) return
    if (n>m(2)) then
      o=n-m(2)
      allocate(zeroes(m(1),o) , source=0.d0)
      allocate(buffer(m(1),n))
      buffer(1:m(1),1:m(2)) = array(1:m(1),1:m(2))
      buffer(1:m(1),m(2)+1:) = zeroes
      deallocate(zeroes)
    else
      allocate(buffer(m(1),n) , source=array(1:m(1),1:n))
    end if
    call move_alloc(buffer,array)
    return
  end subroutine resizeDoubleArray2D

end module moduleResizeArrays 
