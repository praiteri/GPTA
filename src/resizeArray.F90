!disclaimer
module moduleResizeArrays 
  use moduleVariables, only: real64
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
    real(real64), allocatable, dimension(:) :: array
    real(real64), allocatable, dimension(:) :: buffer
    real(real64), allocatable, dimension(:) :: zeroes
    integer :: m, o
    if (n==0) then
      deallocate(array)
      return
    end if
    m = size(array)
    if (m==n) return
    if (n>m) then
      o=n-m 
      allocate(zeroes(o) , source=0.0_real64)
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
    real(real64), allocatable, dimension(:,:) :: array
    real(real64), allocatable, dimension(:,:) :: buffer
    real(real64), allocatable, dimension(:,:) :: zeroes
    integer :: m(2), o
    if (n==0) then
      deallocate(array)
      return
    end if
    m = shape(array)
    if (m(2)==n) return
    if (n>m(2)) then
      o=n-m(2)
      allocate(zeroes(m(1),o) , source=0.0_real64)
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
