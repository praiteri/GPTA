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
module moduleDistances
  implicit none

  interface
    function func_vec(dij)
      implicit none
      real(8)                :: func_vec
      real(8), intent(inout) :: dij(3)
    end function func_vec

    subroutine cart2frac(n,cart,frac)
      implicit none
      integer, intent(in) :: n
      real(8), dimension(3,n), intent(in) :: cart
      real(8), dimension(3,n), intent(out) :: frac
    end subroutine cart2frac

    subroutine frac2cart(n,cart,frac)
      implicit none
      integer, intent(in) :: n
      real(8), dimension(3,n), intent(in) :: cart
      real(8), dimension(3,n), intent(out) :: frac
    end subroutine frac2cart

    subroutine cart2fracNINT(n,cart,frac)
      implicit none
      integer, intent(in) :: n
      real(8), dimension(3,n), intent(in) :: cart
      real(8), dimension(3,n), intent(out) :: frac
    end subroutine cart2fracNINT

    subroutine cart2fracNoWrap(n,cart,frac)
      implicit none
      integer, intent(in) :: n
      real(8), dimension(3,n), intent(in) :: cart
      real(8), dimension(3,n), intent(out) :: frac
    end subroutine cart2fracNoWrap

  end interface
  procedure (func_vec), pointer :: computeDistanceSquaredPBC => null ()
  procedure (cart2frac), pointer :: cartesianToFractional => null ()
  procedure (frac2cart), pointer :: fractionalToCartesian => null ()
  procedure (cart2fracNINT), pointer :: cartesianToFractionalNINT => null ()
  procedure (cart2fracNoWrap), pointer :: cartesianToFractionalNoWrap => null ()

contains
  subroutine initialisePBC(lflag)
    use moduleVariables, only : safeDistances
    implicit none
    character(len=*), intent(inout) :: lflag
    
    if (safeDistances) lflag = 'safe'

    computeDistanceSquaredPBC => null ()
    if     (lflag=="none") then
      computeDistanceSquaredPBC => computeDistanceSquaredVacuum
      cartesianToFractional => copyCoordinates
      fractionalToCartesian => copyCoordinates
      cartesianToFractionalNINT => copyCoordinates
      cartesianToFractionalNoWrap => copyCoordinates

    else if (lflag=="ortho") then
      computeDistanceSquaredPBC => computeDistanceSquaredOrthogonal
      cartesianToFractional => cartesianToFractionalOrthogonal
      fractionalToCartesian => fractionalToCartesianTriclinic
      cartesianToFractionalNINT => cartesianToFractionalOrthogonalNINT
      cartesianToFractionalNoWrap => cartesianToFractionalOrthogonalNoWrap
      
    else if (lflag=="tri" .or. lflag=="safe") then

      if (lflag == "safe") then
        computeDistanceSquaredPBC => computeDistanceSquaredSafe
      else
        computeDistanceSquaredPBC => computeDistanceSquaredTriclinic
      end if
      cartesianToFractional => cartesianToFractionalTriclinic
      fractionalToCartesian => fractionalToCartesianTriclinic
      cartesianToFractionalNINT => cartesianToFractionalTriclinicNINT
      cartesianToFractionalNoWrap => cartesianToFractionalTriclinicNoWrap
    end if
    
    return
  end subroutine initialisePBC

  function computeDistanceSquaredSafe(dij) result(distance)
    use moduleSystem , only : frame
    implicit none
    real(8) :: distance
    real(8), dimension(3), intent(inout) :: dij
    real(8), dimension(3) :: sij, lsafe_pbc_vec
    
    integer :: i, j, k
    real(8) :: rtmp
   
    distance = sum(dij*dij)
    lsafe_pbc_vec(1:3) = dij(1:3)
    do i=-1,1
      do j=-1,1
        do k=-1,1
          if (i==0 .and. j==0 .and. k==0) cycle
          sij = dij + i*frame % hmat(:,1) + j*frame % hmat(:,2) + k*frame % hmat(:,3)
          rtmp=sum(sij*sij)
          if (rtmp<distance) then
            distance = rtmp
            lsafe_pbc_vec(1:3) = sij(1:3)
          end if
        enddo
      enddo
    enddo
    dij = lsafe_pbc_vec(1:3) 

  end function computeDistanceSquaredSafe

  function computeDistanceSquaredVacuum(dij) result(distance)
    implicit none
    real(8) :: distance
    real(8), dimension(3), intent(inout) :: dij

    distance = sum(dij*dij)

  end function computeDistanceSquaredVacuum

  function computeDistanceSquaredOrthogonal(dij) result(distance)
    use moduleSystem , only : frame
    implicit none
    real(8) :: distance
    real(8), dimension(3), intent(inout) :: dij
    real(8), dimension(3) :: sij

    sij(1) = dij(1) / frame % cell(1)
    sij(2) = dij(2) / frame % cell(2)
    sij(3) = dij(3) / frame % cell(3)
    sij(1:3) = sij(1:3)-nint(sij(1:3))
    dij(1) = sij(1) * frame % cell(1)
    dij(2) = sij(2) * frame % cell(2)
    dij(3) = sij(3) * frame % cell(3)
    distance = sum(dij*dij)

  end function computeDistanceSquaredOrthogonal

  function computeDistanceSquaredTriclinic(dij) result(distance)
    use moduleSystem , only : frame
    implicit none
    real(8) :: distance
    real(8), dimension(3), intent(inout) :: dij
    real(8), dimension(3) :: sij
    
    ! sij(1) = frame % hinv(1,1)*dij(1) + frame % hinv(1,2)*dij(2) + frame % hinv(1,3)*dij(3)
    ! sij(2) = frame % hinv(2,1)*dij(1) + frame % hinv(2,2)*dij(2) + frame % hinv(2,3)*dij(3)
    ! sij(3) = frame % hinv(3,1)*dij(1) + frame % hinv(3,2)*dij(1) + frame % hinv(3,3)*dij(3)
    ! sij(1:3) = sij(1:3)-nint(sij(1:3))
    ! dij(1) = frame % hmat(1,1)*sij(1) + frame % hmat(1,2)*sij(2) + frame % hmat(1,3)*sij(3)
    ! dij(2) = frame % hmat(2,1)*dij(1) + frame % hmat(2,2)*sij(2) + frame % hmat(2,3)*sij(3)
    ! dij(3) = frame % hmat(3,1)*dij(1) + frame % hmat(3,2)*dij(1) + frame % hmat(3,3)*sij(3)

    sij = matmul(frame % hinv,dij)
    sij(1:3) = sij(1:3) - nint(sij(1:3))
    dij = matmul(frame % hmat,sij)
    distance = sum(dij*dij)

  end function computeDistanceSquaredTriclinic

  subroutine copyCoordinates(n,cart,frac)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(in) :: cart
    real(8), dimension(3,n), intent(out) :: frac

    integer :: i
    do i=1,n
      frac(1:3,i) = cart(1:3,i)
    enddo
  end subroutine copyCoordinates

  subroutine cartesianToFractionalOrthogonal(n,cart,frac)
    use moduleSystem , only : frame
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(in) :: cart
    real(8), dimension(3,n), intent(out) :: frac

    integer :: i
    do i=1,n
      frac(1:3,i) = cart(1:3,i) / frame % cell(1:3)
      frac(1:3,i) = frac(1:3,i) - floor(frac(1:3,i))
    enddo
  end subroutine cartesianToFractionalOrthogonal

  subroutine cartesianToFractionalTriclinic(n,cart,frac) 
    use moduleSystem , only : frame
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(in) :: cart
    real(8), dimension(3,n), intent(out) :: frac

    integer :: i
    do i=1,n
      ! frac(1,i) = frame % hinv(1,1)*cart(1,i) + frame % hinv(1,2)*cart(2,i) + frame % hinv(1,3)*cart(3,i)
      ! frac(2,i) = frame % hinv(2,1)*cart(1,i) + frame % hinv(2,2)*cart(2,i) + frame % hinv(2,3)*cart(3,i)
      ! frac(3,i) = frame % hinv(3,1)*cart(1,i) + frame % hinv(3,2)*cart(2,i) + frame % hinv(3,3)*cart(3,i)
      frac(1:3,i) = matmul(frame % hinv , cart(1:3,i))
      frac(1:3,i) = frac(1:3,i) - floor(frac(1:3,i))
    enddo
  end subroutine cartesianToFractionalTriclinic

  subroutine fractionalToCartesianOrthogonal(n,frac,cart) 
    use moduleSystem , only : frame
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(in) :: frac
    real(8), dimension(3,n), intent(out) :: cart

    integer :: i
    do i=1,n
      cart(1:3,i) = frac(1:3,i) * frame % cell(1:3)
    enddo
  end subroutine fractionalToCartesianOrthogonal

  subroutine fractionalToCartesianTriclinic(n,frac,cart) 
    use moduleSystem , only : frame
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(in) :: frac
    real(8), dimension(3,n), intent(out) :: cart

    integer :: i
    do i=1,n
      ! cart(1,i) = frame % hmat(1,1)*frac(1,i) + frame % hmat(1,2)*frac(2,i) + frame % hmat(1,3)*frac(3,i)
      ! cart(2,i) = frame % hmat(2,1)*frac(1,i) + frame % hmat(2,2)*frac(2,i) + frame % hmat(2,3)*frac(3,i)
      ! cart(3,i) = frame % hmat(3,1)*frac(1,i) + frame % hmat(3,2)*frac(2,i) + frame % hmat(3,3)*frac(3,i)
      cart(1:3,i) = matmul(frame % hmat , frac(1:3,i))
    enddo 
  end subroutine fractionalToCartesianTriclinic

  subroutine cartesianToFractionalOrthogonalNINT(n,cart,frac) 
    use moduleSystem , only : frame
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(in) :: cart
    real(8), dimension(3,n), intent(out) :: frac

    integer :: i
    do i=1,n
      frac(1:3,i) = cart(1:3,i) / frame % cell(1:3)
      frac(1:3,i) = frac(1:3,i) - nint(frac(1:3,i))
    enddo
  end subroutine cartesianToFractionalOrthogonalNINT

  subroutine cartesianToFractionalTriclinicNINT(n,cart,frac) 
    use moduleSystem , only : frame
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(in) :: cart
    real(8), dimension(3,n), intent(out) :: frac

    integer :: i
    do i=1,n
      ! frac(1,i) = frame % hinv(1,1)*cart(1,i) + frame % hinv(1,2)*cart(2,i) + frame % hinv(1,3)*cart(3,i)
      ! frac(2,i) = frame % hinv(2,1)*cart(1,i) + frame % hinv(2,2)*cart(2,i) + frame % hinv(2,3)*cart(3,i)
      ! frac(3,i) = frame % hinv(3,1)*cart(1,i) + frame % hinv(3,2)*cart(2,i) + frame % hinv(3,3)*cart(3,i)
      frac(1:3,i) = matmul(frame % hinv , cart(1:3,i))
      frac(1:3,i) = frac(1:3,i) - nint(frac(1:3,i))
    enddo
  end subroutine cartesianToFractionalTriclinicNINT


  subroutine cartesianToFractionalOrthogonalNoWrap(n,cart,frac) 
    use moduleSystem , only : frame
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(in) :: cart
    real(8), dimension(3,n), intent(out) :: frac

    integer :: i
    do i=1,n
      frac(1:3,i) = cart(1:3,i) / frame % cell(1:3)
    enddo
  end subroutine cartesianToFractionalOrthogonalNoWrap

  subroutine cartesianToFractionalTriclinicNoWrap(n,cart,frac) 
    use moduleSystem , only : frame
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(in) :: cart
    real(8), dimension(3,n), intent(out) :: frac

    integer :: i
    do i=1,n
      ! frac(1,i) = frame % hinv(1,1)*cart(1,i) + frame % hinv(1,2)*cart(2,i) + frame % hinv(1,3)*cart(3,i)
      ! frac(2,i) = frame % hinv(2,1)*cart(1,i) + frame % hinv(2,2)*cart(2,i) + frame % hinv(2,3)*cart(3,i)
      ! frac(3,i) = frame % hinv(3,1)*cart(1,i) + frame % hinv(3,2)*cart(2,i) + frame % hinv(3,3)*cart(3,i)
      frac(1:3,i) = matmul(frame % hinv , cart(1:3,i))
    enddo
  end subroutine cartesianToFractionalTriclinicNoWrap

  subroutine wrapCoordinates(n,pos)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(3,n), intent(inout) :: pos
    real(8), allocatable, dimension(:,:) :: ptmp

    allocate(ptmp(3,n))

    call cartesianToFractional(n,pos,ptmp)
    call fractionalToCartesian(n,ptmp,pos)

  end subroutine wrapCoordinates
 
end module moduleDistances
