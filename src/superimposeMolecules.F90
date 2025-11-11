! Rotate "mol" onto "ref"
!
!  natoms = 4
!  allocate(mol(3,natoms))
!  allocate(ref(3,natoms))
!  mol(1:3,1) = [27.828_real64,  21.953_real64,  25.702_real64]
!  mol(1:3,2) = [27.134_real64,  20.906_real64,  25.345_real64]
!  mol(1:3,3) = [27.679_real64,  22.463_real64,  26.872_real64]
!  mol(1:3,4) = [28.536_real64,  22.546_real64,  24.796_real64]
!  ref(1:3,1) = [25.830_real64,  24.186_real64,  24.597_real64]
!  ref(1:3,2) = [26.866_real64,  24.817_real64,  24.984_real64]
!  ref(1:3,3) = [25.013_real64,  24.744_real64,  23.792_real64]
!  ref(1:3,4) = [25.618_real64,  23.010_real64,  25.015_real64]
!
!  allocate(new(3,natoms))
!  call Superimpose(natoms, mol, ref, rmsd, rotmat)
!  do i=1,natoms
!    new(1:3,i) = matmul(rotmat,mol(1:3,i))
!  enddo
!
!  write(0,*) 3*natoms
!  write(0,*)
!  do i=1,natoms
!    write(0,*) "C", mol(1:3,i)
!  end do
!  do i=1,natoms
!    write(0,*) "O", ref(1:3,i)
!  end do
!  do i=1,natoms
!    write(0,*) "S", new(1:3,i)
!  end do

module moduleSuperimposeMolecules

#define PREC 8
  private
  public :: Superimpose, CenterCoords, computeCenter
  
contains

  subroutine Superimpose(nlen, coords1, coords2, rmsd, rot, weight)
    implicit none 
    
    integer, intent(in)  :: nlen
    real(PREC), intent(inout) :: coords1(3,nlen)
    real(PREC), intent(inout) :: coords2(3,nlen)
    real(PREC), intent(out) :: rmsd
    real(PREC), intent(out) :: rot(9)
    real(PREC), intent(in), optional :: weight(nlen)
    
    integer :: ireturn
    real(PREC) :: shift(3)
    real(PREC) :: A(9), E0
    
    if (present(weight)) then
      ! center the structures -- if precentered you can omit this step
      ! call CenterCoords(nlen, coords1, weight)
      call CenterCoords(shift, nlen, coords2, weight)
      
      ! calculate the (weighted) inner product of two structures
      call InnerProduct(nlen, coords1, coords2, A, E0, weight)
    else
      ! center the structures -- if precentered you can omit this step
      ! call CenterCoords(nlen, coords1)
      call CenterCoords(shift, nlen, coords2)
      
      ! calculate the inner product of two structures
      call InnerProduct(nlen, coords1, coords2, A, E0)
    endif
    
    ! calculate the RMSD & rotational matrix
    call FastCalcRMSDAndRotation(A, -1, nlen, E0, rmsd, rot, ireturn)
  
  end subroutine Superimpose


! Fotran version of Quaternion Characteristic Polynomial (QCP) 
! Copyright (c) 2016 Naoto Hori
!
! The original code written in c is distributed at http://theobald.brandeis.edu/qcp
!
! Following is the original notice by the original authors.
!
!/*******************************************************************************
! *  -/_|:|_|_\- 
! *
! *  File:           qcprot.c
! *  Version:        1.4
! *
! *  Function:       Rapid calculation of the least-squares rotation using a 
! *                  quaternion-based characteristic polynomial and 
! *                  a cofactor matrix
! *
! *  Author(s):      Douglas L. Theobald
! *                  Department of Biochemistry
! *                  MS 009
! *                  Brandeis University
! *                  415 South St
! *                  Waltham, MA  02453
! *                  USA
! *
! *                  dtheobald@brandeis.edu
! *                  
! *                  Pu Liu
! *                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
! *                  665 Stockton Drive
! *                  Exton, PA  19341
! *                  USA
! *
! *                  pliu24@its.jnj.com
! * 
! *
! *    If you use this QCP rotation calculation method in a publication, please
! *    reference:
! *
! *      Douglas L. Theobald (2005)
! *      "Rapid calculation of RMSD using a quaternion-based characteristic
! *      polynomial."
! *      Acta Crystallographica A 61(4):478-480.
! *
! *      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
! *      "Fast determination of the optimal rotational matrix for macromolecular 
! *      superpositions."
! *      Journal of Computational Chemistry 31(7):1561-1563.
! *
! *
! *  Copyright (c) 2009-2013 Pu Liu and Douglas L. Theobald
! *  All rights reserved.
! *
! *  Redistribution and use in source and binary forms, with or without modification, are permitted
! *  provided that the following conditions are met:
! *
! *  * Redistributions of source code must retain the above copyright notice, this list of
! *    conditions and the following disclaimer.
! *  * Redistributions in binary form must reproduce the above copyright notice, this list
! *    of conditions and the following disclaimer in the documentation and/or other materials
! *    provided with the distribution.
! *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to
! *    endorse or promote products derived from this software without specific prior written
! *    permission.
! *
! *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
! *
! *  Source:         started anew.
! *
! *  Change History:
! *    2009/04/13      Started source
! *    2010/03/28      Modified FastCalcRMSDAndRotation() to handle tiny qsqr
! *                    If trying all rows of the adjoint still gives too small
! *                    qsqr, then just return identity matrix. (DLT)
! *    2010/06/30      Fixed prob in assigning A[9] = 0 in InnerProduct()
! *                    invalid mem access
! *    2011/02/21      Made CenterCoords use weights
! *    2011/05/02      Finally changed CenterCoords declaration in qcprot.h
! *                    Also changed some functions to static
! *    2011/07/08      put in fabs() to fix taking sqrt of small neg numbers, fp error
! *    2012/07/26      minor changes to comments and main.c, more info (v.1.4)
! *  
! ******************************************************************************/

  subroutine FastCalcRMSDAndRotation(A, minScore, nlen, E0, rmsd, rot, ireturn)
    implicit none 
    
    real(PREC), intent(in)  :: A(9)
    integer,    intent(in)  :: minScore
    integer,    intent(in)  :: nlen
    real(PREC), intent(in)  :: E0
    real(PREC), intent(out) :: rmsd
    real(PREC), intent(out) :: rot(9)
    integer,    intent(out) :: ireturn
    
    real(PREC) :: Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz
    real(PREC) :: Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2, &
    SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2, &
    SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy, &
    SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy
    real(PREC) :: C(3)
    integer :: i
    real(PREC) :: mxEigenV
    real(PREC) :: oldg
    real(PREC) :: b, aa, delta, rms, qsqr
    real(PREC) :: q1, q2, q3, q4, normq
    real(PREC) :: a11, a12, a13, a14, a21, a22, a23, a24
    real(PREC) :: a31, a32, a33, a34, a41, a42, a43, a44
    real(PREC) :: a2, x2, y2, z2
    real(PREC) :: xy, az, zx, ay, yz, ax
    real(PREC) :: a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132
    real(PREC), parameter :: evecprec = 1.0e-6
    real(PREC), parameter :: evalprec = 1.0e-11
    real(PREC) :: a1324_1423, a1224_1422, a1223_1322, a1124_1421, a1123_1321, a1122_1221
    
    Sxx = A(1)
    Sxy = A(2)
    Sxz = A(3)
    Syx = A(4)
    Syy = A(5)
    Syz = A(6)
    Szx = A(7)
    Szy = A(8)
    Szz = A(9)
    
    Sxx2 = Sxx * Sxx
    Syy2 = Syy * Syy
    Szz2 = Szz * Szz
    
    Sxy2 = Sxy * Sxy
    Syz2 = Syz * Syz
    Sxz2 = Sxz * Sxz
    
    Syx2 = Syx * Syx
    Szy2 = Szy * Szy
    Szx2 = Szx * Szx
    
    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2
    
    C(3) = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C(2) =  8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz)
    
    SxzpSzx = Sxz + Szx
    SyzpSzy = Syz + Szy
    SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy
    SxzmSzx = Sxz - Szx
    SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy
    SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2
    
    C(1) = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 &
    + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) &
    + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz)) &
    + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz)) &
    + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz)) &
    + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz))
    
    ! Newton-Raphson
    mxEigenV = E0
    do i = 1, 50
      oldg = mxEigenV
      x2 = mxEigenV*mxEigenV
      b = (x2 + C(3))*mxEigenV
      aa = b + C(2)
      delta = ((aa*mxEigenV + C(1))/(2.0*x2*mxEigenV + b + aa))
      mxEigenV = mxEigenV - delta
      ! write(*,*) "diff[",i,"]:",mxEigenV-oldg, evalprec*mxEigenV, mxEigenV
      if (abs(mxEigenV - oldg) < abs(evalprec*mxEigenV)) then
        exit
      endif
    enddo
    
    if (i == 50) then
      write(*,*) "More than",i,"iterations needed!"
    endif
    
    ! the abs() is to guard against extremely small, but *negative* numbers due to floating point error
    rms = sqrt(abs(2.0 * (E0 - mxEigenV)/real(nlen,kind=PREC)))
    rmsd = rms
    !write(*,*)  rms, E0, 2.0 * (E0 - mxEigenV)/len
    
    if ((minScore > 0) .and. (rms < minScore)) then
      ireturn = -1 ! Don't bother with rotation. 
      return
    endif
    
    a11 = SxxpSyy + Szz-mxEigenV
    a12 = SyzmSzy
    a13 = - SxzmSzx
    a14 = SxymSyx
    a21 = SyzmSzy
    a22 = SxxmSyy - Szz-mxEigenV
    a23 = SxypSyx
    a24= SxzpSzx
    a31 = a13
    a32 = a23
    a33 = Syy-Sxx-Szz - mxEigenV
    a34 = SyzpSzy
    a41 = a14
    a42 = a24
    a43 = a34
    a44 = Szz - SxxpSyy - mxEigenV
    a3344_4334 = a33 * a44 - a43 * a34
    a3244_4234 = a32 * a44-a42*a34
    a3243_4233 = a32 * a43 - a42 * a33
    a3143_4133 = a31 * a43-a41*a33
    a3144_4134 = a31 * a44 - a41 * a34
    a3142_4132 = a31 * a42-a41*a32
    q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133
    q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132
    q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132
    
    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4
    
    !The following code tries to calculate another column in the adjoint matrix when the norm of the 
    !current column is too small.
    !Usually this block will never be activated.  To be absolutely safe this should be
    !uncommented, but it is most likely unnecessary.
    
    if (qsqr < evecprec) then
      q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233
      q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133
      q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132
      q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132
      qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4
      
      if (qsqr < evecprec) then
        a1324_1423 = a13 * a24 - a14 * a23
        a1224_1422 = a12 * a24 - a14 * a22
        a1223_1322 = a12 * a23 - a13 * a22
        a1124_1421 = a11 * a24 - a14 * a21
        a1123_1321 = a11 * a23 - a13 * a21
        a1122_1221 = a11 * a22 - a12 * a21
        
        q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322
        q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321
        q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221
        q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4
        
        if (qsqr < evecprec) then
          q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322
          q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321
          q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221
          q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221
          qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4
          
          if (qsqr < evecprec) then
            ! if qsqr is still too small, return the identity matrix.
            rot(:) = 0.0e0
            rot(1) = 1.0e0
            rot(5) = 1.0e0
            rot(9) = 1.0e0
            
            ireturn = 0
            return
          endif
        endif
      endif
    endif
    
    normq = sqrt(qsqr)
    q1 = q1 / normq
    q2 = q2 / normq
    q3 = q3 / normq
    q4 = q4 / normq
    
    a2 = q1 * q1
    x2 = q2 * q2
    y2 = q3 * q3
    z2 = q4 * q4
    
    xy = q2 * q3
    az = q1 * q4
    zx = q4 * q2
    ay = q1 * q3
    yz = q3 * q4
    ax = q1 * q2
    
    rot(1) = a2 + x2 - y2 - z2
    rot(2) = 2 * (xy + az)
    rot(3) = 2 * (zx - ay)
    rot(4) = 2 * (xy - az)
    rot(5) = a2 - x2 + y2 - z2
    rot(6) = 2 * (yz + ax)
    rot(7) = 2 * (zx + ay)
    rot(8) = 2 * (yz - ax)
    rot(9) = a2 - x2 - y2 + z2
    
    ireturn = 1
  end subroutine FastCalcRMSDAndRotation

  subroutine InnerProduct(nlen, coords1, coords2, A, E0, weight)
    implicit none 
    
    integer,    intent(in) :: nlen
    real(PREC), intent(in) :: coords1(3,nlen)
    real(PREC), intent(in) :: coords2(3,nlen)
    real(PREC), intent(out) :: A(9)
    real(PREC), intent(out) :: E0
    real(PREC), intent(in), optional :: weight(nlen)
    
    integer :: i
    real(PREC) :: G
    
    A(:) = 0.0e0
    G = 0.0e0
    
    if (present(weight)) then
      do i = 1, nlen
        G = G + dot_product(coords1(1:3,i), coords1(1:3,i)) * weight(i) &
        + dot_product(coords2(1:3,i), coords2(1:3,i)) * weight(i)
        
        A(1:3) = A(1:3) + coords1(1,i) * coords2(1:3,i)
        A(4:6) = A(4:6) + coords1(2,i) * coords2(1:3,i)
        A(7:9) = A(7:9) + coords1(3,i) * coords2(1:3,i)
      enddo
    else
      do i = 1, nlen
        G = G + dot_product( coords1(1:3,i), coords1(1:3,i) ) &
        + dot_product( coords2(1:3,i), coords2(1:3,i) )
        
        A(1:3) = A(1:3) + coords1(1,i) * coords2(1:3,i)
        A(4:6) = A(4:6) + coords1(2,i) * coords2(1:3,i)
        A(7:9) = A(7:9) + coords1(3,i) * coords2(1:3,i)
      enddo
    endif
    
    E0 = G * 0.5
  end subroutine InnerProduct
  
  subroutine CenterCoords(s, nlen, coords, weight)
    
    implicit none 
    
    real(PREC), intent(out) :: s(3)
    integer, intent(in) :: nlen
    real(PREC), intent(inout) :: coords(3,nlen)
    real(PREC), intent(in), optional :: weight(nlen)
    
    integer :: i
    real(PREC) :: wsum
    
    s(:) = 0.0
    
    if (present(weight)) then
      wsum = 0.0
      do i = 1, nlen
        s(:) = s(:) + weight(i) * coords(:,i)
        wsum = wsum + weight(i)
      enddo
      
      s(:) = s(:) / wsum
    else
      do i = 1, nlen
        s(:) = s(:) + coords(:,i)
      enddo
      
      s(:) = s(:) / real(nlen, kind=PREC)
    endif
    
    do i = 1, nlen
      coords(:,i) = coords(:,i) - s(:)
    enddo
  end subroutine CenterCoords
  
  function computeCenter(nlen, coords, weight) result(s)
    
    implicit none 
    
    real(PREC) :: s(3)
    integer, intent(in) :: nlen
    real(PREC), intent(in) :: coords(3,nlen)
    real(PREC), intent(in), optional :: weight(nlen)
    
    integer :: i
    real(PREC) :: wsum
    
    s(:) = 0.0
    
    if (present(weight)) then
      wsum = 0.0
      do i = 1, nlen
        s(:) = s(:) + weight(i) * coords(:,i)
        wsum = wsum + weight(i)
      enddo
      
      s(:) = s(:) / wsum
    else
      do i = 1, nlen
        s(:) = s(:) + coords(:,i)
      enddo
      
      s(:) = s(:) / real(nlen, kind=PREC)
    endif
    
  end function computeCenter

  ! subroutine InnerProductNoWeight(nlen, coords1, coords2, A, E0)
    
  !   implicit none 
    
  !   integer,    intent(in) :: nlen
  !   real(PREC), intent(in) :: coords1(3,nlen)
  !   real(PREC), intent(in) :: coords2(3,nlen)
  !   real(PREC), intent(out) :: A(9)
  !   real(PREC), intent(out) :: E0
    
  !   integer :: i
  !   real(PREC) :: G
    
  !   A(:) = 0.0e0
  !   G = 0.0e0
    
  !   do i = 1, nlen
  !     G = G + dot_product( coords1(1:3,i), coords1(1:3,i) ) &
  !     + dot_product( coords2(1:3,i), coords2(1:3,i) )
      
  !     A(1:3) = A(1:3) + coords1(1,i) * coords2(1:3,i)
  !     A(4:6) = A(4:6) + coords1(2,i) * coords2(1:3,i)
  !     A(7:9) = A(7:9) + coords1(3,i) * coords2(1:3,i)
  !   enddo
    
  !   E0 = G * 0.5
  ! end subroutine InnerProductNoWeight
  
  ! subroutine CenterCoordsNoWeight(nlen, coords)
    
  !   implicit none 
    
  !   integer, intent(in) :: nlen
  !   real(PREC), intent(inout) :: coords(3,nlen)
    
  !   integer :: i
  !   real(PREC) :: s(3)
    
  !   s(:) = 0.0
    
  !   do i = 1, nlen
  !     s(:) = s(:) + coords(:,i)
  !   enddo
    
  !   s(:) = s(:) / real(nlen, kind=PREC)
    
  !   do i = 1, nlen
  !     coords(:,i) = coords(:,i) - s(:)
  !   enddo
  ! end subroutine CenterCoordsNoWeight

end module moduleSuperimposeMolecules
