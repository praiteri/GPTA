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
module moduleElements 
  use moduleVariables
  
  implicit none

  integer, parameter :: nelement=107
  type element
    character(len=cp) :: lab
    real(8) :: mass
    real(8) :: rcov
    real(8) :: rvdw
    real(8) :: rion
    real(8) :: chg
  end type
  type(element), target, dimension(0:nelement) :: atom

  character(cp), allocatable, dimension(:) :: ions_list

contains

  subroutine defineIonicSpecies(labels)
    implicit none
    character(*), dimension(:), intent(in) :: labels
    character(cp), allocatable, dimension(:) :: tmp_labels
    integer :: n, m

    if (allocated(ions_list)) then
      m = size(ions_list)
      allocate(tmp_labels(m))
      tmp_labels = ions_list
      deallocate(ions_list)
    else 
      m = 0
    end if

    n = size(labels)
    allocate(ions_list(m+n))
    if (m>0) then
      ions_list(1:m) = tmp_labels
      ions_list(m+1:m+n) = labels
    else
      ions_list = labels
    end if

    return
  end subroutine defineIonicSpecies

  function getAtomicNumber(label) result(id)
    implicit none
    character(*), intent(in) :: label
    integer :: id

    character(cp) :: element
    character(len=1), external :: uppercase
    integer :: i, ln

    element = ''
    do i=1,cp
      if (ichar(label(i:i)) >=48 .and. ichar(label(i:i))<=57 ) exit
      if (i==1) then
        element(i:i) = uppercase(label(i:i))
      else
        element(i:i) = label(i:i)
      end if
    enddo
    ln=len_trim(element)
    id=0
100 continue
    do i=1,nelement
      if ( len_trim(atom(i) % lab) /= ln ) cycle
      if(element(1:ln)==atom(i) % lab(1:ln)) then
        id = i
        return
      end if
    enddo

! This should ensure atoms to have letters after the species name
! for exaple CA or CB can be identified as carbon atoms
    if(id<=0 .and. ln>1) then
      ln = ln-1
      goto 100
    end if

    return
  end function getAtomicNumber

  function getElementMass(label) result(mass)
    implicit none
    character(cp), intent(in) :: label
    real(8) :: mass

    character(cp) :: element
    integer :: id
    integer :: i, ln

    mass = -1.d0
    element = ''
    do i=1,cp
      if (ichar(label(i:i)) >=48 .and. ichar(label(i:i))<=57 ) exit
      element(i:i) = label(i:i)
    enddo
    ln = len_trim(element)
    id = -1
100 do i=0,nelement
      if ( len_trim(atom(i) % lab) /= ln ) cycle
      if(element(1:ln)==atom(i) % lab(1:ln)) then
        id = i
        mass = atom(id) % mass
        return
      end if
    enddo

! This should ensure that atoms to have letters after the species name
! for exaple CD or CB can be identified as carbon atoms
    if (id<0 .and. ln>1) then
      ln = ln-1
      goto 100
    ! else
    !   call message(1,1,0,"Warning: unknown element "//label)
    end if

    return
  end function getElementMass

  function massToElement(mass) result(label)
    implicit none
    real(8), intent(in) :: mass
    character(len=cp) :: label

    integer :: i

    do i=1,nelement
      if(abs(mass-atom(i) % mass)<0.1)exit
    enddo
    if(i>nelement)i=0
    label = trim(atom(i) % lab)
    
    return
  end function massToElement

  function getMaximumBondLength(l1,l2,rscale) result (rbond)
    use moduleSystem , only : sbond_l, sbond_d, sbond_n
    implicit none
    character(len=cp) :: l1, l2
    real(8), intent(in) :: rscale
    real(8) :: rbond, r1, r2
    integer :: id, jd, i

    integer, dimension(4) :: h_allowed

    rbond = 0.0d0

    do i=1,sbond_n
      if ( (l1==sbond_l(1,i) .and. l2==sbond_l(2,i)) .or. (l1==sbond_l(2,i) .and. l2==sbond_l(1,i)) ) then
        rbond = sbond_d(i)
        return
      end if
    enddo

! ions
    if (allocated(ions_list)) then
      if (any(ions_list==l1) .or. any(ions_list==l2)) then
        rbond = -1.d0
        return
      end if
    end if

! Atoms that can bond to an H
! C, N, O, S
!    h_allowed = (/6, 7, 8, 9, 13, 14, 15, 16, 17, 35, 53/)
    h_allowed = (/6, 7, 8, 16/)

    id = getAtomicNumber(l1)
    jd = getAtomicNumber(l2)

    if (id<1) then
      r1 = 0.5
    else
      r1 = atom(id) % rcov
    end if

    if (jd<1) then
      r2 = 0.5
    else
      r2 = atom(jd) % rcov
    end if

    if (r1<0.001 .or. r2<0.001) return

    rbond = (r1 + r2) * rscale

! Exclude H-H bonds
    if (id==1) then
      if (all(h_allowed/=jd)) rbond = 0.0d0
    else if (jd==1) then
      if (all(h_allowed/=id)) rbond = 0.0d0
! Exclude O-O bonds
    else if (id==8 .and. jd==8) then
      rbond = 0.0d0
! O - Cl
    else if ( (id==8 .and. jd==17) .or. (id==17 .and. jd==8)  )then
      rbond = 0.0d0
! N - Cl
    else if ( (id==7 .and. jd==17) .or. (id==17 .and. jd==7)  )then
      rbond = 0.0d0
    end if

    return
  end function getMaximumBondLength

  function getMaximumOverlap(l1,l2,rscale) result (rdist)
    implicit none
    character(len=cp), intent(in) :: l1, l2
    real(8), optional, intent(in) :: rscale
    real(8) :: rdist, r1, r2, rr
    integer :: id, jd

    if (present(rscale)) then
      rr = rscale
    else
      rr = 0.5d0
    end if

    id = getAtomicNumber(l1)
    jd = getAtomicNumber(l2)

    if (id<1) then
      r1 = 0.5
    else
      r1 = atom(id) % rcov
    end if

    if (jd<1) then
      r2 = 0.5
    else
      r2 = atom(jd) % rcov
    end if

    rdist = (r1 + r2) * rr

    return
  end function getMaximumOverlap

  subroutine initialisePeriodicTable()
    implicit none

    type(element), pointer, dimension(:) :: a

    a => atom
!                       lab     mass       rcov      rvdw      rion      chg
    a(   0 ) = element("XX  " ,   0.00d0 , 0.00d0 , 1.00d0 , 1.00d0 ,  0.0d0)
    a(   1 ) = element("H   " ,   1.01d0 , 0.37d0 , 1.08d0 ,-0.24d0 ,  0.0d0)
    a(   2 ) = element("He  " ,   4.00d0 , 0.00d0 , 1.00d0 , 0.00d0 ,  0.0d0)
    a(   3 ) = element("Li  " ,   6.94d0 , 0.68d0 , 1.80d0 , 0.90d0 ,  1.0d0)
    a(   4 ) = element("Be  " ,   9.01d0 , 0.35d0 , 0.52d0 , 0.59d0 ,  0.0d0)
    a(   5 ) = element("B   " ,  10.81d0 , 1.05d0 , 1.70d0 , 0.25d0 ,  0.0d0)
    a(   6 ) = element("C   " ,  12.01d0 , 0.77d0 , 1.53d0 , 0.30d0 ,  0.0d0)
    a(   7 ) = element("N   " ,  14.01d0 , 0.75d0 , 1.48d0 , 1.32d0 ,  0.0d0)
    a(   8 ) = element("O   " ,  16.00d0 , 0.73d0 , 1.36d0 , 1.26d0 , -2.0d0)
    a(   9 ) = element("F   " ,  19.00d0 , 0.90d0 , 1.30d0 , 1.19d0 , -1.0d0)
    a(  10 ) = element("Ne  " ,  20.18d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)
    a(  11 ) = element("Na  " ,  22.99d0 , 0.97d0 , 2.30d0 , 1.16d0 ,  1.0d0)
    a(  12 ) = element("Mg  " ,  24.31d0 , 1.10d0 , 1.64d0 , 0.86d0 ,  2.0d0)
    a(  13 ) = element("Al  " ,  26.98d0 , 1.35d0 , 2.05d0 , 0.68d0 ,  3.0d0)
    a(  14 ) = element("Si  " ,  28.09d0 , 1.20d0 , 2.10d0 , 0.40d0 ,  4.0d0)
    a(  15 ) = element("P   " ,  30.97d0 , 1.05d0 , 1.75d0 , 0.31d0 ,  0.0d0)
    a(  16 ) = element("S   " ,  32.07d0 , 1.02d0 , 1.70d0 , 1.70d0 ,  0.0d0)
    a(  17 ) = element("Cl  " ,  35.45d0 , 0.99d0 , 1.65d0 , 1.67d0 , -1.0d0)
    a(  18 ) = element("Ar  " ,  39.95d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)
    a(  19 ) = element("K   " ,  39.10d0 , 1.33d0 , 2.80d0 , 1.52d0 ,  1.0d0)
    a(  20 ) = element("Ca  " ,  40.08d0 , 0.99d0 , 2.75d0 , 1.14d0 ,  2.0d0)
    a(  21 ) = element("Sc  " ,  44.96d0 , 1.44d0 , 2.15d0 , 0.89d0 ,  0.0d0)
    a(  22 ) = element("Ti  " ,  47.88d0 , 1.47d0 , 2.19d0 , 0.74d0 ,  2.0d0)
    a(  23 ) = element("V   " ,  50.94d0 , 1.33d0 , 1.99d0 , 0.68d0 ,  0.0d0)
    a(  24 ) = element("Cr  " ,  52.00d0 , 1.35d0 , 2.01d0 , 0.76d0 ,  0.0d0)
    a(  25 ) = element("Mn  " ,  54.94d0 , 1.35d0 , 2.01d0 , 0.67d0 ,  2.0d0)
    a(  26 ) = element("Fe  " ,  55.85d0 , 1.34d0 , 2.00d0 , 0.69d0 ,  3.0d0)
    a(  27 ) = element("Co  " ,  58.93d0 , 1.33d0 , 1.99d0 , 0.79d0 ,  0.0d0)
    a(  28 ) = element("Ni  " ,  58.69d0 , 1.50d0 , 1.81d0 , 0.83d0 ,  2.0d0)
    a(  29 ) = element("Cu  " ,  63.55d0 , 1.52d0 , 1.54d0 , 0.87d0 ,  0.0d0)
    a(  30 ) = element("Zn  " ,  65.39d0 , 1.45d0 , 2.16d0 , 0.88d0 ,  2.0d0)
    a(  31 ) = element("Ga  " ,  69.72d0 , 1.22d0 , 1.82d0 , 0.76d0 ,  0.0d0)
    a(  32 ) = element("Ge  " ,  72.61d0 , 1.17d0 , 1.75d0 , 0.67d0 ,  0.0d0)
    a(  33 ) = element("As  " ,  74.92d0 , 1.21d0 , 2.20d0 , 0.72d0 ,  0.0d0)
    a(  34 ) = element("Se  " ,  78.96d0 , 1.22d0 , 2.00d0 , 1.84d0 ,  0.0d0)
    a(  35 ) = element("Br  " ,  79.90d0 , 1.21d0 , 1.80d0 , 1.82d0 , -1.0d0)
    a(  36 ) = element("Kr  " ,  83.80d0 , 1.89d0 , 2.82d0 , 0.00d0 ,  0.0d0)
    a(  37 ) = element("Rb  " ,  85.47d0 , 1.47d0 , 2.19d0 , 1.66d0 ,  1.0d0)
    a(  38 ) = element("Sr  " ,  87.62d0 , 1.12d0 , 1.67d0 , 1.32d0 ,  2.0d0)
    a(  39 ) = element("Y   " ,  88.91d0 , 0.00d0 , 2.66d0 , 1.04d0 ,  3.0d0)
    a(  40 ) = element("Zr  " ,  91.22d0 , 1.56d0 , 2.33d0 , 0.86d0 ,  4.0d0)
    a(  41 ) = element("Nb  " ,  92.91d0 , 1.48d0 , 2.21d0 , 0.78d0 ,  0.0d0)
    a(  42 ) = element("Mo  " ,  95.94d0 , 1.47d0 , 2.19d0 , 0.73d0 ,  2.0d0)
    a(  43 ) = element("Tc  " ,  98.00d0 , 1.35d0 , 2.01d0 , 0.70d0 ,  0.0d0)
    a(  44 ) = element("Ru  " , 101.07d0 , 1.40d0 , 2.09d0 , 0.70d0 ,  0.0d0)
    a(  45 ) = element("Rh  " , 102.91d0 , 1.45d0 , 2.16d0 , 0.69d0 ,  0.0d0)
    a(  46 ) = element("Pd  " , 106.42d0 , 1.50d0 , 2.24d0 , 0.76d0 ,  0.0d0)
    a(  47 ) = element("Ag  " , 107.87d0 , 1.59d0 , 2.37d0 , 1.29d0 ,  0.0d0)
    a(  48 ) = element("Cd  " , 112.41d0 , 1.69d0 , 2.52d0 , 1.09d0 ,  0.0d0)
    a(  49 ) = element("In  " , 114.82d0 , 1.63d0 , 2.43d0 , 0.94d0 ,  0.0d0)
    a(  50 ) = element("Sn  " , 118.71d0 , 1.46d0 , 2.18d0 , 0.83d0 ,  0.0d0)
    a(  51 ) = element("Sb  " , 121.75d0 , 1.46d0 , 2.17d0 , 0.74d0 ,  0.0d0)
    a(  52 ) = element("Te  " , 127.60d0 , 1.47d0 , 2.20d0 , 2.07d0 ,  0.0d0)
    a(  53 ) = element("I   " , 126.91d0 , 2.00d0 , 2.05d0 , 2.06d0 , -1.0d0)
    a(  54 ) = element("Xe  " , 131.29d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)
    a(  55 ) = element("Cs  " , 132.91d0 , 1.67d0 , 2.49d0 , 1.81d0 ,  0.0d0)
    a(  56 ) = element("Ba  " , 137.33d0 , 1.34d0 , 2.00d0 , 1.49d0 ,  2.0d0)
    a(  57 ) = element("La  " , 138.91d0 , 1.87d0 , 2.79d0 , 1.17d0 ,  0.0d0)
    a(  58 ) = element("Ce  " , 140.12d0 , 1.83d0 , 2.73d0 , 1.01d0 ,  0.0d0)
    a(  59 ) = element("Pr  " , 140.91d0 , 1.82d0 , 2.72d0 , 1.13d0 ,  0.0d0)
    a(  60 ) = element("Nd  " , 144.24d0 , 1.81d0 , 2.70d0 , 1.12d0 ,  0.0d0)
    a(  61 ) = element("Pm  " , 145.00d0 , 1.80d0 , 2.69d0 , 1.11d0 ,  0.0d0)
    a(  62 ) = element("Sm  " , 150.36d0 , 1.80d0 , 2.69d0 , 1.10d0 ,  0.0d0)
    a(  63 ) = element("Eu  " , 151.97d0 , 1.99d0 , 2.97d0 , 1.09d0 ,  0.0d0)
    a(  64 ) = element("Gd  " , 157.25d0 , 1.79d0 , 2.67d0 , 1.08d0 ,  0.0d0)
    a(  65 ) = element("Tb  " , 158.93d0 , 1.76d0 , 2.63d0 , 1.06d0 ,  0.0d0)
    a(  66 ) = element("Dy  " , 162.50d0 , 1.75d0 , 2.61d0 , 1.05d0 ,  0.0d0)
    a(  67 ) = element("Ho  " , 164.93d0 , 1.74d0 , 2.60d0 , 1.04d0 ,  0.0d0)
    a(  68 ) = element("Er  " , 167.26d0 , 1.73d0 , 2.58d0 , 1.03d0 ,  0.0d0)
    a(  69 ) = element("Tm  " , 168.93d0 , 1.72d0 , 2.57d0 , 1.02d0 ,  0.0d0)
    a(  70 ) = element("Yb  " , 173.04d0 , 1.94d0 , 2.90d0 , 1.01d0 ,  0.0d0)
    a(  71 ) = element("Lu  " , 174.97d0 , 1.72d0 , 2.57d0 , 1.00d0 ,  0.0d0)
    a(  72 ) = element("Hf  " , 178.49d0 , 1.57d0 , 2.34d0 , 0.85d0 ,  0.0d0)
    a(  73 ) = element("Ta  " , 180.95d0 , 1.43d0 , 2.13d0 , 0.78d0 ,  0.0d0)
    a(  74 ) = element("W   " , 183.85d0 , 1.37d0 , 2.04d0 , 0.65d0 ,  0.0d0)
    a(  75 ) = element("Re  " , 186.21d0 , 1.35d0 , 2.01d0 , 0.67d0 ,  0.0d0)
    a(  76 ) = element("Os  " , 190.20d0 , 1.37d0 , 2.04d0 , 0.53d0 ,  0.0d0)
    a(  77 ) = element("Ir  " , 192.22d0 , 1.32d0 , 1.97d0 , 0.71d0 ,  0.0d0)
    a(  78 ) = element("Pt  " , 195.08d0 , 1.50d0 , 1.97d0 , 0.71d0 ,  0.0d0)
    a(  79 ) = element("Au  " , 196.97d0 , 1.50d0 , 1.85d0 , 0.71d0 ,  0.0d0)
    a(  80 ) = element("Hg  " , 200.59d0 , 1.70d0 , 1.90d0 , 1.16d0 ,  0.0d0)
    a(  81 ) = element("Tl  " , 204.38d0 , 1.55d0 , 2.31d0 , 1.02d0 ,  0.0d0)
    a(  82 ) = element("Pb  " , 207.20d0 , 1.54d0 , 2.30d0 , 1.33d0 ,  2.0d0)
    a(  83 ) = element("Bi  " , 208.98d0 , 1.54d0 , 2.30d0 , 1.17d0 ,  1.0d0)
    a(  84 ) = element("Po  " , 209.00d0 , 1.68d0 , 2.51d0 , 0.81d0 ,  0.0d0)
    a(  85 ) = element("At  " , 210.00d0 , 0.00d0 , 0.00d0 , 0.76d0 ,  0.0d0)
    a(  86 ) = element("Rn  " , 222.00d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)
    a(  87 ) = element("Fr  " , 223.00d0 , 0.00d0 , 0.00d0 , 1.94d0 ,  0.0d0)
    a(  88 ) = element("Ra  " , 226.03d0 , 1.90d0 , 2.84d0 , 1.62d0 ,  0.0d0)
    a(  89 ) = element("Ac  " , 227.03d0 , 1.88d0 , 2.81d0 , 1.26d0 ,  0.0d0)
    a(  90 ) = element("Th  " , 232.04d0 , 1.79d0 , 2.67d0 , 0.00d0 ,  0.0d0)
    a(  91 ) = element("Pa  " , 231.04d0 , 1.61d0 , 2.40d0 , 0.92d0 ,  0.0d0)
    a(  92 ) = element("U   " , 238.03d0 , 1.58d0 , 2.36d0 , 1.03d0 ,  0.0d0)
    ! a(  93 ) = element("Np  " , 237.05d0 , 1.55d0 , 2.31d0 , 0.85d0 ,  0.0d0)
    ! a(  94 ) = element("Pu  " , 244.00d0 , 1.53d0 , 2.28d0 , 0.85d0 ,  0.0d0)
    ! a(  95 ) = element("Am  " , 243.00d0 , 1.51d0 , 2.25d0 , 1.35d0 ,  0.0d0)
    ! a(  96 ) = element("Cm  " , 247.00d0 , 0.00d0 , 0.00d0 , 0.99d0 ,  0.0d0)
    ! a(  97 ) = element("Bk  " , 247.00d0 , 0.00d0 , 0.00d0 , 0.97d0 ,  0.0d0)
    ! a(  98 ) = element("Cf  " , 251.00d0 , 0.00d0 , 0.00d0 , 0.96d0 ,  0.0d0)
    ! a(  99 ) = element("Es  " , 252.00d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)
    ! a( 100 ) = element("Fm  " , 257.00d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)
    ! a( 101 ) = element("Md  " , 258.00d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)
    ! a( 102 ) = element("No  " , 259.00d0 , 0.00d0 , 1.24d0 , 1.24d0 ,  0.0d0)
    ! a( 103 ) = element("Lr  " , 260.00d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)
    ! a( 104 ) = element("Rf  " , 261.00d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)
    ! a( 105 ) = element("Ha  " , 260.00d0 , 0.00d0 , 0.00d0 , 0.00d0 ,  0.0d0)

! Deuterium
    a( 106 ) = element("D   " ,   2.01d0 , 0.37d0 , 1.08d0 ,-0.24d0 ,  0.0d0)
! Tritium
    a( 107 ) = element("T   " ,   3.01d0 , 0.37d0 , 1.08d0 ,-0.24d0 ,  0.0d0)

    return
  end subroutine initialisePeriodicTable

end module moduleElements 
