!disclaimer
module moduleElements 
  use moduleVariables
  
  implicit none

  integer, parameter :: nelement=107
  type element
    character(len=cp) :: lab
    real(real64) :: mass
    real(real64) :: rcov
    real(real64) :: rvdw
    real(real64) :: rion
    real(real64) :: chg
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
      call lowercase(label,element)
      element(1:1) = uppercase(label(1:1))
      ! if (i==1) then
      !   element(i:i) = uppercase(label(i:i))
      ! else
      !   element(i:i) = label(i:i)
      ! end if
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
! for exaple HA or HB can be identified as hydrogen atoms
    if(id<=0 .and. ln>1) then
      ln = ln-1
      goto 100
    end if

    return
  end function getAtomicNumber

  function getElement(label) result(elem)
    implicit none
    character(*), intent(in) :: label
    character(2) :: elem
    integer :: id

    character(cp) :: element
    character(len=1), external :: uppercase
    integer :: i, ln

    element = ''
    do i=1,cp
      if (ichar(label(i:i)) >=48 .and. ichar(label(i:i))<=57 ) exit
      call lowercase(label,element)
      element(1:1) = uppercase(label(1:1))
      ! if (i==1) then
      !   element(i:i) = uppercase(label(i:i))
      ! else
      !   element(i:i) = label(i:i)
      ! end if
    enddo
    ln=len_trim(element)
    id=0
100 continue
    do i=1,nelement
      if ( len_trim(atom(i) % lab) /= ln ) cycle
      if(element(1:ln)==atom(i) % lab(1:ln)) then
        elem = atom(i) % lab(1:2)
        return
      end if
    enddo

! This should ensure atoms to have letters after the species name
! for exaple CA or CB can be identified as hydrogen atoms
    if(id<=0 .and. ln>1) then
      ln = ln-1
      goto 100
    end if

    return
  end function getElement

  function getElementMass(label) result(mass)
    implicit none
    character(cp), intent(in) :: label
    real(real64) :: mass

    character(cp) :: element
    integer :: id
    integer :: i, ln

    mass = 0.0_real64
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
    real(real64), intent(in) :: mass
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
    real(real64), intent(in) :: rscale
    real(real64) :: rbond, r1, r2
    integer :: id, jd, i

    integer, dimension(5) :: h_allowed

    rbond = 0.0_real64

    do i=1,sbond_n
      if ( (l1==sbond_l(1,i) .and. l2==sbond_l(2,i)) .or. (l1==sbond_l(2,i) .and. l2==sbond_l(1,i)) ) then
        rbond = sbond_d(i)
        return
      end if
    enddo

! ions
    if (allocated(ions_list)) then
      if (any(ions_list==l1) .or. any(ions_list==l2)) then
        rbond = -1.0_real64
        return
      end if
    end if

! Atoms that can bond to an H
! C, N, O, S, Cu
!    h_allowed = (/6, 7, 8, 9, 13, 14, 15, 16, 17, 35, 53/)
    h_allowed = (/6, 7, 8, 16, 29/)

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
      if (all(h_allowed/=jd)) rbond = 0.0_real64
    else if (jd==1) then
      if (all(h_allowed/=id)) rbond = 0.0_real64
! Exclude O-O bonds
    else if (id==8 .and. jd==8) then
      rbond = 0.0_real64
! O - Cl
    else if ( (id==8 .and. jd==17) .or. (id==17 .and. jd==8)  )then
      rbond = 0.0_real64
! N - Cl
    else if ( (id==7 .and. jd==17) .or. (id==17 .and. jd==7)  )then
      rbond = 0.0_real64
    end if

    return
  end function getMaximumBondLength

  function getMaximumOverlap(l1,l2,rscale) result (rdist)
    implicit none
    character(len=cp), intent(in) :: l1, l2
    real(real64), optional, intent(in) :: rscale
    real(real64) :: rdist, r1, r2, rr
    integer :: id, jd

    if (present(rscale)) then
      rr = rscale
    else
      rr = 0.5_real64
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
    a(   0 ) = element("XX  " ,   0.00_real64 , 0.00_real64 , 1.00_real64 , 1.00_real64 ,  0.0_real64)
    a(   1 ) = element("H   " ,   1.01_real64 , 0.37_real64 , 1.08_real64 ,-0.24_real64 ,  0.0_real64)
    a(   2 ) = element("He  " ,   4.00_real64 , 0.00_real64 , 1.00_real64 , 0.00_real64 ,  0.0_real64)
    a(   3 ) = element("Li  " ,   6.94_real64 , 0.68_real64 , 1.80_real64 , 0.90_real64 ,  1.0_real64)
    a(   4 ) = element("Be  " ,   9.01_real64 , 0.35_real64 , 0.52_real64 , 0.59_real64 ,  0.0_real64)
    a(   5 ) = element("B   " ,  10.81_real64 , 1.05_real64 , 1.70_real64 , 0.25_real64 ,  0.0_real64)
    a(   6 ) = element("C   " ,  12.01_real64 , 0.77_real64 , 1.53_real64 , 0.30_real64 ,  0.0_real64)
    a(   7 ) = element("N   " ,  14.01_real64 , 0.75_real64 , 1.48_real64 , 1.32_real64 ,  0.0_real64)
    a(   8 ) = element("O   " ,  16.00_real64 , 0.73_real64 , 1.36_real64 , 1.26_real64 , -2.0_real64)
    a(   9 ) = element("F   " ,  19.00_real64 , 0.90_real64 , 1.30_real64 , 1.19_real64 , -1.0_real64)
    a(  10 ) = element("Ne  " ,  20.18_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)
    a(  11 ) = element("Na  " ,  22.99_real64 , 0.97_real64 , 2.30_real64 , 1.16_real64 ,  1.0_real64)
    a(  12 ) = element("Mg  " ,  24.31_real64 , 0.10_real64 , 1.64_real64 , 0.86_real64 ,  2.0_real64)
!    a(  12 ) = element("Mg  " ,  24.31_real64 , 1.10_real64 , 1.64_real64 , 0.86_real64 ,  2.0_real64)
    a(  13 ) = element("Al  " ,  26.98_real64 , 1.35_real64 , 2.05_real64 , 0.68_real64 ,  3.0_real64)
    a(  14 ) = element("Si  " ,  28.09_real64 , 1.20_real64 , 2.10_real64 , 0.40_real64 ,  4.0_real64)
    a(  15 ) = element("P   " ,  30.97_real64 , 1.05_real64 , 1.75_real64 , 0.31_real64 ,  0.0_real64)
    a(  16 ) = element("S   " ,  32.07_real64 , 1.02_real64 , 1.70_real64 , 1.70_real64 ,  0.0_real64)
    a(  17 ) = element("Cl  " ,  35.45_real64 , 0.99_real64 , 1.65_real64 , 1.67_real64 , -1.0_real64)
    a(  18 ) = element("Ar  " ,  39.95_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)
    a(  19 ) = element("K   " ,  39.10_real64 , 1.33_real64 , 2.80_real64 , 1.52_real64 ,  1.0_real64)
    a(  20 ) = element("Ca  " ,  40.08_real64 , 0.99_real64 , 2.75_real64 , 1.14_real64 ,  2.0_real64)
    a(  21 ) = element("Sc  " ,  44.96_real64 , 1.44_real64 , 2.15_real64 , 0.89_real64 ,  0.0_real64)
    a(  22 ) = element("Ti  " ,  47.88_real64 , 1.47_real64 , 2.19_real64 , 0.74_real64 ,  2.0_real64)
    a(  23 ) = element("V   " ,  50.94_real64 , 1.33_real64 , 1.99_real64 , 0.68_real64 ,  0.0_real64)
    a(  24 ) = element("Cr  " ,  52.00_real64 , 1.35_real64 , 2.01_real64 , 0.76_real64 ,  0.0_real64)
    a(  25 ) = element("Mn  " ,  54.94_real64 , 1.35_real64 , 2.01_real64 , 0.67_real64 ,  2.0_real64)
    a(  26 ) = element("Fe  " ,  55.85_real64 , 1.34_real64 , 2.00_real64 , 0.69_real64 ,  3.0_real64)
    a(  27 ) = element("Co  " ,  58.93_real64 , 0.10_real64 , 1.99_real64 , 0.79_real64 ,  0.0_real64)
    a(  28 ) = element("Ni  " ,  58.69_real64 , 1.50_real64 , 1.81_real64 , 0.83_real64 ,  2.0_real64)
    a(  29 ) = element("Cu  " ,  63.55_real64 , 1.52_real64 , 1.54_real64 , 0.87_real64 ,  1.5_real64)
    a(  30 ) = element("Zn  " ,  65.39_real64 , 1.45_real64 , 2.16_real64 , 0.88_real64 ,  2.0_real64)
    a(  31 ) = element("Ga  " ,  69.72_real64 , 1.22_real64 , 1.82_real64 , 0.76_real64 ,  0.0_real64)
    a(  32 ) = element("Ge  " ,  72.61_real64 , 1.17_real64 , 1.75_real64 , 0.67_real64 ,  0.0_real64)
    a(  33 ) = element("As  " ,  74.92_real64 , 1.21_real64 , 2.20_real64 , 0.72_real64 ,  0.0_real64)
    a(  34 ) = element("Se  " ,  78.96_real64 , 1.22_real64 , 2.00_real64 , 1.84_real64 ,  0.0_real64)
    a(  35 ) = element("Br  " ,  79.90_real64 , 1.21_real64 , 1.80_real64 , 1.82_real64 , -1.0_real64)
    a(  36 ) = element("Kr  " ,  83.80_real64 , 1.89_real64 , 2.82_real64 , 0.00_real64 ,  0.0_real64)
    a(  37 ) = element("Rb  " ,  85.47_real64 , 0.17_real64 , 2.19_real64 , 1.66_real64 ,  1.0_real64)
!    a(  37 ) = element("Rb  " ,  85.47_real64 , 1.47_real64 , 2.19_real64 , 1.66_real64 ,  1.0_real64)
    a(  38 ) = element("Sr  " ,  87.62_real64 , 1.12_real64 , 1.67_real64 , 1.32_real64 ,  2.0_real64)
    a(  39 ) = element("Y   " ,  88.91_real64 , 0.00_real64 , 2.66_real64 , 1.04_real64 ,  3.0_real64)
    a(  40 ) = element("Zr  " ,  91.22_real64 , 1.56_real64 , 2.33_real64 , 0.86_real64 ,  4.0_real64)
    a(  41 ) = element("Nb  " ,  92.91_real64 , 1.48_real64 , 2.21_real64 , 0.78_real64 ,  0.0_real64)
    a(  42 ) = element("Mo  " ,  95.94_real64 , 1.47_real64 , 2.19_real64 , 0.73_real64 ,  2.0_real64)
    a(  43 ) = element("Tc  " ,  98.00_real64 , 1.35_real64 , 2.01_real64 , 0.70_real64 ,  0.0_real64)
    a(  44 ) = element("Ru  " , 101.07_real64 , 1.40_real64 , 2.09_real64 , 0.70_real64 ,  0.0_real64)
    a(  45 ) = element("Rh  " , 102.91_real64 , 1.45_real64 , 2.16_real64 , 0.69_real64 ,  0.0_real64)
    a(  46 ) = element("Pd  " , 106.42_real64 , 1.50_real64 , 2.24_real64 , 0.76_real64 ,  0.0_real64)
    a(  47 ) = element("Ag  " , 107.87_real64 , 1.59_real64 , 2.37_real64 , 1.29_real64 ,  0.0_real64)
    a(  48 ) = element("Cd  " , 112.41_real64 , 1.69_real64 , 2.52_real64 , 1.09_real64 ,  0.0_real64)
    a(  49 ) = element("In  " , 114.82_real64 , 1.63_real64 , 2.43_real64 , 0.94_real64 ,  0.0_real64)
    a(  50 ) = element("Sn  " , 118.71_real64 , 1.46_real64 , 2.18_real64 , 0.83_real64 ,  0.0_real64)
    a(  51 ) = element("Sb  " , 121.75_real64 , 1.46_real64 , 2.17_real64 , 0.74_real64 ,  0.0_real64)
    a(  52 ) = element("Te  " , 127.60_real64 , 1.47_real64 , 2.20_real64 , 2.07_real64 ,  0.0_real64)
    a(  53 ) = element("I   " , 126.91_real64 , 2.00_real64 , 2.05_real64 , 2.06_real64 , -1.0_real64)
    a(  54 ) = element("Xe  " , 131.29_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)
    a(  55 ) = element("Cs  " , 132.91_real64 , 1.67_real64 , 2.49_real64 , 1.81_real64 ,  0.0_real64)
    a(  56 ) = element("Ba  " , 137.33_real64 , 1.34_real64 , 2.00_real64 , 1.49_real64 ,  2.0_real64)
    a(  57 ) = element("La  " , 138.91_real64 , 1.87_real64 , 2.79_real64 , 1.17_real64 ,  0.0_real64)
    a(  58 ) = element("Ce  " , 140.12_real64 , 1.83_real64 , 2.73_real64 , 1.01_real64 ,  0.0_real64)
    a(  59 ) = element("Pr  " , 140.91_real64 , 1.82_real64 , 2.72_real64 , 1.13_real64 ,  0.0_real64)
    a(  60 ) = element("Nd  " , 144.24_real64 , 1.81_real64 , 2.70_real64 , 1.12_real64 ,  0.0_real64)
    a(  61 ) = element("Pm  " , 145.00_real64 , 1.80_real64 , 2.69_real64 , 1.11_real64 ,  0.0_real64)
    a(  62 ) = element("Sm  " , 150.36_real64 , 1.80_real64 , 2.69_real64 , 1.10_real64 ,  0.0_real64)
    a(  63 ) = element("Eu  " , 151.97_real64 , 1.99_real64 , 2.97_real64 , 1.09_real64 ,  0.0_real64)
    a(  64 ) = element("Gd  " , 157.25_real64 , 1.79_real64 , 2.67_real64 , 1.08_real64 ,  0.0_real64)
    a(  65 ) = element("Tb  " , 158.93_real64 , 1.76_real64 , 2.63_real64 , 1.06_real64 ,  0.0_real64)
    a(  66 ) = element("Dy  " , 162.50_real64 , 1.75_real64 , 2.61_real64 , 1.05_real64 ,  0.0_real64)
    a(  67 ) = element("Ho  " , 164.93_real64 , 1.74_real64 , 2.60_real64 , 1.04_real64 ,  0.0_real64)
    a(  68 ) = element("Er  " , 167.26_real64 , 1.73_real64 , 2.58_real64 , 1.03_real64 ,  0.0_real64)
    a(  69 ) = element("Tm  " , 168.93_real64 , 1.72_real64 , 2.57_real64 , 1.02_real64 ,  0.0_real64)
    a(  70 ) = element("Yb  " , 173.04_real64 , 1.94_real64 , 2.90_real64 , 1.01_real64 ,  0.0_real64)
    ! a(  71 ) = element("Lu  " , 174.97_real64 , 1.72_real64 , 2.57_real64 , 1.00_real64 ,  0.0_real64)
    a(  71 ) = element("Lu  " , 174.97_real64 , 0.72_real64 , 2.57_real64 , 1.00_real64 ,  0.0_real64)
    a(  72 ) = element("Hf  " , 178.49_real64 , 1.57_real64 , 2.34_real64 , 0.85_real64 ,  0.0_real64)
    a(  73 ) = element("Ta  " , 180.95_real64 , 1.43_real64 , 2.13_real64 , 0.78_real64 ,  0.0_real64)
    a(  74 ) = element("W   " , 183.85_real64 , 0.10_real64 , 2.04_real64 , 0.65_real64 ,  0.0_real64)
    a(  75 ) = element("Re  " , 186.21_real64 , 1.35_real64 , 2.01_real64 , 0.67_real64 ,  0.0_real64)
    a(  76 ) = element("Os  " , 190.20_real64 , 1.37_real64 , 2.04_real64 , 0.53_real64 ,  0.0_real64)
    a(  77 ) = element("Ir  " , 192.22_real64 , 0.32_real64 , 1.97_real64 , 0.71_real64 ,  0.0_real64)
    ! a(  77 ) = element("Ir  " , 192.22_real64 , 1.32_real64 , 1.97_real64 , 0.71_real64 ,  0.0_real64)
    a(  78 ) = element("Pt  " , 195.08_real64 , 1.50_real64 , 1.97_real64 , 0.71_real64 ,  0.0_real64)
    a(  79 ) = element("Au  " , 196.97_real64 , 1.50_real64 , 1.85_real64 , 0.71_real64 ,  0.0_real64)
    a(  80 ) = element("Hg  " , 200.59_real64 , 1.70_real64 , 1.90_real64 , 1.16_real64 ,  0.0_real64)
    a(  81 ) = element("Tl  " , 204.38_real64 , 1.55_real64 , 2.31_real64 , 1.02_real64 ,  0.0_real64)
    a(  82 ) = element("Pb  " , 207.20_real64 , 1.54_real64 , 2.30_real64 , 1.33_real64 ,  2.0_real64)
    a(  83 ) = element("Bi  " , 208.98_real64 , 1.54_real64 , 2.30_real64 , 1.17_real64 ,  1.0_real64)
    a(  84 ) = element("Po  " , 209.00_real64 , 1.68_real64 , 2.51_real64 , 0.81_real64 ,  0.0_real64)
    a(  85 ) = element("At  " , 210.00_real64 , 0.00_real64 , 0.00_real64 , 0.76_real64 ,  0.0_real64)
    a(  86 ) = element("Rn  " , 222.00_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)
    a(  87 ) = element("Fr  " , 223.00_real64 , 0.00_real64 , 0.00_real64 , 1.94_real64 ,  0.0_real64)
    a(  88 ) = element("Ra  " , 226.03_real64 , 1.90_real64 , 2.84_real64 , 1.62_real64 ,  0.0_real64)
    a(  89 ) = element("Ac  " , 227.03_real64 , 1.88_real64 , 2.81_real64 , 1.26_real64 ,  0.0_real64)
    a(  90 ) = element("Th  " , 232.04_real64 , 1.79_real64 , 2.67_real64 , 0.00_real64 ,  0.0_real64)
    a(  91 ) = element("Pa  " , 231.04_real64 , 1.61_real64 , 2.40_real64 , 0.92_real64 ,  0.0_real64)
    a(  92 ) = element("U   " , 238.03_real64 , 1.58_real64 , 2.36_real64 , 1.03_real64 ,  0.0_real64)
    a(  93 ) = element("Ep  " ,   1.01_real64 , 0.37_real64 , 1.08_real64 ,-0.24_real64 ,  0.0_real64)
    ! a(  93 ) = element("Np  " , 237.05_real64 , 1.55_real64 , 2.31_real64 , 0.85_real64 ,  0.0_real64)
    ! a(  94 ) = element("Pu  " , 244.00_real64 , 1.53_real64 , 2.28_real64 , 0.85_real64 ,  0.0_real64)
    ! a(  95 ) = element("Am  " , 243.00_real64 , 1.51_real64 , 2.25_real64 , 1.35_real64 ,  0.0_real64)
    ! a(  96 ) = element("Cm  " , 247.00_real64 , 0.00_real64 , 0.00_real64 , 0.99_real64 ,  0.0_real64)
    ! a(  97 ) = element("Bk  " , 247.00_real64 , 0.00_real64 , 0.00_real64 , 0.97_real64 ,  0.0_real64)
    ! a(  98 ) = element("Cf  " , 251.00_real64 , 0.00_real64 , 0.00_real64 , 0.96_real64 ,  0.0_real64)
    ! a(  99 ) = element("Es  " , 252.00_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)
    ! a( 100 ) = element("Fm  " , 257.00_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)
    ! a( 101 ) = element("Md  " , 258.00_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)
    ! a( 102 ) = element("No  " , 259.00_real64 , 0.00_real64 , 1.24_real64 , 1.24_real64 ,  0.0_real64)
    ! a( 103 ) = element("Lr  " , 260.00_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)
    ! a( 104 ) = element("Rf  " , 261.00_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)
    ! a( 105 ) = element("Ha  " , 260.00_real64 , 0.00_real64 , 0.00_real64 , 0.00_real64 ,  0.0_real64)

! Deuterium
    a( 106 ) = element("D   " ,   2.01_real64 , 0.37_real64 , 1.08_real64 ,-0.24_real64 ,  0.0_real64)
! Tritium
    a( 107 ) = element("T   " ,   3.01_real64 , 0.37_real64 , 1.08_real64 ,-0.24_real64 ,  0.0_real64)

    return
  end subroutine initialisePeriodicTable

end module moduleElements 
