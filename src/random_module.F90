! Mersenne Twister random number generator, taken from 
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
! M. Matsumoto and T. Nishimura, 
! "Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator"
! ACM Trans. on Modeling and Computer Simulation Vol. 8, No. 1, January pp.3-30 (1998)
module moduleRandomNumbers
use, intrinsic :: iso_fortran_env, only: real64
implicit none
!
! period parameters
!
integer, parameter :: n = 624, n1 = n+1, m = 397, mata = -1727483681
!
! constant vector a
!
integer, parameter :: umask = -2147483646
!
! most significant w-r bits
!
integer, parameter :: lmask =  2147483647
!
! least significant r bits
! tempering parameters
!
integer, parameter :: tmaskb= -1658038656, tmaskc= -272236544
!
! the array for the state vector
! mti==n+1 means mt[n] is not initialized
!
integer, save      :: mt(0:n-1), mti = n1
!
! Counter to recover the past sequence of numbers
! I added this to have reproducibility of a restart run
! Note howeverm that it regenerates the same sequence of random number starting
! from the beginning... hence it may take a while if the previous run was quite long :)
! Maybe an option to toggle this off should be implemented...
!
integer, save      :: ncount_rand=0

private
public :: init_genrand, grnd, gauss_rand, ncount_rand

contains

function gauss_rand() result(gr)
  implicit none
  real(real64) :: r1, r2, rr, fac, gr
  real(real64), save :: gset
  integer, save :: iset=0

  if (iset==0) then 
100 r1=2.0_real64*grnd()-1.0_real64
    r2=2.0_real64*grnd()-1.0_real64
    rr=r1*r1+r2*r2
    if (rr.ge.1.0_real64.or.rr.eq.0.0_real64) then 
      goto 100
    endif
    fac=(sqrt(-2.0_real64*log(rr)/rr))
    gr=r1*fac
    gset=r2*fac
    iset=1
  else
    gr=gset
    iset=0
  endif

  return
end function gauss_rand

subroutine init_genrand(seed)
! this initialization is based upon the multiplier given on p.106 of the
! 3rd edition of knuth, the art of computer programming vol. 2.

! this version assumes that integer overflow does not cause a crash.

integer, intent(in)  :: seed

integer  :: latest

mt(0) = seed
latest = seed
do mti = 1, n-1
  latest = ieor( latest, ishft( latest, -30 ) )
  latest = latest * 1812433253 + mti
  mt(mti) = latest
end do

return
end subroutine init_genrand
!***********************************************************************

function grnd() result(fn_val)
double precision :: fn_val

integer, save :: mag01(0:1) = (/ 0, mata /)
!                        mag01(x) = x * mata for x=0,1
integer       :: kk, y

! these statement functions have been replaced with separate functions
! tshftu(y) = ishft(y,-11)
! tshfts(y) = ishft(y,7)
! tshftt(y) = ishft(y,15)
! tshftl(y) = ishft(y,-18)

! counter to recover the past sequence of numbers
ncount_rand=ncount_rand+1

if (mti >= n) then
!                       generate n words at one time
  if (mti == n+1) then
!                            if init_genrand() has not been called,
    call init_genrand(4357)
!                              a default initial seed is used
  end if
  
  do  kk = 0, n-m-1
    y = ior(iand(mt(kk),umask), iand(mt(kk+1),lmask))
    mt(kk) = ieor(ieor(mt(kk+m), ishft(y,-1)),mag01(iand(y,1)))
  end do
  do  kk = n-m, n-2
    y = ior(iand(mt(kk),umask), iand(mt(kk+1),lmask))
  mt(kk) = ieor(ieor(mt(kk+(m-n)), ishft(y,-1)),mag01(iand(y,1)))
  end do
  y = ior(iand(mt(n-1),umask), iand(mt(0),lmask))
  mt(n-1) = ieor(ieor(mt(m-1), ishft(y,-1)),mag01(iand(y,1)))
  mti = 0
end if

y = mt(mti)
mti = mti + 1
y = ieor(y, tshftu(y))
y = ieor(y, iand(tshfts(y),tmaskb))
y = ieor(y, iand(tshftt(y),tmaskc))
y = ieor(y, tshftl(y))

if (y < 0) then
  fn_val = (dble(y) + 2.0_real64**32) / (2.0_real64**32 - 1.0_real64)
else
  fn_val = dble(y) / (2.0_real64**32 - 1.0_real64)
end if

return
end function grnd


function tshftu(y) result(fn_val)
integer, intent(in) :: y
integer             :: fn_val

fn_val = ishft(y,-11)
return
end function tshftu


function tshfts(y) result(fn_val)
integer, intent(in) :: y
integer             :: fn_val

fn_val = ishft(y,7)
return
end function tshfts


function tshftt(y) result(fn_val)
integer, intent(in) :: y
integer             :: fn_val

fn_val = ishft(y,15)
return
end function tshftt


function tshftl(y) result(fn_val)
integer, intent(in) :: y
integer             :: fn_val

fn_val = ishft(y,-18)
return
end function tshftl

end module moduleRandomNumbers
