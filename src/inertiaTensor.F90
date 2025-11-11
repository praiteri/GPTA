!disclaimer
subroutine computeInertiaTensor(n, eigenValues, nat, lSelect)
  use moduleVariables, only: real64
  use moduleSystem
  use moduleElements
  implicit none
  integer, intent(inout) :: n
  real(real64), dimension(*), intent(out) :: eigenValues
  integer, intent(in) :: nat
  integer, dimension(:), intent(in) :: lSelect

  real(real64) :: inertiaTensor(3,3)
  real(real64) :: eigenVectors(3,3)
  real(real64) :: totmass, rmass, xcom(3)
  real(real64) :: dx(3)
  integer :: i, idx

  ! set the numner of properties computed
  if (n<0) then
    n = 3
    return
  end if

  xcom = 0.0_real64
  totmass = 0.0_real64
  do idx=1,nat
    i = lSelect(idx)
    rmass = getElementMass(frame % lab(i))
    totmass = totmass + rmass
    xcom = xcom + rmass * frame % pos(:,i) 
  enddo
  xcom = xcom / totmass

  inertiaTensor = 0.0_real64
  do idx=1,nat
    i = lSelect(idx)
    rmass = getElementMass(frame % lab(i))
    dx(1:3) = frame % pos(1:3,i) - xcom(1:3)
    inertiaTensor(1,1) = inertiaTensor(1,1) + rmass * (dx(2)**2 + dx(3)**2)
    inertiaTensor(2,2) = inertiaTensor(2,2) + rmass * (dx(1)**2 + dx(3)**2)
    inertiaTensor(3,3) = inertiaTensor(3,3) + rmass * (dx(1)**2 + dx(2)**2)
    inertiaTensor(1,2) = inertiaTensor(1,2) - rmass * (dx(1)    * dx(2))
    inertiaTensor(1,3) = inertiaTensor(1,3) - rmass * (dx(1)    * dx(3))
    inertiaTensor(2,3) = inertiaTensor(2,3) - rmass * (dx(2)    * dx(3))
  enddo

  inertiaTensor(1,1) = inertiaTensor(1,1) / totmass
  inertiaTensor(2,2) = inertiaTensor(2,2) / totmass
  inertiaTensor(3,3) = inertiaTensor(3,3) / totmass
  inertiaTensor(1,2) = inertiaTensor(1,2) / totmass
  inertiaTensor(1,3) = inertiaTensor(1,3) / totmass
  inertiaTensor(2,3) = inertiaTensor(2,3) / totmass
  inertiaTensor(2,1) = inertiaTensor(1,2)
  inertiaTensor(3,1) = inertiaTensor(1,3)
  inertiaTensor(3,2) = inertiaTensor(2,3)

  ! inertiaTensor(1:9) = reshape(inertiaTensor, [9])

  call jacobi(3, inertiaTensor, eigenvalues, eigenVectors)
  call eigsrt(3, eigenvalues, eigenVectors)

end subroutine computeInertiaTensor

subroutine jacobi(n,a,d,v)
  use moduleVariables, only: real64
  use moduleMessages
  implicit none
  integer, intent(in) :: n
  real(real64), dimension(n), intent(out) :: d
  real(real64), dimension(n,n), intent(inout) :: a
  real(real64), dimension(n,n), intent(out) :: v

! Computes all eigenvalues and eigenvectors of a real symmetric n x n matrix a.
! On output, elements of a above the diagonal are destroyed.
! d is a vector of length n that returns the eigenvalues of a.
! v is an n x n matrix whose columns contain, on output, the normalized eigenvectors of a.
  integer :: i,iip,iq
  real(real64) :: c,g,h,s,sm,t,tau,theta,tresh
  real(real64), dimension(size(d)) :: b,z
  logical, dimension(size(d),size(d)) :: upper_triangle
  integer, parameter :: maxrot=50

  upper_triangle=.false.
  do iip=1,size(d)
    do iq=iip+1,size(d)
      upper_triangle(iip,iq)=.true.
    enddo
  enddo

! Initialize v to the identity matrix. initialize b and d to the diagonal of a.
  v(:,:)=0.0_real64
  do iip=1,size(d)
    v(iip,iip)=1.0_real64
  enddo

  b(:)=get_diag(a(:,:))
  d(:)=b(:)
  z(:)=0.0

! Diagonalisation
  do i=1,maxrot

! The normal return, which relies on quadratic convergence to machine underflow.
! Sum off-diagonal elements.
    sm=sum(abs(a),mask=upper_triangle)
    if (sm == 0.0) return

    tresh=merge(0.2_real64*sm/n**2,0.0_real64, i < 4 )
    do iip=1,n-1
      do iq=iip+1,n
        g=100.0_real64*abs(a(iip,iq))
! After four sweeps, skip the rotation if the off-diagonal element is small.
        if ((i > 4) .and. (abs(d(iip))+g == abs(d(iip))) .and. (abs(d(iq))+g == abs(d(iq)))) then
          a(iip,iq)=0.0
        else if (abs(a(iip,iq)) > tresh) then
          h=d(iq)-d(iip)
        if (abs(h)+g == abs(h)) then
          t=a(iip,iq)/h
        else
          theta=0.5_real64*h/a(iip,iq)
          t=1.0_real64/(abs(theta)+sqrt(1.0_real64+theta**2))
          if (theta < 0.0) t=-t
        end if
        c=1.0_real64/sqrt(1+t**2)
        s=t*c
        if (abs(s) < 1.0e-10_real64) s = 0.0_real64
        tau=s/(1.0_real64+c)
        if (abs(tau) < 1.0e-10_real64) tau = 0.0_real64
        h=t*a(iip,iq)
        z(iip)=z(iip)-h
        z(iq)=z(iq)+h
        d(iip)=d(iip)-h
        d(iq)=d(iq)+h
        a(iip,iq)=0.0
        call jrotate(a(1:iip-1,iip),a(1:iip-1,iq))
        call jrotate(a(iip,iip+1:iq-1),a(iip+1:iq-1,iq))
        call jrotate(a(iip,iq+1:n),a(iq,iq+1:n))
        call jrotate(v(:,iip),v(:,iq))
        end if
      end do
    end do
    b(:)=b(:)+z(:)
    d(:)=b(:)
    z(:)=0.0
  end do
  ! if (i>maxrot) call message(0,0,0,"JACOBI | Matrix diagonalisation didn't converge")
  call message(-1,"JACOBI | Matrix diagonalisation didn't converge")
contains

  subroutine jrotate(a1,a2)
    real(real64), dimension(:), intent(inout) :: a1,a2
    real(real64), dimension(size(a1)) :: wk1
    wk1(:)=a1(:)
    a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
    a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
  end subroutine jrotate

  function get_diag(mat)
    real(real64), dimension(:,:), intent(in) :: mat
    real(real64), dimension(size(mat,1)) :: get_diag
    integer :: j
    do j=1,size(mat,1)
      get_diag(j)=mat(j,j)
    end do
  end function get_diag

end subroutine jacobi

subroutine eigsrt(n,d,v)
  use moduleVariables, only : real64
  implicit none
  integer, intent(in) :: n
  real(real64), dimension(n), intent(inout) :: d
  real(real64), dimension(n,n), intent(inout) :: v
! given the eigenvalues d and eigenvectors v as output from jacobi (11.1) or tqli (11.3),
! this routine sorts the eigenvalues into descending order,
! and rearranges the columns of v correspondingly. the method is straight insertion.
  integer :: i,j
  integer :: ii(1)

  do i=1,n-1
    ii=minloc(d(i:n))
    j=ii(1)+i-1
    if (j /= i) then
      call swap_r(d(i),d(j))
      call swap_v(v(:,i),v(:,j))
    end if
  end do
  return

contains

  subroutine swap_r(a,b)
    implicit none
    real(real64), intent(inout) :: a,b
    real(real64) :: dum
    dum=a
    a=b
    b=dum
    return
  end subroutine swap_r

  subroutine swap_v(a,b)
    implicit none
    real(real64), dimension(:), intent(inout) :: a,b
    real(real64), dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
    return
  end subroutine swap_v

end subroutine eigsrt
