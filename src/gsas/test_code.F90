! Copyright (c) 2021, Paolo Raiteri, Curtin University.
! All rights reserved.
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 3
! of the License, or (at your option) any later version.
!  
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. No claim 
! is made that this program is free from errors and no liability 
! will be accepted for any loss or damage that may result. The 
! user is responsible for checking the validity of their results.
!  
! See the GNU General Public License for more details.
! 
! The GNU GPL can also be found at http://www.gnu.org
!
program test
  implicit none
  CHARACTER*20 :: SPG
  INTEGER*4     :: LAUE,SGINV,SGLATT,SGUNIQ,SGNOPS,IERR,SGNCEN
  REAL*4        :: SGMTRX(48,3,3),SGTRNS(48,3),SGGEN(24)
  REAL*4        :: RT(5,4,25),CEN(3,4)
  INTEGER*4     :: JRT(3,5,24)
  INTEGER*4     :: SGPOL

  integer :: i, j, k, l, nop
  real(8) :: rot(3,4,48), dij(3), sij(3)

!       SPG    Input  20 Characters containing the space group symbol
!       LAUE   Output The Laue Group no. where
!                 1=1BAR, 2=2/M, 3=MMM, 4=4/M, 5=4/MM, 6=R3R, 7=R3MR,
!                 8=3, 9=3M1, 10=31M, 11=6/M, 12=6/MMM, 13=M3 AND 14=M3M
!       SGUNIQ  Output Unique axis in monoclinic space groups
!                 = 4 on error exits; = -1 for rhombahedral in hexagonal setting
!       SGINV   Output 1Bar flag  (0/1) for (acentric/centric)
!       SGLATT  Output Lattice centering no.
!                 1=P, 2=A, 3=B, 4=C, 5=I, 6=F AND 7=R
!       SGNOPS  Output The no. of matrices generated
!       SGPOL   Output The polar axis flag 
!                 1=x, 2=y, 3=x y, 4=z, 5=x z, 6=y z, 7=xyz, 8=111
!       JRT     Output The NSYM (3,5,NSYM) matrices
!       CEN     Output The lattice centering vectors
!       SGNCEN  Output The no. of lattice centering vectors
!       RT      Scratch array of 500 words needed by sgroup
!       IER     Error flag no.

  SPG="F m m 2"
  SPG="I 4 2 2"
  SPG="P 3 m 1"
  SPG="P 21/c"
  call getarg(1,SPG)

  CALL SGROUPNP(SPG,LAUE,SGUNIQ,SGINV,SGLATT,SGNOPS,SGPOL,JRT,CEN,SGNCEN,RT,IERR)
  if (ierr /= 0) stop "Unknown space group : "//trim(SPG)
  write(0,*) "SPG    " , "SPG    " , SPG   
  write(0,*) "LAUE   " , "LAUENO " , LAUE
  write(0,*) "SGUNIQ " , "NAXIS  " , SGUNIQ
  write(0,*) "SGINV  " , "NCENT  " , SGINV
  write(0,*) "SGLATT " , "LCENT  " , SGLATT
  write(0,*) "SGNOPS " , "NSYM   " , SGNOPS
  write(0,*) "SGPOL  " , "NPOLL  " , SGPOL
  DO K=1,SGNOPS
  write(0,*) "JRT    " , JRT(:,:,k)
  ENDDO
  write(0,*) "CEN    " , "CEN    " , CEN(:,1:SGNCEN)
  write(0,*) "SGNCEN " , "NCV    " , SGNCEN
  write(0,*) "IERR   " , "IERR   " , IERR

  nop = 0
  
  DO L=0,SGINV
    DO K=1,SGNOPS
      nop = nop+1
      DO J=1,3
        DO I=1,3
          rot(i,j,nop) = real(JRT(I,J,K),8) * (1-2*l)
        END DO
        rot(j,4,nop) = real(JRT(j,4,K),8)/12.d0
      END DO
    END DO
  END DO
  
  do k=1,nop
    write(0,*)k
    write(0,*)"Rotation matrix"
    write(0,'(3f4.0)')rot(:,1,k)
    write(0,'(3f4.0)')rot(:,2,k)
    write(0,'(3f4.0)')rot(:,3,k)
    write(0,*)"Translation vector"
    write(0,*)rot(:,4,k)
  enddo
  write(0,*)
  write(0,*)"Centering operators"
  do k=1,sgncen
    write(0,*)k
    write(0,*)CEN(:,K)
  enddo

  stop
end program test
