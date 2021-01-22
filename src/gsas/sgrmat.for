      SUBROUTINE SGRMAT(IOP,RT,N,A,B,C,D,E,F,G,H,O)

!Purpose:      S.R. to create a 3,3 matrix from 9 scalers

!Written by:       Allen C. Larson and Robert B. Von Dreele
!                  MS-H805
!                  Los Alamos National Laboratory
!                  Los Alamos, NM  87545

!                  Copyright, 1984-2000, The Regents of the University of California.

!      This software was produced under a U.S. Government contract (W-7405-ENG-36) by the Los Alamos
!      National Laboratory, which is operated by the University of California for the U.S. Department of
!      Energy. The U.S. Government is licensed to use, reproduce, and distribute this software. Permission
!      is granted to the public to copy and use this software without charge, provided that this notice
!      and any statement of authorship are reproduced on all copies. Neither the Government nor the
!      University makes any warranty, express or implied, or assumes any liability or responsibility for
!      the use of this software.

C       THIS PROGRAM WAS DEVELOPED FOR
C                    THE DIVISION OF CHEMISTRY
C                               OF
C               THE NATIONAL RESEARCH COUNCIL OF CANADA
C                               BY
C       ALLEN C. LARSON, 14 CERRADO LOOP, SANTA FE, NM 87505, USA

!Calling sequence parameters:

      INTEGER*4     IOP                 !Matrix generator count
      REAL*4        RT(5,4,25)          !Matrix to be generated
      INTEGER*4     N                   !Number of the matrix to be generated
      REAL*4        A,B,C,D,E,F,G,H,O   !Matrix terms

!Local varaibles:

!Code:

      RT(1,1,N) = A
      RT(1,2,N) = B
      RT(1,3,N) = C
      RT(1,4,N) = 0.0
      RT(2,1,N) = D
      RT(2,2,N) = E
      RT(2,3,N) = F
      RT(2,4,N) = 0.0
      RT(3,1,N) = G
      RT(3,2,N) = H
      RT(3,3,N) = O
      RT(3,4,N) = 0.0
      RT(4,1,N) = 0.0
      RT(4,2,N) = 0.0
      RT(4,3,N) = 0.0
      RT(4,4,N) = 1.0
      RT(5,1,N) = 81*(2*RT(1,1,N)+3*RT(1,2,N)+4*RT(1,3,N))
     1  +9*(2*RT(2,1,N)+3*RT(2,2,N)+4*RT(2,3,N))
     1  +2*RT(3,1,N)+3*RT(3,2,N)+4*RT(3,3,N)
      RT(5,2,N) = 0.0                                          !Clear the translation info
      RT(5,3,N) = IOP
      RT(5,4,N) = 20.

      RETURN
      END
