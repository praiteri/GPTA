      SUBROUTINE SGMTML(X,I,J,K)

!Purpose:      Form product of operators to generate the full group

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

C       This program was developed for
C                    The Division of Chemistry
C                               of
C               The National Research Council of Canada
C                               by
C       Allen C. Larson, 14 Cerrado Loop, Santa Fe, NM 87505, USA

      REAL*4        X(5,4,25)           
      INTEGER*4     I                   !Input matrix number
      INTEGER*4     J                   !Input matrix number
      INTEGER*4     K                   !Output matrix number

      DO L=1,4
        DO M=1,4
          X(L,M,K) = 0.0
          DO N=1,4
            X(L,M,K) = X(L,M,K)+X(L,N,J)*X(N,M,I)
          END DO
        END DO
      END DO

      X(1,4,K) = MOD(NINT((7.0+X(1,4,K))*12)/12.0,1.0)            !Force the translations to be in the cell
      X(2,4,K) = MOD(NINT((7.0+X(2,4,K))*12)/12.0,1.0)            !Also reset them to the value nearest to n/12
      X(3,4,K) = MOD(NINT((7.0+X(3,4,K))*12)/12.0,1.0)

      X(5,1,K) = 81*(2*X(1,1,K)+3*X(1,2,K)+4*X(1,3,K))            !Calculate a matrix flag
     1  +9*(2*X(2,1,K)+3*X(2,2,K)+4*X(2,3,K))
     1  +2*X(3,1,K)+3*X(3,2,K)+4*X(3,3,K)
      X(5,2,K) = 1728*X(1,4,K)+144*X(2,4,K)+12*X(3,4,K)            !Calculate the translation flag
      X(5,2,K) = NINT(X(5,2,K))                              !These should be whole numbers
      X(5,3,K) = IEOR(NINT(X(5,3,J)),NINT(X(5,3,I)))                  !Note the generator matrix number
      X(5,4,K) = 0.0

      RETURN
      END
