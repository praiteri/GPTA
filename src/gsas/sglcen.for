      SUBROUTINE SGLCEN(LCENT,CEN,NCV)

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

      INTEGER*4     LCENT               !Lattice centering type flag
      REAL*4        CEN(3,4)            !List of lattice centering vectors
      INTEGER*4     NCV                 !Number of lattcie centering vectors

      REAL*4        CENV(3,6)           
      INTEGER*4     NCVT(7)             

      DATA NCVT/1,2,2,2,2,4,3/
      DATA CENV/  0,0.5,0.5,  0.5,0,0.5,  0.5,0.5,0,  0.5,0.5,0.5,
     1  0.3333333,0.6666667,0.6666667,  0.6666667,0.3333333,0.3333333/

      NCV = NCVT(LCENT)
      CEN(1,1) = 0.0
      CEN(2,1) = 0.0
      CEN(3,1) = 0.0
      IF ( NCV.GT.1 ) THEN
        J = LCENT-1
        IF ( LCENT.EQ.6 ) J=1
        IF ( LCENT.EQ.7 ) J=5
        DO I=2,NCV                                          !Copy the lattice centering vectors
          CEN(1,I) = CENV(1,J)
          CEN(2,I) = CENV(2,J)
          CEN(3,I) = CENV(3,J)
          J = J+1
        END DO
      END IF
      RETURN
      END
