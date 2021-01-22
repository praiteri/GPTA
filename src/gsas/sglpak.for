      SUBROUTINE SGLPAK(L,IER)

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
C       ALLEN C. LARSON, P.O.BOX 5898, SANTA FE, NM 87502,USA

      DIMENSION     L(4)                

      IF ( L(2).LT.12 ) IER=4
      IF ( L(2).GT.17 ) IER=4
      L(1) = L(2)
      L(2) = 18
      RETURN
      END
