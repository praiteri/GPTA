! ! Copyright (c) 2021, Paolo Raiteri, Curtin University.
! ! All rights reserved.
! ! 
! ! This program is free software; you can redistribute it and/or modify it 
! ! under the terms of the GNU General Public License as published by the 
! ! Free Software Foundation; either version 3 of the License, or 
! ! (at your option) any later version.
! !  
! ! Redistribution and use in source and binary forms, with or without 
! ! modification, are permitted provided that the following conditions are met:
! ! 
! ! * Redistributions of source code must retain the above copyright notice, 
! !   this list of conditions and the following disclaimer.
! ! * Redistributions in binary form must reproduce the above copyright notice, 
! !   this list of conditions and the following disclaimer in the documentation 
! !   and/or other materials provided with the distribution.
! ! * Neither the name of the <ORGANIZATION> nor the names of its contributors 
! !   may be used to endorse or promote products derived from this software 
! !   without specific prior written permission.
! ! 
! ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! ! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! ! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! ! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! ! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! ! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! ! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! ! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! ! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! ! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ! 
module moduleSymmetry 

  implicit none

  integer, parameter :: nspg=230
  character(len=16), dimension(nspg) :: group_name

  integer(4) :: sgnops, sgncen, sgnops_tot
  real(4) :: cen(3,4)
  real(8) :: rot(3,4,48)

contains  

  function get_spg_name(id) result(name)
    
    implicit none
    integer, intent(in) :: id
    character(16) :: name
    
    call init_spg()
    if (id <1 .or. id>nspg) then
      name = "XXXX"
    else
      name = group_name(id)
    end if
    
    return
  end function get_spg_name
  
  subroutine init_spg()
    implicit none

! Triclinic
    group_name(  1) = "P 1"
    group_name(  2) = "P -1"

! Monoclinic
    group_name(  3) = "P 1 2 1"
    group_name(  4) = "P 1 21 1"
    group_name(  5) = "C 1 2 1"
    group_name(  6) = "P 1 m 1"
    group_name(  7) = "P 1 c 1"
    group_name(  8) = "C 1 m 1"
    group_name(  9) = "C 1 c 1"
    group_name( 10) = "P 1 2/m 1"
    group_name( 11) = "P 1 21/m 1"
    group_name( 12) = "C 1 2/m 1"
    group_name( 13) = "P 1 2/c 1"
    group_name( 14) = "P 1 21/c 1"
    group_name( 15) = "C 1 2/c 1"

! Orthorhombic
    group_name( 16) = "P 2 2 2"
    group_name( 17) = "P 2 2 21"
    group_name( 18) = "P 21 21 2"
    group_name( 19) = "P 21 21 2 1"
    group_name( 20) = "C 2 2 21"
    group_name( 21) = "C 2 2 2"
    group_name( 22) = "F 2 2 2"
    group_name( 23) = "I 2 2 2"
    group_name( 24) = "I 21 21 2 1"
    group_name( 25) = "P m m 2"
    group_name( 26) = "P m c 21"
    group_name( 27) = "P c c 2"
    group_name( 28) = "P m a 2"
    group_name( 29) = "P c a 21"
    group_name( 30) = "P n c 2"
    group_name( 31) = "P m n 21"
    group_name( 32) = "P b a 2"
    group_name( 33) = "P n a 21"
    group_name( 34) = "P n n 2"
    group_name( 35) = "C m m 2"
    group_name( 36) = "C m c 21"
    group_name( 37) = "C c c 2"
    group_name( 38) = "A m m 2"
    group_name( 39) = "A b m 2"
    group_name( 40) = "A m a 2"
    group_name( 41) = "A b a 2"
    group_name( 42) = "F m m 2"
    group_name( 43) = "F d d 2"
    group_name( 44) = "I m m 2"
    group_name( 45) = "I b a 2"
    group_name( 46) = "I m a 2"
    group_name( 47) = "P m m m"
    group_name( 48) = "P n n n"
    group_name( 49) = "P c c m"
    group_name( 50) = "P b a n"
    group_name( 51) = "P m m a"
    group_name( 52) = "P n n a"
    group_name( 53) = "P m n a"
    group_name( 54) = "P c c a"
    group_name( 55) = "P b a m"
    group_name( 56) = "P c c n"
    group_name( 57) = "P b c m"
    group_name( 58) = "P n n m"
    group_name( 59) = "P m m n"
    group_name( 60) = "P b c n"
    group_name( 61) = "P b c a"
    group_name( 62) = "P n m a"
    group_name( 63) = "C m c m"
    group_name( 64) = "C m c a"
    group_name( 65) = "C m m m"
    group_name( 66) = "C c c m"
    group_name( 67) = "C m m a"
    group_name( 68) = "C c c a"
    group_name( 69) = "F m m m"
    group_name( 70) = "F d d d"
    group_name( 71) = "I m m m"
    group_name( 72) = "I b a m"
    group_name( 73) = "I b c a"
    group_name( 74) = "I m m a"

!  Tetragonal
    group_name( 75) = "P 4"
    group_name( 76) = "P 41"
    group_name( 77) = "P 42"
    group_name( 78) = "P 43"
    group_name( 79) = "I 4"
    group_name( 80) = "I 41"
    group_name( 81) = "P -4"
    group_name( 82) = "I -4"
    group_name( 83) = "P 4/m"
    group_name( 84) = "P 42/m"
    group_name( 85) = "P 4/n"
    group_name( 86) = "P 42/n"
    group_name( 87) = "I 4/m"
    group_name( 88) = "I 41/a"
    group_name( 89) = "P 4 2 2"
    group_name( 90) = "P 4 21 2"
    group_name( 91) = "P 41 2 2"
    group_name( 92) = "P 41 21 2"
    group_name( 93) = "P 42 2 2"
    group_name( 94) = "P 42 21 2"
    group_name( 95) = "P 43 2 2"
    group_name( 96) = "P 43 21 2"
    group_name( 97) = "I 4 2 2"
    group_name( 98) = "I 41 2 2"
    group_name( 99) = "P 4 m m"
    group_name(100) = "P 4 b m"
    group_name(101) = "P 42 c m"
    group_name(102) = "P 42 n m"
    group_name(103) = "P 4 c c"
    group_name(104) = "P 4 n c"
    group_name(105) = "P 42 m c"
    group_name(106) = "P 42 b c"
    group_name(107) = "I 4 m m"
    group_name(108) = "I 4 c m"
    group_name(109) = "I 41 m d"
    group_name(110) = "I 41 c d"
    group_name(111) = "P -4 2 m"
    group_name(112) = "P -4 2 c"
    group_name(113) = "P -4 21 m"
    group_name(114) = "P -4 21 c"
    group_name(115) = "P -4 m 2"
    group_name(116) = "P -4 c 2"
    group_name(117) = "P -4 b 2"
    group_name(118) = "P -4 n 2"
    group_name(119) = "I -4 m 2"
    group_name(120) = "I -4 c 2"
    group_name(121) = "I -4 2 m"
    group_name(122) = "I -4 2 d"
    group_name(123) = "P 4/m m m"
    group_name(124) = "P 4/m c c"
    group_name(125) = "P 4/n b m"
    group_name(126) = "P 4/n n c"
    group_name(127) = "P 4/m b m"
    group_name(128) = "P 4/m n c"
    group_name(129) = "P 4/n m m"
    group_name(130) = "P 4/n c c"
    group_name(131) = "P 42/m m c"
    group_name(132) = "P 42/m c m"
    group_name(133) = "P 42/n b c"
    group_name(134) = "P 42/n n m"
    group_name(135) = "P 42/m b c"
    group_name(136) = "P 42/m n m"
    group_name(137) = "P 42/n m c"
    group_name(138) = "P 42/n c m"
    group_name(139) = "I 4/m m m"
    group_name(140) = "I 4/m c m"
    group_name(141) = "I 41/a m d"
    group_name(142) = "I 41/a c d"

! Trigonal
    group_name(143) = "P 3"
    group_name(144) = "P 31"
    group_name(145) = "P 32"
    group_name(146) = "R 3"
    group_name(147) = "P -3"
    group_name(148) = "R -3"
    group_name(149) = "P 3 1 2"
    group_name(150) = "P 3 2 1"
    group_name(151) = "P 31 1 2"
    group_name(152) = "P 31 2 1"
    group_name(153) = "P 32 1 2"
    group_name(154) = "P 32 2 1"
    group_name(155) = "R 3 2"
    group_name(156) = "P 3 m 1"
    group_name(157) = "P 3 1 m"
    group_name(158) = "P 3 c 1"
    group_name(159) = "P 3 1 c"
    group_name(160) = "R 3 m"
    group_name(161) = "R 3 c"
    group_name(162) = "P -3 1 m"
    group_name(163) = "P -3 1 c"
    group_name(164) = "P -3 m 1"
    group_name(165) = "P -3 c 1"
    group_name(166) = "R -3 m"
    group_name(167) = "R -3 c"

  ! Hexagonal
    group_name(168) = "P 6"
    group_name(169) = "P 61"
    group_name(170) = "P 65"
    group_name(171) = "P 62"
    group_name(172) = "P 64"
    group_name(173) = "P 63"
    group_name(174) = "P -6"
    group_name(175) = "P 6/m"
    group_name(176) = "P 63/m"
    group_name(177) = "P 6 2 2"
    group_name(178) = "P 61 2 2"
    group_name(179) = "P 65 2 2"
    group_name(180) = "P 62 2 2"
    group_name(181) = "P 64 2 2"
    group_name(182) = "P 63 2 2"
    group_name(183) = "P 6 m m"
    group_name(184) = "P 6 c c"
    group_name(185) = "P 63 c m"
    group_name(186) = "P 63 m c"
    group_name(187) = "P -6 m 2"
    group_name(188) = "P -6 c 2"
    group_name(189) = "P -6 2 m"
    group_name(190) = "P -6 2 c"
    group_name(191) = "P 6/m m m"
    group_name(192) = "P 6/m c c"
    group_name(193) = "P 63/m c m"
    group_name(194) = "P 63/m m c"
! Cubic
    group_name(195) = "P 2 3"
    group_name(196) = "F 2 3"
    group_name(197) = "I 2 3"
    group_name(198) = "P 21 3"
    group_name(199) = "I 21 3"
    group_name(200) = "P m -3"
    group_name(201) = "P n -3"
    group_name(202) = "F m -3"
    group_name(203) = "F d -3"
    group_name(204) = "I m -3"
    group_name(205) = "P a -3"
    group_name(206) = "I a -3"
    group_name(207) = "P 4 3 2"
    group_name(208) = "P 42 3 2"
    group_name(209) = "F 4 3 2"
    group_name(210) = "F 41 3 2"
    group_name(211) = "I 4 3 2"
    group_name(212) = "P 43 3 2"
    group_name(213) = "P 41 3 2"
    group_name(214) = "I 41 3 2"
    group_name(215) = "P -4 3 m"
    group_name(216) = "F -4 3 m"
    group_name(217) = "I -4 3 m"
    group_name(218) = "P -4 3 n"
    group_name(219) = "F -4 3 c"
    group_name(220) = "I -4 3 d"
    group_name(221) = "P m -3 m"
    group_name(222) = "P n -3 n"
    group_name(223) = "P m -3 n"
    group_name(224) = "P n -3 m"
    group_name(225) = "F m -3 m"
    group_name(226) = "F m -3 c"
    group_name(227) = "F d -3 m"
    group_name(228) = "F d -3 c"
    group_name(229) = "I m -3 m"
    group_name(230) = "I a -3 d"
    return
  end subroutine init_spg

end module moduleSymmetry 
