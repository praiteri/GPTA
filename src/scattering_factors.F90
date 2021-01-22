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
module scattering_factors

  type sfactor

     integer  :: atomicnumber
     character(len=2)  :: elementsymbol
     double precision, dimension(4)  :: a_anis ,b_anis
     double precision c_anis
     double precision isot

  end type sfactor


!!p ------------------------ Scattering factors library -----------------------
!!P http://www.ruppweb.org/Xray/comp/scatfac.html
!!p Cromer-Mann coefficients from :
!!p
!!p P. J. Brown, A. G. Fox, E. N. Maslen, M. A. O'Keefe and B. T. M. Willis
!!p International Tables for Crystallography (2006). Vol. C, ch. 6.1, pp. 554-595
!!p Chapter 6.1. Intensity of diffracted intensities
!!p doi: 10.1107/97809553602060000600
!!p 


!!p Dummy atom
 type (SFactor), parameter :: DM = SFactor( 0,"DM",(/  0.00000d0 ,  0.00000d0 ,  0.00000d0 , 0.00000d0/) , & 
                                                   (/  0.00000d0 ,  0.00000d0 ,  0.00000d0 ,    0.00000d0/) ,  0.00000d0 ,  0.0d0)
 type (SFactor), parameter ::  H = SFactor( 1," H",(/  0.48992d0 ,  0.26200d0 ,  0.19677d0 , 0.04988d0/) , & 
                                                   (/ 20.65930d0 ,  7.74039d0 ,  49.55190d0 ,   2.20159d0/) ,  0.00130d0 ,  1.0d0)
 type (SFactor), parameter :: HE = SFactor( 2,"HE",(/  0.87340d0 ,  0.63090d0 ,  0.31120d0 , 0.17800d0/) , & 
                                                   (/  9.10370d0 ,  3.35680d0 ,  22.92760d0 ,   0.98210d0/) ,  0.00640d0 ,  2.0d0)
 type (SFactor), parameter :: LI = SFactor( 3,"LI",(/  1.12820d0 ,  0.75080d0 ,  0.61750d0 , 0.46530d0/) , & 
                                                   (/  3.95460d0 ,  1.05240d0 ,  85.39050d0 , 168.26100d0/) ,  0.03770d0 ,  3.0d0)
 type (SFactor), parameter :: BE = SFactor( 4,"BE",(/  1.59190d0 ,  1.12780d0 ,  0.53910d0 , 0.70290d0/) , & 
                                                   (/ 43.64270d0 ,  1.86230d0 , 103.48300d0 ,   0.54200d0/) ,  0.03850d0 ,  4.0d0)
 type (SFactor), parameter ::  B = SFactor( 5," B",(/  2.05450d0 ,  1.33260d0 ,  1.09790d0 , 0.70680d0/) , & 
                                                   (/ 23.21850d0 ,  1.02100d0 ,  60.34980d0 ,   0.14030d0/) ,  0.00000d0 ,  5.0d0)
 type (SFactor), parameter ::  C = SFactor( 6," C",(/  2.31000d0 ,  1.02000d0 ,  1.58860d0 , 0.86500d0/) , & 
                                                   (/ 20.84390d0 , 10.20750d0 ,   0.56870d0 ,  51.65120d0/) ,  0.21560d0 ,  6.0d0)
 type (SFactor), parameter ::  N = SFactor( 7," N",(/ 12.21260d0 ,  3.13220d0 ,  2.01250d0 , 1.16630d0/) , & 
                                                   (/  0.00570d0 ,  9.89330d0 ,  28.99750d0 ,   0.58260d0/) ,  0.00000d0 ,  7.0d0)
 type (SFactor), parameter ::  O = SFactor( 8," O",(/  3.04850d0 ,  2.28680d0 ,  1.54630d0 , 0.86700d0/) , & 
                                                   (/ 13.27710d0 ,  5.70110d0 ,   0.32390d0 ,  32.90890d0/) ,  0.25080d0 ,  8.0d0)
 type (SFactor), parameter ::  F = SFactor( 9," F",(/  3.53920d0 ,  2.64120d0 ,  1.51700d0 , 1.02430d0/) , & 
                                                   (/ 10.28250d0 ,  4.29440d0 ,   0.26150d0 ,  26.14760d0/) ,  0.27760d0 ,  9.0d0)
 type (SFactor), parameter :: NE = SFactor(10,"NE",(/  3.95530d0 ,  3.11250d0 ,  1.45460d0 , 1.12510d0/) , & 
                                                   (/  8.40420d0 ,  3.42620d0 ,   0.23060d0 ,  21.71840d0/) ,  0.35150d0 , 10.0d0)
 type (SFactor), parameter :: NA = SFactor(11,"NA",(/  4.76260d0 ,  3.17360d0 ,  1.26740d0 , 1.11280d0/) , & 
                                                   (/  3.28500d0 ,  8.84220d0 ,   0.31360d0 , 129.42400d0/) ,  0.67600d0 , 11.0d0)
 type (SFactor), parameter :: MG = SFactor(12,"MG",(/  5.42040d0 ,  2.17350d0 ,  1.22690d0 , 2.30730d0/) , & 
                                                   (/  2.82750d0 , 79.26110d0 ,   0.38080d0 ,   7.19370d0/) ,  0.85840d0 , 12.0d0)
 type (SFactor), parameter :: AL = SFactor(13,"AL",(/  6.42020d0 ,  1.90020d0 ,  1.59360d0 , 1.96460d0/) , & 
                                                   (/  3.03870d0 ,  0.74260d0 ,  31.54720d0 ,  85.08860d0/) ,  1.11510d0 , 13.0d0)
 type (SFactor), parameter :: SI = SFactor(14,"SI",(/  6.29150d0 ,  3.03530d0 ,  1.98910d0 , 1.54100d0/) , & 
                                                   (/  2.43860d0 , 32.33370d0 ,   0.67850d0 ,  81.69370d0/) ,  1.14070d0 , 14.0d0)
 type (SFactor), parameter ::  P = SFactor(15," P",(/  6.43450d0 ,  4.17910d0 ,  1.78000d0 , 1.49080d0/) , & 
                                                   (/  1.90670d0 , 27.15700d0 ,   0.52600d0 ,  68.16450d0/) ,  1.11490d0 , 15.0d0)
 type (SFactor), parameter ::  S = SFactor(16," S",(/  6.90530d0 ,  5.20340d0 ,  1.43790d0 , 1.58630d0/) , & 
                                                   (/  1.46790d0 , 22.21510d0 ,   0.25360d0 ,  56.17200d0/) ,  0.86690d0 , 16.0d0)
 type (SFactor), parameter :: CL = SFactor(17,"CL",(/ 11.46040d0 ,  7.19640d0 ,  6.25560d0 , 1.64550d0/) , & 
                                                   (/  0.01040d0 ,  1.16620d0 ,  18.51940d0 ,  47.77840d0/) ,  0.00000d0 , 17.0d0)
 type (SFactor), parameter :: AR = SFactor(18,"AR",(/  7.48450d0 ,  6.77230d0 ,  0.65390d0 , 1.64420d0/) , & 
                                                   (/  0.90720d0 , 14.84070d0 ,  43.89830d0 ,  33.39290d0/) ,  1.44450d0 , 18.0d0)
 type (SFactor), parameter ::  K = SFactor(19," K",(/  8.21860d0 ,  7.43980d0 ,  1.05190d0 , 0.86590d0/) , & 
                                                   (/ 12.79490d0 ,  0.77480d0 , 213.18700d0 ,  41.68410d0/) ,  1.42280d0 , 19.0d0)
 type (SFactor), parameter :: CA = SFactor(20,"CA",(/  8.62660d0 ,  7.38730d0 ,  1.58990d0 , 1.02110d0/) , & 
                                                   (/ 10.44210d0 ,  0.65990d0 ,  85.74840d0 , 178.43700d0/) ,  1.37510d0 , 20.0d0)
 type (SFactor), parameter :: SC = SFactor(21,"SC",(/  9.18900d0 ,  7.36790d0 ,  1.64090d0 , 1.46800d0/) , & 
                                                   (/  9.02130d0 ,  0.57290d0 , 136.10800d0 ,  51.35310d0/) ,  1.33290d0 , 21.0d0)
 type (SFactor), parameter :: TI = SFactor(22,"TI",(/  9.75950d0 ,  7.35580d0 ,  1.69910d0 , 1.90210d0/) , & 
                                                   (/  7.85080d0 ,  0.50000d0 ,  35.63380d0 , 116.10500d0/) ,  1.28070d0 , 22.0d0)
 type (SFactor), parameter ::  V = SFactor(23," V",(/ 10.29710d0 ,  7.35110d0 ,  2.07030d0 , 2.05710d0/) , & 
                                                   (/  6.86570d0 ,  0.43850d0 ,  26.89380d0 , 102.47800d0/) ,  1.21990d0 , 23.0d0)
 type (SFactor), parameter :: CR = SFactor(24,"CR",(/ 10.64060d0 ,  7.35370d0 ,  3.32400d0 , 1.49220d0/) , & 
                                                   (/  6.10380d0 ,  0.39200d0 ,  20.26260d0 ,  98.73990d0/) ,  1.18320d0 , 24.0d0)
 type (SFactor), parameter :: MN = SFactor(25,"MN",(/ 11.28190d0 ,  7.35730d0 ,  3.01930d0 , 2.24410d0/) , & 
                                                   (/  5.34090d0 ,  0.34320d0 ,  17.86740d0 ,  83.75430d0/) ,  1.08960d0 , 25.0d0)
 type (SFactor), parameter :: FE = SFactor(26,"FE",(/ 11.76950d0 ,  7.35730d0 ,  3.52220d0 , 2.30450d0/) , & 
                                                   (/  4.76110d0 ,  0.30720d0 ,  15.35350d0 ,  76.88050d0/) ,  1.03690d0 , 26.0d0)
 type (SFactor), parameter :: CO = SFactor(27,"CO",(/ 12.28410d0 ,  7.34090d0 ,  4.00340d0 , 2.34880d0/) , & 
                                                   (/  4.27910d0 ,  0.27840d0 ,  13.53590d0 ,  71.16920d0/) ,  1.01180d0 , 27.0d0)
 type (SFactor), parameter :: NI = SFactor(28,"NI",(/ 12.83760d0 ,  7.29200d0 ,  4.44380d0 , 2.38000d0/) , & 
                                                   (/  3.87850d0 ,  0.25650d0 ,  12.17630d0 ,  66.34210d0/) ,  1.03410d0 , 28.0d0)
 type (SFactor), parameter :: CU = SFactor(29,"CU",(/ 13.33800d0 ,  7.16760d0 ,  5.61580d0 , 1.67350d0/) , & 
                                                   (/  3.58280d0 ,  0.24700d0 ,  11.39660d0 ,  64.81260d0/) ,  1.19100d0 , 29.0d0)
 type (SFactor), parameter :: ZN = SFactor(30,"ZN",(/ 14.07430d0 ,  7.03180d0 ,  5.16520d0 , 2.41000d0/) , & 
                                                   (/  3.26550d0 ,  0.23330d0 ,  10.31630d0 ,  58.70970d0/) ,  1.30410d0 , 30.0d0)
 type (SFactor), parameter :: GA = SFactor(31,"GA",(/ 15.23540d0 ,  6.70060d0 ,  4.35910d0 , 2.96230d0/) , & 
                                                   (/  3.06690d0 ,  0.24120d0 ,  10.78050d0 ,  61.41350d0/) ,  1.71890d0 , 31.0d0)
 type (SFactor), parameter :: GE = SFactor(32,"GE",(/ 16.08160d0 ,  6.37470d0 ,  3.70680d0 , 3.68300d0/) , & 
                                                   (/  2.85090d0 ,  0.25160d0 ,  11.44680d0 ,  54.76250d0/) ,  2.13130d0 , 32.0d0)
 type (SFactor), parameter :: AS = SFactor(33,"AS",(/ 16.67230d0 ,  6.07010d0 ,  3.43130d0 , 4.27790d0/) , & 
                                                   (/  2.63450d0 ,  0.26470d0 ,  12.94790d0 ,  47.79720d0/) ,  2.53100d0 , 33.0d0)
 type (SFactor), parameter :: SE = SFactor(34,"SE",(/ 17.00060d0 ,  5.81960d0 ,  3.97310d0 , 4.35430d0/) , & 
                                                   (/  2.40980d0 ,  0.27260d0 ,  15.23720d0 ,  43.81630d0/) ,  2.84090d0 , 34.0d0)
 type (SFactor), parameter :: BR = SFactor(35,"BR",(/ 17.17890d0 ,  5.23580d0 ,  5.63770d0 , 3.98510d0/) , & 
                                                   (/  2.17230d0 , 16.57960d0 ,   0.26090d0 ,  41.43280d0/) ,  2.95570d0 , 35.0d0)
 type (SFactor), parameter :: KR = SFactor(36,"KR",(/ 17.35550d0 ,  6.72860d0 ,  5.54930d0 , 3.53750d0/) , & 
                                                   (/  1.93840d0 , 16.56230d0 ,   0.22610d0 ,  39.39720d0/) ,  2.82500d0 , 36.0d0)
 type (SFactor), parameter :: RB = SFactor(37,"RB",(/ 17.17840d0 ,  9.64350d0 ,  5.13990d0 , 1.52920d0/) , & 
                                                   (/  1.78880d0 , 17.31510d0 ,   0.27480d0 , 164.93400d0/) ,  3.48730d0 , 37.0d0)
 type (SFactor), parameter :: SR = SFactor(38,"SR",(/ 17.56630d0 ,  9.81840d0 ,  5.42200d0 , 2.66940d0/) , & 
                                                   (/  1.55640d0 , 14.09880d0 ,   0.16640d0 , 132.37600d0/) ,  2.50640d0 , 38.0d0)
 type (SFactor), parameter ::  Y = SFactor(39," Y",(/ 17.77600d0 , 10.29460d0 ,  5.72629d0 , 3.26588d0/) , & 
                                                   (/  1.40290d0 , 12.80060d0 ,   0.12560d0 , 104.35400d0/) ,  1.91213d0 , 39.0d0)
 type (SFactor), parameter :: ZR = SFactor(40,"ZR",(/ 17.87650d0 , 10.94800d0 ,  5.41732d0 , 3.65721d0/) , & 
                                                   (/  1.27618d0 , 11.91600d0 ,   0.11762d0 ,  87.66270d0/) ,  2.06929d0 , 40.0d0)
 type (SFactor), parameter :: NB = SFactor(41,"NB",(/ 17.61420d0 , 12.01440d0 ,  4.04183d0 , 3.53346d0/) , & 
                                                   (/  1.18865d0 , 11.76600d0 ,   0.20478d0 ,  69.79570d0/) ,  3.75591d0 , 41.0d0)
 type (SFactor), parameter :: MO = SFactor(42,"MO",(/  3.70250d0 , 17.23560d0 , 12.88760d0 , 3.74290d0/) , & 
                                                   (/  0.27720d0 ,  1.09580d0 ,  11.00400d0 ,  61.65840d0/) ,  4.38750d0 , 42.0d0)
 type (SFactor), parameter :: TC = SFactor(43,"TC",(/ 19.13010d0 , 11.09480d0 ,  4.64901d0 , 2.71263d0/) , & 
                                                   (/  0.86413d0 ,  8.14487d0 ,  21.57070d0 ,  86.84720d0/) ,  5.40428d0 , 43.0d0)
 type (SFactor), parameter :: RU = SFactor(44,"RU",(/ 19.26740d0 , 12.91820d0 ,  4.86337d0 , 1.56756d0/) , & 
                                                   (/  0.80852d0 ,  8.43467d0 ,  24.79970d0 ,  94.29280d0/) ,  5.37874d0 , 44.0d0)
 type (SFactor), parameter :: RH = SFactor(45,"RH",(/ 19.29570d0 , 14.35010d0 ,  4.73425d0 , 1.28918d0/) , & 
                                                   (/  0.75154d0 ,  8.21758d0 ,  25.87490d0 ,  98.60620d0/) ,  5.32800d0 , 45.0d0)
 type (SFactor), parameter :: PD = SFactor(46,"PD",(/ 19.33190d0 , 15.50170d0 ,  5.29537d0 , 0.60584d0/) , & 
                                                   (/  0.69866d0 ,  7.98929d0 ,  25.20520d0 ,  76.89860d0/) ,  5.26593d0 , 46.0d0)
 type (SFactor), parameter :: AG = SFactor(47,"AG",(/ 19.28080d0 , 16.68850d0 ,  4.80450d0 , 1.04630d0/) , & 
                                                   (/  0.64460d0 ,  7.47260d0 ,  24.66050d0 ,  99.81560d0/) ,  5.17900d0 , 47.0d0)
 type (SFactor), parameter :: CD = SFactor(48,"CD",(/ 19.22140d0 , 17.64440d0 ,  4.46100d0 , 1.60290d0/) , & 
                                                   (/  0.59460d0 ,  6.90890d0 ,  24.70080d0 ,  87.48250d0/) ,  5.06940d0 , 48.0d0)
 type (SFactor), parameter :: IN = SFactor(49,"IN",(/ 19.16240d0 , 18.55960d0 ,  4.29480d0 , 2.03960d0/) , & 
                                                   (/  0.54760d0 ,  6.37760d0 ,  25.84990d0 ,  92.80290d0/) ,  4.93910d0 , 49.0d0)
 type (SFactor), parameter :: SN = SFactor(50,"SN",(/ 19.18890d0 , 19.10050d0 ,  4.45850d0 , 2.46630d0/) , & 
                                                   (/  5.83030d0 ,  0.50310d0 ,  26.89090d0 ,  83.95710d0/) ,  4.78210d0 , 50.0d0)
 type (SFactor), parameter :: SB = SFactor(51,"SB",(/ 19.64180d0 , 19.04550d0 ,  5.03710d0 , 2.68270d0/) , & 
                                                   (/  5.30340d0 ,  0.46070d0 ,  27.90740d0 ,  75.28250d0/) ,  4.59090d0 , 51.0d0)
 type (SFactor), parameter :: TE = SFactor(52,"TE",(/ 19.96440d0 , 19.01380d0 ,  6.14487d0 , 2.52390d0/) , & 
                                                   (/  4.81742d0 ,  0.42089d0 ,  28.52840d0 ,  70.84030d0/) ,  4.35200d0 , 52.0d0)
 type (SFactor), parameter ::  I = SFactor(53," I",(/ 20.14720d0 , 18.99490d0 ,  7.51380d0 , 2.27350d0/) , & 
                                                   (/  4.34700d0 ,  0.38140d0 ,  27.76600d0 ,  66.87760d0/) ,  4.07120d0 , 53.0d0)
 type (SFactor), parameter :: XE = SFactor(54,"XE",(/ 20.29330d0 , 19.02980d0 ,  8.97670d0 , 1.99000d0/) , & 
                                                   (/  3.92820d0 ,  0.34400d0 ,  26.46590d0 ,  64.26580d0/) ,  3.71180d0 , 54.0d0)
 type (SFactor), parameter :: CS = SFactor(55,"CS",(/ 20.38920d0 , 19.10620d0 , 10.66200d0 , 1.49530d0/) , & 
                                                   (/  3.56900d0 ,  0.31070d0 ,  24.38790d0 , 213.90400d0/) ,  3.33520d0 , 55.0d0)
 type (SFactor), parameter :: BA = SFactor(56,"BA",(/ 20.33610d0 , 19.29700d0 , 10.88800d0 , 2.69590d0/) , & 
                                                   (/  3.21600d0 ,  0.27560d0 ,  20.20730d0 , 167.20200d0/) ,  2.77310d0 , 56.0d0)
 type (SFactor), parameter :: LA = SFactor(57,"LA",(/ 20.57800d0 , 19.59900d0 , 11.37270d0 , 3.28719d0/) , & 
                                                   (/  2.94817d0 ,  0.24447d0 ,  18.77260d0 , 133.12400d0/) ,  2.14678d0 , 57.0d0)
 type (SFactor), parameter :: CE = SFactor(58,"CE",(/ 21.16710d0 , 19.76950d0 , 11.85130d0 , 3.33049d0/) , & 
                                                   (/  2.81219d0 ,  0.22684d0 ,  17.60830d0 , 127.11300d0/) ,  1.86264d0 , 58.0d0)
 type (SFactor), parameter :: PR = SFactor(59,"PR",(/ 22.04400d0 , 19.66970d0 , 12.38560d0 , 2.82428d0/) , & 
                                                   (/  2.77393d0 ,  0.22209d0 ,  16.76690d0 , 143.64400d0/) ,  2.05830d0 , 59.0d0)
 type (SFactor), parameter :: ND = SFactor(60,"ND",(/ 22.68450d0 , 19.68470d0 , 12.77400d0 , 2.85137d0/) , & 
                                                   (/  2.66248d0 ,  0.21063d0 ,  15.88500d0 , 137.90300d0/) ,  1.98486d0 , 60.0d0)
 type (SFactor), parameter :: PM = SFactor(61,"PM",(/ 23.34050d0 , 19.60950d0 , 13.12350d0 , 2.87516d0/) , & 
                                                   (/  2.56270d0 ,  0.20209d0 ,  15.10090d0 , 132.72100d0/) ,  2.02876d0 , 61.0d0)
 type (SFactor), parameter :: SM = SFactor(62,"SM",(/ 24.00420d0 , 19.42580d0 , 13.43960d0 , 2.89604d0/) , & 
                                                   (/  2.47274d0 ,  0.19645d0 ,  14.39960d0 , 128.00700d0/) ,  2.20963d0 , 62.0d0)
 type (SFactor), parameter :: EU = SFactor(63,"EU",(/ 24.62740d0 , 19.08860d0 , 13.76030d0 , 2.92270d0/) , & 
                                                   (/  2.38790d0 ,  0.19420d0 ,  13.75460d0 , 123.17400d0/) ,  2.57450d0 , 63.0d0)
 type (SFactor), parameter :: GD = SFactor(64,"GD",(/ 25.07090d0 , 19.07980d0 , 13.85180d0 , 3.54545d0/) , & 
                                                   (/  2.25341d0 ,  0.18195d0 ,  12.93310d0 , 101.39800d0/) ,  2.41960d0 , 64.0d0)
 type (SFactor), parameter :: TB = SFactor(65,"TB",(/ 25.89760d0 , 18.21850d0 , 14.31670d0 , 2.95354d0/) , & 
                                                   (/  2.24256d0 ,  0.19614d0 ,  12.66480d0 , 115.36200d0/) ,  3.58324d0 , 65.0d0)
 type (SFactor), parameter :: DY = SFactor(66,"DY",(/ 26.50700d0 , 17.63830d0 , 14.55960d0 , 2.96577d0/) , & 
                                                   (/  2.18020d0 ,  0.20217d0 ,  12.18990d0 , 111.87400d0/) ,  4.29728d0 , 66.0d0)
 type (SFactor), parameter :: HO = SFactor(67,"HO",(/ 26.90490d0 , 17.29400d0 , 14.55830d0 , 3.63837d0/) , & 
                                                   (/  2.07051d0 ,  0.19794d0 ,  11.44070d0 ,  92.65660d0/) ,  4.56796d0 , 67.0d0)
 type (SFactor), parameter :: ER = SFactor(68,"ER",(/ 27.65630d0 , 16.42850d0 , 14.97790d0 , 2.98233d0/) , & 
                                                   (/  2.07356d0 ,  0.22354d0 ,  11.36040d0 , 105.70300d0/) ,  5.92046d0 , 68.0d0)
 type (SFactor), parameter :: TM = SFactor(69,"TM",(/ 28.18190d0 , 15.88510d0 , 15.15420d0 , 2.98706d0/) , & 
                                                   (/  2.02859d0 ,  0.23885d0 ,  10.99750d0 , 102.96100d0/) ,  6.75621d0 , 69.0d0)
 type (SFactor), parameter :: YB = SFactor(70,"YB",(/ 28.66410d0 , 15.43450d0 , 15.30870d0 , 2.98963d0/) , & 
                                                   (/  1.98890d0 ,  0.25712d0 ,  10.66470d0 , 100.41700d0/) ,  7.56672d0 , 70.0d0)
 type (SFactor), parameter :: LU = SFactor(71,"LU",(/ 28.94760d0 , 15.22080d0 , 15.10000d0 , 3.71601d0/) , & 
                                                   (/  1.90182d0 ,  9.98519d0 ,   0.26103d0 ,  84.32980d0/) ,  7.97628d0 , 71.0d0)
 type (SFactor), parameter :: HF = SFactor(72,"HF",(/ 29.14400d0 , 15.17260d0 , 14.75860d0 , 4.30013d0/) , & 
                                                   (/  1.83262d0 ,  9.59990d0 ,   0.27512d0 ,  72.02900d0/) ,  8.58154d0 , 72.0d0)
 type (SFactor), parameter :: TA = SFactor(73,"TA",(/ 29.20240d0 , 15.22930d0 , 14.51350d0 , 4.76492d0/) , & 
                                                   (/  1.77333d0 ,  9.37046d0 ,   0.29598d0 ,  63.36440d0/) ,  9.24354d0 , 73.0d0)
 type (SFactor), parameter ::  W = SFactor(74," W",(/ 29.08180d0 , 15.43000d0 , 14.43270d0 , 5.11982d0/) , & 
                                                   (/  1.72029d0 ,  9.22590d0 ,   0.32170d0 ,  57.05600d0/) ,  9.88750d0 , 74.0d0)
 type (SFactor), parameter :: RE = SFactor(75,"RE",(/ 28.76210d0 , 15.71890d0 , 14.55640d0 , 5.44174d0/) , & 
                                                   (/  1.67191d0 ,  9.09227d0 ,   0.35050d0 ,  52.08610d0/) , 10.47200d0 , 75.0d0)
 type (SFactor), parameter :: OS = SFactor(76,"OS",(/ 28.18940d0 , 16.15500d0 , 14.93050d0 , 5.67589d0/) , & 
                                                   (/  1.62903d0 ,  8.97948d0 ,   0.38266d0 ,  48.16470d0/) , 11.00050d0 , 76.0d0)
 type (SFactor), parameter :: IR = SFactor(77,"IR",(/ 27.30490d0 , 16.72960d0 , 15.61150d0 , 5.83377d0/) , & 
                                                   (/  1.59279d0 ,  8.86553d0 ,   0.41792d0 ,  45.00110d0/) , 11.47220d0 , 77.0d0)
 type (SFactor), parameter :: PT = SFactor(78,"PT",(/ 27.00590d0 , 17.76390d0 , 15.71310d0 , 5.78370d0/) , & 
                                                   (/  1.51293d0 ,  8.81174d0 ,   0.42459d0 ,  38.61030d0/) , 11.68830d0 , 78.0d0)
 type (SFactor), parameter :: AU = SFactor(79,"AU",(/ 16.88190d0 , 18.59130d0 , 25.55820d0 , 5.86000d0/) , & 
                                                   (/  0.46110d0 ,  8.62160d0 ,   1.48260d0 ,  36.39560d0/) , 12.06580d0 , 79.0d0)
 type (SFactor), parameter :: HG = SFactor(80,"HG",(/ 20.68090d0 , 19.04170d0 , 21.65750d0 , 5.96760d0/) , & 
                                                   (/  0.54500d0 ,  8.44840d0 ,   1.57290d0 ,  38.32460d0/) , 12.60890d0 , 80.0d0)
 type (SFactor), parameter :: TL = SFactor(81,"TL",(/ 27.54460d0 , 19.15840d0 , 15.53800d0 , 5.52593d0/) , & 
                                                   (/  0.65515d0 ,  8.70751d0 ,   1.96347d0 ,  45.81490d0/) , 13.17460d0 , 81.0d0)
 type (SFactor), parameter :: PB = SFactor(82,"PB",(/ 31.06170d0 , 13.06370d0 , 18.44200d0 , 5.96960d0/) , & 
                                                   (/  0.69020d0 ,  2.35760d0 ,   8.61800d0 ,  47.25790d0/) , 13.41180d0 , 82.0d0)
 type (SFactor), parameter :: BI = SFactor(83,"BI",(/ 33.36890d0 , 12.95100d0 , 16.58770d0 , 6.46920d0/) , & 
                                                   (/  0.70400d0 ,  2.92380d0 ,   8.79370d0 ,  48.00930d0/) , 13.57820d0 , 83.0d0)
 type (SFactor), parameter :: PO = SFactor(84,"PO",(/ 34.67260d0 , 15.47330d0 , 13.11380d0 , 7.02588d0/) , & 
                                                   (/  0.70100d0 ,  3.55078d0 ,   9.55642d0 ,  47.00450d0/) , 13.67700d0 , 84.0d0)
 type (SFactor), parameter :: AT = SFactor(85,"AT",(/ 35.31630d0 , 19.02110d0 ,  9.49887d0 , 7.42518d0/) , & 
                                                   (/  0.68587d0 ,  3.97458d0 ,  11.38240d0 ,  45.47150d0/) , 13.71080d0 , 85.0d0)
 type (SFactor), parameter :: RN = SFactor(86,"RN",(/ 35.56310d0 , 21.28160d0 ,  8.00370d0 , 7.44330d0/) , & 
                                                   (/  0.66310d0 ,  4.06910d0 ,  14.04220d0 ,  44.24730d0/) , 13.69050d0 , 86.0d0)
 type (SFactor), parameter :: FR = SFactor(87,"FR",(/ 35.92990d0 , 23.05470d0 , 12.14390d0 , 2.11253d0/) , & 
                                                   (/  0.64645d0 ,  4.17619d0 ,  23.10520d0 , 150.64500d0/) , 13.72470d0 , 87.0d0)
 type (SFactor), parameter :: RA = SFactor(88,"RA",(/ 35.76300d0 , 22.90640d0 , 12.47390d0 , 3.21097d0/) , & 
                                                   (/  0.61634d0 ,  3.87135d0 ,  19.98870d0 , 142.32500d0/) , 13.62110d0 , 88.0d0)
 type (SFactor), parameter :: AC = SFactor(89,"AC",(/ 35.65970d0 , 23.10320d0 , 12.59770d0 , 4.08655d0/) , & 
                                                   (/  0.58909d0 ,  3.65155d0 ,  18.59900d0 , 117.02000d0/) , 13.52660d0 , 89.0d0)
 type (SFactor), parameter :: TH = SFactor(90,"TH",(/ 35.56450d0 , 23.42190d0 , 12.74730d0 , 4.80703d0/) , & 
                                                   (/  0.56336d0 ,  3.46204d0 ,  17.83090d0 ,  99.17220d0/) , 13.43140d0 , 90.0d0)
 type (SFactor), parameter :: PA = SFactor(91,"PA",(/ 35.88470d0 , 23.29480d0 , 14.18910d0 , 4.17287d0/) , & 
                                                   (/  0.54775d0 ,  3.41519d0 ,  16.92350d0 , 105.25100d0/) , 13.42870d0 , 91.0d0)
 type (SFactor), parameter ::  U = SFactor(92," U",(/ 36.02280d0 , 23.41280d0 , 14.94910d0 , 4.18800d0/) , & 
                                                   (/  0.52930d0 ,  3.32530d0 ,  16.09270d0 , 100.61300d0/) , 13.39660d0 , 92.0d0)
 type (SFactor), parameter :: NP = SFactor(93,"NP",(/ 36.18740d0 , 23.59640d0 , 15.64020d0 , 4.18550d0/) , & 
                                                   (/  0.51193d0 ,  3.25396d0 ,  15.36220d0 ,  97.49080d0/) , 13.35730d0 , 93.0d0)
 type (SFactor), parameter :: PU = SFactor(94,"PU",(/ 36.52540d0 , 23.80830d0 , 16.77070d0 , 3.47947d0/) , & 
                                                   (/  0.49938d0 ,  3.26371d0 ,  14.94550d0 , 105.98000d0/) , 13.38120d0 , 94.0d0)
 type (SFactor), parameter :: AM = SFactor(95,"AM",(/ 36.67060d0 , 24.09920d0 , 17.34150d0 , 3.49331d0/) , & 
                                                   (/  0.48363d0 ,  3.20647d0 ,  14.31360d0 , 102.27300d0/) , 13.35920d0 , 95.0d0)
 type (SFactor), parameter :: CM = SFactor(96,"CM",(/ 36.64880d0 , 24.40960d0 , 17.39900d0 , 4.21665d0/) , & 
                                                   (/  0.46515d0 ,  3.08997d0 ,  13.43460d0 ,  88.48340d0/) , 13.28870d0 , 96.0d0)
 type (SFactor), parameter :: BK = SFactor(97,"BK",(/ 36.78810d0 , 24.77360d0 , 17.89190d0 , 4.23284d0/) , & 
                                                   (/  0.45102d0 ,  3.04619d0 ,  12.89460d0 ,  86.00300d0/) , 13.27540d0 , 97.0d0)
 type (SFactor), parameter :: CF = SFactor(98,"CF",(/ 36.91850d0 , 25.19950d0 , 18.33170d0 , 4.24391d0/) , & 
                                                   (/  0.43753d0 ,  3.00775d0 ,  12.40440d0 ,  83.78810d0/) , 13.26740d0 , 98.0d0)

 integer, parameter  :: SF_Library_dim = 98

 type(SFactor), parameter, dimension(0:SF_Library_dim)  :: SF_Library = &
     (/DM,H,HE,LI,BE,B,C,N,O,F,NE,NA,MG,AL,SI,P,S,CL,AR,K,CA,SC,TI,V,CR,MN,FE,CO,NI,&
       CU,ZN,GA,GE,AS,SE,BR,KR,RB,SR,Y,ZR,NB,MO,TC,RU,RH,PD,AG,CD,IN,SN, &
       SB,TE,I,XE,CS,BA,LA,CE,PR,ND,PM,SM,EU,GD,TB,DY,HO,ER,TM,YB,LU,HF,TA,W,RE,OS,&
       IR,PT,AU,HG,TL,PB,BI,PO,AT,RN,FR,RA,AC,TH,PA,U,NP,PU,AM,CM,BK,CF/)

end module scattering_factors
