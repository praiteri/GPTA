!disclaimer
module scattering_factors
  use moduleVariables, only: real64
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
 type (SFactor), parameter :: DM = SFactor( 0,"DM",(/  0.00000_real64 ,  0.00000_real64 ,  0.00000_real64 , 0.00000_real64/) , & 
                                                   (/  0.00000_real64 ,  0.00000_real64 ,  0.00000_real64 ,    0.00000_real64/) ,  0.00000_real64 ,  0.0_real64)
 type (SFactor), parameter ::  H = SFactor( 1," H",(/  0.48992_real64 ,  0.26200_real64 ,  0.19677_real64 , 0.04988_real64/) , & 
                                                   (/ 20.65930_real64 ,  7.74039_real64 ,  49.55190_real64 ,   2.20159_real64/) ,  0.00130_real64 ,  1.0_real64)
 type (SFactor), parameter :: HE = SFactor( 2,"HE",(/  0.87340_real64 ,  0.63090_real64 ,  0.31120_real64 , 0.17800_real64/) , & 
                                                   (/  9.10370_real64 ,  3.35680_real64 ,  22.92760_real64 ,   0.98210_real64/) ,  0.00640_real64 ,  2.0_real64)
 type (SFactor), parameter :: LI = SFactor( 3,"LI",(/  1.12820_real64 ,  0.75080_real64 ,  0.61750_real64 , 0.46530_real64/) , & 
                                                   (/  3.95460_real64 ,  1.05240_real64 ,  85.39050_real64 , 168.26100_real64/) ,  0.03770_real64 ,  3.0_real64)
 type (SFactor), parameter :: BE = SFactor( 4,"BE",(/  1.59190_real64 ,  1.12780_real64 ,  0.53910_real64 , 0.70290_real64/) , & 
                                                   (/ 43.64270_real64 ,  1.86230_real64 , 103.48300_real64 ,   0.54200_real64/) ,  0.03850_real64 ,  4.0_real64)
 type (SFactor), parameter ::  B = SFactor( 5," B",(/  2.05450_real64 ,  1.33260_real64 ,  1.09790_real64 , 0.70680_real64/) , & 
                                                   (/ 23.21850_real64 ,  1.02100_real64 ,  60.34980_real64 ,   0.14030_real64/) ,  0.00000_real64 ,  5.0_real64)
 type (SFactor), parameter ::  C = SFactor( 6," C",(/  2.31000_real64 ,  1.02000_real64 ,  1.58860_real64 , 0.86500_real64/) , & 
                                                   (/ 20.84390_real64 , 10.20750_real64 ,   0.56870_real64 ,  51.65120_real64/) ,  0.21560_real64 ,  6.0_real64)
 type (SFactor), parameter ::  N = SFactor( 7," N",(/ 12.21260_real64 ,  3.13220_real64 ,  2.01250_real64 , 1.16630_real64/) , & 
                                                   (/  0.00570_real64 ,  9.89330_real64 ,  28.99750_real64 ,   0.58260_real64/) ,  0.00000_real64 ,  7.0_real64)
 type (SFactor), parameter ::  O = SFactor( 8," O",(/  3.04850_real64 ,  2.28680_real64 ,  1.54630_real64 , 0.86700_real64/) , & 
                                                   (/ 13.27710_real64 ,  5.70110_real64 ,   0.32390_real64 ,  32.90890_real64/) ,  0.25080_real64 ,  8.0_real64)
 type (SFactor), parameter ::  F = SFactor( 9," F",(/  3.53920_real64 ,  2.64120_real64 ,  1.51700_real64 , 1.02430_real64/) , & 
                                                   (/ 10.28250_real64 ,  4.29440_real64 ,   0.26150_real64 ,  26.14760_real64/) ,  0.27760_real64 ,  9.0_real64)
 type (SFactor), parameter :: NE = SFactor(10,"NE",(/  3.95530_real64 ,  3.11250_real64 ,  1.45460_real64 , 1.12510_real64/) , & 
                                                   (/  8.40420_real64 ,  3.42620_real64 ,   0.23060_real64 ,  21.71840_real64/) ,  0.35150_real64 , 10.0_real64)
 type (SFactor), parameter :: NA = SFactor(11,"NA",(/  4.76260_real64 ,  3.17360_real64 ,  1.26740_real64 , 1.11280_real64/) , & 
                                                   (/  3.28500_real64 ,  8.84220_real64 ,   0.31360_real64 , 129.42400_real64/) ,  0.67600_real64 , 11.0_real64)
 type (SFactor), parameter :: MG = SFactor(12,"MG",(/  5.42040_real64 ,  2.17350_real64 ,  1.22690_real64 , 2.30730_real64/) , & 
                                                   (/  2.82750_real64 , 79.26110_real64 ,   0.38080_real64 ,   7.19370_real64/) ,  0.85840_real64 , 12.0_real64)
 type (SFactor), parameter :: AL = SFactor(13,"AL",(/  6.42020_real64 ,  1.90020_real64 ,  1.59360_real64 , 1.96460_real64/) , & 
                                                   (/  3.03870_real64 ,  0.74260_real64 ,  31.54720_real64 ,  85.08860_real64/) ,  1.11510_real64 , 13.0_real64)
 type (SFactor), parameter :: SI = SFactor(14,"SI",(/  6.29150_real64 ,  3.03530_real64 ,  1.98910_real64 , 1.54100_real64/) , & 
                                                   (/  2.43860_real64 , 32.33370_real64 ,   0.67850_real64 ,  81.69370_real64/) ,  1.14070_real64 , 14.0_real64)
 type (SFactor), parameter ::  P = SFactor(15," P",(/  6.43450_real64 ,  4.17910_real64 ,  1.78000_real64 , 1.49080_real64/) , & 
                                                   (/  1.90670_real64 , 27.15700_real64 ,   0.52600_real64 ,  68.16450_real64/) ,  1.11490_real64 , 15.0_real64)
 type (SFactor), parameter ::  S = SFactor(16," S",(/  6.90530_real64 ,  5.20340_real64 ,  1.43790_real64 , 1.58630_real64/) , & 
                                                   (/  1.46790_real64 , 22.21510_real64 ,   0.25360_real64 ,  56.17200_real64/) ,  0.86690_real64 , 16.0_real64)
 type (SFactor), parameter :: CL = SFactor(17,"CL",(/ 11.46040_real64 ,  7.19640_real64 ,  6.25560_real64 , 1.64550_real64/) , & 
                                                   (/  0.01040_real64 ,  1.16620_real64 ,  18.51940_real64 ,  47.77840_real64/) ,  0.00000_real64 , 17.0_real64)
 type (SFactor), parameter :: AR = SFactor(18,"AR",(/  7.48450_real64 ,  6.77230_real64 ,  0.65390_real64 , 1.64420_real64/) , & 
                                                   (/  0.90720_real64 , 14.84070_real64 ,  43.89830_real64 ,  33.39290_real64/) ,  1.44450_real64 , 18.0_real64)
 type (SFactor), parameter ::  K = SFactor(19," K",(/  8.21860_real64 ,  7.43980_real64 ,  1.05190_real64 , 0.86590_real64/) , & 
                                                   (/ 12.79490_real64 ,  0.77480_real64 , 213.18700_real64 ,  41.68410_real64/) ,  1.42280_real64 , 19.0_real64)
 type (SFactor), parameter :: CA = SFactor(20,"CA",(/  8.62660_real64 ,  7.38730_real64 ,  1.58990_real64 , 1.02110_real64/) , & 
                                                   (/ 10.44210_real64 ,  0.65990_real64 ,  85.74840_real64 , 178.43700_real64/) ,  1.37510_real64 , 20.0_real64)
 type (SFactor), parameter :: SC = SFactor(21,"SC",(/  9.18900_real64 ,  7.36790_real64 ,  1.64090_real64 , 1.46800_real64/) , & 
                                                   (/  9.02130_real64 ,  0.57290_real64 , 136.10800_real64 ,  51.35310_real64/) ,  1.33290_real64 , 21.0_real64)
 type (SFactor), parameter :: TI = SFactor(22,"TI",(/  9.75950_real64 ,  7.35580_real64 ,  1.69910_real64 , 1.90210_real64/) , & 
                                                   (/  7.85080_real64 ,  0.50000_real64 ,  35.63380_real64 , 116.10500_real64/) ,  1.28070_real64 , 22.0_real64)
 type (SFactor), parameter ::  V = SFactor(23," V",(/ 10.29710_real64 ,  7.35110_real64 ,  2.07030_real64 , 2.05710_real64/) , & 
                                                   (/  6.86570_real64 ,  0.43850_real64 ,  26.89380_real64 , 102.47800_real64/) ,  1.21990_real64 , 23.0_real64)
 type (SFactor), parameter :: CR = SFactor(24,"CR",(/ 10.64060_real64 ,  7.35370_real64 ,  3.32400_real64 , 1.49220_real64/) , & 
                                                   (/  6.10380_real64 ,  0.39200_real64 ,  20.26260_real64 ,  98.73990_real64/) ,  1.18320_real64 , 24.0_real64)
 type (SFactor), parameter :: MN = SFactor(25,"MN",(/ 11.28190_real64 ,  7.35730_real64 ,  3.01930_real64 , 2.24410_real64/) , & 
                                                   (/  5.34090_real64 ,  0.34320_real64 ,  17.86740_real64 ,  83.75430_real64/) ,  1.08960_real64 , 25.0_real64)
 type (SFactor), parameter :: FE = SFactor(26,"FE",(/ 11.76950_real64 ,  7.35730_real64 ,  3.52220_real64 , 2.30450_real64/) , & 
                                                   (/  4.76110_real64 ,  0.30720_real64 ,  15.35350_real64 ,  76.88050_real64/) ,  1.03690_real64 , 26.0_real64)
 type (SFactor), parameter :: CO = SFactor(27,"CO",(/ 12.28410_real64 ,  7.34090_real64 ,  4.00340_real64 , 2.34880_real64/) , & 
                                                   (/  4.27910_real64 ,  0.27840_real64 ,  13.53590_real64 ,  71.16920_real64/) ,  1.01180_real64 , 27.0_real64)
 type (SFactor), parameter :: NI = SFactor(28,"NI",(/ 12.83760_real64 ,  7.29200_real64 ,  4.44380_real64 , 2.38000_real64/) , & 
                                                   (/  3.87850_real64 ,  0.25650_real64 ,  12.17630_real64 ,  66.34210_real64/) ,  1.03410_real64 , 28.0_real64)
 type (SFactor), parameter :: CU = SFactor(29,"CU",(/ 13.33800_real64 ,  7.16760_real64 ,  5.61580_real64 , 1.67350_real64/) , & 
                                                   (/  3.58280_real64 ,  0.24700_real64 ,  11.39660_real64 ,  64.81260_real64/) ,  1.19100_real64 , 29.0_real64)
 type (SFactor), parameter :: ZN = SFactor(30,"ZN",(/ 14.07430_real64 ,  7.03180_real64 ,  5.16520_real64 , 2.41000_real64/) , & 
                                                   (/  3.26550_real64 ,  0.23330_real64 ,  10.31630_real64 ,  58.70970_real64/) ,  1.30410_real64 , 30.0_real64)
 type (SFactor), parameter :: GA = SFactor(31,"GA",(/ 15.23540_real64 ,  6.70060_real64 ,  4.35910_real64 , 2.96230_real64/) , & 
                                                   (/  3.06690_real64 ,  0.24120_real64 ,  10.78050_real64 ,  61.41350_real64/) ,  1.71890_real64 , 31.0_real64)
 type (SFactor), parameter :: GE = SFactor(32,"GE",(/ 16.08160_real64 ,  6.37470_real64 ,  3.70680_real64 , 3.68300_real64/) , & 
                                                   (/  2.85090_real64 ,  0.25160_real64 ,  11.44680_real64 ,  54.76250_real64/) ,  2.13130_real64 , 32.0_real64)
 type (SFactor), parameter :: AS = SFactor(33,"AS",(/ 16.67230_real64 ,  6.07010_real64 ,  3.43130_real64 , 4.27790_real64/) , & 
                                                   (/  2.63450_real64 ,  0.26470_real64 ,  12.94790_real64 ,  47.79720_real64/) ,  2.53100_real64 , 33.0_real64)
 type (SFactor), parameter :: SE = SFactor(34,"SE",(/ 17.00060_real64 ,  5.81960_real64 ,  3.97310_real64 , 4.35430_real64/) , & 
                                                   (/  2.40980_real64 ,  0.27260_real64 ,  15.23720_real64 ,  43.81630_real64/) ,  2.84090_real64 , 34.0_real64)
 type (SFactor), parameter :: BR = SFactor(35,"BR",(/ 17.17890_real64 ,  5.23580_real64 ,  5.63770_real64 , 3.98510_real64/) , & 
                                                   (/  2.17230_real64 , 16.57960_real64 ,   0.26090_real64 ,  41.43280_real64/) ,  2.95570_real64 , 35.0_real64)
 type (SFactor), parameter :: KR = SFactor(36,"KR",(/ 17.35550_real64 ,  6.72860_real64 ,  5.54930_real64 , 3.53750_real64/) , & 
                                                   (/  1.93840_real64 , 16.56230_real64 ,   0.22610_real64 ,  39.39720_real64/) ,  2.82500_real64 , 36.0_real64)
 type (SFactor), parameter :: RB = SFactor(37,"RB",(/ 17.17840_real64 ,  9.64350_real64 ,  5.13990_real64 , 1.52920_real64/) , & 
                                                   (/  1.78880_real64 , 17.31510_real64 ,   0.27480_real64 , 164.93400_real64/) ,  3.48730_real64 , 37.0_real64)
 type (SFactor), parameter :: SR = SFactor(38,"SR",(/ 17.56630_real64 ,  9.81840_real64 ,  5.42200_real64 , 2.66940_real64/) , & 
                                                   (/  1.55640_real64 , 14.09880_real64 ,   0.16640_real64 , 132.37600_real64/) ,  2.50640_real64 , 38.0_real64)
 type (SFactor), parameter ::  Y = SFactor(39," Y",(/ 17.77600_real64 , 10.29460_real64 ,  5.72629_real64 , 3.26588_real64/) , & 
                                                   (/  1.40290_real64 , 12.80060_real64 ,   0.12560_real64 , 104.35400_real64/) ,  1.91213_real64 , 39.0_real64)
 type (SFactor), parameter :: ZR = SFactor(40,"ZR",(/ 17.87650_real64 , 10.94800_real64 ,  5.41732_real64 , 3.65721_real64/) , & 
                                                   (/  1.27618_real64 , 11.91600_real64 ,   0.11762_real64 ,  87.66270_real64/) ,  2.06929_real64 , 40.0_real64)
 type (SFactor), parameter :: NB = SFactor(41,"NB",(/ 17.61420_real64 , 12.01440_real64 ,  4.04183_real64 , 3.53346_real64/) , & 
                                                   (/  1.18865_real64 , 11.76600_real64 ,   0.20478_real64 ,  69.79570_real64/) ,  3.75591_real64 , 41.0_real64)
 type (SFactor), parameter :: MO = SFactor(42,"MO",(/  3.70250_real64 , 17.23560_real64 , 12.88760_real64 , 3.74290_real64/) , & 
                                                   (/  0.27720_real64 ,  1.09580_real64 ,  11.00400_real64 ,  61.65840_real64/) ,  4.38750_real64 , 42.0_real64)
 type (SFactor), parameter :: TC = SFactor(43,"TC",(/ 19.13010_real64 , 11.09480_real64 ,  4.64901_real64 , 2.71263_real64/) , & 
                                                   (/  0.86413_real64 ,  8.14487_real64 ,  21.57070_real64 ,  86.84720_real64/) ,  5.40428_real64 , 43.0_real64)
 type (SFactor), parameter :: RU = SFactor(44,"RU",(/ 19.26740_real64 , 12.91820_real64 ,  4.86337_real64 , 1.56756_real64/) , & 
                                                   (/  0.80852_real64 ,  8.43467_real64 ,  24.79970_real64 ,  94.29280_real64/) ,  5.37874_real64 , 44.0_real64)
 type (SFactor), parameter :: RH = SFactor(45,"RH",(/ 19.29570_real64 , 14.35010_real64 ,  4.73425_real64 , 1.28918_real64/) , & 
                                                   (/  0.75154_real64 ,  8.21758_real64 ,  25.87490_real64 ,  98.60620_real64/) ,  5.32800_real64 , 45.0_real64)
 type (SFactor), parameter :: PD = SFactor(46,"PD",(/ 19.33190_real64 , 15.50170_real64 ,  5.29537_real64 , 0.60584_real64/) , & 
                                                   (/  0.69866_real64 ,  7.98929_real64 ,  25.20520_real64 ,  76.89860_real64/) ,  5.26593_real64 , 46.0_real64)
 type (SFactor), parameter :: AG = SFactor(47,"AG",(/ 19.28080_real64 , 16.68850_real64 ,  4.80450_real64 , 1.04630_real64/) , & 
                                                   (/  0.64460_real64 ,  7.47260_real64 ,  24.66050_real64 ,  99.81560_real64/) ,  5.17900_real64 , 47.0_real64)
 type (SFactor), parameter :: CD = SFactor(48,"CD",(/ 19.22140_real64 , 17.64440_real64 ,  4.46100_real64 , 1.60290_real64/) , & 
                                                   (/  0.59460_real64 ,  6.90890_real64 ,  24.70080_real64 ,  87.48250_real64/) ,  5.06940_real64 , 48.0_real64)
 type (SFactor), parameter :: IN = SFactor(49,"IN",(/ 19.16240_real64 , 18.55960_real64 ,  4.29480_real64 , 2.03960_real64/) , & 
                                                   (/  0.54760_real64 ,  6.37760_real64 ,  25.84990_real64 ,  92.80290_real64/) ,  4.93910_real64 , 49.0_real64)
 type (SFactor), parameter :: SN = SFactor(50,"SN",(/ 19.18890_real64 , 19.10050_real64 ,  4.45850_real64 , 2.46630_real64/) , & 
                                                   (/  5.83030_real64 ,  0.50310_real64 ,  26.89090_real64 ,  83.95710_real64/) ,  4.78210_real64 , 50.0_real64)
 type (SFactor), parameter :: SB = SFactor(51,"SB",(/ 19.64180_real64 , 19.04550_real64 ,  5.03710_real64 , 2.68270_real64/) , & 
                                                   (/  5.30340_real64 ,  0.46070_real64 ,  27.90740_real64 ,  75.28250_real64/) ,  4.59090_real64 , 51.0_real64)
 type (SFactor), parameter :: TE = SFactor(52,"TE",(/ 19.96440_real64 , 19.01380_real64 ,  6.14487_real64 , 2.52390_real64/) , & 
                                                   (/  4.81742_real64 ,  0.42089_real64 ,  28.52840_real64 ,  70.84030_real64/) ,  4.35200_real64 , 52.0_real64)
 type (SFactor), parameter ::  I = SFactor(53," I",(/ 20.14720_real64 , 18.99490_real64 ,  7.51380_real64 , 2.27350_real64/) , & 
                                                   (/  4.34700_real64 ,  0.38140_real64 ,  27.76600_real64 ,  66.87760_real64/) ,  4.07120_real64 , 53.0_real64)
 type (SFactor), parameter :: XE = SFactor(54,"XE",(/ 20.29330_real64 , 19.02980_real64 ,  8.97670_real64 , 1.99000_real64/) , & 
                                                   (/  3.92820_real64 ,  0.34400_real64 ,  26.46590_real64 ,  64.26580_real64/) ,  3.71180_real64 , 54.0_real64)
 type (SFactor), parameter :: CS = SFactor(55,"CS",(/ 20.38920_real64 , 19.10620_real64 , 10.66200_real64 , 1.49530_real64/) , & 
                                                   (/  3.56900_real64 ,  0.31070_real64 ,  24.38790_real64 , 213.90400_real64/) ,  3.33520_real64 , 55.0_real64)
 type (SFactor), parameter :: BA = SFactor(56,"BA",(/ 20.33610_real64 , 19.29700_real64 , 10.88800_real64 , 2.69590_real64/) , & 
                                                   (/  3.21600_real64 ,  0.27560_real64 ,  20.20730_real64 , 167.20200_real64/) ,  2.77310_real64 , 56.0_real64)
 type (SFactor), parameter :: LA = SFactor(57,"LA",(/ 20.57800_real64 , 19.59900_real64 , 11.37270_real64 , 3.28719_real64/) , & 
                                                   (/  2.94817_real64 ,  0.24447_real64 ,  18.77260_real64 , 133.12400_real64/) ,  2.14678_real64 , 57.0_real64)
 type (SFactor), parameter :: CE = SFactor(58,"CE",(/ 21.16710_real64 , 19.76950_real64 , 11.85130_real64 , 3.33049_real64/) , & 
                                                   (/  2.81219_real64 ,  0.22684_real64 ,  17.60830_real64 , 127.11300_real64/) ,  1.86264_real64 , 58.0_real64)
 type (SFactor), parameter :: PR = SFactor(59,"PR",(/ 22.04400_real64 , 19.66970_real64 , 12.38560_real64 , 2.82428_real64/) , & 
                                                   (/  2.77393_real64 ,  0.22209_real64 ,  16.76690_real64 , 143.64400_real64/) ,  2.05830_real64 , 59.0_real64)
 type (SFactor), parameter :: ND = SFactor(60,"ND",(/ 22.68450_real64 , 19.68470_real64 , 12.77400_real64 , 2.85137_real64/) , & 
                                                   (/  2.66248_real64 ,  0.21063_real64 ,  15.88500_real64 , 137.90300_real64/) ,  1.98486_real64 , 60.0_real64)
 type (SFactor), parameter :: PM = SFactor(61,"PM",(/ 23.34050_real64 , 19.60950_real64 , 13.12350_real64 , 2.87516_real64/) , & 
                                                   (/  2.56270_real64 ,  0.20209_real64 ,  15.10090_real64 , 132.72100_real64/) ,  2.02876_real64 , 61.0_real64)
 type (SFactor), parameter :: SM = SFactor(62,"SM",(/ 24.00420_real64 , 19.42580_real64 , 13.43960_real64 , 2.89604_real64/) , & 
                                                   (/  2.47274_real64 ,  0.19645_real64 ,  14.39960_real64 , 128.00700_real64/) ,  2.20963_real64 , 62.0_real64)
 type (SFactor), parameter :: EU = SFactor(63,"EU",(/ 24.62740_real64 , 19.08860_real64 , 13.76030_real64 , 2.92270_real64/) , & 
                                                   (/  2.38790_real64 ,  0.19420_real64 ,  13.75460_real64 , 123.17400_real64/) ,  2.57450_real64 , 63.0_real64)
 type (SFactor), parameter :: GD = SFactor(64,"GD",(/ 25.07090_real64 , 19.07980_real64 , 13.85180_real64 , 3.54545_real64/) , & 
                                                   (/  2.25341_real64 ,  0.18195_real64 ,  12.93310_real64 , 101.39800_real64/) ,  2.41960_real64 , 64.0_real64)
 type (SFactor), parameter :: TB = SFactor(65,"TB",(/ 25.89760_real64 , 18.21850_real64 , 14.31670_real64 , 2.95354_real64/) , & 
                                                   (/  2.24256_real64 ,  0.19614_real64 ,  12.66480_real64 , 115.36200_real64/) ,  3.58324_real64 , 65.0_real64)
 type (SFactor), parameter :: DY = SFactor(66,"DY",(/ 26.50700_real64 , 17.63830_real64 , 14.55960_real64 , 2.96577_real64/) , & 
                                                   (/  2.18020_real64 ,  0.20217_real64 ,  12.18990_real64 , 111.87400_real64/) ,  4.29728_real64 , 66.0_real64)
 type (SFactor), parameter :: HO = SFactor(67,"HO",(/ 26.90490_real64 , 17.29400_real64 , 14.55830_real64 , 3.63837_real64/) , & 
                                                   (/  2.07051_real64 ,  0.19794_real64 ,  11.44070_real64 ,  92.65660_real64/) ,  4.56796_real64 , 67.0_real64)
 type (SFactor), parameter :: ER = SFactor(68,"ER",(/ 27.65630_real64 , 16.42850_real64 , 14.97790_real64 , 2.98233_real64/) , & 
                                                   (/  2.07356_real64 ,  0.22354_real64 ,  11.36040_real64 , 105.70300_real64/) ,  5.92046_real64 , 68.0_real64)
 type (SFactor), parameter :: TM = SFactor(69,"TM",(/ 28.18190_real64 , 15.88510_real64 , 15.15420_real64 , 2.98706_real64/) , & 
                                                   (/  2.02859_real64 ,  0.23885_real64 ,  10.99750_real64 , 102.96100_real64/) ,  6.75621_real64 , 69.0_real64)
 type (SFactor), parameter :: YB = SFactor(70,"YB",(/ 28.66410_real64 , 15.43450_real64 , 15.30870_real64 , 2.98963_real64/) , & 
                                                   (/  1.98890_real64 ,  0.25712_real64 ,  10.66470_real64 , 100.41700_real64/) ,  7.56672_real64 , 70.0_real64)
 type (SFactor), parameter :: LU = SFactor(71,"LU",(/ 28.94760_real64 , 15.22080_real64 , 15.10000_real64 , 3.71601_real64/) , & 
                                                   (/  1.90182_real64 ,  9.98519_real64 ,   0.26103_real64 ,  84.32980_real64/) ,  7.97628_real64 , 71.0_real64)
 type (SFactor), parameter :: HF = SFactor(72,"HF",(/ 29.14400_real64 , 15.17260_real64 , 14.75860_real64 , 4.30013_real64/) , & 
                                                   (/  1.83262_real64 ,  9.59990_real64 ,   0.27512_real64 ,  72.02900_real64/) ,  8.58154_real64 , 72.0_real64)
 type (SFactor), parameter :: TA = SFactor(73,"TA",(/ 29.20240_real64 , 15.22930_real64 , 14.51350_real64 , 4.76492_real64/) , & 
                                                   (/  1.77333_real64 ,  9.37046_real64 ,   0.29598_real64 ,  63.36440_real64/) ,  9.24354_real64 , 73.0_real64)
 type (SFactor), parameter ::  W = SFactor(74," W",(/ 29.08180_real64 , 15.43000_real64 , 14.43270_real64 , 5.11982_real64/) , & 
                                                   (/  1.72029_real64 ,  9.22590_real64 ,   0.32170_real64 ,  57.05600_real64/) ,  9.88750_real64 , 74.0_real64)
 type (SFactor), parameter :: RE = SFactor(75,"RE",(/ 28.76210_real64 , 15.71890_real64 , 14.55640_real64 , 5.44174_real64/) , & 
                                                   (/  1.67191_real64 ,  9.09227_real64 ,   0.35050_real64 ,  52.08610_real64/) , 10.47200_real64 , 75.0_real64)
 type (SFactor), parameter :: OS = SFactor(76,"OS",(/ 28.18940_real64 , 16.15500_real64 , 14.93050_real64 , 5.67589_real64/) , & 
                                                   (/  1.62903_real64 ,  8.97948_real64 ,   0.38266_real64 ,  48.16470_real64/) , 11.00050_real64 , 76.0_real64)
 type (SFactor), parameter :: IR = SFactor(77,"IR",(/ 27.30490_real64 , 16.72960_real64 , 15.61150_real64 , 5.83377_real64/) , & 
                                                   (/  1.59279_real64 ,  8.86553_real64 ,   0.41792_real64 ,  45.00110_real64/) , 11.47220_real64 , 77.0_real64)
 type (SFactor), parameter :: PT = SFactor(78,"PT",(/ 27.00590_real64 , 17.76390_real64 , 15.71310_real64 , 5.78370_real64/) , & 
                                                   (/  1.51293_real64 ,  8.81174_real64 ,   0.42459_real64 ,  38.61030_real64/) , 11.68830_real64 , 78.0_real64)
 type (SFactor), parameter :: AU = SFactor(79,"AU",(/ 16.88190_real64 , 18.59130_real64 , 25.55820_real64 , 5.86000_real64/) , & 
                                                   (/  0.46110_real64 ,  8.62160_real64 ,   1.48260_real64 ,  36.39560_real64/) , 12.06580_real64 , 79.0_real64)
 type (SFactor), parameter :: HG = SFactor(80,"HG",(/ 20.68090_real64 , 19.04170_real64 , 21.65750_real64 , 5.96760_real64/) , & 
                                                   (/  0.54500_real64 ,  8.44840_real64 ,   1.57290_real64 ,  38.32460_real64/) , 12.60890_real64 , 80.0_real64)
 type (SFactor), parameter :: TL = SFactor(81,"TL",(/ 27.54460_real64 , 19.15840_real64 , 15.53800_real64 , 5.52593_real64/) , & 
                                                   (/  0.65515_real64 ,  8.70751_real64 ,   1.96347_real64 ,  45.81490_real64/) , 13.17460_real64 , 81.0_real64)
 type (SFactor), parameter :: PB = SFactor(82,"PB",(/ 31.06170_real64 , 13.06370_real64 , 18.44200_real64 , 5.96960_real64/) , & 
                                                   (/  0.69020_real64 ,  2.35760_real64 ,   8.61800_real64 ,  47.25790_real64/) , 13.41180_real64 , 82.0_real64)
 type (SFactor), parameter :: BI = SFactor(83,"BI",(/ 33.36890_real64 , 12.95100_real64 , 16.58770_real64 , 6.46920_real64/) , & 
                                                   (/  0.70400_real64 ,  2.92380_real64 ,   8.79370_real64 ,  48.00930_real64/) , 13.57820_real64 , 83.0_real64)
 type (SFactor), parameter :: PO = SFactor(84,"PO",(/ 34.67260_real64 , 15.47330_real64 , 13.11380_real64 , 7.02588_real64/) , & 
                                                   (/  0.70100_real64 ,  3.55078_real64 ,   9.55642_real64 ,  47.00450_real64/) , 13.67700_real64 , 84.0_real64)
 type (SFactor), parameter :: AT = SFactor(85,"AT",(/ 35.31630_real64 , 19.02110_real64 ,  9.49887_real64 , 7.42518_real64/) , & 
                                                   (/  0.68587_real64 ,  3.97458_real64 ,  11.38240_real64 ,  45.47150_real64/) , 13.71080_real64 , 85.0_real64)
 type (SFactor), parameter :: RN = SFactor(86,"RN",(/ 35.56310_real64 , 21.28160_real64 ,  8.00370_real64 , 7.44330_real64/) , & 
                                                   (/  0.66310_real64 ,  4.06910_real64 ,  14.04220_real64 ,  44.24730_real64/) , 13.69050_real64 , 86.0_real64)
 type (SFactor), parameter :: FR = SFactor(87,"FR",(/ 35.92990_real64 , 23.05470_real64 , 12.14390_real64 , 2.11253_real64/) , & 
                                                   (/  0.64645_real64 ,  4.17619_real64 ,  23.10520_real64 , 150.64500_real64/) , 13.72470_real64 , 87.0_real64)
 type (SFactor), parameter :: RA = SFactor(88,"RA",(/ 35.76300_real64 , 22.90640_real64 , 12.47390_real64 , 3.21097_real64/) , & 
                                                   (/  0.61634_real64 ,  3.87135_real64 ,  19.98870_real64 , 142.32500_real64/) , 13.62110_real64 , 88.0_real64)
 type (SFactor), parameter :: AC = SFactor(89,"AC",(/ 35.65970_real64 , 23.10320_real64 , 12.59770_real64 , 4.08655_real64/) , & 
                                                   (/  0.58909_real64 ,  3.65155_real64 ,  18.59900_real64 , 117.02000_real64/) , 13.52660_real64 , 89.0_real64)
 type (SFactor), parameter :: TH = SFactor(90,"TH",(/ 35.56450_real64 , 23.42190_real64 , 12.74730_real64 , 4.80703_real64/) , & 
                                                   (/  0.56336_real64 ,  3.46204_real64 ,  17.83090_real64 ,  99.17220_real64/) , 13.43140_real64 , 90.0_real64)
 type (SFactor), parameter :: PA = SFactor(91,"PA",(/ 35.88470_real64 , 23.29480_real64 , 14.18910_real64 , 4.17287_real64/) , & 
                                                   (/  0.54775_real64 ,  3.41519_real64 ,  16.92350_real64 , 105.25100_real64/) , 13.42870_real64 , 91.0_real64)
 type (SFactor), parameter ::  U = SFactor(92," U",(/ 36.02280_real64 , 23.41280_real64 , 14.94910_real64 , 4.18800_real64/) , & 
                                                   (/  0.52930_real64 ,  3.32530_real64 ,  16.09270_real64 , 100.61300_real64/) , 13.39660_real64 , 92.0_real64)
 type (SFactor), parameter :: NP = SFactor(93,"NP",(/ 36.18740_real64 , 23.59640_real64 , 15.64020_real64 , 4.18550_real64/) , & 
                                                   (/  0.51193_real64 ,  3.25396_real64 ,  15.36220_real64 ,  97.49080_real64/) , 13.35730_real64 , 93.0_real64)
 type (SFactor), parameter :: PU = SFactor(94,"PU",(/ 36.52540_real64 , 23.80830_real64 , 16.77070_real64 , 3.47947_real64/) , & 
                                                   (/  0.49938_real64 ,  3.26371_real64 ,  14.94550_real64 , 105.98000_real64/) , 13.38120_real64 , 94.0_real64)
 type (SFactor), parameter :: AM = SFactor(95,"AM",(/ 36.67060_real64 , 24.09920_real64 , 17.34150_real64 , 3.49331_real64/) , & 
                                                   (/  0.48363_real64 ,  3.20647_real64 ,  14.31360_real64 , 102.27300_real64/) , 13.35920_real64 , 95.0_real64)
 type (SFactor), parameter :: CM = SFactor(96,"CM",(/ 36.64880_real64 , 24.40960_real64 , 17.39900_real64 , 4.21665_real64/) , & 
                                                   (/  0.46515_real64 ,  3.08997_real64 ,  13.43460_real64 ,  88.48340_real64/) , 13.28870_real64 , 96.0_real64)
 type (SFactor), parameter :: BK = SFactor(97,"BK",(/ 36.78810_real64 , 24.77360_real64 , 17.89190_real64 , 4.23284_real64/) , & 
                                                   (/  0.45102_real64 ,  3.04619_real64 ,  12.89460_real64 ,  86.00300_real64/) , 13.27540_real64 , 97.0_real64)
 type (SFactor), parameter :: CF = SFactor(98,"CF",(/ 36.91850_real64 , 25.19950_real64 , 18.33170_real64 , 4.24391_real64/) , & 
                                                   (/  0.43753_real64 ,  3.00775_real64 ,  12.40440_real64 ,  83.78810_real64/) , 13.26740_real64 , 98.0_real64)

 integer, parameter  :: SF_Library_dim = 98

 type(SFactor), parameter, dimension(0:SF_Library_dim)  :: SF_Library = &
     (/DM,H,HE,LI,BE,B,C,N,O,F,NE,NA,MG,AL,SI,P,S,CL,AR,K,CA,SC,TI,V,CR,MN,FE,CO,NI,&
       CU,ZN,GA,GE,AS,SE,BR,KR,RB,SR,Y,ZR,NB,MO,TC,RU,RH,PD,AG,CD,IN,SN, &
       SB,TE,I,XE,CS,BA,LA,CE,PR,ND,PM,SM,EU,GD,TB,DY,HO,ER,TM,YB,LU,HF,TA,W,RE,OS,&
       IR,PT,AU,HG,TL,PB,BI,PO,AT,RN,FR,RA,AC,TH,PA,U,NP,PU,AM,CM,BK,CF/)

end module scattering_factors
