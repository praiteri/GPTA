!disclaimer
module moduleVariables
#ifdef GPTA_MPI
    use mpi
#endif
#if GPTA_OMP
    use omp_lib
#endif
  use, intrinsic :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_INT, c_f_pointer, C_NULL_CHAR
  use, intrinsic :: iso_fortran_env, only: real64

  implicit none

#ifdef GPTA_MPI
  integer :: ierr_mpi
#endif

  integer :: io = 6
  ! character(len=3) :: adv_char = "no "
  integer :: nProgress = 1000
  integer :: me 
  integer :: readingCPU = 0
  logical :: masterReadingOnly = .false.
  integer :: numberOfMpiProcesses 
  integer :: numberOfActiveProcesses
  integer :: MPI_Working_Group, MPI_world
  integer :: MPI_Working_Comm

  integer :: ompNumThreads = 0
  logical :: ldebug = .false.
  
  ! Variables tha can be redefined on the command line
  real(real64) :: distanceScaling = 1.15_real64
  integer :: randomNumberSeed = 0
  logical :: safeDistances = .false.
  logical :: forceVerletList = .false.

  ! String parsing
  integer, parameter :: STRLEN   = 1000
  integer, parameter :: MAXWORDS = 1000  ! Max words on a line
  integer, parameter :: MAXFILES = 100   ! Max words on a line

  ! Precision
  integer, parameter :: cp       = 4  ! character type for atoms' names
  integer, parameter :: fp       = 6  ! file type

  ! Pi and friends
  real(real64), parameter :: pi       = 3.1415926535898_real64
  real(real64), parameter :: pih      = 1.5707963267949_real64
  real(real64), parameter :: twopi    = 6.2831853071796_real64
  real(real64), parameter :: sqrpi    = 1.772453850905515_real64
  real(real64), parameter :: rbohr    = 0.529177249_real64

  real(real64), dimension(3,3) :: identityMatrix = reshape([1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 1.0_real64], [3,3])

  ! The data type located in libxdrfile for XTC files
  type, public, bind(C) :: xdrfile
    type(C_PTR) :: fp, xdr
    character(kind=C_CHAR) :: mode
    integer(C_INT) :: buf1, buf1size, buf2, buf2size
  end type xdrfile

  ! File type
  type fileTypeDef
    integer                :: funit
    character(len=STRLEN)  :: fname
    character(len=fp)      :: ftype
    character(len=11)      :: fform
    character(len=6)       :: fpos
    logical                :: first_access
    type(xdrfile), pointer :: xd
  end type fileTypeDef
  
  type frameTypeDef
    integer :: nframe
    integer :: natoms
    real(real64) :: origin(3) = 0.0_real64
    real(real64) :: cell(6)
    real(real64) :: hmat(3,3)
    real(real64) :: hinv(3,3)
    real(real64) :: volume
    character(cp), allocatable, dimension(:) :: lab
    character(2), allocatable, dimension(:) :: element
    real(real64), allocatable, dimension(:) :: mass
    real(real64), allocatable, dimension(:,:) :: pos
    real(real64), allocatable, dimension(:,:) :: frac
    real(real64), allocatable, dimension(:) :: chg
    ! real(real64), allocatable, dimension(:,:) :: vel
    real(real64), allocatable, dimension(:,:) :: forces

  end type frameTypeDef
  logical :: resetFrameLabels =.true.
  logical :: resetFrameCharges = .true.
  logical :: resetFrameElements = .true.

  type moleculeTypeDef
    integer                                  :: ID
    character(fp)                            :: resname
    integer                                  :: numberOfAtoms
    integer, allocatable, dimension(:)       :: listOfAtoms
    character(cp), allocatable, dimension(:) :: listOfLabels
  real(real64), dimension(3)                 :: centreOfMass
    logical                                  :: brokenBonds = .false.
  end type

  integer, parameter :: maxVariables = 200
  integer, parameter :: actionNameLength = 20
  integer, parameter :: maximimNumberOfActions = 100
  integer, parameter :: maximumActionStringLength = 1000
  character(len=maximumActionStringLength), dimension(0:maximimNumberOfActions) :: actionDetails
  
!!!!!!!!!!! Molecular properties
  abstract interface
    subroutine molecularProperty(numberOfLocalMolecules,list,property)
      use, intrinsic :: iso_fortran_env, only: real64
      integer, intent(in) :: numberOfLocalMolecules
      integer, dimension(:), intent(in) :: list
      real(real64), dimension(:), intent(out) :: property
    end subroutine molecularProperty

    subroutine systemProperty(numberOfProperties, property, numberOfSelectedAtoms, selectionList)
      use, intrinsic :: iso_fortran_env, only: real64
      integer, intent(inout) :: numberOfProperties
      real(real64), dimension(*), intent(out) :: property
      integer, intent(in), optional :: numberOfSelectedAtoms
      integer, dimension(:), intent(in), optional :: selectionList
    end subroutine systemProperty
  end interface

  type :: molecularPropertyProcedure
    character(len=STRLEN) :: command
    character(len=STRLEN) :: name
    procedure (molecularProperty), pointer, nopass :: property => null ()
  end type molecularPropertyProcedure

  type :: systemPropertyProcedure
    character(len=STRLEN) :: name
    procedure (systemProperty), pointer, nopass :: property => null ()
  end type systemPropertyProcedure

  type actionTypeDef
    character(len=STRLEN) :: name
    character(len=maximumActionStringLength) :: actionDetails
    logical :: actionInitialisation = .true.
    logical :: firstAction = .true.
    logical :: requiresNeighboursList = .false.
    logical :: requiresNeighboursListDouble = .false.
    logical :: requiresNeighboursListUpdates = .true.
    real(real64) :: cutoffNeighboursList

    real(real64) :: timeTally = 0.0_real64
    integer :: tallyExecutions
    
    integer,               dimension(maxVariables) :: integerVariables = 0
    real(real64),          dimension(maxVariables) :: doubleVariables  = 0.0_real64
    logical,               dimension(maxVariables) :: logicalVariables = .false.
    character(len=STRLEN), dimension(maxVariables) :: stringVariables  = ''
    character(cp),         dimension(maxVariables) :: labelVariables   = ''

    type(frameTypeDef)                         :: localFrame
    character(cp), allocatable, dimension(:)   :: localLabels
    real(real64),  allocatable, dimension(:,:) :: localPositions
    real(real64),  allocatable, dimension(:)   :: localCharges
    integer,       allocatable, dimension(:)   :: localIndices

    integer                                     :: numberOfBins
    real(real64), dimension(2)                  :: distributionLimits
    integer, allocatable, dimension(:)          :: arrayCounter

    real(real64), allocatable, dimension(:)     :: array1D
    real(real64), allocatable, dimension(:,:)   :: array2D
    real(real64), allocatable, dimension(:,:,:) :: array3D
    integer, dimension(2)                       :: numberOfBins2D

    integer, allocatable, dimension(:)     :: Iarray1D
    integer, allocatable, dimension(:,:)   :: Iarray2D
    integer, allocatable, dimension(:,:,:) :: Iarray3D

    logical :: updateAtomsSelection = .true.
    logical, allocatable, dimension(:,:) :: isSelected
    integer, allocatable, dimension(:,:) :: idxSelection
    integer, allocatable, dimension(:,:) :: idxToSelection

    type(fileTypeDef) :: inputFile
    type(fileTypeDef) :: outputFile
    type(fileTypeDef) :: auxiliaryFile

    type(molecularPropertyProcedure) :: molprop

    type(systemPropertyProcedure) :: sysprop
    
  end type actionTypeDef

  interface 
    function cross_product(a,b) result(c)
      use, intrinsic :: iso_fortran_env, only: real64
      real(real64), dimension (3), intent(in)  :: a
      real(real64), dimension (3), intent(in)  :: b
      real(real64), dimension (3) :: c
    end function cross_product
  end interface

  interface
    subroutine dumpXYZ(iounit,nat,pos,lab,hmat)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none
      integer, intent(in) :: iounit
      integer, intent(in) :: nat
      real(real64), dimension(3,nat), intent(in) :: pos
      character(4), dimension(nat), optional, intent(in) :: lab
      real(real64), dimension(3,3), optional, intent(in) :: hmat
    end subroutine
  end interface

contains

  function timing() result(t1)
    implicit none
    real(real64) :: t1
#ifdef GPTA_MPI
    t1 = MPI_wtime()
#elif GPTA_OMP
    t1 = OMP_GET_WTIME()
#else
    call cpu_time(t1)
#endif
    return
  end function timing

end module moduleVariables
