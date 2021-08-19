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
module moduleVariables
#ifdef GPTA_MPI
    use mpi
#endif
#if GPTA_OMP
    use omp_lib
#endif
  use, intrinsic :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_INT, c_f_pointer, C_NULL_CHAR

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
  real(8) :: distanceScaling = 1.15d0
  integer :: randomNumberSeed = 12345
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
  real(8), parameter :: pi       = 3.1415926535898d0
  real(8), parameter :: pih      = 1.5707963267949d0
  real(8), parameter :: twopi    = 6.2831853071796d0
  real(8), parameter :: sqrpi    = 1.772453850905515d0
  real(8), parameter :: rbohr    = 0.529177249d0

  real(8), dimension(3,3) :: identityMatrix = reshape([1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0], [3,3])

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
    real(8) :: origin(3) = 0.d0
    real(8) :: cell(6)
    real(8) :: hmat(3,3)
    real(8) :: hinv(3,3)
    real(8) :: volume
    character(cp), allocatable, dimension(:) :: lab
    character(2), allocatable, dimension(:) :: element
    real(8), allocatable, dimension(:) :: mass
    real(8), allocatable, dimension(:,:) :: pos
    real(8), allocatable, dimension(:,:) :: frac
    real(8), allocatable, dimension(:) :: chg
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
    real(8), dimension(3)                    :: centreOfMass
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
      integer, intent(in) :: numberOfLocalMolecules
      integer, dimension(:), intent(in) :: list
      real(8), dimension(:), intent(out) :: property
    end subroutine molecularProperty

    subroutine systemProperty(numberOfProperties, property, numberOfSelectedAtoms, selectionList)
      integer, intent(inout) :: numberOfProperties
      real(8), dimension(*), intent(out) :: property
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
    real(8) :: cutoffNeighboursList

    real(8) :: timeTally = 0.d0
    integer :: tallyExecutions
    
    integer,               dimension(maxVariables) :: integerVariables = 0
    real(8),               dimension(maxVariables) :: doubleVariables  = 0.d0
    logical,               dimension(maxVariables) :: logicalVariables = .false.
    character(len=STRLEN), dimension(maxVariables) :: stringVariables  = ''
    character(cp),         dimension(maxVariables) :: labelVariables   = ''

    type(frameTypeDef)                         :: localFrame
    character(cp), allocatable, dimension(:)   :: localLabels
    real(8),       allocatable, dimension(:,:) :: localPositions
    real(8),       allocatable, dimension(:)   :: localCharges
    integer,       allocatable, dimension(:)   :: localIndices

    integer                                :: numberOfBins
    real(8), dimension(2)                  :: distributionLimits
    integer, allocatable, dimension(:)     :: arrayCounter
    real(8), allocatable, dimension(:)     :: array1D
    integer, dimension(2)                  :: numberOfBins2D
    real(8), allocatable, dimension(:,:)   :: array2D
    real(8), allocatable, dimension(:,:,:) :: array3D
    integer, allocatable, dimension(:,:,:) :: Iarray3D

    logical :: updateAtomsSelection = .true.
    logical, allocatable, dimension(:,:) :: isSelected
    integer, allocatable, dimension(:,:) :: idxSelection

    type(fileTypeDef) :: inputFile
    type(fileTypeDef) :: outputFile
    type(fileTypeDef) :: auxiliaryFile

    type(molecularPropertyProcedure) :: molprop

    type(systemPropertyProcedure) :: sysprop
    
  end type actionTypeDef

  interface 
    function cross_product(a,b) result(c)
      real(8), dimension (3), intent(in)  :: a
      real(8), dimension (3), intent(in)  :: b
      real(8), dimension (3) :: c
    end function cross_product
  end interface

  interface
    subroutine dumpXYZ(iounit,nat,pos,lab,hmat)
      implicit none
      integer, intent(in) :: iounit
      integer, intent(in) :: nat
      real(8), dimension(3,nat), intent(in) :: pos
      character(4), dimension(nat), optional, intent(in) :: lab
      real(8), dimension(3,3), optional, intent(in) :: hmat
    end subroutine
  end interface

contains

  function timing() result(t1)
    implicit none
    real(8) :: t1
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
