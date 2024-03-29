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
module moduleProperties
  use moduleSystem
  use moduleResizeArrays
  use moduleVariables
  implicit none

  public :: workData, associateProperties
  private

  type :: propertyProcedure
    procedure (initialisePropertyRoutine), pointer, nopass :: initialise => null ()
    procedure (computePropertyRoutine), pointer, nopass :: compute => null ()
    procedure (dumpPropertyRoutine), pointer, nopass :: dump => null ()
    procedure (extractPropertyRoutine), pointer, nopass :: extract => null ()
  end type propertyProcedure
  type (propertyProcedure), save :: workData
  
  type propertyType
    character(len=20)                    :: type
    integer                              :: counter
    integer                              :: numberOfBins
    integer, dimension(2)                :: numberOfBins2D
    real(8), allocatable                 :: lowerLimit
    real(8), allocatable                 :: upperLimit
    real(8),              dimension(2)   :: lowerLimits
    real(8),              dimension(2)   :: upperLimits

    real(8)                              :: average
    real(8)                              :: average2
    real(8), allocatable, dimension(:)   :: dist1D
    real(8), allocatable, dimension(:,:) :: dist2D
    integer, allocatable, dimension(:)   :: tally

    integer                              :: numberOfValues
    integer                              :: iounit

  end type propertyType

  integer :: numberOfProperties = 0
  type(propertyType), dimension(1000) :: property

contains

  subroutine initialisePropertyRoutine(ID, whatToDo, lowerLimits, upperLimits, numberOfBins, numberOfValues, iounit)
    integer,               intent(inout)        :: ID
    character(len=*),      intent(in)           :: whatToDo
    real(8), dimension(:), intent(in), optional :: lowerLimits  
    real(8), dimension(:), intent(in), optional :: upperLimits  
    integer, dimension(:), intent(in), optional :: numberOfBins
    integer,               intent(in), optional :: numberOfValues
    integer,               intent(in), optional :: iounit

    if (ID <= 0) then
      numberOfProperties = numberOfProperties + 1
      ID = numberOfProperties
    end if

    property(ID) % counter = 0
    property(ID) % type = trim(whatToDo)

    if (whatToDo == "dump" .or. whatToDo == "dumpSeq" ) then
      property(ID) % iounit = iounit
      
    else if (whatToDo == "store") then
      allocate(property(ID) % dist1D(property(ID) % numberOfBins), source=0.d0)

    else if (whatToDo == "average") then
        property(ID) % average = 0.d0
        property(ID) % average2 = 0.d0

    else if (whatToDo == "multiAverage") then
      property % numberOfValues = numberOfValues
      allocate(property(ID) % dist1D(2 * property(ID) % numberOfValues), source=0.d0)

    else if (whatToDo == "histogram") then
      allocate(property(ID) % dist1D(numberOfBins(1)), source=0.d0)
      property(ID) % numberOfBins = numberOfBins(1)
      property(ID) % lowerLimit   = lowerLimits(1)
      property(ID) % upperLimit   = upperLimits(1)

    else if (whatToDo == "dist2D" .or. whatToDo == "prob2D") then
      allocate(property(ID) % dist2D(numberOfBins(1),numberOfBins(2)), source=0.d0)
      property(ID) % numberOfBins2D   = numberOfBins(1:2)
      property(ID) % lowerLimits(1:2) = lowerLimits(1:2)
      property(ID) % upperLimits(1:2) = upperLimits(1:2)

    else if (whatToDo == "avgDist") then
      allocate(property(ID) % dist1D(numberOfBins(1)), source=0.d0)
      allocate(property(ID) % tally(numberOfBins(1)), source=0)
      property(ID) % numberOfBins = numberOfBins(1)
      property(ID) % lowerLimit   = lowerLimits(1)
      property(ID) % upperLimit   = upperLimits(1)

    end if

  end subroutine initialisePropertyRoutine

  subroutine computePropertyRoutine(ID, numberOfValues, xValues, yValues)
    integer,               intent(in)           :: ID
    integer,               intent(in)           :: numberOfValues
    real(8), dimension(:), intent(in), optional :: xValues
    real(8), dimension(:), intent(in), optional :: yValues

    integer :: i, idx, jdx
    integer :: nbin, mbin
    real(8) :: dx, xmin, xmax
    real(8) :: dy, ymin, ymax

    if (property(ID) % type == "dump") then
      property(ID) % counter = property(ID) % counter + 1
      write(property(ID) % iounit,'(i6,1x)',advance='no') property(ID) % counter
      do idx=1,numberOfValues
        write(property(ID) % iounit,'(f20.10)',advance='no') xValues(idx)
      end do
      write(property(ID) % iounit,*)
      call flush(property(ID) % iounit)

    else if (property(ID) % type == "dumpSeq") then
      property(ID) % counter = property(ID) % counter + 1
      write(property(ID) % iounit,'("#",2x,i0)') property(ID) % counter
      do idx=1,numberOfValues
        write(property(ID) % iounit,'(f20.10)') xValues(idx)
      end do
      write(property(ID) % iounit,*)
      write(property(ID) % iounit,*)
      call flush(property(ID) % iounit)

    else if (property(ID) % type == "store") then
      idx = property(ID) % counter 
      if (idx+numberOfValues > property(ID) % numberOfBins) then
        property(ID) % numberOfBins = property(ID) % numberOfBins + 1000
        call resizeArray(property(ID) % dist1D , property(ID) % numberOfBins)
      endif
      do i=1,numberOfValues
        property(ID) % dist1D(idx+i) = xValues(i)
      end do
      property(ID) % counter = property(ID) % counter + numberOfValues
      
    else if (property(ID) % type == "average") then
        do i=1,numberOfValues
          property(ID) % average = property(ID) % average + xValues(i)
          property(ID) % average2 = property(ID) % average2 + xValues(i)**2
        end do
        property(ID) % counter = property(ID) % counter + numberOfValues

    else if (property(ID) % type == "multiAverage") then
        idx = 0
        do i=1,numberOfValues
          property(ID) % dist1D(idx+1) = property(ID) % dist1D(idx+1) + xValues(i)
          property(ID) % dist1D(idx+2) = property(ID) % dist1D(idx+2) + xValues(i)**2
          idx = idx + 2
        end do
        property(ID) % counter = property(ID) % counter + 1

    else if (property(ID) % type == "histogram") then
      xmin = property(ID) % lowerLimit
      xmax = property(ID) % upperLimit
      nbin = property(ID) % numberOfBins
      dx = (xmax-xmin) / dble(nbin)
      do i=1,numberOfValues
        idx = int( (xValues(i) - xmin) / dx) + 1
        if (idx <= 0 .or. idx > nbin) cycle
        property(ID) % dist1D(idx) = property(ID) % dist1D(idx) + 1
      end do
      property(ID) % counter = property(ID) % counter + 1

    else if (property(ID) % type == "dist2D" .or. property(ID) % type == "prob2D") then
      xmin = property(ID) % lowerLimits(1)
      xmax = property(ID) % upperLimits(1)
      nbin = property(ID) % numberOfBins2D(1)
      dx = (xmax-xmin) / dble(nbin)
      ymin = property(ID) % lowerLimits(2)
      ymax = property(ID) % upperLimits(2)
      mbin = property(ID) % numberOfBins2D(2)
      dy = (ymax-ymin) / dble(mbin)
      do i=1,numberOfValues
        idx = int( (xValues(i) - xmin) / dx) + 1
        jdx = int( (yValues(i) - ymin) / dy) + 1
        if (idx <= 0 .or. idx > nbin) cycle
        if (jdx <= 0 .or. jdx > mbin) cycle
        property(ID) % dist2D(idx,jdx) = property(ID) % dist2D(idx,jdx) + 1
      end do
      property(ID) % counter = property(ID) % counter + 1

    else if (property(ID) % type == "avgDist") then
      xmin = property(ID) % lowerLimit
      xmax = property(ID) % upperLimit
      nbin = property(ID) % numberOfBins
      dx = (xmax-xmin) / dble(nbin)
      do i=1,numberOfValues
        idx = int( (yValues(i) - ymin) / dx) + 1
        if (idx <= 0 .or. idx > nbin) cycle
        property(ID) % tally(idx) = property(ID) % tally(idx) + 1
        property(ID) % dist1D(idx) = property(ID) % dist1D(idx) + xValues(i)
      end do
      property(ID) % counter = property(ID) % counter + 1

    end if

  end subroutine computePropertyRoutine

#ifdef GPTA_MPI
  subroutine communicateProperty(ID)
    implicit none
    integer, intent(in) :: ID
    integer :: nn, i, ierror
    integer, allocatable, dimension(:) :: ncounts, displ
    real(8), allocatable, dimension(:) :: allValues
    
    if (property(ID) % type == "store") then
      ! gather the number of values that each processor has stored
      allocate(ncounts(numberOfMpiProcesses))
      call MPI_GATHER(property(ID) % counter, 1, MPI_INT, ncounts, 1, MPI_INT, 0, MPI_Working_Comm, ierror )

      ! define the starting point for where the received data will be stored
      allocate(displ(numberOfMpiProcesses), source=0)
      do i=2,numberOfMpiProcesses
        displ(i) = ncounts(i-1) + displ(i-1)
      enddo
      
      ! gather all values
      allocate(allValues(sum(ncounts)))
      call MPI_GATHERV(property(ID) % dist1D, property(ID) % counter, MPI_DOUBLE, &
                       allValues, ncounts, displ, MPI_DOUBLE, 0, MPI_Working_Comm, ierror )
      if (me == 0) call move_alloc(allValues, property(ID) % dist1D)

    else
      
      if (me == 0) then
        call MPI_reduce(MPI_IN_PLACE, property(ID) % counter, 1, MPI_INT, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        if (property(ID) % type == "average") then
          call MPI_reduce(MPI_IN_PLACE, property(ID) % average , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)
          call MPI_reduce(MPI_IN_PLACE, property(ID) % average2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        else if (property(ID) % type == "multiAverage") then
          call MPI_reduce(MPI_IN_PLACE, property(ID) % dist1D , 2 * property(ID) % numberOfValues, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        else if (property(ID) % type == "histogram") then
          call MPI_reduce(MPI_IN_PLACE, property(ID) % dist1D , property(ID) % numberOfBins, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        else if (property(ID) % type == "dist2D" .or. property(ID) % type == "prob2D") then
          nn = product(property(ID) % numberOfBins2D)
          call MPI_reduce(MPI_IN_PLACE, property(ID) % dist2D , nn, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        else if (property(ID) % type == "avgDist") then
          call MPI_reduce(MPI_IN_PLACE, property(ID) % dist1D , property(ID) % numberOfBins, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)
          call MPI_reduce(MPI_IN_PLACE, property(ID) % tally  , property(ID) % numberOfBins, MPI_INT   , MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        end if

      else
        call MPI_reduce(property(ID) % counter, property(ID) % counter, 1, MPI_INT, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)
       
        if (property(ID) % type == "average") then
          call MPI_reduce(property(ID) % average , property(ID) % average , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)
          call MPI_reduce(property(ID) % average2, property(ID) % average2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        else if (property(ID) % type == "multiAverage") then
          call MPI_reduce(property(ID) % dist1D, property(ID) % dist1D , 2 * property(ID) % numberOfValues, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        else if (property(ID) % type == "histogram") then
          call MPI_reduce(property(ID) % dist1D, property(ID) % dist1D, property(ID) % numberOfBins, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        else if (property(ID) % type == "dist2D" .or. property(ID) % type == "prob2D") then
          nn = product(property(ID) % numberOfBins2D)
          call MPI_reduce(property(ID) % dist2D, property(ID) % dist2D, nn, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        else if (property(ID) % type == "avgDist") then
          call MPI_reduce(property(ID) % dist1D, property(ID) % dist1D, property(ID) % numberOfBins, MPI_DOUBLE, MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)
          call MPI_reduce(property(ID) % tally , property(ID) % tally , property(ID) % numberOfBins, MPI_INT   , MPI_SUM, 0, MPI_Working_Comm, ierr_mpi)

        end if
      end if
    end if
  end subroutine communicateProperty
#endif

  subroutine dumpPropertyRoutine(ID, iounit, lowerLimits, upperLimits, normalisation, columns)
    integer, intent(in) :: ID
    integer, intent(in) :: iounit
    real(8), dimension(:), intent(in), optional :: lowerLimits  
    real(8), dimension(:), intent(in), optional :: upperLimits  
    character(len=*),      intent(in), optional :: normalisation  
    integer,               intent(in), optional :: columns 

    integer :: i, j, idx, jdx, kdx
    integer :: nbin, mbin
    real(8) :: dx, xmin, xmax, dx2
    real(8) :: dy, ymin, ymax
    real(8) :: norm1, norm0, stdev
    integer :: ncols

    integer :: normalisationType
    real(8) :: sigma, rtmp, w
    real(8), allocatable, dimension(:) :: tmpArray

#ifdef GPTA_MPI
    call communicateProperty(ID)
    if (me /=0) return
#endif

    if (present(lowerLimits)) then
      if (property(ID) % type == "histogram") then
        property(ID) % lowerLimit   = lowerLimits(1)
      else if (property(ID) % type == "dist2D" .or. property(ID) % type == "prob2D") then
        property(ID) % lowerLimits(1:2) = lowerLimits(1:2)
      else if (property(ID) % type == "avgDist") then
        property(ID) % lowerLimit   = lowerLimits(1)
      end if
    end if
    
    if (present(upperLimits)) then
      if (property(ID) % type == "histogram") then
        property(ID) % upperLimit   = upperLimits(1)
      else if (property(ID) % type == "dist2D" .or. property(ID) % type == "prob2D") then
        property(ID) % upperLimits(1:2) = upperLimits(1:2)
      else if (property(ID) % type == "avgDist") then
        property(ID) % upperLimit   = upperLimits(1)
      end if
    end if
    
    ! CHANGE HERE
    if (present(normalisation)) then
      if (index(normalisation,"none") > 0) then
        normalisationType = 1
      else if (index(normalisation,"probability") > 0) then
        normalisationType = 2
      else if (index(normalisation,"custom") > 0) then
        normalisationType = 3
        idx = index(normalisation,"=") + 1
        read(normalisation(idx:),*) rtmp
      else if (index(normalisation,"kernel") > 0) then
        normalisationType = 4
        idx = index(normalisation,"=") + 1
        read(normalisation(idx:),*) sigma
      else
        normalisationType = 1
      end if
    end if
    
    if (present(columns)) then
      ncols = columns
    else
      ncols = 1
    end if
    
    if (property(ID) % type == "store") then
      do idx=1,property(ID) % counter
        write(iounit,'(f20.10)',advance='no') property(ID) % dist1D(idx)
        if (mod(idx,ncols)==0) write(iounit,*)
      end do
      property(ID) % counter = 0
    
    else if (property(ID) % type == "average") then
      norm0 = property(ID) % average  / property(ID) % counter
      norm1 = property(ID) % average2 / property(ID) % counter
      stdev = norm1 - norm0*norm0
      if (stdev < 0.d0) then
        stdev = 0.d0
      else
        stdev = sqrt(stdev)
      end if
      write(iounit,'(f20.10,1x,f20.10)') norm0, stdev
    
    else if (property(ID) % type == "multiAverage") then
      idx = 0
      do i=1,property(ID) % numberOfValues
        norm0 = property(ID) % dist1D(idx+1) / property(ID) % counter
        norm1 = property(ID) % dist1D(idx+2) / property(ID) % counter
        stdev = norm1 - norm0*norm0
        if (stdev < 0.d0) then
          stdev = 0.d0
        else
          stdev = sqrt(stdev)
        end if
        write(iounit,'(3(f15.5,1x))') norm0, stdev, stdev/sqrt( real(property(ID) % numberOfValues) )
        idx = idx + 2
      end do

    else if (property(ID) % type == "histogram" .or. property(ID) % type == "histo") then
      xmin = property(ID) % lowerLimit
      xmax = property(ID) % upperLimit
      nbin = property(ID) % numberOfBins
      dx = (xmax-xmin) / dble(nbin)

      if ( normalisationType == 1 ) then
        norm1 = 1.d0
      else if ( normalisationType == 2 ) then
        norm1 = integrate(property(ID) % numberOfBins, property(ID) % dist1D, dx)
      else if ( normalisationType == 3 ) then
        norm1 = rtmp
      else if ( normalisationType == 4 ) then
        allocate(tmpArray(property(ID) % numberOfBins), source=0.d0)
  
        dx2 = -0.5d0 * dx**2 / sigma**2
        idx = nint(5.d0*sigma / dx)

        do i=1,property(ID) % numberOfBins
          jdx=max(1,i-idx)
          kdx=min(property(ID) % numberOfBins,i+idx)
          do j=jdx,kdx
            w = exp( dx2 * (i-j)**2) * property(ID) % dist1D(i)
            tmpArray(j) = tmpArray(j) + w
          end do
        end do

        property(ID) % dist1D = tmpArray
        norm1 = integrate(property(ID) % numberOfBins, property(ID) % dist1D, dx)
      end if

      do idx=1,property(ID) % numberOfBins
        write(iounit,'(2(e20.10,1x))',advance='no' ) xmin + dx*(idx-0.5)
        write(iounit,'(2(e20.10))'   ,advance='yes') property(ID) % dist1D(idx) / norm1
      enddo
        
    else if (property(ID) % type == "prob2D") then
      xmin = property(ID) % lowerLimits(1)
      xmax = property(ID) % upperLimits(1)
      nbin = property(ID) % numberOfBins2D(1)
      dx = (xmax-xmin) / dble(nbin)
      ymin = property(ID) % lowerLimits(2)
      ymax = property(ID) % upperLimits(2)
      mbin = property(ID) % numberOfBins2D(2)
      dy = (ymax-ymin) / dble(mbin)
      norm1 = sum(property(ID) % dist2D) * dx * dy 
      do jdx=1,property(ID) % numberOfBins2D(2)
        norm0 = integrate(property(ID) % numberOfBins2D(1), property(ID) % dist2D(:,jdx) , dx)
        do idx=1,property(ID) % numberOfBins2D(1)
          write(iounit,'(e20.10,1x)',advance='no' ) xmin + dx*(idx-0.5)
        write(iounit,'(e20.10,1x)',advance='no' ) ymin + dy*(jdx-0.5)
          write(iounit,'(e20.10,1x)',advance='no' ) property(ID) % dist2D(idx,jdx) / norm1
          if (norm0 > epsilon(1.d0)) then
            write(iounit,'(e20.10)',advance='yes') property(ID) % dist2D(idx,jdx) / norm0
          else 
            write(iounit,'(e20.10)',advance='yes') - epsilon(1.d0)
          end if
        enddo
        write(iounit,*)
      enddo
        
    else if (property(ID) % type == "dist2D") then
      xmin = property(ID) % lowerLimits(1)
      xmax = property(ID) % upperLimits(1)
      nbin = property(ID) % numberOfBins2D(1)
      dx = (xmax-xmin) / dble(nbin)
      
      ymin = property(ID) % lowerLimits(2)
      ymax = property(ID) % upperLimits(2)
      mbin = property(ID) % numberOfBins2D(2)
      dy = (ymax-ymin) / dble(mbin)
      
      norm1 = sum(property(ID) % dist2D) * dx * dy
      do jdx=1,property(ID) % numberOfBins2D(2)
        
        norm0 = integrate(property(ID) % numberOfBins2D(1), property(ID) % dist2D(:,jdx) , dx)
        do idx=1,property(ID) % numberOfBins2D(1)
          write(iounit,'(e20.10,1x)',advance='no' ) xmin + dx*(idx-0.5)
          write(iounit,'(e20.10,1x)',advance='no' ) ymin + dy*(jdx-0.5)
          write(iounit,'(e20.10,1x)',advance='no') property(ID) % dist2D(idx,jdx) / property(ID) % counter / norm1 
          if (norm0 > epsilon(1.d0)) then
            write(iounit,'(e20.10,1x)',advance='yes') property(ID) % dist2D(idx,jdx) / norm0
          else 
            write(iounit,'(e20.10,1x)',advance='yes') -tiny(1.d0)
          end if
        enddo
        write(iounit,*)
      enddo
        
    else if (property(ID) % type == "avgDist") then
      xmin = property(ID) % lowerLimit
      xmax = property(ID) % upperLimit
      nbin = property(ID) % numberOfBins
      dx = (xmax-xmin) / dble(nbin)
      do idx=1,property(ID) % numberOfBins
        if (property(ID) % tally(idx) == 0) then
          write(iounit,'(2(e20.10,1x),i0)') xmin + dx*(idx-0.5) &
          , property(ID) % dist1D(idx) &
          , property(ID) % tally(idx)
        else
          write(iounit,'(2(e20.10,1x),i0)') xmin + dx*(idx-0.5) &
          , property(ID) % dist1D(idx) / property(ID) % tally(idx) &
          , property(ID) % tally(idx) 
        end if
      end do
        
    end if

    call flush(iounit)
    
  end subroutine dumpPropertyRoutine
    
  subroutine extractPropertyRoutine(ID, values, tally, numberOfCounts)
    implicit none
    integer, intent(in) :: ID
    real(8), allocatable, dimension(:) :: values
    integer, allocatable, dimension(:), optional :: tally
    integer, optional :: numberOfCounts
    
    integer :: nbin, mbin, i, j, k, idx
      
#ifdef GPTA_MPI
    call communicateProperty(ID)
    if (me /= 0) return
#endif

    if (present(numberOfCounts)) numberOfCounts = property(ID) % counter

    if (property(ID) % type == "store") then
      allocate(values(property(ID) % counter), source=property(ID) % dist1D(1:numberOfCounts))

    else if (property(ID) % type == "average") then
      allocate(values(1))
      values(1) = property(ID) % average  / property(ID) % counter

    else if (property(ID) % type == "multiAverage") then
      allocate(values(property(ID) % numberOfValues))
      idx = 0
      do i=1,property(ID) % numberOfValues
        values(i) = property(ID) % dist1D(idx+1) / property(ID) % counter
        idx = idx + 2
      end do

    else if (property(ID) % type == "histogram") then
      nbin = property(ID) % numberOfBins
      allocate(values(nbin), source=property(ID) % dist1D)

    else if (property(ID) % type == "dist2D") then
      nbin = property(ID) % numberOfBins2D(1)
      mbin = property(ID) % numberOfBins2D(2)
      allocate(values(nbin*mbin))
      k=0
      do j=1,mbin
        do i=1,nbin
          k = k + 1
          values(k) = property(ID) % dist2D(i,j)
        end do
      end do

    else if (property(ID) % type == "avgDist") then
      nbin = property(ID) % numberOfBins
      allocate(values(nbin), source=property(ID) % dist1D)
      allocate(tally(nbin), source=property(ID) % tally)

    end if  
  end subroutine extractPropertyRoutine

  subroutine associateProperties()
    implicit none
    workData % initialise => initialisePropertyRoutine
    workData % compute => computePropertyRoutine
    workData % dump => dumpPropertyRoutine
    workData % extract => extractPropertyRoutine
  end subroutine associateProperties

  function integrate(n,f,dx) result(int)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in), dimension(n) :: f
    real(8), intent(in) :: dx
    real(8) :: int
    integer :: i

    int = 0.5d0+f(1)
    do i=2,n-1
      int = int + f(i)
    end do
    int = int + 0.5d0*f(n)
    int = int * dx

  end function integrate

end module moduleProperties
