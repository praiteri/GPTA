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
module moduleRDF

  use moduleVariables
  use moduleSystem
  use moduleFiles
  use moduleStrings

  implicit none

  public :: computeRadialPairDistribution, computeRadialPairDistributionHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  
  integer, pointer :: numberOfBins
  integer, pointer :: tallyExecutions
  real(8), pointer :: rcut
  real(8), pointer :: averageVolume
  real(8), pointer, dimension(:) :: dist
  logical, pointer :: soluteRDF
  logical, pointer :: logscale

contains

  subroutine computeRadialPairDistributionHelp()
    use moduleMessages
    implicit none
    call message(0,"This action computes the Radial Pair Distribution function")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb traj.dcd --gofr +s Ca O1,O2 +out gofr.out +nbin 100")
  end subroutine computeRadialPairDistributionHelp

  subroutine computeRadialPairDistribution(a)
    use moduleMessages
    implicit none
    type(actionTypeDef), target :: a

    call associatePointers(a)
    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      if (firstAction) then         
        call dumpScreenInfo()

        if (resetFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if
        
        call selectAtoms(2,actionCommand,a)

        ! create a list of the atoms' indices for each group
        call createSelectionList(a,2)

        ! Check groups have atoms
        if (count(a % isSelected(:,1)) == 0) call message(-1,"--gofr - no atoms selected in the first group")
        if (count(a % isSelected(:,2)) == 0) call message(-1,"--gofr - no atoms selected in the second group")

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else

        if (a % updateAtomsSelection) then
          call selectAtoms(2,actionCommand,a)
          call createSelectionList(a,2)
        end if

      end if

      if (soluteRDF) then
        if (logscale) then
          call computeActionSoluteLogScale(a)
        else
          call computeActionSolute(a)
        endif
      else
        call computeAction(a)
      end if

    end if

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) then
      if (logscale) then
        call finaliseActionLogScale(a)
      else
        call finaliseAction(a)
      end if
    endif

  end subroutine computeRadialPairDistribution

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    outputFile           => a % outputFile
    numberOfBins         => a % numberOfBins
    rcut                 => a % doubleVariables(1)
    averageVolume        => a % doubleVariables(2)
    soluteRDF            => a % logicalVariables(1)
    logscale             => a % logicalVariables(2)
    dist                 => a % array1D
    
  end subroutine  

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    a % actionInitialisation = .false.
    a % cutoffNeighboursList = 6.0d0

    call assignFlagValue(actionCommand,"+solute",soluteRDF,.false.)
    call assignFlagValue(actionCommand,"+log",logscale,.false.)
    
    if (soluteRDF) then
      a % requiresNeighboursList = .false.
      a % requiresNeighboursListUpdates = .false.
    else
      a % requiresNeighboursList = .true.
      a % requiresNeighboursListUpdates = .true.
    end if
    a % requiresNeighboursListDouble = .false.

    call assignFlagValue(actionCommand,"+rcut",rcut,a % cutoffNeighboursList)
    call assignFlagValue(actionCommand,"+nbin",numberOfBins,100)
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'gofr.out')
    allocate(a % array1D(numberOfBins) , source=0.d0)

    a % cutoffNeighboursList = rcut
    tallyExecutions = 0 

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    use moduleMessages 
    implicit none
    call message(0,"Settings for the radial pair distribution function")
    call message(0,"...Cutoff distance",r=rcut)
    call message(0,"...Number of bins",i=numberOfBins)
    call message(0,"...Output file",str=outputFile % fname)

  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleVariables
    use moduleSystem
    implicit none
    type(actionTypeDef), target :: a

    integer :: iatm, ineigh, jatm, idx
    integer :: imol, jmol
    real(8) :: dr, dn
    real(8), allocatable, dimension(:) :: local_dist

    dr = rcut / dble(numberOfBins)
    if (computeNeighboursListDouble) then
      dn = 0.5d0
    else
      dn = 1.0d0
    end if
    averageVolume = averageVolume + frame % volume

    allocate(local_dist(numberOfBins) , source=0.d0)

    if (numberOfMolecules > 0) then
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (iatm,ineigh,jatm,idx)
!$OMP DO REDUCTION(+:local_dist)
      do iatm=1,frame % natoms
        imol = atomToMoleculeIndex(iatm)
        do ineigh=1,nneigh(iatm)
          if (rneigh(ineigh,iatm) > rcut) cycle
          jatm=lneigh(ineigh,iatm)
          jmol = atomToMoleculeIndex(jatm)
          if (imol == jmol) cycle
          idx = int(rneigh(ineigh,iatm)/dr) + 1
          if ( (a % isSelected(iatm,1) .and. a % isSelected(jatm,2))) local_dist(idx) = local_dist(idx) + dn
          if ( (a % isSelected(iatm,2) .and. a % isSelected(jatm,1))) local_dist(idx) = local_dist(idx) + dn
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
    else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (iatm,ineigh,jatm,idx)
!$OMP DO REDUCTION(+:local_dist)
      do iatm=1,frame % natoms
        do ineigh=1,nneigh(iatm)
          if (rneigh(ineigh,iatm) > rcut) cycle
          jatm=lneigh(ineigh,iatm)
          idx = int(rneigh(ineigh,iatm)/dr) + 1
          if ( (a % isSelected(iatm,1) .and. a % isSelected(jatm,2))) local_dist(idx) = local_dist(idx) + dn
          if ( (a % isSelected(iatm,2) .and. a % isSelected(jatm,1))) local_dist(idx) = local_dist(idx) + dn
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
    end if

    dist = dist + local_dist
    deallocate(local_dist)

  end subroutine computeAction

  subroutine computeActionSolute(a)
    use moduleVariables
    use moduleSystem 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a

    integer :: i, j, iatm, jatm, idx, nsel1, nsel2
    real(8) :: dr, dij(3), r2, dist2
    integer, allocatable, dimension(:) :: local_dist

    dr = rcut / dble(numberOfBins)
    r2 = rcut**2
    averageVolume = averageVolume + frame % volume

    nsel1 = count(a % isSelected(:,1))
    nsel2 = count(a % isSelected(:,2))

    allocate(local_dist(numberOfBins) , source=0)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,j,iatm,jatm,dij,dist2,idx)
!$OMP DO REDUCTION(+:local_dist)
    do i=1,nsel1
      iatm = a % idxSelection(i,1)
      ! write(0,*)tallyExecutions,iatm
      do j=1,nsel2
        jatm = a % idxSelection(j,2)

        dij = frame % pos(:,jatm) - frame % pos(:,iatm)
        dist2 = computeDistanceSquaredPBC(dij)

        if (dist2 > r2) cycle

        idx = int(sqrt(dist2)/dr) + 1
        local_dist(idx) = local_dist(idx) + 1
      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    dist = dist + real(local_dist,8)

  end subroutine computeActionSolute

  subroutine finaliseAction(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: i
    integer :: nsel1, nsel2
    real(8) :: dr, rtmp, rtmp1, rtmp2, integral, nn, numberDensity
    real(8) :: pi43, avg

#ifdef GPTA_MPI
    call actionCommunication()
    if (me /= 0) return
#endif

    call initialiseFile(outputFile,outputFile % fname)

    dr = rcut / dble(numberOfBins)
    avg = averageVolume / dble(tallyExecutions)
    
    nsel1 = count(a % isSelected(:,1))
    nsel2 = count(a % isSelected(:,2))
    numberDensity = nsel2 / avg
    
    rtmp = nsel1 * dble(tallyExecutions)
    dist = dist / rtmp 

    write(outputFile % funit,"(a)") "# Distance | g(r) | Coordination Number | KB integral"
    integral = 0.d0
    pi43 = 4.0d0/3.0d0*pi
    do i=2,numberOfBins
      rtmp1 = dr*(i)
      rtmp2 = dr*(i-1)
      nn = pi43*(rtmp1**3-rtmp2**3)*numberDensity
      integral = integral + dist(i)
      write(outputFile % funit,'(4(f15.5,1x))') &
        0.5d0*(rtmp1+rtmp2),                    & ! distance
        dist(i)/nn,                             & ! g(r)
        integral,                               & ! coordination number
        integral/numberDensity - pi43*rtmp1**3    ! Kirkwood-Buff integral
    enddo
    close(outputFile % funit)

  end subroutine finaliseAction

  subroutine computeActionSoluteLogScale(a)
    use moduleVariables
    use moduleSystem 
    use moduleDistances
    implicit none
    type(actionTypeDef), target :: a

    integer :: i, j, iatm, jatm, idx, nsel1, nsel2
    real(8) :: dr, dij(3), r2, dist2, dd
    integer, allocatable, dimension(:) :: local_dist

    dr = log(rcut) / dble(numberOfBins)
    r2 = rcut**2
    averageVolume = averageVolume + frame % volume

    nsel1 = count(a % isSelected(:,1))
    nsel2 = count(a % isSelected(:,2))

    allocate(local_dist(numberOfBins) , source=0)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (i,j,iatm,jatm,dij,dist2,idx)
!$OMP DO REDUCTION(+:local_dist)
    do i=1,nsel1
      iatm = a % idxSelection(i,1)
      do j=1,nsel2
        jatm = a % idxSelection(j,2)

        dij = frame % pos(:,jatm) - frame % pos(:,iatm)
        dist2 = computeDistanceSquaredPBC(dij)

        if (dist2 > r2) cycle

        dd = sqrt(dist2)
        idx = int(log(dd)/dr) + 1
        local_dist(idx) = local_dist(idx) + 1
      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    dist = dist + real(local_dist,8)

  end subroutine computeActionSoluteLogScale

  subroutine finaliseActionLogScale(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: i
    integer :: nsel1, nsel2
    real(8) :: dr, rtmp, rtmp1, rtmp2, integral, nn, numberDensity
    real(8) :: pi43, avg

#ifdef GPTA_MPI
    call actionCommunication()
    if (me /= 0) return
#endif

    call initialiseFile(outputFile,outputFile % fname)

    dr = log(rcut) / dble(numberOfBins)
    avg = averageVolume / dble(tallyExecutions)
    
    nsel1 = count(a % isSelected(:,1))
    nsel2 = count(a % isSelected(:,2))
    numberDensity = nsel2 / avg
    
    rtmp = nsel1 * dble(tallyExecutions)
    dist = dist / rtmp 
    
    write(outputFile % funit,"(a)") "# Distance | g(r) | Coordination Number | KB integral"
    integral = 0.d0
    pi43 = 4.0d0/3.0d0*pi
    do i=2,numberOfBins
      rtmp1 = exp(dr*(i))
      rtmp2 = exp(dr*(i-1))
      nn = pi43*(rtmp1**3-rtmp2**3)*numberDensity
      integral = integral + dist(i)
      write(outputFile % funit,'(4(f10.5,1x))') &
        0.5d0*(rtmp1+rtmp2),                    & ! distance
        dist(i)/nn,                             & ! g(r)
        integral,                               & ! coordination number
        integral/numberDensity - pi43*rtmp1**3    ! Kirkwood-Buff integral
    enddo
    close(outputFile % funit)

  end subroutine finaliseActionLogScale

#ifdef GPTA_MPI
  subroutine actionCommunication()
    implicit none
    if (me == 0) then
      call MPI_reduce(MPI_IN_PLACE,averageVolume  ,1           ,MPI_DOUBLE,MPI_SUM,0,MPI_Working_Comm,ierr_mpi)
      call MPI_reduce(MPI_IN_PLACE,tallyExecutions,1           ,MPI_INT   ,MPI_SUM,0,MPI_Working_Comm,ierr_mpi)
      call MPI_reduce(MPI_IN_PLACE,dist           ,numberOfBins,MPI_DOUBLE,MPI_SUM,0,MPI_Working_Comm,ierr_mpi)
    else
      call MPI_reduce(averageVolume  ,averageVolume  ,1           ,MPI_DOUBLE,MPI_SUM,0,MPI_Working_Comm,ierr_mpi)
      call MPI_reduce(tallyExecutions,tallyExecutions,1           ,MPI_INT   ,MPI_SUM,0,MPI_Working_Comm,ierr_mpi)
      call MPI_reduce(dist           ,dist           ,numberOfBins,MPI_DOUBLE,MPI_SUM,0,MPI_Working_Comm,ierr_mpi)
    end if
  end subroutine actionCommunication
#endif

end module moduleRDF
