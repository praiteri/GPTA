!disclaimer
module moduleDensityMap2D
  use moduleVariables
  use moduleStrings
  use moduleMessages 
  use moduleSystem 
  use moduleFiles
  use moduleDistances
  use moduleProperties

  public :: computeDensityMap2D, computeDensityMap2DHelp
  private

  character(:), pointer :: actionCommand

  type(fileTypeDef), pointer :: outputFile
  
  integer, pointer :: tallyExecutions
  integer, pointer :: numberOfLocalAtoms
  integer, pointer :: numberOfReplicas
  integer, pointer, dimension(:) :: numberOfBins
  integer, pointer, dimension(:) :: hkl

  real(real64), pointer :: thickness
  real(real64), pointer, dimension(:) :: origin
  real(real64), pointer, dimension(:) :: normalUnitVector
  real(real64), pointer, dimension(:) :: delta

  logical, pointer :: inputCell
  real(real64), pointer, dimension(:) :: surfaceCell

  integer, pointer :: ID

contains

  subroutine computeDensityMap2DHelp()
    implicit none
    call message(0,"This action computes the 2D density map projected on an arbitrary plane defined by its Miller indices.")
    call message(0,"The location of the plane in the cell is chosen by defining the origin of the axes of the new cell.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --dmap2D +s O2 +hkl 1,0,0 +origin 1,0,0 +nbin 50,50 +out dmap2D.out")
    call message(0,"  gpta.x --i coord.pdb --dmap2D +s O2 +hkl 1,1,0 +origin 46.12,0,0  +nbin 50,50 +thick 1.")
    call message(0,"  gpta.x --i coord.pdb --dmap2D +s O2 +hkl 1,1,0 +origin 0,0,0 +thick 1. +fill 1")
    call message(0,"  gpta.x --i coord.pdb --dmap2D +s O2 +cell 46.12,0,0,-30.7036,44.9762,0,0,0,86.4160")
  end subroutine computeDensityMap2DHelp

  subroutine computeDensityMap2D(a)
    implicit none
    type(actionTypeDef), target :: a
    
    integer, dimension(2) :: itmp
    
    real(real64), allocatable, dimension(:,:) :: cartesianCoord
    real(real64), allocatable, dimension(:,:) :: fractionalCoord
    
    integer :: iatm, jatm, ntmp, mtmp, i, j, k
    real(real64) :: rtmp(3), vol
    real(real64) :: dist, sij(3), dij(3)
    logical, allocatable, dimension(:) :: ldel
    real(real64), dimension(3,3) :: hnew, hinv

    integer :: nrepl
    real(real64), allocatable, dimension(:,:) :: vrepl
    character(len=STRLEN) :: stringCell

    ! Associate variables
    actionCommand          => a % actionDetails
    tallyExecutions        => a % tallyExecutions
    
    outputFile             => a % outputFile

    numberOfBins(1:2)      => a % integerVariables(1:2)
    hkl(1:3)               => a % integerVariables(3:5)
    numberOfLocalAtoms     => a % integerVariables(6)
    numberOfReplicas       => a % integerVariables(7)
    ID                     => a % integerVariables(8)

    origin(1:3)            => a % doubleVariables(1:3)
    normalUnitVector(1:3)  => a % doubleVariables(4:6)
    thickness              => a % doubleVariables(7)
    delta                  => a % doubleVariables(8:9)
    
    inputCell              => a % logicalVariables(1)
    surfaceCell(1:9)       => a % doubleVariables(10:18)

    if (a % actionInitialisation ) then
      a % actionInitialisation  = .false.
      a % requiresNeighboursList = .false.

      call assignFlagValue(actionCommand,"+out",outputFile % fname,'dmap2D.out')

      call assignFlagValue(actionCommand,"+hkl",hkl,[0,0,1])
      normalUnitVector = dble(hkl) 
      normalUnitVector = normalUnitVector / sqrt( sum( normalUnitVector * normalUnitVector ) )

      call assignFlagValue(actionCommand,"+origin",origin,[0.0_real64,0.0_real64,0.0_real64])

      call assignFlagValue(actionCommand,"+thick",thickness,1.0_real64)

      call assignFlagValue(actionCommand,"+nbin",itmp,[100,100])
      numberOfBins(1:2) = itmp

      call assignFlagValue(actionCommand,"+fill",numberOfReplicas,0)

      call assignFlagValue(actionCommand,"+cell",inputCell,.false.)
      if (inputCell) then
        call assignFlagValue(actionCommand,"+cell ", stringCell, "NONE")
        call readCellFreeFormat(stringCell, hnew)
        surfaceCell = reshape(hnew, [9])
      end if

      call workData % initialise(ID, "dist2D", numberOfBins=numberOfBins, lowerLimits=[0.0_real64,0.0_real64], upperLimits=[1.0_real64,1.0_real64])

      allocate(a % array2D(itmp(1),itmp(2)) , source=0.0_real64)
      tallyExecutions = 0
      return
    end if

    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      ! Select atoms
      if (a % updateAtomsSelection) then

        if (.not. inputCell) then
          ! looking for two vectors in the surface plane
          call defineSurfaceVectors(hkl, frame % hmat, hnew)

          call getInverseCellMatrix(hnew,hinv,vol)
          if (vol<0.0) then
            rtmp = hnew(1:3,2)
            hnew(1:3,2) = hnew(1:3,1)
            hnew(1:3,1) = rtmp 
          end if

          surfaceCell = reshape(hnew, [9])
        end if

        ! grid spacing in fractional coordinates
        delta = 1.0_real64 / numberOfBins

        ! write actions info
        call dumpActionInfo()

        ! select atoms
        call selectAtoms(1,actionCommand,a)
        call createSelectionList(a,1)
        a % updateAtomsSelection=.false.

        numberOfLocalAtoms = count(a % isSelected(:,1))
      end if
      
      ! replicate the cell to fill the surface unit cell
      nrepl = 0
      allocate(vrepl(3, (2*numberOfReplicas + 1)**3 ))
      do i=-numberOfReplicas,numberOfReplicas
        do j=-numberOfReplicas,numberOfReplicas
          do k=-numberOfReplicas,numberOfReplicas
            nrepl = nrepl + 1
            vrepl(1:3,nrepl) = i*frame % hmat(1:3,1) + j*frame % hmat(1:3,2) + k*frame % hmat(1:3,3) - origin(1:3)
          end do
        end do
      end do

      ! Local array with the coordinates of only the selected atoms
      allocate(cartesianCoord(3,numberOfLocalAtoms * nrepl))
      ntmp = 0
      do i=1,nrepl
        do iatm=1,numberOfLocalAtoms
          ntmp = ntmp + 1
          cartesianCoord(1:3,ntmp) = frame % pos(1:3,a % idxSelection(iatm,1)) + vrepl(1:3,i)
        end do
      end do

      ! Rotate the selected atoms with the surface unit cell
      hnew = reshape(surfaceCell, [3,3])

      call makeUpperTriangularCell(hnew, cartesianCoord, ntmp)
      call getInverseCellMatrix(hnew, hinv, vol)

      allocate(fractionalCoord(3,ntmp))
      do iatm=1,ntmp
        dij(1:3) = cartesianCoord(1:3,iatm)
        sij = matmul(hinv,dij)
        sij(1:3) = sij(1:3) - floor(sij(1:3))
        dij = matmul(hnew,sij)
        fractionalCoord(1:3,iatm) = sij(1:3)
        cartesianCoord(1:3,iatm) = dij(1:3)
      end do

     ! Delete atoms if they overlap
      if (nrepl > 1) then
        allocate(ldel(ntmp), source=.false.)
        do iatm=1,ntmp
          if (ldel(iatm)) cycle
          do jatm=iatm+1,ntmp
            if (ldel(jatm)) cycle
            dij(1:3) = cartesianCoord(1:3,iatm) - cartesianCoord(1:3,jatm)
            sij = matmul(hinv,dij)
            sij(1:3) = sij(1:3) - nint(sij(1:3))
            dij = matmul(hnew,sij)
            dist = sqrt( sum(dij*dij) )
            if (dist < 0.1_real64) ldel(jatm) = .true.
          end do
        end do
        mtmp = ntmp
        ntmp = 0
        do iatm=1,mtmp
          if (ldel(iatm)) cycle
          ntmp = ntmp + 1
          cartesianCoord(1:3,ntmp) = cartesianCoord(1:3,iatm)
          fractionalCoord(1:3,ntmp) = fractionalCoord(1:3,iatm)
        end do
        deallocate(ldel)
      end if

      ! extract the atoms that are within the desired thickness
      ntmp = 0
      do iatm=1,numberOfLocalAtoms * nrepl
        if (cartesianCoord(3,iatm) > thickness .or. cartesianCoord(3,iatm) < -0.0_real64) cycle
        ntmp = ntmp + 1
        fractionalCoord(:,ntmp) = fractionalCoord(:,iatm)
      end do

      ! write(0,*)ntmp,mtmp
      ! allocate(label(ntmp), source="Ar  ")
      ! call dumpPDB(123,ntmp,cartesianCoord,label,hnew)

      call workData % compute(ID, numberOfValues=ntmp, xValues=fractionalCoord(1,:), yValues=fractionalCoord(2,:))

    end if ! frames read succefsully

    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) then
      call finaliseAction()
    end if

  end subroutine ComputeDensityMap2D

  subroutine finaliseAction()
    implicit none
    real(real64) :: da(2), db(2), pp(2), vol
    real(real64), dimension(3,3) :: hnew
    real(real64), allocatable, dimension(:) :: dataArray
    integer :: i, j, k
    integer :: numberOfCounts

    call workData % extract(ID, dataArray, numberOfCounts=numberOfCounts)
    if (me /= readingCPU) return

    hnew = reshape(surfaceCell, [3,3])
    call makeUpperTriangularCell(hnew, origin, 1)
    da = hnew(1:2,1) / numberOfBins(1)
    db = hnew(1:2,2) / numberOfBins(2)
    vol = thickness * abs(da(1)*db(2) - da(2)*db(1))
 
    dataArray = dataArray / numberOfCounts
    dataArray = dataArray !/ vol

    call initialiseFile(outputFile, outputFile % fname)
    write(outputFile % funit,'("# X [angstrom] | Y [angstrom] | Density [atom/angstrom^3]")')
    k=0
    do j=1,numberOfBins(2)
      do i=1,numberOfBins(1)
        pp = i*da + j*db
        k = k + 1
        write(outputFile % funit,'(3(1x,f12.5))') pp, dataArray(k) ! dmap(i,j) 
      end do
      write(outputFile % funit,*)
    end do
    close(outputFile % funit)

  end subroutine finaliseAction

  subroutine dumpActionInfo()
    implicit none
    
    call message(0,"Compute 2D density map in a plane")
    call message(0,"...Miller indices of the plane",iv=hkl)
    call message(0,"...Surface cell origin",rv=origin)
    call message(0,"...Surface vector A",rv=surfaceCell(1:3))
    call message(0,"...Surface vector B",rv=surfaceCell(4:6))
    call message(0,"...Slice thickness",r=thickness)
    call message(0,"...Number of replicas to fill surface cell",i=numberOfReplicas)
    call message(0,"...Number of bins",iv=numberOfBins)
    
  end subroutine dumpActionInfo

end module moduleDensityMap2D
