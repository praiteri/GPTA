!disclaimer
module moduleSteinhardt
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleFiles
  use moduleProperties
  use moduleMessages
  use moduleDistances
  use moduleElements  

  implicit none

  public :: computeSteinhardt, computeSteinhardtHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions
  
  integer, pointer :: ID, QX, dumpEvery
  integer, pointer :: numberOfBins
  real(8), pointer :: cutoffRadius

  logical, pointer :: QXdump       
  logical, pointer :: QXaverage    
  logical, pointer :: QXrange      
  logical, pointer :: QXhistogram  
  logical, pointer :: QXprobability
  real(8), dimension(:), pointer :: distributionLimits

contains

  subroutine computeSteinhardtHelp()
    use moduleMessages
    implicit none
    call message(0,"This action computes the Steinhardt'sorder parameter")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i geopt.pdb t.0.dcd --stein +q6 +rcut 3.2 +range 0.5,5 +s O O")
    call message(0,"  gpta.x --i geopt.pdb t.0.dcd --stein +q3 +rcut 3.2 +range 0.5,5 +s O O")
    call message(0,"  gpta.x --i geopt.pdb t.0.dcd --stein +q6 +rcut 3.2 +s O O +prob 0,1 +nbin 100 +ndump 100")
  end subroutine computeSteinhardtHelp

  subroutine computeSteinhardt(a)
    implicit none
    type(actionTypeDef), target :: a

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      ! Atoms' selection
      if (firstAction) then
        ! dump info about the action on the screen
        call dumpScreenInfo()

        ! select atoms
        call selectAtoms(2,actionCommand,a)
        call createSelectionList(a,2)

        ! Throw a warning for unused flags
        call checkUsedFlags(actionCommand)
        firstAction = .false.

      end if

      if (QX==3) then
        call computeAction_q3(a)
      else if (QX==6) then
        call computeAction_q6(a)
      end if

      if (dumpEvery>0 .and. mod(tallyExecutions,dumpEvery)==0) then
        call finaliseAction(.false.)
      end if
  
    end if

    if (endOfCoordinatesFiles) then
      call finaliseAction(.true.)
    end if 

  end subroutine computeSteinhardt

  subroutine initialiseAction(a)

    implicit none
    type(actionTypeDef), target :: a
    logical :: lflag(6) = .false.

    ! Local pointers
    actionCommand           => a % actionDetails
    firstAction             => a % firstAction
    tallyExecutions         => a % tallyExecutions
    outputFile              => a % outputFile

    ID                      => a % integerVariables(1)
    QX                      => a % integerVariables(2)
    numberOfBins            => a % integerVariables(3)
    dumpEvery               => a % integerVariables(4)

    cutoffRadius            => a % doubleVariables(1)
    distributionLimits(1:2) => a % doubleVariables(2:3)

    QXdump                  => a % logicalVariables(1)
    QXaverage               => a % logicalVariables(2)
    QXrange                 => a % logicalVariables(3)
    QXhistogram             => a % logicalVariables(4)
    QXprobability           => a % logicalVariables(5)
    
    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .true.
    a % requiresNeighboursListDouble = .true.
    a % cutoffNeighboursList = 3.2d0

    call assignFlagValue(actionCommand,"+q3 ",lflag(3),.false.)
    call assignFlagValue(actionCommand,"+q6 ",lflag(6),.false.)
    if (count(lflag)==0) then
      call message(-1,"Steinhardt order parameter is missing either +q3/+q6")
    else if (count(lflag)==0) then
      call message(-1,"Steinhardt order parameter only one order can be used +q3/+q6")
    end if
    if (lflag(3)) QX = 3
    if (lflag(6)) QX = 6

    ! get output file name from the command line, if present
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'stein.out')
    call initialiseFile(outputFile, outputFile % fname)

    ! to dump distribution every ndump frames
    call assignFlagValue(actionCommand,"+ndump",dumpEvery,-1)

    ! get cutoff radius from the command line, if present
    call assignFlagValue(actionCommand,"+rcut",cutoffRadius,3.2d0)
    
    call assignFlagValue(actionCommand,"+dump " ,a % logicalVariables(1), .false.) ! Dump all
    call assignFlagValue(actionCommand,"+avg "  ,a % logicalVariables(2), .false.) ! Average
    call assignFlagValue(actionCommand,"+range ",a % logicalVariables(3), .false.) ! Count num in a range
    call assignFlagValue(actionCommand,"+histo ",a % logicalVariables(4), .false.) ! Histogram - non-normalised
    call assignFlagValue(actionCommand,"+prob " ,a % logicalVariables(5), .false.) ! Probability distribution - normalised to 1

    if (count(a % logicalVariables(1:5)) > 1) then
      call message(-1,"Only one between +dump/+avg/+thres/+histo/+prob can be used")
    end if

    if (QXdump) then 
      call workData % initialise(ID, "dump", iounit=outputFile % funit)
      
    else if (QXaverage) then
      call workData % initialise(ID, "average")

    else if (QXrange) then
      call assignFlagValue(actionCommand,"+range", distributionLimits(1:2),[0.d0,0.5d0])
      call workData % initialise(ID, "range", iounit=outputFile % funit, \
        lowerLimits=[distributionLimits(1)], upperLimits=[distributionLimits(2)])

    else if (QXhistogram .or. QXprobability) then
      call assignFlagValue(actionCommand,"+nbin ",numberOfBins,100)
      if (QXhistogram) then
        call assignFlagValue(actionCommand,"+histo", distributionLimits(1:2),[0.d0,0.5d0])
      else
        call assignFlagValue(actionCommand,"+prob", distributionLimits(1:2),[0.d0,0.5d0])
      end if
      call workData % initialise(ID, "histogram", \
      numberOfBins=[numberOfBins], lowerLimits=[distributionLimits(1)], upperLimits=[distributionLimits(2)])
    end if

    tallyExecutions = 0

  end subroutine initialiseAction

  subroutine finaliseAction(closeFile)
    implicit none
    logical :: closeFile 

    if (QXhistogram) then
      call workData % dump(ID, outputFile % funit, normalisation='none')
    else if (QXprobability) then
      call workData % dump(ID, outputFile % funit, normalisation='probability', \
                           lowerLimits=[distributionLimits(1)], upperLimits=[distributionLimits(2)])
    else
      call workData % dump(ID, outputFile % funit)
    end if

    if (closeFile) close(outputFile % funit)

  end subroutine finaliseAction

  subroutine dumpScreenInfo()
    implicit none
    call message(0, "Computing Steinhardt's order parameter")
    call message(0, "...Output file",str=outputFile % fname)
    call message(0, "...Cutoff radius",r=cutoffRadius)
    call message(0, "...Number of bins",i=numberOfBins)

  end subroutine dumpScreenInfo

  subroutine computeAction_q6(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: nsel

    integer :: i, ineigh, m, ii
    integer :: iatm, jatm, idx, jdx
    integer :: nn
    real(8) :: dist2, dij(3)

    integer, parameter :: nElem=6

    complex(8), dimension(-nElem:nElem) :: q6
    complex(8), allocatable, dimension(:,:) :: atmQ6
    real(8), allocatable, dimension(:) :: normQ6, nbond

    real(8) :: normaliz(0:nElem)
    real(8), dimension(0:nElem) :: fact1, fact2, coeff_poly
    real(8) :: rsoft = 5.5d0
    real(8) :: frr, pol_ass, u, nnqk
    complex(8) :: com1, comm1(0:nElem)
    complex(8) :: yy(-nElem:nElem)

    real(8), allocatable, dimension(:) :: dotp
    real(8) :: dotik

    fact1 = [720.0, 120.0, 24.0, 6.0, 2.0, 1.0, 1.0]
    fact2 = [720.,5040.,40320.,362880.,3628800.,39916800.,479001600.]
    coeff_poly = [-0.3125,0.0,6.5625,0.0,-19.6875,0.0,14.4375]
    nsel = count(a % isSelected(:,1))
    allocate(atmQ6(-nElem:nElem,nsel))
    allocate(normQ6(nsel))
    allocate(dotp(nsel))
    allocate(nbond(nsel))

    do m = 0, nElem
      normaliz(m)=sqrt((7.0*fact1(m))/(4.0*pi*fact2(m)))*(-1.d0)**m
    end do

    nbond=0
    
!$omp parallel do default(shared) &
!$omp private(i, iatm, nn, q6, ineigh, jatm, dij, dist2, u, frr) &
!$omp private(com1, m, pol_ass, yy, nnqk) &
!$omp schedule(dynamic, 8)
    do i=1,nsel
      iatm = a % idxSelection(i,1)
      nn=0
      q6 = cmplx(0.0d0, 0.0d0)
      
      do ineigh=1,nneigh(iatm)
        jatm = lneigh(ineigh,iatm)
        if ( .not. a % isSelected(jatm,2) ) cycle
        
        ! Early exit to avoid unnecessary computations
        if ( rneigh(ineigh,iatm) > cutoffRadius ) cycle
        nn=nn+1

        ! Cache rneigh value to avoid repeated lookups
        u = rneigh(ineigh,iatm)
        dij(1:3) = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
        dist2 = computeDistanceSquaredPBC(dij)
        dij = dij / u

        ! Precompute frr
        u = u / rsoft
        frr = (1.0d0-u**8)/(1.0d0-u**16)

        nbond(i) = nbond(i) + frr
        com1 = cmplx(dij(1), dij(2))

        ! Unroll inner loops where possible
        do m=0,nElem
          pol_ass = deriv_poly_6(coeff_poly, m, dij(3))
          yy(m) = normaliz(m) * pol_ass * com1**m
          yy(-m) = (-1.0d0)**m * conjg(yy(m))
        end do
        
        do m=-nElem,nElem
          q6(m) = q6(m) + frr * yy(m)
        enddo

      end do 

      if (nn == 0) cycle
      
      q6 = q6 / nbond(i)
      nnqk = q6(0) * q6(0)
      
      do m=1,nElem
        nnqk = nnqk + 2.0d0 * (REAL(q6(m))**2 + DIMAG(q6(m))**2)
      enddo
      
      atmQ6(:,i) = q6
      normq6(i) = dsqrt(nnqk)
    end do
!$omp end parallel do

    dotp = 0.d0
!$omp parallel do &
!$omp private(idx, iatm, nn, ineigh, jatm, jdx, u, frr, dotik, m) &
!$omp shared(nsel, a, nneigh, lneigh, cutoffRadius, rneigh, rsoft, &
!$omp        atmQ6, normq6, nbond, dotp),&
!$omp schedule(dynamic, 8)
    do idx=1,nsel
      iatm = a % idxSelection(idx,1)
      nn=0
      do ineigh=1,nneigh(iatm)
        jatm = lneigh(ineigh,iatm)
        if ( .not. a % isSelected(jatm,2) ) cycle
        if ( rneigh(ineigh,iatm) > cutoffRadius ) cycle
        nn=nn+1

        jdx = a % idxToSelection(jatm,2)
        u = rneigh(ineigh,iatm) / rsoft
        u = u **8
        frr = (1.0d0-u) / (1.0d0-u**2)

        dotik = 0.
        !$omp simd reduction(+:dotik)
        do m = -nElem, nElem
          dotik = dotik + atmQ6(m, idx) * conjg(atmQ6(m, jdx))
        end do
        !$omp end simd

        dotp(idx) = dotp(idx) + frr * dotik / normq6(idx) / normq6(jdx)
      end do
      if (nn==0) cycle
      dotp(idx) = dotp(idx) / nbond(idx)
    end do
!$omp end parallel do

    call workData % compute(ID, numberOfValues=nsel, xValues=dotp)

  end subroutine computeAction_q6

  subroutine computeAction_q3(a)
    implicit none
    type(actionTypeDef), target :: a

    integer :: nsel

    integer :: i, ineigh, m, ii
    integer :: iatm, jatm, idx, jdx
    integer :: nn
    real(8) :: dist2, dij(3)

    integer, parameter :: nElem=3

    complex(8), dimension(-nElem:nElem) :: q3
    complex(8), allocatable, dimension(:,:) :: atmQ3
    real(8), allocatable, dimension(:) :: normQ3, nbond

    real(8) :: normaliz(0:nElem)
    real(8), dimension(0:nElem) :: fact1, fact2, coeff_poly
    real(8) :: rsoft = 5.5d0
    real(8) :: frr, pol_ass, u, nnqk
    complex(8) :: com1, comm1(0:nElem)
    complex(8) :: yy(-nElem:nElem)

    real(8), allocatable, dimension(:) :: dotp
    real(8) :: dotik

    fact1 = [6.0, 2.0,   1.0,   1.0]
    fact2 = [6.0, 24.0,  120.0, 720.]
    coeff_poly = [0.d0, -1.5d0, 0.d0, 2.5d0]
    nsel = count(a % isSelected(:,1))
    allocate(atmq3(-nElem:nElem,nsel))
    allocate(normq3(nsel))
    allocate(dotp(nsel))
    allocate(nbond(nsel))

    do m = 0, nElem
      normaliz(m)=sqrt((7.0*fact1(m))/(4.0*pi*fact2(m)))*(-1.d0)**m
    end do

    nbond=0
    do i=1,nsel
      iatm = a % idxSelection(i,1)

      nn=0
      q3 = cmplx(0.0 ,0.0)
      do ineigh=1,nneigh(iatm)
        jatm = lneigh(ineigh,iatm)
        if ( .not. a % isSelected(jatm,2) ) cycle
        if ( rneigh(ineigh,iatm) > cutoffRadius ) cycle
        nn=nn+1

        dij(1:3) = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
        dist2 = computeDistanceSquaredPBC(dij)
        dij = dij / rneigh(ineigh,iatm)

        u = rneigh(ineigh,iatm) / rsoft
        frr = (1.0d0-u**8)/(1.0d0-u**16)

        nbond(i) = nbond(i) + frr
        com1=cmplx( dij(1) , dij(2) )

        do m=0,nElem
          pol_ass=deriv_poly_3(coeff_poly,m,dij(3))
          yy(m)=normaliz(m)*pol_ass*com1**m
          yy(-m)=(-1.0d0)**m*conjg(yy(m))
       end do
       do m= -nElem,nElem
          q3(m) = q3(m) + frr*yy(m)
       enddo

      end do 
      if (nn==0) cycle
      q3 = q3/nbond(i)
      nnqk = q3(0)*q3(0)
      do m = 1,nElem
         nnqk = nnqk + 2.d0*(REAL(q3(m))**2 + DIMAG(q3(m))**2)
      enddo
      
      atmq3(:,i) = q3
      normq3(i) = dsqrt(nnqk)
    end do
      
    dotp = 0.d0
    do idx=1,nsel
      iatm = a % idxSelection(idx,1)

      nn=0
      do ineigh=1,nneigh(iatm)
        jatm = lneigh(ineigh,iatm)
        if ( .not. a % isSelected(jatm,2) ) cycle
        if ( rneigh(ineigh,iatm) > cutoffRadius ) cycle
        nn=nn+1

        jdx = a % idxToSelection(jatm,2)
        u = rneigh(ineigh,iatm) / rsoft
        frr = (1.0d0-u**8) / (1.0d0-u**16)

        dotik = 0.
        do m = -nElem,nElem
          dotik = dotik + atmq3(m,idx)*conjg(atmq3(m,jdx))
        enddo
        dotp(idx) = dotp(idx) + frr * dotik / normq3(idx) / normq3(jdx)
      end do
      if (nn==0) cycle
      dotp(idx) = dotp(idx) / nbond(idx)
    end do

    call workData % compute(ID, numberOfValues=nsel, xValues=dotp)

  end subroutine computeAction_q3

  FUNCTION deriv_poly_3(coeff_poly,order,x) result(res)
    integer, parameter                   :: dbl=kind(1.d0)
    real(dbl), dimension(0:3), intent(in):: coeff_poly
    integer, intent(in)                  :: order
    real(dbl)                            :: x,xi
    real(dbl)                            :: res
    integer                              :: i,j,fact

    res = 0.0_dbl
    xi  = 1.0_dbl
    do i=order,ubound(coeff_poly,1)
       ! fact= i!/(i-order)!
       fact=1
       do j=(i-order+1),i
          fact=fact*j
       end do
       res=res+coeff_poly(i)*dble(fact)*xi
       xi=xi*x
    end do

    return
  END FUNCTION deriv_poly_3

  FUNCTION deriv_poly_6(coeff_poly,order,x) result(res)
    integer, parameter                   :: dbl=kind(1.d0)
    real(dbl), dimension(0:6), intent(in):: coeff_poly
    integer, intent(in)                  :: order
    real(dbl)                            :: x,xi
    real(dbl)                            :: res
    integer                              :: i,j,fact

    res = 0.0_dbl
    xi  = 1.0_dbl
    do i=order,ubound(coeff_poly,1)
       ! fact= i!/(i-order)!
       fact=1
       do j=(i-order+1),i
          fact=fact*j
       end do
       res=res+coeff_poly(i)*dble(fact)*xi
       xi=xi*x
    end do

    return
    END FUNCTION deriv_poly_6

end module moduleSteinhardt
