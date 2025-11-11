!disclaimer
module moduleMolecularProperties
  use moduleVariables
  use moduleFiles
  use moduleProperties
  
  implicit none

  public :: computeMolecularProperties, computeMolecularPropertiesHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction

  type(fileTypeDef), pointer :: outputFile
  integer, pointer :: tallyExecutions

  character(len=STRLEN), pointer :: targetMolecule

  integer, pointer :: numberOfLocalMolecules
  integer, pointer, dimension(:) :: localIndices
  
  integer, pointer :: dumpUnit, ID
  
  integer, pointer :: numberOfBins
  integer, dimension(:), pointer :: numberOfBins2D
  real(real64), dimension(:), pointer :: distributionLimits

  logical, pointer :: dumpProperty    ! Output dump all molecule on the same line
  logical, pointer :: dumpPropertySeq ! Output dump sequentially
  logical, pointer :: averageProperty ! Average
  logical, pointer :: histoProperty   ! Histogram
  logical, pointer :: distProperty    ! Distribution (normalised)
  logical, pointer :: distZProperty   ! Distribution of Averages (z)
  logical, pointer :: distZProperty2D
  logical, pointer :: degrees

  integer, pointer, dimension(:) :: torsionIndices
  integer, pointer, dimension(:) :: normalIndices
  integer, pointer, dimension(:) :: angleIndices
  integer, pointer, dimension(:) :: bondIndices

contains

  subroutine computeMolecularPropertiesHelp()
    use moduleMessages
    implicit none
    call message(0,"This action computes a variety of molecular properties of the system.")
    call message(0,"The properties that can currently be computed are")
    call message(0,"  the molecular dipole.")
    call message(0,"  the polar angle of the molecular dipole.")
    call message(0,"  the polar angle for the vector normal to the plane defined by three atoms in the molecule.")
    call message(0,"  the angle between three atoms in the molecule.")
    call message(0,"  the torsional angle between four atoms in the molecule.")
    call message(0,"  the radius of gyration.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i c.pdb --top --molprop +id M1 +bond 1,2 +avg")
    call message(0,"  gpta.x --i c.pdb --top --molprop +id M1 +angle 1,2,3 +avg")
    call message(0,"  gpta.x --i c.pdb --top --molprop +id M1 +torsion 1,2,3,4 +histo -pi,pi +nbin 180")
    call message(0,"  gpta.x --i c.pdb --top --molprop +id M1 +torsion 1,2,3,4 +dist -pi,pi +nbin 180")
    call message(0,"  gpta.x --i c.pdb --top --molprop +id M1 +dipole +dist 0,4 +nbin 50 +out dipole.out")
    call message(0,"  gpta.x --i c.pdb --top --molprop +id M1 +dipoleZ +distZ +nbinZ 100")
    call message(0,"  gpta.x --i c.pdb --top --molprop +id M1 +normalZ 1,2,3 +dist2D 0,pi +nbinZ 100 +nbin 90")
    call message(0,"  gpta.x --i c.pdb --top --molprop +id M1 +dipole +dump +out dipole.out")
  end subroutine computeMolecularPropertiesHelp
  
  subroutine computeMolecularProperties(a)
    use moduleSystem 
    use moduleDistances
    use moduleMessages
    implicit none
    type(actionTypeDef), target :: a

    real(real64), allocatable, dimension(:) :: localProperty
    real(real64), allocatable, dimension(:,:) :: cartesianCoord
    real(real64), allocatable, dimension(:,:) :: fractionalCoord
    integer :: ilocal, imol

    call associatePointers(a)
    
    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if
    ! Normal processing of the frame
    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      if (firstAction) then

        if (resetFrameLabels) then
          a % updateAtomsSelection = .false.
        else
          a % updateAtomsSelection = .true.
        end if

        if (numberOfMolecules == 0) call runInternalAction("topology","NULL")

        call getLocalIndices(targetMolecule, a, numberOfLocalMolecules)
        localIndices => a % localIndices
  
        call dumpScreenInfo(a)

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      else
        
        if (a % updateAtomsSelection) then
          call getLocalIndices(targetMolecule, a, numberOfLocalMolecules)
          localIndices => a % localIndices
        end if

      end if

      call associatePointers(a)

      ! call computeMoleculesCOM()
      allocate(cartesianCoord(3,numberOfLocalMolecules))       
      allocate(fractionalCoord(3,numberOfLocalMolecules))       
      do ilocal=1,numberOfLocalMolecules
        imol = localIndices(ilocal)
        cartesianCoord(1:3,ilocal) = listOfMolecules(imol) % centreOfMass(1:3)      
      enddo
      call cartesianToFractional(numberOfLocalMolecules, cartesianCoord, fractionalCoord)

      allocate(localProperty(numberOfLocalMolecules), source=0.0_real64)
      call a % molprop % property(numberOfLocalMolecules, localIndices, localProperty)
      
      if ( trim(a % molprop % command) == "+angle" .or. trim(a % molprop % command) == "+torsion") then
        if (degrees) localProperty = localProperty / pi * 180
      end if
      call workData % compute(ID, numberOfValues=numberOfLocalMolecules, xValues=localProperty, yValues=fractionalCoord(3,:))

    end if
    
    ! Normal processing of the frame - finalise calculation and write output
    if (endOfCoordinatesFiles) call finaliseAction()

  end subroutine computeMolecularProperties

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    actionCommand           => a % actionDetails
    firstAction             => a % firstAction
    tallyExecutions         => a % tallyExecutions
    outputFile              => a % outputFile

    targetMolecule          => a % stringVariables(1)
    
    numberOfBins            => a % numberOfBins
    numberOfBins2D(1:2)     => a % numberOfBins2D

    dumpProperty            => a % logicalVariables(1) ! Output dump all molecule on the same line
    dumpPropertySeq         => a % logicalVariables(2) ! Output dump sequentially
    averageProperty         => a % logicalVariables(3) ! Average
    histoProperty           => a % logicalVariables(4) ! Histogram
    distProperty            => a % logicalVariables(5) ! Distribution (normalised)
    distZProperty           => a % logicalVariables(6) ! Distribution of Averages (z)
    distZProperty2D         => a % logicalVariables(7) ! 2D map
    degrees                 => a % logicalVariables(8) ! angles in degrees

    numberOfLocalMolecules  => a % integerVariables(1)
    localIndices            => a % localIndices

    dumpUnit                => a % outputFile % funit

    ID                      => a % integerVariables(2)

    torsionIndices(1:4)     => a % integerVariables(3:6)
    normalIndices(1:3)      => a % integerVariables(7:9)
    angleIndices(1:3)       => a % integerVariables(10:12)
    bondIndices(1:2)       => a % integerVariables(13:14)

    distributionLimits(1:2) => a % doubleVariables(1:2)

  end subroutine 

  subroutine initialiseAction(a)
    use moduleStrings
    use moduleMessages
    use moduleSystem
    implicit none
    type(actionTypeDef), target :: a

    integer :: n

    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .false.
    a % cutoffNeighboursList = 3.0_real64
    
    call assignFlagValue(actionCommand,"+id",targetMolecule,"NULL")
    
    call assignFlagValue(actionCommand,"+out",outputFile % fname,'molprop.out')
    call initialiseFile(outputFile,outputFile % fname)
    
    n=0
    if (flagExists(actionCommand,"+dipole ")) then
      n = n + 1
      a % molprop % command  = "+dipole"
      a % molprop % name     = "Molecular dipoles"
      a % molprop % property => computeMolecularDipoles
    end if
    
    if (flagExists(actionCommand,"+dipoleZ ")) then
      n = n + 1
      a % molprop % command  = "+dipoleZ"
      a % molprop % name     = "Molecular dipoles (polar angle)"
      a % molprop % property => computeMolecularDipolesProjection
    end if
    
    if (flagExists(actionCommand,"+normalZ ")) then
      n = n + 1
      a % molprop % command  = "+normalZ"
      a % molprop % name     =  "Normal vector (polar angle)"
      a % molprop % property => computeNormalProjection
      call assignFlagValue(actionCommand,"+normal",normalIndices(1:3),[0,0,0])
    end if
    
    if (flagExists(actionCommand,"+torsion ")) then
      n = n + 1
      a % molprop % command  = "+torsion"
      a % molprop % name     = "Torsional angle"
      a % molprop % property => computeTorsion
      call assignFlagValue(actionCommand,"+torsion",torsionIndices(1:4),[0,0,0,0])
    end if
    
    if (flagExists(actionCommand,"+angle ")) then
      n = n + 1
      a % molprop % command  = "+angle"
      a % molprop % name     =  "Covalent angle"
      a % molprop % property => computeAngle
      call assignFlagValue(actionCommand,"+angle",angleIndices(1:3),[0,0,0])
    end if
  
    if (flagExists(actionCommand,"+bond ")) then
      n = n + 1
      a % molprop % command  = "+bond"
      a % molprop % name     =  "Covalent bond length"
      a % molprop % property => computeBond
      call assignFlagValue(actionCommand,"+bond",bondIndices(1:2),[0,0])
    end if
  
    if (flagExists(actionCommand,"+rgyr" )) then
      n = n + 1
      a % molprop % command  = "+rgyr"
      a % molprop % name     = "Radius Of Gyration"
      a % molprop % property => computeRadiusOfGyration
    end if
      
    if (n == 0) call message(-1,"--molprop | no property specified")
    if (n /= 1) call message(-1,"--molprop | only one property per action can be specified")
    
    ! Extract to do with the properties
    call assignFlagValue(actionCommand,"+dumpMol ",dumpProperty   , .false.) ! Output dump
    call assignFlagValue(actionCommand,"+dumpSeq ",dumpPropertySeq, .false.) ! Dump sequentialy

    call assignFlagValue(actionCommand,"+avg "    ,averageProperty, .false.) ! Average
    call assignFlagValue(actionCommand,"+histo "  ,histoProperty  , .false.) ! Histogram
    call assignFlagValue(actionCommand,"+dist "   ,distProperty   , .false.) ! Distribution (normalised)
    call assignFlagValue(actionCommand,"+distZ "  ,distZProperty  , .false.) ! Distribution of Averages (z)
    call assignFlagValue(actionCommand,"+dist2D " ,distZProperty2D, .false.) ! Distribution of Averages (z)

    call assignFlagValue(actionCommand,"+deg " ,degrees, .false.) ! Distribution of Averages (z)

    if (dumpProperty) then 
      call workData % initialise(ID, "dump", iounit=outputFile % funit)
      
    else if (dumpPropertySeq) then 
      call workData % initialise(ID, "dumpSeq", iounit=outputFile % funit)
      
    else if (averageProperty) then
      allocate(a % array1D(2))
      call workData % initialise(ID, "average")

    else if (distProperty) then
      call assignFlagValue(actionCommand,"+nbin ",numberOfBins,100)
      call assignFlagValue(actionCommand,"+dist ", distributionLimits(1:2),[1.0_real64,2.0_real64])
      call workData % initialise(ID, "histogram", numberOfBins=[numberOfBins], lowerLimits=[distributionLimits(1)], upperLimits=[distributionLimits(2)])

    else if (histoProperty) then
      call assignFlagValue(actionCommand,"+nbin ",numberOfBins,100)
      call assignFlagValue(actionCommand,"+histo ", distributionLimits(1:2),[1.0_real64,2.0_real64])
      call workData % initialise(ID, "histogram", numberOfBins=[numberOfBins], lowerLimits=[distributionLimits(1)], upperLimits=[distributionLimits(2)])

    else if (distZProperty) then
      call assignFlagValue(actionCommand,"+nbin ",numberOfBins,100)
      call workData % initialise(ID, "avgDist ", numberOfBins=[numberOfBins], lowerLimits=[0.0_real64], upperLimits=[1.0_real64])


    else if (distZProperty2D) then
      call assignFlagValue(actionCommand,"+dist2D ",distributionLimits(1:2),[0.0_real64,0.0_real64])
      call assignFlagValue(actionCommand,"+nbin ",numberOfBins2D(1),100)
      call assignFlagValue(actionCommand,"+nbinZ ",numberOfBins2D(2),nint(frame % hmat(3,3)))
      call workData % initialise(ID, "dist2D ", numberOfBins=[numberOfBins2D], lowerLimits=[distributionLimits(1),0.0_real64], upperLimits=[distributionLimits(2),1.0_real64])

    else 
      dumpPropertySeq = .true.
      call workData % initialise(ID, "dumpSeq", iounit=outputFile % funit)

    end if

    tallyExecutions = 0 

  end subroutine initialiseAction

  subroutine dumpScreenInfo(a)
    use moduleMessages 
    implicit none
    type(actionTypeDef), target :: a

    call message(0,"Computing molecular properties")
    call message(0,"...Output file",str=outputFile % fname)
    
    call message(0,"...Molecule's type selected",str=targetMolecule)
    call message(0,"...Number of molecules selected",i=numberOfLocalMolecules)
    
    call message(0,"...Computing -> "//trim(a % molprop % name))
    if (a % molprop % command == "+torsion") call message(0,"......Atom indices used",iv=torsionIndices(1:4))
    if (a % molprop % command == "+normal")  call message(0,"......Atom indices used",iv=normalIndices(1:3))
    if (a % molprop % command == "+angle")  call message(0,"......Atom indices used",iv=angleIndices(1:3))
    if (a % molprop % command == "+bond")  call message(0,"......Atom indices used",iv=bondIndices(1:2))
    
    if (averageProperty) then
      call message(0,"...Computing the average property")
    else if (dumpProperty) then
      call message(0,"...Writing the property of each molecule")
    else if (dumpPropertySeq) then
      call message(0,"...Writing the property of each molecule sequentially")
    else if (distProperty) then
      call message(0,"...Computing the normalised distribution")
      call message(0,"......Distribution limits",rv=distributionLimits(1:2))
      call message(0,"......Number of bins",i=numberOfBins)
    else if (histoProperty) then
      call message(0,"...Computing the histogram")
      call message(0,"......Distribution limits",rv=distributionLimits(1:2))
      call message(0,"......Number of bins",i=numberOfBins)
    else if (distZProperty) then
      call message(0,"...Computing the distribution of averages along z")
      call message(0,"......Number of bins",i=numberOfBins)
    else if (distZProperty2D) then
      call message(0,"...Computing the distribution of the property along z")
      call message(0,"......Distribution limits",rv=distributionLimits)
      call message(0,"......Number of bins for the property",i=numberOfBins2D(1))
      call message(0,"......Number of bins for z",i=numberOfBins2D(2))
    end if
    call message(2)

  end subroutine dumpScreenInfo

  subroutine finaliseAction()
    use moduleSystem
    use moduleMessages
    implicit none
    real(real64), allocatable, dimension(:) :: values
    character(len=100) :: str

    if ( trim(outputFile % fname) == "stdout") then
      call message(0,"Molecular properties")
      call workData % extract(ID, values)
      if (averageProperty) call message(0,"...Average value",r=values(1))
      
    else
      call message(0,"# Molecular properties")
      if (distZProperty2D) then
         call workData % dump(ID, outputFile % funit, upperLimits=[distributionLimits(2), frame % hmat(3,3)])
      else if (distZProperty) then
        call workData % dump(ID, outputFile % funit, upperLimits=[frame % hmat(3,3)])
      else if (histoProperty) then
        call workData % dump(ID, outputFile % funit, normalisation="none")
      else if (distProperty) then
        call workData % dump(ID, outputFile % funit, normalisation="probability")
      else
        call workData % dump(ID, outputFile % funit, normalisation="none")
      end if
      close(outputFile % funit)

    end if
    call message(2)

  end subroutine finaliseAction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Molecule's properties start here !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! here the dipole is defined as going from the negative to the positive charge (iatm)
  subroutine computeMolecularDipoles(numberOfLocalMolecules,list,property) 
    use moduleSystem 
    implicit none
    integer, intent(in) :: numberOfLocalMolecules
    integer, dimension(:), intent(in) :: list
    real(real64), dimension(:), intent(out) :: property

    integer :: ilocal, imol, iatm, jatm, i
    real(real64) :: dipoleVector(3), dij(3), rtmp

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (ilocal, imol, iatm, jatm, i, dipoleVector, dij, rtmp)
!$OMP DO
    do ilocal=1,numberOfLocalMolecules
      imol = list(ilocal)
      dipoleVector = 0.0_real64
      iatm = listOfMolecules(imol) % listOfAtoms(1)
      do i=2,listOfMolecules(imol) % numberOfAtoms
        jatm = listOfMolecules(imol) % listOfAtoms(i)
        dij = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
        dipoleVector(1:3) = dipoleVector(1:3) + frame % chg(jatm) * dij
      enddo
      rtmp = sqrt(sum(dipoleVector*dipoleVector))
      property(ilocal) = rtmp * 4.8032047_real64
    enddo
!$OMP END DO
!$OMP END PARALLEL

  end subroutine computeMolecularDipoles

  subroutine computeMolecularDipolesProjection(numberOfLocalMolecules,list,property) 
    use moduleSystem 
    implicit none
    integer, intent(in) :: numberOfLocalMolecules
    integer, dimension(:), intent(in) :: list
    real(real64), dimension(:), intent(out) :: property

    integer :: ilocal, imol, iatm, jatm, i
    real(real64) :: dipoleVector(3), dij(3), rtmp
    
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (ilocal, imol, iatm, jatm, i, dipoleVector, dij, rtmp)
!$OMP DO
    do ilocal=1,numberOfLocalMolecules
      imol = list(ilocal)
      dipoleVector = 0.0_real64
      iatm = listOfMolecules(imol) % listOfAtoms(1)
      do i=2,listOfMolecules(imol) % numberOfAtoms
        jatm = listOfMolecules(imol) % listOfAtoms(i)
        dij = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
        dipoleVector(1:3) = dipoleVector(1:3) + frame % chg(jatm) * dij
      enddo
      ! z component of the molecular dipole
      rtmp = sqrt(sum(dipoleVector*dipoleVector))
      property(ilocal) = dipoleVector(3) * 4.8032047_real64
    enddo
!$OMP END DO
!$OMP END PARALLEL

  end subroutine computeMolecularDipolesProjection

  subroutine computeTorsion(numberOfLocalMolecules,list,property) 
      use moduleActions, only : actionTypeDef
      use moduleSystem 
      implicit none
      integer, intent(in) :: numberOfLocalMolecules
      integer, dimension(:), intent(in) :: list
      real(real64), dimension(:), intent(out) :: property

      ! real(real64), external :: computeTorsionalAngle
      real(real64) :: torsionPositions(3,4)

      integer :: ilocal, i, j, imol

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (ilocal, imol, i, j, torsionPositions)
!$OMP DO
     do ilocal=1,numberOfLocalMolecules
        imol = list(ilocal)
        do i=1,4
          j = torsionIndices(i)
          torsionPositions(1:3,i) = frame % pos(1:3,listOfMolecules(imol) % listOfAtoms(j))
        enddo
        property(ilocal) = computeTorsionalAngle(torsionPositions)
      enddo
!$OMP END DO
!$OMP END PARALLEL

  end subroutine computeTorsion

  function computeTorsionalAngle(pos) result(phi)
    use moduleVariables
    implicit none
    real(real64), dimension(3,4), intent(in) :: pos
    real(real64), dimension(3) :: d12, d23, d34, nvec1, nvec2, nvec3
    real(real64) :: cosphi, sinphi, phi

    d12(1:3) = pos(1:3,2) - pos(1:3,1)
    d23(1:3) = pos(1:3,3) - pos(1:3,2)
    d34(1:3) = pos(1:3,4) - pos(1:3,3)

    ! Get the normal vectors
    nvec1 = cross_product(d12,d23)
    nvec2 = cross_product(d23,d34)

    ! Calculate the dihedral angle
    nvec3 = cross_product(nvec1,nvec2)

    cosphi = dot_product(nvec1,nvec2)
    sinphi = dot_product(nvec3,d23) / sqrt(sum(d23*d23))
    phi    = atan2(sinphi,cosphi)
    if (phi < 0.0_real64 ) phi = phi + twopi

  end function computeTorsionalAngle

  subroutine computeNormalProjection(numberOfLocalMolecules,list,property) 
      use moduleActions, only : actionTypeDef
      use moduleSystem 
      use moduleMessages
      implicit none
      integer, intent(in) :: numberOfLocalMolecules
      integer, dimension(:), intent(in) :: list
      real(real64), dimension(:), intent(out) :: property

      integer :: ilocal, imol, iatm, jatm, katm
      real(real64), dimension(3) :: dij, dik, normalVector

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (ilocal, imol, iatm, jatm, katm, dij, dik, normalVector)
!$OMP DO
      do ilocal=1,numberOfLocalMolecules
        imol = list(ilocal)
      
        iatm = listOfMolecules(imol) % listOfAtoms(normalIndices(1))
        jatm = listOfMolecules(imol) % listOfAtoms(normalIndices(2))
        katm = listOfMolecules(imol) % listOfAtoms(normalIndices(3))
        dij = frame % pos(1:3,iatm) - frame % pos(1:3,jatm)
        dik = frame % pos(1:3,iatm) - frame % pos(1:3,katm)
        normalVector(1:3) = cross_product(dij,dik)
      
        normalVector(1:3) = cross_product(dij,dik)
        if (sqrt(sum(normalVector*normalVector)) < 1.0e-4_real64) &
          call message(-1,"Cannot compute normal, bonds are parallel")
          
        property(ilocal) = normalVector(3) / sqrt(sum(normalVector*normalVector))
      
      enddo
!$OMP END DO
!$OMP END PARALLEL

  end subroutine computeNormalProjection

  subroutine computeAngle(numberOfLocalMolecules,list,property) 
      use moduleActions, only : actionTypeDef
      use moduleSystem 
      use moduleMessages
      implicit none
      integer, intent(in) :: numberOfLocalMolecules
      integer, dimension(:), intent(in) :: list
      real(real64), dimension(:), intent(out) :: property

      integer :: ilocal, imol, iatm, jatm, katm
      real(real64), dimension(3) :: dij, dik
      real(real64) :: dprod

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (ilocal, imol, iatm, jatm, katm, dij, dik, dprod)
!$OMP DO
      do ilocal=1,numberOfLocalMolecules
        imol = list(ilocal)
      
        jatm = listOfMolecules(imol) % listOfAtoms(angleIndices(1))
        iatm = listOfMolecules(imol) % listOfAtoms(angleIndices(2))
        katm = listOfMolecules(imol) % listOfAtoms(angleIndices(3))
        dij = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
        dik = frame % pos(1:3,katm) - frame % pos(1:3,iatm)
        
        dij = dij / sqrt(sum(dij*dij))
        dik = dik / sqrt(sum(dik*dik))

        dprod = dot_product(dij,dik)
        property(ilocal) = acos(dprod)
      
      enddo
!$OMP END DO
!$OMP END PARALLEL

  end subroutine computeAngle

  subroutine computeBond(numberOfLocalMolecules,list,property) 
    use moduleActions, only : actionTypeDef
    use moduleSystem 
    use moduleMessages
    implicit none
    integer, intent(in) :: numberOfLocalMolecules
    integer, dimension(:), intent(in) :: list
    real(real64), dimension(:), intent(out) :: property

    integer :: ilocal, imol, iatm, jatm, katm
    real(real64), dimension(3) :: dij, dik
    real(real64) :: dprod

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (ilocal, imol, iatm, jatm, katm, dij, dik, dprod)
!$OMP DO
    do ilocal=1,numberOfLocalMolecules
      imol = list(ilocal)
    
      jatm = listOfMolecules(imol) % listOfAtoms(bondIndices(1))
      iatm = listOfMolecules(imol) % listOfAtoms(bondIndices(2))
      dij = frame % pos(1:3,jatm) - frame % pos(1:3,iatm)
      property(ilocal) = sqrt(sum(dij*dij))
    
    enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine computeBond

subroutine computeRadiusOfGyration(numberOfLocalMolecules,list,property) 
    use moduleActions, only : actionTypeDef
    use moduleSystem 
    use moduleElements
    use moduleMessages
    implicit none
    integer, intent(in) :: numberOfLocalMolecules
    integer, dimension(:), intent(in) :: list
    real(real64), dimension(:), intent(out) :: property
    real(real64) :: xcom(3)
    integer :: ilocal, imol, iatm, n, idx
    real(real64) :: rgyr
    real(real64) :: dx(3)
    
    logical, save :: firstTimeIn = .true.
    real(real64), save, allocatable, dimension(:) :: rmass
    real(real64), save :: moleculeMass
    
    if (firstTimeIn) then
      firstTimeIn = .false.
      
      imol = list(1)
      n = listOfMolecules(imol) % numberOfAtoms
      allocate(rmass(n))
      
      moleculeMass = 0.0_real64
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        iatm = listOfMolecules(imol) % listOfAtoms(idx)
        rmass(idx) = getElementMass(frame % lab(iatm))
        moleculeMass = moleculeMass + rmass(idx)
      enddo
    end if

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (ilocal, xcom, imol, idx, rgyr, dx)
!$OMP DO
    do ilocal=1,numberOfLocalMolecules
      xcom = 0.0_real64
      imol = list(ilocal)
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        iatm = listOfMolecules(imol) % listOfAtoms(idx)
        xcom = xcom + rmass(idx) * frame % pos(:,iatm) 
      enddo
      xcom = xcom / moleculeMass
      
      rgyr = 0.0_real64
      do idx=1,listOfMolecules(imol) % numberOfAtoms
        iatm = listOfMolecules(imol) % listOfAtoms(idx)
        dx(1:3) = frame % pos(1:3,iatm) - xcom(1:3)
        rgyr = rgyr + rmass(idx) * sum(dx*dx)
      enddo
      property(ilocal) = sqrt(rgyr / moleculeMass)
    enddo
!$OMP END DO
!$OMP END PARALLEL
    
  end subroutine computeRadiusOfGyration
  
end module moduleMolecularProperties

subroutine getLocalIndices(targetMolecule, a, numberOfLocalMolecules)
  use moduleSystem
  use moduleResizeArrays
  use moduleStrings
  use moduleVariables, only : actionTypeDef
  implicit none
  type(actionTypeDef), target :: a

  character(len=*), intent(in) :: targetMolecule
  integer, intent(out) :: numberOfLocalMolecules
 
  integer :: imol
  integer :: numberOfWords
  character(len=STRLEN), dimension(100) :: listOfWords

  if (allocated(a % localIndices)) deallocate(a % localIndices)

  call parse(targetMolecule,",",listOfWords,numberOfWords)
  allocate(a % localIndices(numberOfMolecules))
  
  if (targetMolecule=="NULL") then

    numberOfLocalMolecules = numberOfMolecules
    do imol=1,numberOfMolecules
      a % localIndices(imol) = imol
    enddo
    
  else
    
    numberOfLocalMolecules = 0
    do imol=1,numberOfMolecules
      if (all(listOfWords(1:numberOfWords) /= listOfMolecules(imol) % resname)) cycle
      numberOfLocalMolecules = numberOfLocalMolecules + 1
      a % localIndices(numberOfLocalMolecules) = imol
    enddo
    call resizeArray(a % localIndices, numberOfLocalMolecules)

  end if

end subroutine getLocalIndices

