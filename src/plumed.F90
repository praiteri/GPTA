!disclaimer
module modulePlumedInterface
  use moduleVariables
  implicit none

  public :: plumedInterface
  private

#ifdef GPTA_PLUMED
  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  character(len=STRLEN), pointer :: plumedInputFile
  character(len=STRLEN), pointer :: plumedOutputFile

  integer, pointer :: tallyExecutions

  real(real64) :: energyUnits, lengthUnits, timeUnits
#endif

contains

subroutine initialiseAction(a)
  use moduleStrings
  use moduleSystem
  implicit none
  type(actionTypeDef), target :: a

#ifdef GPTA_PLUMED
  a % actionInitialisation = .false.
  a % requiresNeighboursList = .false.
  a % requiresNeighboursListUpdates = .false.
  a % requiresNeighboursListDouble = .false.
  a % cutoffNeighboursList = 1.0_real64

  call assignFlagValue(actionCommand,"+f",plumedInputFile,'NULL')
  call assignFlagValue(actionCommand,"+out",plumedOutputFile,'plumed_gpta.out')

  call dumpScreenInfo()
  tallyExecutions = 0 
#else
  call message(-1,"Recompile with the interface to PLUMED 2.x")
#endif
  end subroutine initialiseAction

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

#ifdef GPTA_PLUMED
    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    plumedInputFile      => a % stringVariables(1)
    plumedOutputFile     => a % stringVariables(2)
#endif
  end subroutine associatePointers

#ifdef GPTA_PLUMED
  subroutine dumpScreenInfo()
    use moduleMessages 
    implicit none
    call message(0,"Creating PLUMED interface")
  end subroutine dumpScreenInfo
#endif

  subroutine plumedInterface(a)
    use moduleVariables
    use moduleElements
    use moduleSystem
    use moduleMessages
    implicit none
    type(actionTypeDef), target :: a

#ifdef GPTA_PLUMED

    integer :: plumedVariable
    real(real64), pointer :: mass(:), force(:,:), virial(:,:)

    integer :: i
    real(real64) :: etot

    call associatePointers(a)
    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      if (firstAction) then
        energyUnits = 96.48530749925792_real64
        lengthUnits = 0.1_real64
        timeUnits = 1.0_real64

        call plumed_f_installed(plumedVariable)
        if(plumedVariable<=0) call message(-1,"PLUMED 2 not available")
        
        if (plumedInputFile=="NULL") call message(-1,"Plumed input file must be specified with +f")
        
        call message(1,"creating plumed from GPTA")
        call plumed_f_gcreate()
        call plumed_f_gcmd("setRealPrecision"//char(0),8)
        call plumed_f_gcmd("setMDEnergyUnits"//char(0),energyUnits)
        call plumed_f_gcmd("setMDLengthUnits"//char(0),lengthUnits)
        call plumed_f_gcmd("setMDTimeUnits"//char(0),timeUnits)
        call plumed_f_gcmd("setNatoms"//char(0),frame % natoms)
        call plumed_f_gcmd("setPlumedDat"//char(0),trim(plumedInputFile)//char(0))
        call plumed_f_gcmd("setLogFile"//char(0),trim(plumedOutputFile)//char(0))
        call plumed_f_gcmd("setMDEngine"//char(0),"driver");
        call plumed_f_gcmd("setTimestep"//char(0),1.0_real64);
        call plumed_f_gcmd("init"//char(0),0);

        call dumpScreenInfo()

        call checkUsedFlags(actionCommand)
        firstAction = .false.

      end if

      allocate(mass(frame % natoms))
      allocate(force(3,frame % natoms), source=0.0_real64)
      allocate(virial(3,3), source=0.0_real64)
      etot=0.0_real64
      do i=1,frame % natoms
        mass(i) = getElementMass(frame % lab(i))
      enddo
      CALL plumed_f_gcmd("setStep"//char(0),tallyExecutions)
      CALL plumed_f_gcmd("setMasses"//char(0),mass)
      CALL plumed_f_gcmd("setForces"//char(0),force)
      CALL plumed_f_gcmd("setVirial"//char(0),virial)
      CALL plumed_f_gcmd("setEnergy"//char(0),etot)
      CALL plumed_f_gcmd("setPositions"//char(0),frame % pos)
      CALL plumed_f_gcmd("setBox"//char(0),frame % hmat)
      CALL plumed_f_gcmd("calc"//char(0),0)
      deallocate(mass,force,virial)

    end if

#endif
  end subroutine plumedInterface

end module modulePlumedInterface
