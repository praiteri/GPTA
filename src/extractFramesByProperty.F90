!disclaimer
module moduleExtractFramesByProperty
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleFiles
  use moduleMessages
  use moduleProperties
  use moduleResizeArrays
  
  implicit none

  public :: extractFramesByProperty, extractFramesByPropertyHelp
  private

  character(:), pointer :: actionCommand
  logical, pointer :: firstAction
  type(fileTypeDef), pointer :: inputfile
  type(fileTypeDef), pointer :: outputfile
  integer, pointer :: tallyExecutions
  real(real64), pointer, dimension(:) :: limits

  integer, pointer :: nProperties
  real(real64), pointer, dimension(:,:) :: properties
  logical, pointer :: header 

contains

  subroutine extractFramesByPropertyHelp()
    implicit none
    call message(0,"This action extracts frames based on properties read from a file.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --framesByProperty +f colvar_gpta +out out.1.dcd \\")
    call message(0,"     +min *,1.15,-1.0,0.0 +max *,1.30,-0.9,0.2")
  end subroutine extractFramesByPropertyHelp


  subroutine extractFramesByProperty(a)
    implicit none
    type(actionTypeDef), target :: a

    call associatePointers(a)

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then
      tallyExecutions = tallyExecutions + 1

      ! Atoms' selection
      if (firstAction) then
        ! write info 
        call dumpScreenInfo()
 
        ! Throw a warning for unused flags
        call checkUsedFlags(actionCommand)
        firstAction = .false.          
      end if

      call computeAction(a)
    end if

    if (endOfCoordinatesFiles) then
      ! call finaliseAction(a)
    end if 

  end subroutine extractFramesByProperty

  subroutine associatePointers(a)
    implicit none
    type(actionTypeDef), target :: a

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    tallyExecutions      => a % tallyExecutions
    inputfile            => a % inputfile
    outputfile           => a % outputfile
    limits(1:)           => a % doubleVariables(1:)
    nProperties          => a % integerVariables(1)
    properties           => a % array2D
    header               => a % logicalVariables(1)

  end subroutine associatePointers

  subroutine initialiseAction(a)
    implicit none
    type(actionTypeDef), target :: a
    character(STRLEN) :: flagString
    character(STRLEN), dimension(100) :: str
    integer :: i, nval
    
    a % actionInitialisation = .false.
    a % requiresNeighboursList = .true.
    a % requiresNeighboursListUpdates = .false.
    a % requiresNeighboursListDouble = .true.
    a % cutoffNeighboursList = 3.0_real64

    call assignFlagValue(actionCommand,"+f",inputfile % fname,'properties')
    call assignFlagValue(actionCommand,"+out",outputfile % fname,'out.pdb')
    call initialiseFile(outputFile, outputFile % fname)

    block
      integer :: io, ierr, nl, i
      character(len=200) :: line, words(10)
      open(newunit=io, file= a % inputfile % fname, status='old')

      nProperties = 0
      ! count the number of lines
      nl = 0
      do
        read(io,'(a200)',iostat=ierr)line
        if (index(line,"#") > 0) cycle ! skip comments
        if (len_trim(line) == 0) cycle ! skip empty lines
        if (ierr/=0) exit ! exit if EOF or EOL
        nl = nl + 1
        if (nProperties == 0) then
          do nProperties=1,20
            read(line,*,iostat=ierr)words(1:nProperties)
            if (ierr/=0) exit
          end do
          nProperties = nProperties - 1
        end if
      end do
      rewind(io)

      ! count the number of properties in the line
      allocate(a % array2D(nProperties,nl))

      i=0
      do while (i<nl)
        read(io,'(a200)',iostat=ierr)line
        if (index(line,"#") > 0) cycle ! skip comments
        if (len_trim(line) == 0) cycle ! skip empty lines
        i = i + 1
        read(line,*,iostat=ierr)a % array2D(:,i)
        if (ierr/=0) then
          call message(-1,line)
        endif
      end do
      close(io)
    end block

    call extractFlag(actionCommand,"+min",flagString)
    call parse(flagString, ",", str, nval)
    if (nval /= nProperties) call message(-1,"Number of values /= from number of columns",iv=[nval,nProperties])
    do i=1,nProperties
      if (str(i) == "*" .or. str(i) == "INF") then
        limits(i) = minval(a % array2D(i,:)) - 1.0_real64
      else
        read(str(i),*)limits(i)
      end if
    end do

    call extractFlag(actionCommand,"+max",flagString)
    call parse(flagString, ",", str, nval)
    if (nval /= nProperties) call message(-1,"Number of values /= from number of columns",iv=[nval,nProperties])
    do i=1,nProperties
      if (str(i) == "*" .or. str(i) == "INF") then
        limits(i+nProperties) = maxval(a % array2D(i,:)) + 1.0_real64
      else
        read(str(i),*)limits(i+nProperties)
      end if
    end do
    
    header = .true.
    tallyExecutions = 0 

  end subroutine initialiseAction

  subroutine finaliseAction(a)
    implicit none
    type(actionTypeDef), target :: a

  end subroutine finaliseAction

  subroutine dumpScreenInfo()
    implicit none
    integer :: i
    character(len=5) :: str
    call message(0,"Extract frames by properties")
    call message(0,"...Properties file",str=inputfile % fname)
    call message(0,"...Number of properties in file",i=nProperties)
    do i=1,nProperties
      write(str,'(i0)')i
      call message(0,"......Limits for property #"//trim(str),rv=[limits(i),limits(i+nProperties)])
    end do
    call message(0,"...Output file",str=outputfile % fname)
    call message(2)
  end subroutine dumpScreenInfo

  subroutine computeAction(a)
    use moduleSystem

    implicit none
    type(actionTypeDef), target :: a

    integer :: i
    logical, allocatable, dimension(:) :: dumpFrame

    allocate(dumpFrame(nProperties), source=.false.)
    do i=1,nProperties
      if ( properties(i,tallyExecutions) > limits(i) .and. &
           properties(i,tallyExecutions) < limits(i+nProperties) ) dumpFrame(i) = .true.
    end do
    
    if (count(dumpFrame) == nProperties) then
      call dumpCoordinates(outputFile % ftype, &
                           outputFile % funit, &
                           frame % natoms, &
                           frame % pos, &
                           frame % lab, &
                           hmat = frame % hmat, &
                           header = header)
      header = .false.
    end if
 
  end subroutine computeAction

end module moduleExtractFramesByProperty
