!disclaimer
subroutine parseActionLine()
  use moduleActions, only : numberOfActions, actionType, actionDetails
  use moduleMessages
  use moduleHelp

  implicit none

  integer :: iarg
  integer :: nargs
  character(len=200) :: word

  ! Number of command line arguments
  nargs=iargc()

  ! Count the commands
  numberOfActions=0
  iarg=0
  do while(iarg<nargs)
    iarg = iarg + 1
    call getarg(iarg,word)
    if (word(1:2)=="--") then
      numberOfActions=numberOfActions+1
      actionType(numberOfActions) = trim(word)
      actionDetails(numberOfActions) = ""
    else 
      if (word(1:1)=="-1") then
        if ( ((ichar(word(2:2)) >= 48) .and. (ichar(word(2:2)) <= 57)) .or. (ichar(word(2:2)) == 46) )then
          actionDetails(numberOfActions) = trim(actionDetails(numberOfActions))//" "//trim(word)
        else
          write(0,*)ichar(word(2:2)),word(2:2)
          call message(-1,"Unknown command",str=word)
        end if
      else
        actionDetails(numberOfActions) = trim(actionDetails(numberOfActions))//" "//trim(word)
      end if  
    end if
  enddo

  if (any(actionType == "--help")) call help()

  do iarg=1,numberOfActions
    if ( index(actionDetails(iarg),"+help") > 0) call commandHelp(actionType(iarg))
  end do

end subroutine parseActionLine
