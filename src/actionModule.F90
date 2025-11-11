!disclaimer
module moduleActions
  use moduleVariables
  implicit none

  integer :: numberOfActions
  character(len=actionNameLength), dimension(maximimNumberOfActions) :: actionType

  type(actionTypeDef), allocatable, dimension(:), target :: action

  interface
    subroutine command(action)
      use moduleVariables
      implicit none
      type(actionTypeDef), target :: action
    end subroutine command
  end interface
  ! procedure (command), pointer :: work => null ()
  procedure (command), pointer :: extraWork => null ()
  
  type :: allActionsProcedure
    procedure (command), pointer, nopass :: work => null ()
  end type allActionsProcedure
  type(allActionsProcedure) :: allActions(100)

end module moduleActions
