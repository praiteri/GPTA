!disclaimer
module moduleMessages 

  implicit none

  private
  public :: message

  interface message_int
    subroutine message(istop,comment,icol,i,r,iv,rv,str,l)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none
      integer, intent(in)           :: istop
      character(len=*), optional    :: comment
      integer, intent(in), optional :: icol
      integer, intent(in), optional :: i
      real(real64), intent(in), optional :: r
      integer, dimension(:), intent(in), optional :: iv
      real(real64), dimension(:), intent(in), optional :: rv
      character(len=*), intent(in), optional :: str
      logical, intent(in), optional :: l
    end subroutine message      
  end interface message_int

end module moduleMessages 

subroutine message(istop,comment,icol,i,r,iv,rv,str,l)
  use moduleVariables

  implicit none
  integer, intent(in)                         :: istop
  character(len=*)                 , optional :: comment
  integer, intent(in)              , optional :: icol
  integer, intent(in)              , optional :: i
  real(real64), intent(in)              , optional :: r
  integer, dimension(:), intent(in), optional :: iv
  real(real64), dimension(:), intent(in), optional :: rv
  character(len=*),      intent(in), optional :: str
  logical, intent(in)              , optional :: l
  
  integer, parameter :: line_length=110
  integer            :: ndashes
  character(len=30)  :: fmtdash
  character(len=30)  :: fmtdash1
  integer            :: cpos

  ! character(len=line_length), parameter :: fmtstring='(2x,a50," : ",a)'
  character(len=line_length), parameter :: fmti='(i11)'
  character(len=line_length) :: fmt0
  character(len=line_length) :: fstr

  character(len=line_length*2) :: line
  character(len=line_length) :: str0=""
  character(len=line_length) :: str1=""
  character(len=line_length) :: str2=""
  character(len=line_length) :: str3=""

  integer :: idx, n
  ! character(len=line_length), external :: get_format_r
  character(len=256) :: fmt_tmp

  logical, save :: lastDash = .false.
  logical, save :: lastFrame = .false.
  logical :: commentOnly

  if (me /= 0) return

  commentOnly = .false.

!!!!
  ndashes = line_length
  cpos = int(line_length/2)
  write(fstr,'("(a",i0,")")') line_length
  write(fmtdash,'("(",i0,"a1)")') line_length
  write(fmtdash1,'("(",i0,"a1)")') 72
!!!!

  if (istop==-10) then
    write(io,'(a)')"+--------------------------------------------------------------------------------------------------------+"
    write(io,'(a)')"|     General Purpose Trajectory Analyser (GPTA)                                                         |"
    write(io,'(a)')"|     Copyright (C) 2021 - Paolo Raiteri - Curtin University                                             |"
    write(io,'(a)')"|                                                                                                        |"
    write(io,'(a)')"|     This program is free software: you can redistribute it and/or modify it under the terms of the     |"
    write(io,'(a)')"|     GNU General Public License as published by the Free Software Foundation, either version 3 of       |"
    write(io,'(a)')"|     the License, or (at your option) any later version.                                                |"
    write(io,'(a)')"|                                                                                                        |"
    write(io,'(a)')"|     This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;          |"
    write(io,'(a)')"|     without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.          |"
    write(io,'(a)')"|     See the GNU General Public License for more details.                                               |"
    write(io,'(a)')"+--------------------------------------------------------------------------------------------------------+"
    return
  end if  

  if (istop==-1) write(io,*)"--- ERROR DETECTED ---"

  ! delete line
  if (istop == -2) then
    if (io == 6) call deleteLine(io,line_length)
    return
  end if

  ! Avoid writing double lines
  if (lastDash .and. istop == 2) return

  if (istop==1 .or. istop==2 .or. istop==4 .or. lastFrame) then
    lastDash = .true.
  else
    lastDash = .false.
  end if

  if (istop == -3) then
#ifdef GPTA_MPI
    call MPI_Abort(MPI_COMM_WORLD, 9, idx )
#else
    call exit(1)
#endif
  end if

  line = ""
  str0 = ""
  str1 = ""
  str2 = ""
  str3 = ""

  ! write only the dashes
  if (istop==2 .or. istop==4) then
    write(io,fmtdash)("_",idx=1,ndashes)
    return

  else if (istop==3) then
    write(io,fmtdash1)("-",idx=1,72)
    return
    
  else
    line(1:cpos-1)  = trim(comment)
    line(cpos:cpos+1) = ": "
    if (present(icol)) then
      write(str0,'(i0)') icol
      line(1:cpos-1)  = trim(comment)//" "//trim(str0)
    endif

    if (present(i)) then
      fmt0 = '('//trim(fmti)//')'
      write(str1,fmt0) i
      line(cpos+2:) = trim(str1)

    else if (present(r)) then
      fmt0 = "("
      call get_format_r(r, fmt_tmp)
      fmt0 = trim(fmt0)//trim(fmt_tmp)
      fmt0 = trim(fmt0)//")"
      write(str1,fmt0) r
      line(cpos+2:) = trim(str1)

    else if (present(iv)) then
      n = size(iv)
      fmt0 = "("
      do idx=1,n
        fmt0 = trim(fmt0)//trim(fmti)//",4x,"
      enddo
      fmt0 = trim(fmt0)//")"
      write(str1,fmt0) iv
      call compact(str1)
      line(cpos+2:) = trim(str1)

    else if (present(rv)) then
      n = size(rv)
      fmt0 = "('[',"
      do idx=1,n
        call get_format_r(rv(idx), fmt_tmp)
        fmt0 = trim(fmt0)//trim(fmt_tmp)//","
        if (idx<n) fmt0 = trim(fmt0)//"',',2x,"
      enddo
      fmt0 = trim(fmt0)//"']')"
      write(str1,fmt0) rv
      line(cpos+2:) = trim(str1)

    else if (present(str)) then
      write(str1,fstr) trim(str)
      call compact(str1)
      line(cpos+2:) = trim(str1)

    else if (present(l)) then
      write(str1,'(a)') trim(str)
      if (l) then
        line(cpos+9:) = "TRUE"
      else
        line(cpos+8:) = "FALSE"
      end if
        
    else
      commentOnly = .true.
      line(cpos:cpos) = " "
    end if
    
    if (line(cpos:cpos) == ":") then
      n = len_trim(line(1:cpos-1))
      do idx=n+1,cpos-1
        line(idx:idx)="."
      enddo
    end if

    if (index(comment,"---Frames") > 0) then
      write(io,'(a)',advance='no') trim(line)
      lastFrame = .true.
    else
      if (lastFrame) then
        if (io == 6) then
          call deleteLine(io,line_length)
        else
          write(io,fmtdash)("_",idx=1,ndashes)
        end if
      end if
      if (commentOnly) then
        write(io,'(a)') trim(comment)
      else
        write(io,'(a)') trim(line)
      endif
      lastFrame = .false.
    end if
    
    if (istop==1) write(io,fmtdash)("_",idx=1,ndashes)
  end if 
  call flush(io)

  if (istop<0) then
#ifdef GPTA_MPI
    call MPI_Abort(MPI_COMM_WORLD, 9, idx )
#else
    call exit(1)
#endif
  end if 

  return
end subroutine message

subroutine deleteLine(io,line_length)
  implicit none
  integer, intent(in) :: io, line_length
  integer :: idx
  character(len=1) :: back=char(8)

  do idx=1,line_length
    write(io,'(a1)',advance='no')back
  enddo
  write(io,'(100x)',advance='no')
  do idx=1,line_length
    write(io,'(a1)',advance='no')back
  enddo

end subroutine deleteLine

subroutine get_format_r(r, fmt0)
  use moduleVariables, only: real64
  implicit none
  real(real64), intent(in) :: r
  character(len=*), intent(out) :: fmt0
  character(len=5) :: str1, str2

  integer :: nd, ntot
  ntot = 10

  if ( abs(r) < tiny(1.0_real64) ) then
    nd = 0
  else
    nd = int(log10(abs(r)))
  end if

  if (nd < -4 .or. nd >= 5) then
    write(str1,'(i0)') ntot+5
    write(str2,'(i0)') ntot-4
    fmt0 = 'e'//trim(str1)//'.'//trim(str2)//',1x,",",1x'
  else
    fmt0 = "f11.4,2x"
  end if
  
  return
end subroutine