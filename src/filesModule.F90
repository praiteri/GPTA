!disclaimer
module moduleFiles

  implicit none

  interface init_file
    subroutine initialiseFile(file,fname,fposition,fstatus,fformat)
      use moduleVariables
      use moduleMessages 
      implicit none
      type(fileTypeDef), intent(out) :: file
      character(len=*), intent(in) :: fname
      character(len=*), intent(in), optional :: fposition
      character(len=*), intent(in), optional :: fstatus
      character(len=*), intent(in), optional :: fformat
    end subroutine initialiseFile
  end interface init_file
    
end module moduleFiles
  
subroutine initialiseFile(file,fname_input,fposition,fstatus,fformat)
  use moduleVariables
  use moduleMessages 
  use moduleRandomNumbers
  implicit none
  type(fileTypeDef), intent(out) :: file
  character(len=*), intent(in) :: fname_input
  character(len=*), intent(in), optional :: fposition
  character(len=*), intent(in), optional :: fstatus
  character(len=*), intent(in), optional :: fformat
  
  character(len=len(fname_input)) :: fname

  integer :: ilen
  character(fp) :: ftype
  character(len=11) :: fform
  
  character(len=6) :: fpos
  character(len=7) :: fstat
  
  integer :: idx, ierr
  
  file % first_access = .true.
  
  fname = fname_input

  block
    integer :: exitstat, cmdstat
    character(256) :: ff, fnew
    character(500) :: cmd
    character(256) :: cmdmsg, rnd
    
    if (index(trim(fname),"@") > 0) then
      call get_filename(trim(fname),ff)
      write(rnd,'(i5)') int(grnd()*99999) + 1

      fnew = "_tmp-gpta_"//trim(rnd)//"_"//trim(ff)
      call message(0,"Fetching remote file",str=fname)
      call message(0,"...local temporary file",str=fnew)

      cmd = "scp "//trim(fname)//" ./"//trim(fnew)//" > /dev/null 2>&1"

      call execute_command_line(cmd, .true., EXITSTAT=exitstat, CMDSTAT=cmdstat, CMDMSG=cmdmsg)

      if (cmdstat /= 0) then
          write(0,*) 'Command execution error: ', trim(cmdmsg)
      else if (exitstat /= 0) then
          write(0,*) 'Command failed with exit status:', exitstat
      end if
      call message(2)
      fname = fnew
    end if

  end block

  if (present(fposition)) then
    fpos = fposition
  else
    fpos = 'rewind'
  end if
  
  if (present(fstatus)) then
    fstat = fstatus
  else
    fstat = 'unknown'
  end if
  
  ! iounit=iounit+1
  
  ilen = len_trim(fname)
  idx = index(fname,".",back=.true.)
  if (ilen-idx>0) then
    ftype = fname(idx+1:ilen)
  else
    ftype='NULL'
  end if
  
  if (present(fformat)) then
    fform = fformat
  else
    if (ftype=="dcd") then
      fform = "unformatted"
    else
      fform = "formatted"
    end if
  end if
  
  file % fname = trim(fname)
  file % ftype = ftype
  file % fform = fform
  file % fpos  = fpos
  
  ! this is a FORTRAN08 feature
  open(newunit  = file % funit, &
       file     = file % fname, &
       status   = fstat       , &
       form     = file % fform, &
       position = file % fpos , &
       iostat   = ierr        )
  
  if (ierr/=0) then
    call message(-1,"Error opening file",str=trim(fname))
  end if
  
  return
end subroutine initialiseFile
