! Copyright (c) 2021, Paolo Raiteri, Curtin University.
! All rights reserved.
! 
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License as published by the 
! Free Software Foundation; either version 3 of the License, or 
! (at your option) any later version.
!  
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! * Redistributions of source code must retain the above copyright notice, 
!   this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright notice, 
!   this list of conditions and the following disclaimer in the documentation 
!   and/or other materials provided with the distribution.
! * Neither the name of the <ORGANIZATION> nor the names of its contributors 
!   may be used to endorse or promote products derived from this software 
!   without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
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
  
subroutine initialiseFile(file,fname,fposition,fstatus,fformat)
  use moduleVariables
  use moduleMessages 
  implicit none
  type(fileTypeDef), intent(out) :: file
  character(len=*), intent(in) :: fname
  character(len=*), intent(in), optional :: fposition
  character(len=*), intent(in), optional :: fstatus
  character(len=*), intent(in), optional :: fformat
  
  integer :: ilen
  character(fp) :: ftype
  character(len=11) :: fform
  
  character(len=6) :: fpos
  character(len=7) :: fstat
  
  integer :: idx, ierr
  
  file % first_access = .true.
  
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
