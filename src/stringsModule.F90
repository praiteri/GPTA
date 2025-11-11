!disclaimer
module moduleStrings

  use moduleVariables, only : STRLEN
private :: iWord, numberOfWords, listOfWords

  integer :: iWord
  integer :: numberOfWords
  character(len=STRLEN), dimension(100) :: listOfWords

  interface readTokens
    module procedure read_dble_scalar, read_int_scalar, read_dble_vector, read_int_vector, read_string, read_chr_vector
  end interface readTokens

  interface assignFlagValue
    module procedure findScalarInteger, findVectorInteger, findArrayInteger, &
                     findScalarDouble, findVectorDouble, findArrayDouble, &
                     findScalarReal, findVectorReal, &
                     findScalarString, findVectorString, findArrayString,&
                     findLogical
  end interface 

contains

  !**********************************************************************

  subroutine removeCharacter(str,chr)
    implicit none
    character(len=*), intent(inout) :: str
    character(len=1), intent(in) :: chr
    integer :: i
    
    do i=1,len_trim(str)
      if (str(i:i) == chr) str(i:i) = " "
    end do

  end subroutine removeCharacter

  subroutine extractFlag(string,flag,result)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    character(len=*), intent(out) :: result

    integer :: iWord
    integer :: nl

    nl = len_trim(flag)+1
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        result = listOfWords(iWord)(nl:)
        exit
      end if
    end do
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

  end subroutine extractFlag

  subroutine parse(str,delims,args,nargs)
  
  ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
  ! the delimiters contained in the string 'delims'. Preceding a delimiter in
  ! 'str' by a backslash (\) makes this particular instance not a delimiter.
  ! The integer output variable nargs contains the number of arguments found.
  implicit none

  character(len=*) :: str,delims
  character(len=len_trim(str)) :: strwork
  character(len=*),dimension(:) :: args

  integer :: nargs, na, k, i, lenstr
  nargs=0
  if (len_trim(str) == 0) return
  strwork = trim(str)
  call compact(strwork)
  if (strwork(1:1)==delims .or. strwork(1:1)=="%") then
    strwork(1:1)=" "
    call compact(strwork)
  end if
  na=size(args)
  do i=1,na
    args(i)=' '
  end do

  lenstr=len_trim(strwork)
  if(lenstr==0) return
  k=0
  
  do
    if(len_trim(strwork) == 0) exit
    nargs=nargs+1
    call split(strwork,delims,args(nargs))
    call remove_special(args(nargs))
  end do
  
  end subroutine parse
  
  subroutine count_delim(str, delim, count)
    implicit none
    character(len=*), intent(in) :: str
    character(len=1), intent(in) :: delim
    integer, intent(out) :: count

    integer :: i
    count = 0

    do i = 1, len_trim(str)
      if (str(i:i) == delim) then
        count = count + 1
      end if
    end do
  end subroutine count_delim

  subroutine parse2(str,delims,args,nargs)
  
  ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
  ! the delimiters contained in the string 'delims'. Preceding a delimiter in
  ! 'str' by a backslash (\) makes this particular instance not a delimiter.
  ! The integer output variable nargs contains the number of arguments found.
  implicit none

  character(len=*) :: str,delims
  character(len=len_trim(str)) :: strwork
  character(len=*),dimension(:) :: args
  integer :: nargs
  integer :: i

  if (len_trim(str) == 0) return
  strwork = trim(str)
  call compact(strwork)

  call count_delim(strwork, delims, nargs)
  if (nargs == 0) return

  call parse(strwork,delims,args,nargs)
  do i=1,nargs
    args(i) = trim(delims) // trim(args(i))
  end do
 
  end subroutine parse2

  !**********************************************************************

  subroutine readline(nunitr,line,ios)
  
  ! Reads line from unit=nunitr, ignoring blank lines
  ! and deleting comments beginning with "#"
  implicit none
  integer, intent(in) :: nunitr
  character(len=*), intent(inout) :: line
  integer, intent(inout) :: ios
  
  integer :: ipos, jpos

  do
    read(nunitr,'(a)', iostat=ios) line      ! read input line
    if(ios /= 0) return
    line=adjustl(line)
    ipos=index(line,'#')
    jpos=index(line,'@')
    if(ipos == 1 .and. jpos /=2) cycle
    if(ipos > 1) line=line(:ipos-1)
    call compact(line)
    if(len_trim(line) /= 0) exit
  end do

  end subroutine readline


  !**********************************************************************
  
  subroutine remove_special(str)
  
  ! Removes backslash (\) characters
  ! Removes percent ( % ) characters 
  implicit none
  character(len=*):: str
  character(len=1):: ch
  character(len=len_trim(str))::outstr
  
  integer :: i, k, lenstr

  str=adjustl(str)
  lenstr=len_trim(str)
  outstr=' '
  k=0  
  do i=1,lenstr
    ch=str(i:i)
    if(ch == '\') cycle
    if(ch == ' % ') cycle
    k=k+1
    outstr(k:k)= ch
  end do
  
  str=adjustl(outstr)
  
  end subroutine remove_special

  subroutine split(str,delims,before,sep)
  
  ! Routine finds the first instance of a character from 'delims' in the
  ! the string 'str'. The characters before the found delimiter are
  ! output in 'before'. The characters after the found delimiter are
  ! output in 'str'. The optional output character 'sep' contains the
  ! found delimiter. A delimiter in 'str' is treated like an ordinary
  ! character if it is preceded by a backslash (\). If the backslash
  ! character is desired in 'str', then precede it with another backslash.
  
  implicit none
  character(len=*) :: str,delims,before
  character,optional :: sep
  logical :: pres
  character :: ch, cha
  
  integer :: i, k, ipos, ibsl, lenstr, iposa

  pres=present(sep)
  str=adjustl(str)
  call compact(str)
  lenstr=len_trim(str)
  if(lenstr == 0) return        ! string str is empty
  k=0
  ibsl=0                        ! backslash initially inactive
  before=' '
  do i=1,lenstr
     ch=str(i:i)
     if(ibsl == 1) then          ! backslash active
        k=k+1
        before(k:k)=ch
        ibsl=0
        cycle
     end if
     if(ch == '\') then          ! backslash with backslash inactive
        k=k+1
        before(k:k)=ch
        ibsl=1
        cycle
     end if
     ipos = max( index(delims,ch) , index("%",ch) )
     if(ipos == 0) then          ! character is not a delimiter
        k=k+1
        before(k:k)=ch
        cycle
     end if
     if(ch /= ' ') then          ! character is a delimiter that is not a space
        str=str(i+1:)
        if(pres) sep=ch
        exit
     end if
     cha=str(i+1:i+1)            ! character is a space delimiter
     iposa = max( index(delims,cha) , index("%",cha) )
     if(iposa > 0) then          ! next character is a delimiter
        str=str(i+2:)
        if(pres) sep=cha
        exit
     else
        str=str(i+1:)
        if(pres) sep=ch
        exit
     end if
  end do
  if(i >= lenstr) str=''
  str=adjustl(str)              ! remove initial spaces
  return
  
  end subroutine split
  
  subroutine read_dble_scalar(cmd,vec)
    use moduleMessages 
    implicit none
    character(len=*), intent(in) :: cmd
    real(8), intent(out) :: vec
    integer :: nw
    character(len=STRLEN), dimension(50) :: words
    integer :: ios

    call parse(cmd," ",words,nw)
    if (nw/=2) call message(-1,"Wrong number of arguments for +"//trim(cmd))
    if (words(2) == "pi") then
      vec = acos(-1.d0)
    else if (words(2) == "-pi") then
      vec = -acos(-1.d0)
    else if (words(2) == "twopi") then
      vec = 2.d0*acos(-1.d0)
    else if (words(2) == "-twopi") then
      vec = -2.d0*acos(-1.d0)
    else if (words(2) == "pih") then
      vec =  acos(-1.d0) / 2.d0
    else if (words(2) == "-pih") then
      vec = -acos(-1.d0) / 2.d0
    else if (words(2) == "piq") then
      vec =  acos(-1.d0) / 4.d0
    else if (words(2) == "-piq") then
      vec = -acos(-1.d0) / 4.d0
    else
      read(words(2),*,iostat=ios) vec
    end if    
    if (ios /= 0) call message(-1,"Wrong arguments for +"//trim(cmd))

    return
  end subroutine

  subroutine read_int_scalar(cmd,vec)
    use moduleMessages 
    implicit none
    character(len=*), intent(in) :: cmd
    integer, intent(out) :: vec
    integer :: nw
    character(len=STRLEN), dimension(50) :: words
    integer :: ios

    call parse(cmd," ",words,nw)
    if (nw/=2) call message(-1,"Wrong number of arguments for +"//trim(cmd))
    read(words(2),*,iostat=ios) vec
    if (ios /= 0) call message(-1,"Wrong arguments for +"//trim(cmd))

    return
  end subroutine

  subroutine read_string(cmd,vec)
    use moduleMessages 
    implicit none
    character(len=*), intent(in) :: cmd
    character(len=*), intent(out) :: vec
    integer :: nw
    character(len=STRLEN), dimension(50) :: words

    call parse(cmd," ",words,nw)
    if (nw/=2) call message(-1,"Wrong number of arguments for +"//trim(cmd))
    vec = trim(words(2))

    return
  end subroutine

  subroutine read_dble_vector(cmd,vec)
    use moduleMessages 
    implicit none
    character(len=*), intent(in) :: cmd
    real(8), dimension(:), intent(out) :: vec
    integer :: i, nw, nt, n
    character(len=STRLEN), dimension(50) :: words, tokens
    integer :: ios

    call parse(cmd," ",words,nw)
    n=size(vec)
    call parse(words(2),",",tokens,nt)
    if (nt/=n) call message(-1,"Wrong number of arguments for +"//trim(cmd))
    do i=1,n
      if (tokens(i) == "pi") then
        vec(i) = acos(-1.d0)
      else if (tokens(i) == "-pi") then
        vec(i) = -acos(-1.d0)
      else if (tokens(i) == "twopi") then
        vec(i) = 2.d0*acos(-1.d0)
      else if (tokens(i) == "-twopi") then
        vec(i) = -2.d0*acos(-1.d0)
      else
        read(tokens(i),*,iostat=ios) vec(i)
      endif
      if (ios /= 0) call message(-1,"Wrong arguments for +"//trim(cmd))
    enddo

    return
  end subroutine

  subroutine read_int_vector(cmd,vec)
    use moduleMessages 
    implicit none
    character(len=*), intent(in) :: cmd
    integer, dimension(:), intent(out) :: vec
    integer :: i, nw, nt, n
    character(len=STRLEN), dimension(50) :: words, tokens
    integer :: ios

    call parse(cmd," ",words,nw)
    n=size(vec)
    call parse(words(2),",",tokens,nt)
    if (nt/=n) call message(-1,"Wrong number of arguments for +"//trim(cmd))
    do i=1,n
      read(tokens(i),*,iostat=ios)vec(i)
      if (ios /= 0) call message(-1,"Wrong arguments for +"//trim(cmd))
    enddo

    return
  end subroutine

  subroutine read_chr_vector(cmd,vec)
    use moduleMessages 
    implicit none
    character(len=*), intent(in) :: cmd
    character(len=*), dimension(:), intent(out) :: vec
    integer :: i, nw, nt, n
    character(len=STRLEN), dimension(50) :: words, tokens

    call parse(cmd," ",words,nw)
    n=size(vec)
    call parse(words(2),",",tokens,nt)
    if (nt/=n) call message(-1,"Wrong number of arguments for +"//trim(cmd))
    do i=1,n
      read(tokens(i),*)vec(i)
    enddo

    return
  end subroutine

!!!!!!!!

  subroutine findScalarInteger(string,flag,result,default)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    integer, intent(out) :: result
    integer, intent(in) :: default

    result = default
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call readTokens(listOfWords(iWord),result)
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findScalarInteger

  subroutine findVectorInteger(string,flag,result,default)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    integer, dimension(:), intent(out) :: result
    integer, dimension(:), intent(in) :: default

    result = default
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call readTokens(listOfWords(iWord),result)
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findVectorInteger

  subroutine findArrayInteger(string,flag,result)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    integer, allocatable, dimension(:), intent(out) :: result

    character(len=30), dimension(100) :: words, tokens
    integer :: nw, nt

    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call parse(listOfWords(iWord)," ",words,nw)
        call parse(words(2),",",tokens,nt)
        allocate(result(nt))
        call readTokens(listOfWords(iWord),result(1:nt))
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findArrayInteger

  subroutine findScalarDouble(string,flag,result,default)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    real(8), intent(out) :: result
    real(8), intent(in) :: default

    result = default
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call readTokens(listOfWords(iWord),result)
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findScalarDouble

  subroutine findVectorDouble(string,flag,result,default)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    real(8), dimension(:), intent(out) :: result
    real(8), dimension(:), intent(in) :: default

    result = default
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call readTokens(listOfWords(iWord),result)
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findVectorDouble

  subroutine findArrayDouble(string,flag,result)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    real(8), allocatable, dimension(:), intent(out) :: result

    character(len=30), dimension(100) :: words, tokens
    integer :: nw, nt

    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call parse(listOfWords(iWord)," ",words,nw)
        call parse(words(2),",",tokens,nt)
        allocate(result(nt))
        call readTokens(listOfWords(iWord),result(1:nt))
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findArrayDouble

  subroutine findScalarReal(string,flag,result,default)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    real(8), intent(out) :: result
    real(4), intent(in) :: default

    result = dble(default)
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call readTokens(listOfWords(iWord),result)
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findScalarReal

  subroutine findVectorReal(string,flag,result,default)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    real(8), dimension(:), intent(out) :: result
    real(4), dimension(:), intent(in) :: default

    result = dble(default)
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call readTokens(listOfWords(iWord),result)
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findVectorReal

  subroutine findScalarString(string,flag,result,default)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    character(len=*), intent(out):: result
    character(len=*), intent(in) :: default

    result = default
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call readTokens(listOfWords(iWord),result)
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findScalarString

  subroutine findVectorString(string,flag,result,default)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    character(len=*), dimension(:), intent(out) :: result
    character(len=*), intent(in) :: default

    result = default
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call readTokens(listOfWords(iWord),result)
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findVectorString

  subroutine findArrayString(string,flag,result)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    character(len=*), allocatable, dimension(:), intent(out) :: result

    character(len=200), dimension(100) :: words, tokens
    integer :: nw, nt, i

    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call parse(listOfWords(iWord)," ",words,nw)
        call parse(words(2),",",tokens,nt)
        allocate(result(nt))
        do i=1,nt
          result(i) = tokens(i)
        enddo
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

  end subroutine findArrayString
  
  subroutine findLogical(string,flag,result,default)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag
    logical, intent(out):: result
    logical, intent(in) :: default

    result = default
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        result = (.not. default)
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

    return
  end subroutine findLogical

  logical function flagExists(string,flag) 
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag

    flagExists = .false.
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        flagExists = .true.
        exit
      end if
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

  end function flagExists

  subroutine removeFlag(string,flag) 
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag

    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) exit
    enddo
    call concatenateWords(numberOfWords, listOfWords, string, iWord)

  end subroutine removeFlag

  integer function countFlagArguments(string,flag)
    implicit none
    character(len=*), intent(inout) :: string
    character(len=*), intent(in) :: flag

    character(len=200), dimension(100) :: words, tokens
    integer :: nw, nt

    nt = 0
    call parse(string,"+",listOfWords,numberOfWords)
    do iWord=1,numberOfWords
      if (index(listOfWords(iWord),flag(2:))==1) then
        call parse(listOfWords(iWord)," ",words,nw)
        call parse(words(2),",",tokens,nt)
        countFlagArguments = nt
        exit
      end if
    enddo
    
    return
  end function countFlagArguments

  subroutine concatenateWords(n,w,string,idx)
    implicit none
    integer, intent(in) :: n, idx
    character(len=*), dimension(n), intent(in) :: w
    character(len=*), intent(inout) :: string
    character(len=len_trim(string)) :: strsave
    integer :: i

    strsave = string
    call compact(strsave)

    string = ""
    do i=1,n
      if (i==idx) then
        string = trim(string)//" %"//trim(w(i))
      else
        string = trim(string)//" +"//trim(w(i))
      end if
    end do
    call compact(string)

    ! Remove leading + in case the command's first argument wasn't a flag i.e. --o a.pdb
    if (len_trim(string)-len_trim(strsave) == 1) then
      string(1:1) = " "
      call compact(string)
    end if

    do i=1,len_trim(string)
      if (strsave(i:i)=="%" .and. string(i:i)=="+") string(i:i) = "%"
    end do

  end subroutine concatenateWords

  function extractSelectionIndices(string) result(indices)
    character(len=*), intent(in) :: string
    integer, dimension(3) :: indices

    integer :: ntmp
    character(len=10), dimension(3) :: localNumbers

    indices(3) = 1
    call parse(string,":",localNumbers,ntmp)
    if (ntmp==1) then
      read(localNumbers(1),*) indices(1)
      indices(2) = indices(1)
    else if (ntmp==2) then
      read(localNumbers(1),*) indices(1)
      read(localNumbers(2),*) indices(2)
    else if (ntmp==3) then
      read(localNumbers(1),*) indices(1)
      read(localNumbers(2),*) indices(2)
      read(localNumbers(3),*) indices(3)
    end if

  end function extractSelectionIndices

    function extractFrameIndices(string) result(indices)
    character(len=*), intent(in) :: string
    integer, dimension(3) :: indices

    integer :: ntmp
    character(len=10), dimension(3) :: localNumbers

    indices(3) = 1
    call parse(string,":",localNumbers,ntmp)
    if (ntmp==1) then
      read(localNumbers(1),*) indices(1)
      indices(2) = huge(1)
    else if (ntmp==2) then
      read(localNumbers(1),*) indices(1)
      read(localNumbers(2),*) indices(2)
    else if (ntmp==3) then
      read(localNumbers(1),*) indices(1)
      read(localNumbers(2),*) indices(2)
      read(localNumbers(3),*) indices(3)
    end if

  end function extractFrameIndices

end module moduleStrings

  !**********************************************************************
  
  subroutine compact(str)
  
  ! Converts multiple spaces and tabs to single spaces; deletes control characters;
  ! removes initial spaces.
  
  character(len=*):: str
  character(len=1):: ch
  character(len=len_trim(str)):: outstr
  
  str=adjustl(str)
  lenstr=len_trim(str)
  outstr=' '
  isp=0
  k=0
  
  do i=1,lenstr
    ch=str(i:i)
    ich=iachar(ch)
  
    select case(ich)
  
      case(9,32)     ! space or tab character
        if(isp==0) then
          k=k+1
          outstr(k:k)=' '
        end if
        isp=1
  
      case(33:)      ! not a space, quote, or control character
        k=k+1
        outstr(k:k)=ch
        isp=0
  
    end select
  
  end do
  
  str=adjustl(outstr)
  
  end subroutine compact

  subroutine checkUsedFlags(cmd)
    use moduleMessages
    implicit none
    character(len=*), intent(in) :: cmd
    character(len=100) :: str
    integer :: i
    call compact(cmd)
    i = index(cmd,"+")
    if (i==0) return

    call message(2)
    do i=1,len_trim(cmd)
      if (cmd(i:i)=="+") then
        read(cmd(i:),*) str
        call message(0,"###### WARNING ###### unrecognised flag",str=str)
      end if
    end do
    call message(2)

  end subroutine checkUsedFlags
