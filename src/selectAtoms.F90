!disclaimer
subroutine selectAtoms(nsel,cmd,a)
  use moduleVariables
  use moduleSystem 
  use moduleStrings
  use moduleMessages 
  implicit none
  integer, intent(in) :: nsel
  character(len=*), intent(inout) :: cmd
  type(actionTypeDef), intent(inout) :: a

  integer :: numberOfFlags, numberOfWords, numberOfStrings
  character(len=STRLEN), dimension(100) :: listOfWords
  character(len=STRLEN), dimension(100) :: listOfFlags
  
  character(len=STRLEN) :: selectionString
  character(len=STRLEN) :: flagString
  character(len=STRLEN) :: posString
  integer :: selectionID
  
  integer :: ntmp
  integer :: idx, jdx, indices(3)
  integer :: i, j, k, itmp, imol, iatm
  character(len=1) :: chr
  
  character(cp), dimension(100) :: localLabels
  character(len=STRLEN), dimension(100) :: localString
  character(len=STRLEN), dimension(100) :: tmpString
  logical, dimension(3) :: lflag

  real(real64), dimension(2,3) :: posLimits

  selectionString = ''

  if (.not. allocated(a % isSelected)) then
    allocate(a % isSelected(frame % natoms , nsel), source=.false.)
  else
    a % isSelected = .false.
  end if 

  selectionID = 0
  call parse2(cmd,"+",listOfFlags,numberOfFlags)
  do idx=1,numberOfFlags
  end do

  do idx=1,numberOfFlags

    ! select atoms by label
    if ( index(listOfFlags(idx),"+s ") > 0) then
      selectionString = trim(selectionString) // " " // trim(listOfFlags(idx))
      
      call extractFlag(listOfFlags(idx),"+s",flagString)
      call parse(flagString," ",listOfWords,numberOfWords)
      
      do jdx=1,numberOfWords
        selectionID = selectionID + 1
        if (selectionID > nsel) call message(-1,"Too many selection flags or blocks",str=selectionString)

        call parse(listOfWords(jdx),",",localLabels,ntmp)   
        if (ntmp == 0) call message(-1,"Cannot parse selection by label",str=listOfWords(jdx))

        if (ntmp == 1 .and. (localLabels(1) == "all" .or. localLabels(1) == "ALL")) then
          a % isSelected(: , selectionID) = .true.
        else
          do i=1,frame % natoms
            if (any(localLabels(1:ntmp)==frame % lab(i))) a % isSelected(i , selectionID) = .true.
          end do
        end if
      end do

    ! select atoms by index
    else if ( index(listOfFlags(idx),"+i ") > 0) then
      selectionString = trim(selectionString) // " " // trim(listOfFlags(idx))

      call extractFlag(listOfFlags(idx),"+i ",flagString)
      call parse(flagString," ",listOfWords,numberOfWords)

      do jdx=1,numberOfWords
        selectionID = selectionID + 1
        if (selectionID > nsel) call message(-1,"Too many selection flags or blocks",str=selectionString)

        call parse(trim(listOfWords(jdx)),",",localString,ntmp) 
        if (ntmp == 0) call message(-1,"Cannot parse selection by index",str=listOfWords(jdx))

        if (ntmp == 1 .and. (localString(1) == "all" .or. localString(1) == "ALL")) then
          a % isSelected(: , selectionID) = .true.
        else
          do j=1,ntmp
            indices = extractSelectionIndices(localString(j))
            do i=indices(1),indices(2),indices(3)
              a % isSelected(i , selectionID) = .true.
            end do
          end do
        
        end if

      end do

    ! select atoms by molecule
    else if ( index(listOfFlags(idx),"+mol ") > 0) then
      selectionString = trim(selectionString) // " " // trim(listOfFlags(idx))
  
      if (numberOfMolecules == 0) call message(-1,"Cannot perform selection by molecules, use --top")

      call extractFlag(listOfFlags(idx),"+mol ",flagString)
      call parse(flagString," ",listOfWords,numberOfWords)

      do jdx=1,numberOfWords
        selectionID = selectionID + 1
        if (selectionID > nsel) call message(-1,"Too many selection flags or blocks",str=selectionString)

        call parse(listOfWords(jdx),",",localLabels,ntmp)    
        if (ntmp == 0) call message(-1,"Cannot parse selection by label",str=listOfWords(jdx))

        if (ntmp == 1 .and. (localLabels(1) == "all" .or. localLabels(1) == "ALL")) then
          a % isSelected(: , selectionID) = .true.
        else
          do itmp=1,ntmp
            do imol=1,numberOfMolecules
              if (listOfMolecules(imol) % resname == localLabels(itmp)) then
                do iatm=1,listOfMolecules(imol) % numberOfAtoms
                  i = listOfMolecules(imol) % listOfAtoms(iatm)
                  a % isSelected(i , selectionID) = .true.
                end do
              end if
            end do
          end do
        end if
      end do

    else if ( index(listOfFlags(idx),"+pos ") > 0) then
      selectionID = selectionID + 1

      call extractFlag(listOfFlags(idx),"+pos ",flagString)
      call parse(flagString," ",listOfWords,numberOfWords)
      
      do i=1,numberOfWords
        call parse(listOfWords(i),",",localString,numberOfStrings)
        if (numberOfStrings /= 3) call message(-1,"Wrong number of blocks for +pos selection :",numberOfStrings)
        do j=1,numberOfStrings
          posString = trim(localString(j))
          call removeCharacter(posString,"[")
          call removeCharacter(posString,"]")
          call compact(posString)

          call parse(posString,":",tmpString,ntmp)
          if (trim(tmpString(1)) == "inf") then
            posLimits(1,j) = -huge(1.0_real64)
          else
            read(tmpString(1),*) posLimits(1,j) 
          endif

          if (trim(tmpString(2)) == "inf") then
            posLimits(2,j) = huge(1.0_real64)
          else
            read(tmpString(2),*) posLimits(2,j) 
          endif
        end do
        
        do k=1,frame % natoms
          lflag = .false.
          do j=1,3
            if ( frame%pos(j,k) > posLimits(1,j) .and. frame%pos(j,k) < posLimits(2,j) )  lflag(j) = .true.
          end do
          if (count(lflag) == 3) a % isSelected(k , selectionID) = .true.
        end do
      end do

      if (numberOfMolecules > 0) then
        do k=1,frame % natoms
          if ( .not. a % isSelected(k , selectionID) ) cycle
          imol = atomToMoleculeIndex(k)
          do j=1,listOfMolecules(imol) % numberOfAtoms
            i = listOfMolecules(imol) % listOfAtoms(j)
            a % isSelected(i , selectionID) = .true.
          end do
        end do
      end if
    else
      cycle
    end if
  end do

  if (selectionID < nsel) call message(-1,"Not enough selection flags or blocks",str=selectionString)

  ! mark the selection flags as used - not very elegant
  do 
    idx = maxval( [index(cmd,"+s") , index(cmd,"+i") , index(cmd,"+mol")] )
    if (idx == 0) exit
    cmd(idx:idx) = '%'
  end do

  if (a % firstAction) then
    call message(0,"...Atoms selection for "//trim(a % name))
    call message(0,"......Selection command",str=selectionString)
    do i=1,nsel
      write(chr,'(i1)') i
      call message(0,"......Atoms selected in group "//chr,i=count(a % isSelected(:,i)))
    end do
    call message(2)
  end if
  
end subroutine selectAtoms

subroutine getNumberOfAtomSelections(cmd,selectionID)
  use moduleVariables
  use moduleStrings

  implicit none
  character(len=*), intent(inout) :: cmd

  integer :: idx
  integer :: selectionID
  integer :: numberOfFlags
  character(len=STRLEN), dimension(100) :: listOfFlags
  
  selectionID = 0
  call parse2(cmd,"+",listOfFlags,numberOfFlags)
  do idx=1,numberOfFlags

    ! select atoms by label
    if ( index(listOfFlags(idx),"+s ") > 0) selectionID = selectionID + 1
    if ( index(listOfFlags(idx),"+i ") > 0) selectionID = selectionID + 1
    if ( index(listOfFlags(idx),"+mol ") > 0) selectionID = selectionID + 1
    if ( index(listOfFlags(idx),"+pos ") > 0) selectionID = selectionID + 1

  end do
end subroutine getNumberOfAtomSelections