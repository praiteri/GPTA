!disclaimer
subroutine defineVariables(cmd)
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleNeighbours 
  use moduleElements
  use moduleMessages
  use moduleNeighbours

  implicit none

  character(len=STRLEN), intent(in) :: cmd

  integer :: iword, i
  integer :: numberOfWords, nw
  character(len=STRLEN), dimension(100) :: listOfWords, words

  character(len=STRLEN) :: wstr

  call parse(cmd, " ", listOfWords, numberOfWords)

  iword = 0
  do while (iword < numberOfWords)
    iword = iword + 1
    i = len_trim(listOfWords(iword)) + 1
    call lowercase (listOfWords(iword),wstr)
    wstr(i:) = ""
    select case (wstr)
      case default
        call message(-1,"--define - unknown variable",str=wstr)

      case ("seed")
        iword = iword + 1
        read(listOfWords(iword),*) randomNumberSeed

      case ("rscale")
        iword = iword + 1
        read(listOfWords(iword),*) distanceScaling
        
      case ("safedist" , "safepbc")
        safeDistances = .true.

      case ("verlet")
        forceVerletList = .true.

      case ("rcov")

        iword = iword + 1
        call parse(listOfWords(iword),"=",words,nw)
        if (nw /= 2) then
          call message(-1,"--define - cannot read covalent radius",str=listOfWords(iword))
        else
          do i=1,nelement
            if (words(1) == atom(i) % lab) then
              read(words(2),*) atom(i) % rcov
              exit
            end if
          end do
        end if

      case ("bond")
        iword = iword + 1
        call parse(listOfWords(iword),"=",words,nw)
        if (nw /= 2) call message(-1,"Cannot parse special bond",str=listOfWords(iword))
        sbond_n = sbond_n + 1
        read(words(2),*) sbond_d(sbond_n)
        
        wstr = words(1)
        call parse(wstr,"-",words,nw)
        if (nw /= 2) then
          call parse(wstr,",",words,nw)
          if (nw /= 2) call message(-1,"Cannot parse special bond",str=wstr)
        end if
        read(words(1),*) sbond_l(1,sbond_n)
        read(words(2),*) sbond_l(2,sbond_n)

      case ("reax")
        resetFrameLabels = .false.
        resetFrameCharges = .false.
        resetFrameElements= .false.
        
    end select
  end do

end subroutine defineVariables

subroutine defineVariablesHelp()
  use moduleMessages
  implicit none
  call message(0,"This action is executed once before any other command and alter the default values of some global variables.")
  call message(0,"The allowed arguments for this action are:")

  call message(0,"  * seed     -> set the random number generator seed")
  call message(0,"  * rscale   -> set the scaling parameter for the maximum bond lengths")
  call message(0,"  * safepbc  -> use a safe (and slow) algorithm for the calculation of distances and PBC")
  call message(0,"  * safedist -> use a safe (and slow) algorithm for the calculation of distances and PBC")
  call message(0,"  * verlet   -> use the Verlet algorithm for the neighbours' list")
  call message(0,"  * rcov     -> sets the covalent radius of an atom")
  call message(0,"  * bond     -> set the maximum bond length for a pair")
  call message(0,"  * reax     -> assumes a reactive trajectory with different connectivity in each frame (Experimental)")
  call message(0,"Examples:")
  call message(0,"  gpta.x --define seed 12345")
  call message(0,"  gpta.x --define rscale 1.15")
  call message(0,"  gpta.x --define safepbc")
  call message(0,"  gpta.x --define safedist")
  call message(0,"  gpta.x --define verlet")
  call message(0,"  gpta.x --define rcov Ca=0.0")
  call message(0,"  gpta.x --define bond Ca-OW=0.1")
  call message(0,"  gpta.x --define reax")
end subroutine defineVariablesHelp 

subroutine dumpVariablesInfo()
  use moduleVariables
  use moduleSystem
  use moduleStrings
  use moduleNeighbours 
  use moduleMessages
  
  integer :: i

  call message(0,"Global variables settings")
  call message(0,"...Random number seed",i=randomNumberSeed)
  call message(0,"...Scaling factor for covalent bonds",r=distanceScaling)
  if (safeDistances) call message(0,"...Compute distances with loop over neighbour cells")
  if (sbond_n>0) then
    do i=1,sbond_n
      call message(0,"...Special bond length for "//trim(sbond_l(1,i))//"-"//trim(sbond_l(2,i)),r=sbond_d(i))
    end do
  end if
  call message(2)

end subroutine dumpVariablesInfo
