!disclaimer
module moduleXrayAction
  use moduleVariables
  use moduleMessages
  use moduleSystem
  use scattering_factors, only : SF_Library_dim, SF_Library
  implicit none
  
  public :: computeXrayPowder, computeXrayPowderHelp
  private
  
  ! real(real64), parameter  :: acclat   = 0.000001_real64
  integer, parameter  :: npeaks_max = 100000
  
  real(real64) :: d_fwhm      !(in two theta)
  real(real64) :: step        !(in two theta)
  real(real64) :: n_cauchy 
  real(real64) :: lambda
  real(real64) :: min2th
  real(real64) :: max2th
  integer :: kmax

  real(real64), save, allocatable, dimension (:)  :: pattern
  real(real64), save, allocatable, dimension (:)  :: twotheta
  integer, save, allocatable, dimension (:)  :: libndx

  integer :: ipts, npts
  integer :: npeaks

  real(real64), dimension(npeaks_max) :: peak_fr, inte_fr, dhkl_fr
  integer, dimension(3,npeaks_max) :: hkl_fr
  integer, dimension(npeaks_max) :: order

  character(:), pointer :: actionCommand
  type(fileTypeDef), pointer :: outputFile
  logical, pointer :: firstAction

contains

subroutine computeXrayPowderHelp()
    implicit none
    call message(0,"This action computes the X-ray powder diffraction spectrum for the input structure.")
    call message(0,"Examples:")
    call message(0,"  gpta.x --i coord.pdb --xray +lambda 1.54056")
  end subroutine computeXrayPowderHelp

  subroutine computeXrayPowder(a)
    ! use moduleVariables
    use moduleSystem
    use moduleElements
    implicit none
    type(actionTypeDef), target :: a

    integer :: iatm, ndx, id

    if (a % actionInitialisation) then
      call initialiseAction(a)
      return
    end if

    if (frameReadSuccessfully) then

      if (firstAction) then
        
        call dumpScreenInfo()
        
        ndx = 0
        allocate(libndx(frame % natoms))
        
        do iatm=1,frame % natoms
          ndx = ndx + 1
          id = getAtomicNumber( frame % lab(iatm) )
          if (id==0) call message(0,"Cannot assign scattering factor to ",str=frame % lab(iatm))
          libndx(ndx) = id
        enddo
        
        call checkUsedFlags(actionCommand)
        firstAction = .false.
      end if

      call sfact_calc(frame % natoms, frame % hinv, frame % pos)

    end if

    if (endOfCoordinatesFiles) then
      call finaliseAction()
    end if 
  end subroutine computeXrayPowder

  subroutine initialiseAction(a)
    use moduleSystem
    use moduleStrings
    implicit none
    type(actionTypeDef), target :: a
    integer :: i

    a % actionInitialisation = .false.

    ! Local pointers
    actionCommand        => a % actionDetails
    firstAction          => a % firstAction
    outputFile           => a % outputFile

    call assignFlagValue(actionCommand,"+out",outputFile % fname,'xray.out')

    call assignFlagValue(actionCommand,"+min",min2th,5.0_real64)
    call assignFlagValue(actionCommand,"+max",max2th,50.0_real64)
    call assignFlagValue(actionCommand,"+lambda",lambda,1.540560_real64)
    call assignFlagValue(actionCommand,"+kmax",kmax,30)
    
    call assignFlagValue(actionCommand,"+fwhm",d_fwhm,0.1_real64)
    call assignFlagValue(actionCommand,"+d2t",step,0.01_real64)
    call assignFlagValue(actionCommand,"+cauchy",n_cauchy,0.5_real64)

    npts= int((max2th - min2th)/step) + 1
    allocate(twotheta(npts))
    do i=1,npts
      twotheta(i)=min2th+(i-1)*step
    enddo
    allocate(pattern(npts), source=0.0_real64)

  end subroutine initialiseAction

  subroutine dumpScreenInfo()
    use moduleMessages 
    implicit none
    
    call message(0,"Computing Xray Powder Diffraction Pattern")
    call message(0,"...Output file",str=outputFile % fname)
    call message(0,"...Wavelength [nm]",r=lambda)
    call message(0,"...Twotheta range",rv=[min2th,max2th])
    call message(0,"...Twotheta resolution",r=step)
    call message(0,"...Number of kspace vectors",i=kmax)
    call message(0,"...Smoothing parameters")
    call message(0,"......FWHM",r=d_fwhm)
    call message(0,"......Cauchy n parameter",r=n_cauchy)

  end subroutine dumpScreenInfo

  subroutine finaliseAction()
    ! use m_rnkpar, only : rnkpar
    use m_mrgrnk, only : mrgrnk
    use moduleFiles
    implicit none
    integer :: i, ipeaks, ndx
    real(real64) :: rtmp

    integer :: funit

    call initialiseFile(outputFile,outputFile % fname)
    funit = outputFile % funit

    ! call rnkpar(dhkl_fr(1:npeaks) , order(1:npeaks), npeaks)
    call mrgrnk(dhkl_fr(1:npeaks) , order(1:npeaks))
    rtmp=1000.0_real64
    ndx=0
    do ipeaks=npeaks,1,-1
      ipts=order(ipeaks)
      if (rtmp-dhkl_fr(ipts) > 0.000001_real64) then
        rtmp=dhkl_fr(ipts)
        ndx=ndx+1
        if (ndx==1) then
          write(funit,'("# index    d_hkl    twoth  intensity      (h , k , l) ...")')
        else
          write(funit,*)
        endif
        write(funit,'("#",i6,1x,2(f8.3,1x),f10.0,5x,3(i3,1x))',advance='no') &
          ndx,dhkl_fr(ipts),peak_fr(ipts),inte_fr(ipts),hkl_fr(1:3,ipts)
      else
        write(funit,'(5x,3(i3,1x))',advance='no')hkl_fr(1:3,ipts)
      endif
    enddo
    write(funit,*)

    rtmp=maxval(pattern)
    rtmp=1.0_real64/rtmp
    write(funit,'("# Two Theta | Intensity [a.u.]")')
    do i=1,npts
      write(funit,'(f7.3,1x,f12.5)')twotheta(i),rtmp*pattern(i)
    enddo

  end subroutine finaliseAction

  subroutine sfact_calc(nn,hinv,pos)
    implicit none
    integer :: nn
    real(real64), dimension (3,3)  :: hinv
    real(real64), dimension (3,nn) :: pos
    integer :: i, ndx
    integer :: ih, ik, il
    real(real64), dimension (3)  :: sij, g
    real(real64) :: rg, sa, sinth, sin2th, cos2th, twoth
    real(real64) :: ricos, risin
    integer  :: ipeaks, isf_lib
    real(real64) :: ri, rip, s_el, b, p, s, x, rtmp

    npeaks = 0
    peak_fr = 0.0_real64
    inte_fr = 0.0_real64
! Loop over hkl
    do ih=0,kmax
      do ik=-kmax,kmax
        do il=-kmax,kmax
          g(1)=real(ih,8)
          g(2)=real(ik,8)
          g(3)=real(il,8)
          sij=matmul(g,hinv)
          rg=sqrt(dot_product(sij,sij))

          if(all(g==0.0_real64))cycle

! 1/d where d is the distance between the hkl planes (length of the hkl vector)
          sa = 0.5_real64*rg
! sinth = lambda/2d
          sinth=0.5_real64*rg*lambda
          if(abs(sinth)>1.)cycle

! compute twotheta in degrees
          twoth=2.0_real64*asin(sinth)
          cos2th=cos(twoth)
          sin2th=sin(twoth)
          twoth=twoth*360.0_real64/twopi
          if (twoth<=min2th .or. twoth>=max2th) cycle

          ricos=0.0_real64
          risin=0.0_real64
          ndx=0
! Loop over the atoms to calculate the structure factor
          do i=1,nn
            ndx=ndx+1
            isf_lib=libndx(ndx)
            sij=matmul(hinv,pos(:,i))
            sij=sij-int(sij)
            sij=sij-int(2*sij)
            s=twopi*dot_product(g,sij)

! Atomic scattering factors
            s_el = SF_Library(isf_lib)%a_anis(1)*exp(-SF_Library(isf_lib)%b_anis(1)*sa*sa)+ &
                   SF_Library(isf_lib)%a_anis(2)*exp(-SF_Library(isf_lib)%b_anis(2)*sa*sa)+ &
                   SF_Library(isf_lib)%a_anis(3)*exp(-SF_Library(isf_lib)%b_anis(3)*sa*sa)+ &
                   SF_Library(isf_lib)%a_anis(4)*exp(-SF_Library(isf_lib)%b_anis(4)*sa*sa)+ &
                   SF_Library(isf_lib)%c_anis

! Thermal vibration correction
! not sure about the 0.05 and 0.06
             if(SF_Library(isf_lib)%AtomicNumber.eq.1)then
                B=2*TWOPI*TWOPI*0.05_real64
             else
                B=2*TWOPI*TWOPI*0.06_real64
             endif
             s_el = s_el*exp(-B*sa*sa)
  !!!p Real and Imaginary part of the intensity
             ricos=ricos+cos(s)*s_el
             risin=risin+sin(s)*s_el
          enddo
! Scattered intensity
          ri=(ricos*ricos+risin*risin)
! Lorentz polarisation factor - the 1/2 might be redundant
          p=(1.0_real64+cos2th*cos2th)/(sinth*sin2th)
          ! p=0.5_real64*(1.0_real64+cos2th*cos2th)/(sinth*sin2th)
          rip=ri*p
          npeaks=npeaks + 1
          if (npeaks.gt.npeaks_max) call message (0,"Increase NPEAKS_MAX. Exiting")
          peak_fr(npeaks) = twoth
          hkl_fr(1:3,npeaks) = (/ih,ik,il/)
          dhkl_fr(npeaks) = 1.0_real64/rg
          if (ih.eq.0) then
             inte_fr(npeaks) = rip
          elseif(ih.ne.0) then
             inte_fr(npeaks) = 2*rip
          endif
        enddo
      enddo
    enddo

! Smooth it out
    do ipeaks = 1, npeaks
      do ipts = 1, npts
        if ( abs(peak_fr(ipeaks)-twotheta(ipts)).gt.30*d_fwhm) cycle

        x = (twotheta(ipts)-peak_fr(ipeaks))/(0.5_real64*d_fwhm)
        rtmp = 0.0_real64
! prevents underflow in the exponent
        if (-log(2.0_real64)*x**2>-100) rtmp = (1-n_cauchy)*exp(-log(2.0_real64)*x**2)

! compute the contribution of this peak on the specific twoth
        rtmp = n_cauchy/(1+x**2) + rtmp

        pattern(ipts) = pattern(ipts) + inte_fr(ipeaks)*rtmp
      enddo
    enddo

    return
  end subroutine sfact_calc

end module moduleXrayAction
