subroutine create_new_cell(hmat_new)
  use moduleVariables, only: real64, cp, frameTypeDef, moleculeTypeDef, fileTypeDef
  use moduleSystem, only: frame, numberOfMolecules, listOfMolecules
  use moduleMessages
  use moduleDistances
  use moduleNeighbours

  implicit none
  real(real64), intent(in) :: hmat_new(3,3)
  real(real64) :: hmat_work(3,3)
  ! type(frameTypeDef), intent(inout) :: frame

  integer :: h, k, l
  integer :: iatm, imol, idx, jmol, nmax, ii
  integer :: kmax
  real(real64) ::  hinv_new(3,3), rtmp
  real(real64) ::  hmat_fin(3,3), hinv_fin(3,3)
  
  type(moleculeTypeDef), allocatable, dimension(:) :: molecule_new
  real(real64), allocatable, dimension(:,:) :: pos_new
  real(real64), allocatable, dimension(:) :: chg_new
  character(cp), allocatable, dimension(:) :: lab_new
  integer, allocatable, dimension(:) :: atm2mol_new
  real(real64), allocatable, dimension(:,:) :: zmol

  real(real64), dimension(3) :: dij, sij, rdist, xcom
  integer :: nmols_new, natoms_new

  hmat_work = hmat_new

  ! make sure that the three axes form a right-handed cell
  call getInverseCellMatrix(hmat_work,hinv_new,rtmp)
  if (rtmp<0.0) then
    hmat_work(1:3,2)= -hmat_work(1:3,2)
    call getInverseCellMatrix(hmat_work,hinv_new,rtmp)
  endif

  if ( abs((rtmp/frame % volume)-nint((rtmp/frame % volume))) > 1.0e-4_real64 ) &
    call message(-1,"Incommensurate volume for the new cell",r=rtmp/frame % volume)

  ! Look for molecules inside the new cell
  iatm=0
  nmols_new=0
  nmax = nint( rtmp/frame % volume ) + 5
  allocate(molecule_new(numberOfMolecules*nmax))
  allocate(zmol(3,numberOfMolecules*nmax))
  natoms_new=0
  allocate(pos_new(3,  frame % natoms*nmax))
  allocate(lab_new(    frame % natoms*nmax))
  allocate(chg_new(    frame % natoms*nmax))
  allocate(atm2mol_new(frame % natoms*nmax))
  pos_new = 0.d0
  lab_new = "XX"
  chg_new = 0.d0
  atm2mol_new = 0

  if (numberOfMolecules == 0)then
    call runInternalAction("topology","NULL")
  end if

  do imol=1,numberOfMolecules
    xcom=0.0_real64
    do idx=1,listOfMolecules(imol) % numberOfAtoms
      xcom = xcom + frame % pos(1:3,listOfMolecules(imol)%listOfAtoms(idx))
    enddo
    xcom=xcom/listOfMolecules(imol) % numberOfAtoms

    kmax=30
    do h=-kmax,kmax
      do k=-kmax,kmax
        do l=-kmax,kmax
          dij = xcom + h * frame % hmat(:,1) + k * frame % hmat(:,2) + l * frame % hmat(:,3)

          sij(1) = hinv_new(1,1)*dij(1) + hinv_new(1,2)*dij(2) + hinv_new(1,3)*dij(3)
          sij(2) = hinv_new(2,1)*dij(1) + hinv_new(2,2)*dij(2) + hinv_new(2,3)*dij(3)
          sij(3) = hinv_new(3,1)*dij(1) + hinv_new(3,2)*dij(2) + hinv_new(3,3)*dij(3)
          if (sij(1)>=0.0_real64 .and. sij(1)<=1.0_real64 .and. &
              sij(2)>=0.0_real64 .and. sij(2)<=1.0_real64 .and. &
              sij(3)>=0.0_real64 .and. sij(3)<=1.0_real64) then

            ! Remove molecule if it overlaps with something already in the cell
            do jmol=1,nmols_new
              rdist(1:3) = dij(1:3) - zmol(1:3,jmol)
              sij(1) = hinv_new(1,1)*rdist(1) + hinv_new(1,2)*rdist(2) + hinv_new(1,3)*rdist(3)
              sij(2) = hinv_new(2,1)*rdist(1) + hinv_new(2,2)*rdist(2) + hinv_new(2,3)*rdist(3)
              sij(3) = hinv_new(3,1)*rdist(1) + hinv_new(3,2)*rdist(2) + hinv_new(3,3)*rdist(3)
              sij(1:3) = sij(1:3)-nint(sij(1:3))
              rdist(1) = hmat_work(1,1)*sij(1) + hmat_work(1,2)*sij(2) + hmat_work(1,3)*sij(3)
              rdist(2) = hmat_work(2,1)*sij(1) + hmat_work(2,2)*sij(2) + hmat_work(2,3)*sij(3)
              rdist(3) = hmat_work(3,1)*sij(1) + hmat_work(3,2)*sij(2) + hmat_work(3,3)*sij(3)
              rtmp = sqrt(sum(rdist*rdist))
              if (rtmp < 0.01_real64) exit
            enddo
            if (jmol<=nmols_new) cycle

            ! New molecule found
            nmols_new = nmols_new+1
            natoms_new = natoms_new+listOfMolecules(imol) % numberOfAtoms

            ! Add the new molecule's info
            molecule_new(nmols_new) % ID            = listOfMolecules(imol) % ID
            molecule_new(nmols_new) % numberOfAtoms = listOfMolecules(imol) % numberOfAtoms
            allocate(molecule_new(nmols_new) % listOfAtoms(listOfMolecules(imol) % numberOfAtoms))
            zmol(1:3,nmols_new) = dij(1:3)

            ! Add the new atoms' info
            do idx=1,listOfMolecules(imol) % numberOfAtoms
              iatm = iatm + 1
              molecule_new(nmols_new) % listOfAtoms(idx) = iatm
              ii = listOfMolecules(imol) % listOfAtoms(idx)
              lab_new(iatm) = frame % lab(ii)
              pos_new(1:3,iatm) = frame % pos(1:3,ii) + h * frame % hmat(:,1) + k * frame % hmat(:,2) + l * frame % hmat(:,3)
              chg_new(iatm) = frame % chg(ii)
              atm2mol_new(iatm) = nmols_new
            enddo

          endif
        enddo
      enddo
    enddo
  enddo

! Rotate cell, the atoms' positions and molecules' com
  hmat_fin = hmat_work

  call makeUpperTriangularCell(hmat_work,pos_new,natoms_new)
  call makeUpperTriangularCell(hmat_fin,zmol,nmols_new)
!  call straightenCell(hmat_fin,"z",hmat_work)
!  hmat_fin=hmat_work
  call getInverseCellMatrix(hmat_fin,hinv_fin,rtmp)

  call deleteSystemArrays(frame)
  call createSystemArrays(frame, natoms_new ) ! also sets frame % natoms
  frame % hmat = hmat_work
  do iatm=1,natoms_new
    frame % pos(1:3,iatm) = pos_new(1:3,iatm)
  end do
  do iatm=1,natoms_new
    frame % lab(iatm) = lab_new(iatm)
  end do
  do iatm=1,natoms_new
    frame % chg(iatm) = chg_new(iatm)
  end do

  call getInverseCellMatrix(frame % hmat, frame % hinv, frame % volume)
  call hmat2cell (frame % hmat, frame % cell, "DEG")

  call cartesianToFractional(frame % natoms, frame % pos, frame % frac)

  call setUpNeigboursList()
  call updateNeighboursList(.true.)

  if (numberOfMolecules > 0) call runInternalAction("topology","+update +rebuild")
  
  do imol=1,numberOfMolecules
    xcom=0.0_real64
    do idx=1,listOfMolecules(imol) % numberOfAtoms
      xcom = xcom + frame % pos(1:3,listOfMolecules(imol)%listOfAtoms(idx))
    end do
    listOfMolecules(imol) % centreOfMass = xcom / listOfMolecules(imol) % numberOfAtoms
  end do

end subroutine create_new_cell
