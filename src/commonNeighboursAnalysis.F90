!disclaimer
module molecule_analysis
    implicit none

    ! Define a custom type to represent a bond
    type :: bond_type
        integer :: atom1, atom2
    end type bond_type

    ! Define a type to store the graph representation
    type :: graph_type
        integer :: num_atoms
        logical, allocatable :: adjacency(:,:)
    end type graph_type

contains
    subroutine find_longest_chain(bonds, num_bonds, longest_chain, chain_length)
        type(bond_type), intent(in) :: bonds(:)
        integer, intent(in) :: num_bonds
        integer, allocatable, intent(out) :: longest_chain(:)
        integer, intent(out) :: chain_length

        type(graph_type) :: graph
        integer :: i, start_atom
        integer, allocatable :: current_chain(:), best_chain(:)
        logical, allocatable :: visited(:)
        integer :: max_chain_length, num_atoms

        ! Find the number of atoms
        num_atoms = get_num_atoms(bonds, num_bonds)

        ! Create graph representation
        call create_graph(bonds, num_bonds, num_atoms, graph)

        ! Initialize variables for DFS
        allocate(visited(num_atoms))
        allocate(current_chain(num_atoms))
        allocate(best_chain(num_atoms))
        max_chain_length = 0

        ! Try starting from each atom
        do start_atom = 1, num_atoms
            visited = .false.
            current_chain = 0

            ! Perform DFS starting from this atom
            call dfs(start_atom, 1, graph, visited, current_chain, best_chain, max_chain_length)
        end do

        ! Prepare output
        chain_length = max_chain_length
        allocate(longest_chain(chain_length))
        longest_chain = best_chain(1:chain_length)

        ! Cleanup
        deallocate(visited, current_chain, best_chain)
        deallocate(graph%adjacency)

    contains
        recursive subroutine dfs(current, depth, graph, visited, current_chain, best_chain, max_length)
            integer, intent(in) :: current, depth
            type(graph_type), intent(in) :: graph
            logical, intent(inout) :: visited(:)
            integer, intent(inout) :: current_chain(:), best_chain(:), max_length
            integer :: neighbor

            visited(current) = .true.
            current_chain(depth) = current

            ! Update best chain if current is longer
            if (depth > max_length) then
                max_length = depth
                best_chain(1:depth) = current_chain(1:depth)
            end if

            ! Explore neighbors
            do neighbor = 1, graph%num_atoms
                if (graph%adjacency(current, neighbor) .and. (.not. visited(neighbor))) then
                    call dfs(neighbor, depth + 1, graph, visited, current_chain, best_chain, max_length)
                end if
            end do

            visited(current) = .false.
        end subroutine dfs
    end subroutine find_longest_chain

    function get_num_atoms(bonds, num_bonds) result(num_atoms)
        type(bond_type), intent(in) :: bonds(:)
        integer, intent(in) :: num_bonds
        integer :: num_atoms
        integer :: i

        num_atoms = 0
        do i = 1, num_bonds
            num_atoms = max(num_atoms, max(bonds(i)%atom1, bonds(i)%atom2))
        end do
    end function get_num_atoms

    subroutine create_graph(bonds, num_bonds, num_atoms, graph)
        type(bond_type), intent(in) :: bonds(:)
        integer, intent(in) :: num_bonds, num_atoms
        type(graph_type), intent(out) :: graph
        integer :: i

        graph%num_atoms = num_atoms
        allocate(graph%adjacency(num_atoms, num_atoms))
        graph%adjacency = .false.

        do i = 1, num_bonds
            graph%adjacency(bonds(i)%atom1, bonds(i)%atom2) = .true.
            graph%adjacency(bonds(i)%atom2, bonds(i)%atom1) = .true.
        end do
    end subroutine create_graph
end module molecule_analysis

module moduleCommonNeighboursAnalysis
    use moduleVariables
    use moduleSystem
    use moduleStrings
    use moduleFiles
    use moduleProperties
    use moduleMessages
    use moduleDistances
    
    implicit none

    public :: cnaAction, cnaActionHelp
    private

    character(:), pointer :: actionCommand
    logical, pointer :: firstAction
    type(fileTypeDef), pointer :: outputFile
    integer, pointer :: tallyExecutions
    
    real(8), pointer :: cutoffRadius
    integer, pointer, dimension(:) :: ID
    
contains

    subroutine cnaActionHelp()
        use moduleMessages
        implicit none
        call message(0,"This action computes the Common Neighbours Analisys.")
        call message(0,"Examples:")
        call message(0,"  gpta.x --i coord.pdb traj.dcd --cna +s Co +rcut 3.2 +out cna.out")
    end subroutine cnaActionHelp
    
    subroutine cnaAction(a)
        implicit none
        type(actionTypeDef), target :: a
        integer :: ngroups, maxNeigh, nsel1

        call associatePointers(a)
        if (a % actionInitialisation) then
            call initialiseAction(a)
            return
        end if
    
        if (frameReadSuccessfully) then
            tallyExecutions = tallyExecutions + 1
    
            ! Atoms' selection for ngroups = x
            ngroups = 1
            if (firstAction) then
                ! dump info about the action on the screen
                call dumpScreenInfo()
    
                ! select ngroups group(s) of atoms
                call selectAtoms(ngroups,actionCommand,a)
    
                ! update selection if "reactive" trajectory
                ! i.e. if the atom may change name 
                if (resetFrameLabels) then
                    a % updateAtomsSelection = .false.
                else
                    a % updateAtomsSelection = .true.
                end if
    
                ! create a list of the atoms' indices for each group
                call createSelectionList(a,ngroups)

                ! Temporary working arrays
                allocate(a % Iarray1D(frame % natoms), source=0)
                allocate(a % Iarray2D(50, frame % natoms), source=0)

                ! Throw a warning for unused flags
                call checkUsedFlags(actionCommand)
                firstAction = .false.
    
            else
                
                ! Repeat selection for reactive trajectories
                if (a % updateAtomsSelection) then
                    ! select two groups of atoms
                    call selectAtoms(ngroups,actionCommand,a)
                    ! create a list of the atoms' indices for each group
                    call createSelectionList(a,ngroups)
                end if 
                    
            end if
    
            call computeAction(a)
        end if
    
        if (endOfCoordinatesFiles) then
            call finaliseAction(a)
        end if 
    
    end subroutine cnaAction
    
    subroutine associatePointers(a)
        implicit none
        type(actionTypeDef), target :: a

        ! Local pointers
        actionCommand   => a % actionDetails
        firstAction     => a % firstAction
        tallyExecutions => a % tallyExecutions
        outputFile      => a % outputFile

        cutoffRadius    => a % doubleVariables(1)
        ID(1:7)         => a % integerVariables(1:7)
        
    end subroutine associatePointers

    subroutine initialiseAction(a)

        implicit none
        type(actionTypeDef), target :: a
        integer :: i
        character(3) :: strings(8)

        a % actionInitialisation = .false.
        a % requiresNeighboursList = .true.
        a % requiresNeighboursListUpdates = .true.
        a % requiresNeighboursListDouble = .true.
        a % cutoffNeighboursList = 3.2d0

        ! get output file name from the command line, if present
        call assignFlagValue(actionCommand,"+out",outputFile % fname,'cna.out')
        
        ! get cutoff radius from the command line, if present
        call assignFlagValue(actionCommand,"+rcut",cutoffRadius,a % cutoffNeighboursList)
        
        ! get number of bins for the distribution from the command line, if present
        ! call assignFlagValue(actionCommand,"+nbin",numberOfBins,100)
        
        ! UNK, FCC, HCP, BCC, ICO, NORMAL
        ! do i=1,6
        !     call workData % initialise(ID(i), "store")
        ! end do
        ! call workData % initialise(ID_unk, "dump")
        ! call workData % initialise(ID_theta, "average")

        a % cutoffNeighboursList = cutoffRadius
        tallyExecutions = 0

        outputFile % fname = "cna.out"
        call initialiseFile(outputFile, outputFile % fname)
        strings = ["UNK","FCC","HCP","BCC","ICO","NX ","NY ","NZ "]
        write(a % outputFile % funit,"('#',8(a8,1x))")strings

    end subroutine initialiseAction

    subroutine finaliseAction(a)
        implicit none
        type(actionTypeDef), target :: a

        ! real(8), allocatable, dimension(:) :: localDistribution
        ! integer :: i, ntmp, nvalues, nsel1, ires, nf
        ! integer :: nhcp, nfcc
        ! real(real64) :: vector(3)
        ! real(real64), allocatable, dimension(:) :: values
        ! real(real64), allocatable, dimension(:) :: bins, hist

        ! character(3), allocatable, dimension(:) :: strings
        ! real(real64), allocatable, dimension(:) :: results

        ! nf = 8 
        ! allocate(results(nf))
        ! allocate(strings(nf))
        ! strings = ["UNK","FCC","HCP","BCC","ICO","NX ","NY ","NZ "]

        ! ires = 1
        ! nsel1 = count(a % isSelected(:,1))
        ! do i=1,5
        !     call workData % extract(ID(i), localDistribution, numberOfCounts=ntmp)
        !     if (i == 2) nfcc = localDistribution(1)
        !     if (i == 3) nhcp = localDistribution(1)
        !     results(i) = localDistribution(1)
        !     ires = ires + 1
        !     deallocate(localDistribution)
        ! end do
        ! call workData % extract(ID(i), localDistribution, numberOfCounts=ntmp)

        ! allocate(values(ntmp/3))

        ! nvalues = 0
        ! vector = 0.d0
        ! do i=1,ntmp/3
        !     vector = vector + localDistribution(3*i-2:3*i)
        !     if (norm2(localDistribution(3*i-2:3*i)) > 1e-6) then
        !         nvalues = nvalues + 1
        !         values(nvalues) = localDistribution(3*i)
        !     end if
        ! end do
        ! deallocate(localDistribution)
        ! ! write(outputFile % funit,"('Normal vector:',3(f8.3,1x))") vector/nvalues
        ! results(ires:ires+2) = vector/nvalues
        ! ires = ires + 3
        
        ! write(a % outputFile % funit,"(8(a8,1x))")strings(1:ires-1)
        ! write(a % outputFile % funit,"(5(i8,1x),3(f8.3,1x))")int(results(1:5)),results(6:ires-1)

        ! ! write(0,*) 'Number of values:',nvalues
        ! ! allocate(bins(50),hist(50))
        ! ! call compute_histogram(values(1:nvalues), nvalues, 0.d0, 1.0d0, 50, bins, hist)
        ! ! do i=1,50
        ! !     write(124,*) bins(i),hist(i)
        ! ! end do
        close(a % outputFile % funit)

    end subroutine finaliseAction

    subroutine dumpScreenInfo()
        implicit none
        call message(0,"Test interface")
        call message(0,"...Output file",str=outputFile % fname)
        call message(0,"...Cutoff radius",r=cutoffRadius)
        ! call message(0,"...Number of bins",i=numberOfBins)

    end subroutine dumpScreenInfo

    subroutine computeAction(a)
        implicit none
        type(actionTypeDef), target :: a

        integer :: idx, jdx, ineigh
        integer :: iatm, jatm
        integer :: nsel1
        real(8) :: dij(3), dist2, rcut2
        integer :: ctype(5)

        logical :: lbcc, lfcc, lhcp, lico

        character(3) :: phase
        real(real64) :: vector(3), rtmp(1)
        character(3), allocatable, dimension(:) :: environments
        real(real64), allocatable, dimension(:,:) :: normal_vector

        character(3), allocatable, dimension(:) :: strings
        integer, dimension(5) :: results


        rcut2 = cutoffRadius**2

        nsel1 = count(a % isSelected(:,1))
        allocate(environments(nsel1))
        allocate(normal_vector(3,nsel1))

        a % Iarray1D = 0
        a % Iarray2D = 0

        ! loop over over the selected atoms
        do idx=1,nsel1
            iatm = a % idxSelection(idx,1)
            do ineigh=1,nneigh(iatm)
                jatm = lneigh(ineigh,iatm)
                if (a % isSelected(jatm,1)) then
                    dij = frame % pos(:,iatm) - frame % pos(:,jatm)
                    dist2 = computeDistanceSquaredPBC(dij)
                    if (dist2 < rcut2) then
                        a % Iarray1D(iatm) = a % Iarray1D(iatm) + 1
                        a % Iarray2D(a % Iarray1D(iatm),iatm) = jatm
                    end if
                end if
            end do
        end do

        do idx=1,nsel1
            iatm = a % idxSelection(idx,1)
            phase = 'UNK'
            vector = 0.d0
            if (a % Iarray1D(iatm) == 12 .or. a % Iarray1D(iatm) == 14) then
                call indentity_environment(a,iatm,phase,vector)
            end if
            environments(idx) = phase
            normal_vector(:,idx) = vector
            ! call workData % compute(ID(6), numberOfValues=3, xValues=vector)
        end do

        ! vector = sum(normal_vector, dim=2) / norm2(sum(normal_vector, dim=2))
        vector = sum(normal_vector, dim=2) / count(environments=="HCP",dim=1)
        ! write(0,*) 'Average normal vector:',vector

        block
            integer :: i
            strings = ["UNK","FCC","HCP","BCC","ICO"]
            do i=1,5
                results(i) = count(environments==strings(i),dim=1)
                ! call workData % compute(ID(i), numberOfValues=1, xValues=rtmp)
            end do
        end block
        ! write(0,*)ctype
        ! write(a % outputFile % funit,*)ctype
        write(a % outputFile % funit,"(5(i8,1x),3(f8.3,1x))")results(1:5),vector

        ! write(123,*) count(environments=="FCC",dim=1)

    end subroutine computeAction

    subroutine indentity_environment(a,iatm,phase,normal_vector) 
        use molecule_analysis
        implicit none
        type(actionTypeDef), target :: a
        integer, intent(in) :: iatm
        character(3), intent(out) :: phase
        real(real64), intent(out) :: normal_vector(3)

        integer :: jatm, katm
        integer :: ineigh
        integer, allocatable :: common_arr(:)
        integer :: common_count
        integer :: nbonds
        integer :: i, j, k, ii, jj, kk
    
        type(bond_type) :: list_of_bonds(20)
        integer, allocatable :: longest_chain(:)
        integer :: chain_length
        character(3), dimension(20) :: signature

        integer :: indices(6)
        real(real64) :: vtmp(3), dij(3), dik(3), d2

        phase = 'UNK'
        normal_vector = 0.d0

        ! write(0,*) 'Atom ',iatm,' has ',a % Iarray1D(iatm),' neighbours'
        do ineigh=1,a % Iarray1D(iatm)
            jatm = a % Iarray2D(ineigh,iatm)
            ! write(0,*) 'Neighbour ',ineigh,' is atom ',jatm

            ! Common neighbours
            call find_common_elements( a % Iarray2D(:,iatm), a % Iarray2D(:,jatm), &
                                       a % Iarray1D(iatm), a % Iarray1D(jatm), &
                                       common_arr, common_count)
                           
            ! Number of bonds between the common neighbours
            nbonds = 0
            do i=1,common_count
                ii = common_arr(i)
                do j=i+1,common_count
                    jj = common_arr(j)
                    if ( findloc( a % Iarray2D(:,ii), jj, dim=1 ) > 0 ) then
                        nbonds = nbonds + 1
                        if (nbonds > 20) then
                            ! write(0,*) 'Too many bonds between common neighbours'
                            return
                        end if
                        list_of_bonds(nbonds) = bond_type(i,j)
                    end if
                end do
            end do

            ! Maximum chain length
            call find_longest_chain(list_of_bonds(1:nbonds), nbonds, longest_chain, chain_length)

            ! write(0,*) 'Common neighbours: ',common_count
            ! write(0,*) 'List of common neighbours:', common_arr
            ! write(0,*) 'Number of bonds:',nbonds
            ! do i=1,nbonds
            !     write(0,*) 'Bond ',i,': ',list_of_bonds(i)%atom1,'-',list_of_bonds(i)%atom2
            ! end do
            ! write(0,*) 'Longest chain (# atoms): ',chain_length
            ! write(0,*) 'Longest chain (# bonds): ',chain_length-1
            ! write(0,*) 'Atoms in the chain: ',longest_chain

            write(signature(ineigh),"(i1,i1,i1)")common_count,nbonds,chain_length-1
        end do
        call identify_signatures(signature,a % Iarray1D(iatm),phase)

        ! Compute the normal for the HCP structure
        if (phase /= 'HCP') return

        i = 0
        do ineigh=1,a % Iarray1D(iatm)
            if (trim(signature(ineigh)) == "422") then
                i = i + 1
                indices(i) = a % Iarray2D(ineigh,iatm)
            end if
        end do

        ! do i=1,6
        !     write(*,*) 'C',frame % pos(:,indices(i))
        ! end do

        do i = 1, 4
            ii = indices(i)
            do j = i + 1, 5
                jj = indices(j)
                do k = j + 1, 6
                    kk = indices(k)
                    dij = frame % pos(:,jj) - frame % pos(:,ii)
                    d2 = computeDistanceSquaredPBC(dij)
                    dik = frame % pos(:,kk) - frame % pos(:,ii)
                    d2 = computeDistanceSquaredPBC(dik)
                    vtmp = cross_product(dij,dik) 
                    if (vtmp(3) < 0) vtmp = -vtmp
                    normal_vector = normal_vector + vtmp/norm2(vtmp)
                end do    
            end do    
        end do    
        normal_vector = normal_vector/norm2(normal_vector)

    end subroutine indentity_environment

    subroutine identify_signatures(signature,n,phase)
        implicit none
        character(len=3), intent(in) :: signature(:)
        character(len=3), intent(out) :: phase
        integer, intent(in) :: n
        integer :: i, idx
        character(len=3), allocatable :: unique_signatures(:)
        integer, allocatable :: unique_counts(:)
        integer :: count

        allocate(unique_signatures(n))
        allocate(unique_counts(n), source=0)
        unique_signatures(1) = signature(1)
        unique_counts(1) = 1
        count = 1
        do i=2,n
            idx = findloc(unique_signatures, signature(i),dim=1)
            if (idx == 0) then
                count = count + 1
                unique_signatures(count) = signature(i)
                unique_counts(count) = 1
            else
                unique_counts(idx) = unique_counts(idx) + 1
            end if
        end do
        ! do i=1,count
        !     write(0,*) 'Signature ',unique_signatures(i),' appears ',unique_counts(i),' times'
        ! end do
        ! stop

        if (count == 1) then
            if (unique_signatures(1) == "421" .and. unique_counts(1) == 12) then
            !    write(*,*) 'FCC structure identified'
                phase = 'FCC'
            else if (unique_signatures(1) == "422" .and. unique_counts(1) == 12) then
            else if (unique_signatures(1) == "555" .and. unique_counts(1) == 12) then
                ! write(*,*) 'ICO structure identified'
                phase = 'ICO'
            end if
        else if (count == 2) then
            if ( (unique_signatures(1) == "421" .and. unique_counts(1) == 6 .and. &
                  unique_signatures(2) == "422" .and. unique_counts(2) == 6) .or. &
                 (unique_signatures(1) == "422" .and. unique_counts(1) == 6 .and. &
                  unique_signatures(2) == "421" .and. unique_counts(2) == 6) ) then
                ! write(*,*) 'HCP structure identified'
                phase = 'HCP'
            end if
        else
            ! write(*,*) 'Unknown structure'
            phase = 'UNK'
        end if

    end subroutine identify_signatures

    subroutine find_common_elements(arr1, arr2, n1, n2, common_arr, common_count)
        implicit none
        
        ! Arguments
        integer, intent(in) :: n1, n2                ! Sizes of input arrays
        integer, intent(in) :: arr1(n1), arr2(n2)    ! Input arrays
        integer, allocatable, intent(out) :: common_arr(:)  ! Output array
        integer, intent(out) :: common_count         ! Number of common elements
        
        ! Local variables
        integer :: i, j, k
        logical :: found
        integer, allocatable :: temp_common(:)
        
        ! Initialize
        allocate(temp_common(min(n1, n2)))  ! Maximum possible common elements
        common_count = 0
        
        ! Find common elements
        do i = 1, n1
            found = .false.
            ! Check if this element was already found
            do k = 1, common_count
                if (arr1(i) == temp_common(k)) then
                    found = .true.
                    exit
                end if
            end do
            
            ! If not already found, search in arr2
            if (.not. found) then
                do j = 1, n2
                    if (arr1(i) == arr2(j)) then
                        common_count = common_count + 1
                        temp_common(common_count) = arr1(i)
                        exit
                    end if
                end do
            end if
        end do
        
        ! Allocate and fill the output array
        if (common_count > 0) then
            allocate(common_arr(common_count))
            common_arr = temp_common(1:common_count)
        else
            allocate(common_arr(0))  ! Zero-size array if no common elements
        end if
        
        deallocate(temp_common)
        
    end subroutine find_common_elements

end module moduleCommonNeighboursAnalysis
