!disclaimer
subroutine initialiseActions()
  use moduleActions
  use moduleMessages
  use moduleStrings
  use moduleTemplateAction

  use moduleAtomsAttributes
  use moduleDeleteAtoms
  use moduleMolecularTopology
  use moduleReplaceMolecule
  use moduleDumpCoordinates
  use moduleCreateSurface
  use moduleAlignMolecule
  
  use XMLtopology
  
  use moduleExtractSystemProperties
  use moduleRDF
  use moduleModifyCoordinates
  use moduleDensityProfile
  use moduleDensityMap2D
  use moduleDensityMap3D
  use moduleMolecularProperties
  use moduleAddAtoms
  use moduleResidenceTime
  use moduleXrayAction
  use moduleMeanSquareDisplacement
  use moduleExtractClustersAction
  use moduleSolvationShell
  use moduleExtractFramesByProperty
  use moduleSteinhardt
  use moduleCommonNeighboursAnalysis
  use moduleMergeFiles

#ifdef GPTA_PLUMED
  use modulePlumedInterface
#endif

  implicit none
  integer :: i
  allocate(action(0:numberOfActions))
  do i=1,numberOfActions
    action(i) % name = actionType(i)
    action(i) % actionDetails = actionDetails(i)
  enddo

  do i=1,numberOfActions
    
    ! print *, "WWWW",me, trim(actionType(i)),i
    nullify(allActions(i) % work)

    if (actionType(i) == "--test"             ) allActions(i) % work => testAction
        
    if (actionType(i) == "--top"              ) allActions(i) % work => computeTopology
    if (actionType(i) == "--pbc"              ) allActions(i) % work => applyPeriodicboundaryConditions
    if (actionType(i) == "--unwrap"           ) allActions(i) % work => unwrapCoordinates
    if (actionType(i) == "--shift"            ) allActions(i) % work => shiftCoordinates
    if (actionType(i) == "--rescale"          ) allActions(i) % work => rescaleCell
    if (actionType(i) == "--newcell"          ) allActions(i) % work => defineNewCell
    if (actionType(i) == "--fixcom"           ) allActions(i) % work => shiftCOM
    if (actionType(i) == "--repl"             ) allActions(i) % work => replicateCell
    if (actionType(i) == "--mirror"           ) allActions(i) % work => mirrorCell
    if (actionType(i) == "--noclash"          ) allActions(i) % work => removeOverlappingMolecules
    if (actionType(i) == "--delete"           ) allActions(i) % work => deleteAtoms
    if (actionType(i) == "--set"              ) allActions(i) % work => setAtomAttributes
    if (actionType(i) == "--replace"          ) allActions(i) % work => replaceMolecules
    if (actionType(i) == "--add"              ) allActions(i) % work => addAtoms
    if (actionType(i) == "--surface"          ) allActions(i) % work => createSurface
    if (actionType(i) == "--align"            ) allActions(i) % work => alignMolecule
    if (actionType(i) == "--merge"            ) allActions(i) % work => mergeFiles
         
    if (actionType(i) == "--fixcell"          ) allActions(i) % work => fixCell
    if (actionType(i) == "--xml"              ) allActions(i) % work => writeXML

    if (actionType(i) == "--extract "         ) allActions(i) % work => extractSystemProperties
    if (actionType(i) == "--cluster"          ) allActions(i) % work => extractClusters
    if (actionType(i) == "--framesByProperty" ) allActions(i) % work => extractFramesByProperty

    if (actionType(i) == "--gofr"             ) allActions(i) % work => computeRadialPairDistribution
    if (actionType(i) == "--dmap1D"           ) allActions(i) % work => computeDensityProfile
    if (actionType(i) == "--dmap2D"           ) allActions(i) % work => computeDensityMap2D 
    if (actionType(i) == "--dmap3D"           ) allActions(i) % work => computeDensityMap3D 
    if (actionType(i) == "--solvation"        ) allActions(i) % work => solvationShell 
    if (actionType(i) == "--molprop"          ) allActions(i) % work => computeMolecularProperties
    if (actionType(i) == "--restime"          ) allActions(i) % work => computeResidenceTime
    if (actionType(i) == "--xray"             ) allActions(i) % work => computeXrayPowder
    if (actionType(i) == "--msd"              ) allActions(i) % work => computeMSD
    if (actionType(i) == "--stein"            ) allActions(i) % work => computeSteinhardt
    if (actionType(i) == "--cna"              ) allActions(i) % work => cnaAction

#ifdef GPTA_PLUMED
    if (actionType(i) == "--plumed" ) allActions(i) % work => plumedInterface
#endif

    if (actionType(i)(1:3) == "--o" ) allActions(i) % work => writeCoordinates

    ! Check the command exists
    if (.not. associated(allActions(i) % work)) call message(-1,"Unknown Command : "//trim(actionType(i)))

  enddo

  ! call runAllActions()
  ! Run over all the action commands
  block
    integer :: iact
    do iact=1,numberOfActions
      call allActions(iact) % work(action(iact))
    enddo
  end block
  
  return
end subroutine initialiseActions

subroutine runAllActions()
  use moduleActions
  use moduleMessages 
  ! use moduleSystem, only : computeNeighboursList, frameReadSuccessfully, computeNeighboursList, computeNeighboursListOnce
  use moduleNeighbours

  integer :: iact
  real(real64) :: startTime, endTime

  ! Run over the processing commands
  do iact=1,numberOfActions
    ! Compute the neighbours' list before the first action that requires it
    if (action(iact) % requiresNeighboursList .and. computeNeighboursList) call updateNeighboursList(.false.)
    
    startTime = timing()
    call allActions(iact) % work(action(iact))
    endTime = timing()
    action(iact) % timeTally = action(iact) % timeTally + (endTime-startTime)
  enddo
  
end subroutine runAllActions

subroutine runInternalAction(actionName, actionFlags)
    use moduleActions
    use moduleMessages 
    use moduleMolecularTopology
    use moduleModifyCoordinates
    use moduleDumpCoordinates
    character(len=*), intent(in) :: actionName
    character(len=*), intent(in) :: actionFlags
    character(len=maximumActionStringLength) :: flags

    call message(0,"Running required action ",str=actionName)
    call message(1,"...Action flags",str=actionFlags)
    
    flags = actionFlags

    if (actionName == "topology" .or. actionName == "top") then
      action(0) % name = "--top"
      extraWork => computeTopology      
    
    else if (actionName == "pbc") then
      action(0) % name = "--pbc"
      extraWork => applyPeriodicboundaryConditions      
    
    else if (actionName == "dump") then
      action(0) % name = "--o"
      extraWork => writeCoordinates
    
    end if

    action(0) % actionDetails = actionFlags
    action(0) % actionInitialisation = .true.
    action(0) % firstAction = .true.

    call extraWork(action(0)) ! Initialisation
    call extraWork(action(0)) ! First execution
    nullify(extraWork)

end subroutine runInternalAction

subroutine executeOneOffActions()
  use moduleVariables
  use moduleSystem 
  use moduleStrings
  use moduleMessages 
  use moduleActions
  implicit none

  integer :: iact, jcmd
  integer :: i, nw, idx(3)
  character(len=STRLEN) :: w(50)

  first_frame = 1
  last_frame = huge(1)
  stride_frame = 1

  ! Frames to process
  jcmd=0
  do iact=1,numberOfActions

    if (actionType(iact) == "--frames") then

      nlist_frames = 0
      call parse(actionDetails(iact),",",w,nw) 

      if (nw==1) then
        idx = extractFrameIndices(w(1))
        first_frame = idx(1)
        last_frame = idx(2)
        stride_frame = idx(3)

      else
        nlist_frames = nw
        allocate(list_frames(nw))
        do i=1,nw
          read(w(i),*) list_frames(i)
        enddo

      endif

    else if (actionType(iact) == "--frame") then
      read(actionDetails(iact),*) first_frame
      last_frame = first_frame

    else if (actionType(iact) == "--skip") then
      read(actionDetails(iact),*) stride_frame

    else if (actionType(iact) == "--log") then
      call compact(actionDetails(iact))
      open(newunit=io,file=actionDetails(iact),status='unknown')
      ! adv_char = 'yes'
      nProgress = huge(1)

    else if (actionType(iact) == "--last") then
      if (numberOfMpiProcesses >1) call message(-1,"--last only available in serial")
      lastFrameOnly = .true.

    else if (actionType(iact) == "--nt") then
      read(actionDetails(iact),*) ompNumThreads

    else if (actionType(iact) == "--define") then
      call defineVariables(actionDetails(iact))

    else
      jcmd = jcmd + 1
      actionType(jcmd) = actionType(iact)
      actionDetails(jcmd) = actionDetails(iact)
    end if
  enddo
  numberOfActions = jcmd

  return
end subroutine executeOneOffActions

