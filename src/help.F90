!disclaimer
module moduleHelp
  
  implicit none

  public :: help, commandHelp
  private

contains

  subroutine help()
    use moduleVariables
    use moduleMessages
    implicit none
#ifdef GPTA_MPI
    integer :: idx
#endif

    call message(-10)

    call message(0,"GPTA HELP")
    call message(0,'For help about specific commands use the +help flag e.g.')
    call message(1,"gpta.x --gofr +help ")
    call message(0,'List of Available Commands in GPTA')
    
    call message(0,"%%%%% Global")
    call message(0," --nt          :: sets the number of openMP threads to use" )
    call message(0," --define      :: alters the value of some global parameters" )
    
    call message(0,"%%%%% Input/Output")
    call message(0," --i           :: reads the coordinates file(s)" )
    call message(0," --o           :: writes the coordinates file" )
    call message(0," --frame       :: selects one frame for processing" )
    call message(0," --frames      :: selects multiple frames for processing" )
    call message(0," --last        :: selects the last frame for processing" )
    call message(0," --skip        :: sets the stride for processing the frames" )
    call message(0," --log         :: sets the filename for the program output, it replaces the screen output" )

    call message(2)
    call message(0,"%%%%% Actions that change the atomic coordinates")
    call message(0," --pbc         :: applies PBC to the system")
    call message(0," --unwrap      :: removes the PBC wrapping in a trajectory")
    call message(0," --shift       :: translated the atoms")
    call message(0," --rescale     :: rescales the cell")
    call message(0," --fixcom      :: translated the centre of mass")
    call message(0," --repl        :: replicates the cell")
    call message(0," --mirror      :: adds the mirror image of the system")
    call message(0," --noclash     :: removes overlapping molecules")
    call message(0," --delete      :: deletes atoms/molecules in the system")
    call message(0," --set         :: changes the properties of the atoms")
    call message(0," --replace     :: replace molecules in the system")
    call message(0," --add         :: adds atoms/molecules to the system")
    call message(0," --surface     :: creates a surface of a crystal")
    call message(0," --surface     :: aligns the coordinates of a molecule to a reference")
    call message(0," --merge       :: merges more than one coordinates file into one")
    call message(0," --fixcell     :: fixes the cell parameters")
    call message(0," --xml         :: writes the coordinates in XML format")
    call message(0," --cna         :: computes the common neighbours analysis")

    call message(2)
    call message(0,"%%%%% Actions that calculate properties ")
    call message(0," --top         :: calculates of the topology")
    call message(0," --extract     :: computes system-wide properties from the trajectory file")
    call message(0," --cluster     :: extract spherical clusters from the coordinates")

    call message(0," --gofr        :: computes the radial pair distribution function")
    call message(0," --dmap1D      :: computes the 1D density map")
    call message(0," --dmap2D      :: computes the 2D density map")
    call message(0," --dmap3D      :: computes the 3D density map")
    call message(0," --solvation   :: computes the solvent density map around the solute")
    call message(0," --molprop     :: computes molecular properties")
    call message(0," --restime     :: computes the solvent residence time")
    call message(0," --xray        :: computes the x-ray powder diffraction spectrum")
    call message(0," --msd         :: computes the mean square displacement")
    call message(0," --stein       :: computes the Steinhardt's order parameter")

    call message(2)
    call message(0,"%%%%% Special actions ")
#ifdef GPTA_PLUMED
    call message(0, " --plumed      :: conputes properties of the system using PLUMED")
#endif

    call message(1," --test        :: test routine for develpment")

#ifdef GPTA_MPI
    call MPI_Abort(MPI_COMM_WORLD, 9, idx )
#else
    call exit(0)
#endif

  end subroutine help

  subroutine commandHelp(cmd)
    use moduleMessages

    use moduleTemplateAction

    use moduleAtomsAttributes
    use moduleDeleteAtoms
    use moduleMolecularTopology
    use moduleReplaceMolecule
    use moduleDumpCoordinates
    use moduleCreateSurface
    use moduleAlignMolecule
    
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
    character(len=*), intent(in) :: cmd
    call message(-10)

    call message(0,"GPTA HELP for action "//trim(cmd))
    call message(0,"For a more detailed description of the command please refer to the manual.")
    call message(0,"..........................................................................................................")

    if (cmd == "--i" ) then
      call message(0,"This action reads the atomic coordinate from file(s).")
      call message(0,"It takes one or more filenames as argument; it can also be repeated.")
      call message(0,"There is also an alternative form, --iABC, where ABC is the typical file extension identifying the coordinates' type")
      call message(0,"These are the currently supported input coordinates file types")
      call message(0,"     xyz -> Extended XYZ")
      call message(0,"     pdb -> Protein Data Bank file")
      call message(0,"     dcd -> CHARMM DCD file ")
      call message(0,"    dcd2 -> CHARMM DCD file (safe/slow routine)")
      call message(0,"     arc -> TINKER ARC file (with connectivity)")
      call message(0,"    txyz -> TINKER XYZ input file")
      call message(0,"     gin -> GULP input file")
      call message(0,"     gau -> GAUSSIAN output file")
      call message(0,"     gro -> GROMOS input file")
      call message(0,"     lmp -> LAMMPS input file")
      call message(0,"  lmptrj -> LAMMPS custom trajectory file (partially supported)")
      call message(0,"     cif -> Crystallographic Information File ")
      call message(0,"     xtc -> GROMACS XTC trajectory file")
      call message(0,"     trr -> GROMACS TRR trajectory file")
      call message(0,"Examples:")
      call message(0,"  gpta.x --i coord.pdb trajectory.dcd")
      call message(0,"  gpta.x --ixyz coord.data --i traj.pdb")
    end if

    if (cmd == "--o" ) then
      call message(0,"This action write the atomic coordinates to file.")
      call message(0,"It takes one filenames as argument, it can be repeated.")
      call message(0,"There is also an alternative form, --oABC, where ABC is the typical file extension identifying the coordinates' type")
      call message(0,"These are the currently supported input coordinates file types")
      call message(0,"     xyz -> Extended XYZ")
      call message(0,"     pdb -> Protein Data Bank file")
      call message(0,"     pdb -> Protein Data Bank file (with CONECT section)")
      call message(0,"     dcd -> CHARMM DCD file ")
      call message(0,"     lmp -> LAMMPS input file")
      call message(0,"  lmptrj -> LAMMPS custom trajectory file (partially supported)")
      call message(0,"     psf -> CHARMM PSF file")
      call message(0,"     arc -> TINKER ARC file (with connectivity)")
      call message(0,"     gin -> GULP input file")
      call message(0,"    gin2 -> GULP input file (fractional coordinates)")
      call message(0,"     xtc -> GROMACS XTC trajectory file")
      call message(0,"Examples:")
      call message(0,"  gpta.x --i coord.pdb trajectory.dcd")
      call message(0,"  gpta.x --ixyz coord.data --i traj.pdb")
    end if

    if (cmd == "--test"             ) call testActionHelp()

    if (cmd == "--top" ) then
      call message(0,"This action computes the system's topology and groups the atoms in molecules.")
      call message(0,"It is a prerequisity to many actions.")
      call message(0,"Examples:")
      call message(0,"  gpta.x --i coord.pdb --top")
    end if

    if (cmd == "--define"           ) call defineVariablesHelp()

    if (cmd == "--rescale"          ) call rescaleCellHelp()
    if (cmd == "--newcell"          ) call defineNewCellHelp()
    if (cmd == "--shift"            ) call shiftCoordinatesHelp()
    if (cmd == "--fixcom"           ) call shiftCOMHelp()
    if (cmd == "--pbc"              ) call applyPeriodicboundaryConditionsHelp()
    if (cmd == "--unwrap"           ) call unwrapCoordinatesHelp()
    if (cmd == "--repl"             ) call replicateCellHelp()
    if (cmd == "--mirror"           ) call mirrorCellHelp()
    if (cmd == "--noclash"          ) call removeOverlappingMoleculesHelp()
    
    ! if (cmd == "--xml"              ) call writeXMLHelp()
    ! if (cmd == "--frames"           ) call selectFramesHelp()
    ! if (cmd == "--frame"            ) call selectFrameHelp()
    ! if (cmd == "--last"             ) call selectLastFrameHelp()
    ! if (cmd == "--skip"             ) call selectSkipFramesHelp()
    ! if (cmd == "--log"              ) call setLogFileHelp()
    ! if (cmd == "--nt"               ) call setNumberOfThreadsHelp()

    if (cmd == "--delete"           ) call deleteAtomsHelp()
    if (cmd == "--set"              ) call setAtomAttributesHelp()
    if (cmd == "--replace"          ) call replaceMoleculesHelp()
    if (cmd == "--add"              ) call addAtomsHelp()
    if (cmd == "--surface"          ) call createSurfaceHelp()
    if (cmd == "--cluster"          ) call extractClustersHelp()
        
    if (cmd == "--align"            ) call alignMoleculeHelp()
    if (cmd == "--merge"            ) call mergeFilesHelp()
    if (cmd == "--fixcell"          ) call fixCellHelp()
        
    if (cmd == "--extract"          ) call extractSystemPropertiesHelp()
    if (cmd == "--framesByProperty" ) call extractFramesByPropertyHelp()
    if (cmd == "--gofr"             ) call computeRadialPairDistributionHelp()
    if (cmd == "--dmap1D"           ) call computeDensityProfileHelp()
    if (cmd == "--dmap2D"           ) call computeDensityMap2DHelp()
    if (cmd == "--dmap3D"           ) call computeDensityMap3DHelp()
    if (cmd == "--solvation"        ) call solvationShellHelp()
    if (cmd == "--molprop"          ) call computeMolecularPropertiesHelp()
    if (cmd == "--restime"          ) call computeResidenceTimeHelp()
    if (cmd == "--xray"             ) call computeXrayPowderHelp()
    if (cmd == "--msd"              ) call computeMSDHelp()
    if (cmd == "--stein"            ) call computeSteinhardtHelp()
    if (cmd == "--cna"              ) call cnaActionHelp()


#ifdef GPTA_PLUMED
     if (cmd == "--plumed" ) call plumedInterfaceHelp()
#endif

    call message(-3)

  end subroutine

end module moduleHelp
