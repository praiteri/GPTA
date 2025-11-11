cmake_minimum_required(VERSION 3.10)
project(gpta Fortran C)

# Add the source files
set(SOURCES
    src/xdr.F90
    src/xdr/xdrfile.c
    src/xdr/xdrfile.h
    src/xdr/xdrfile_trr.c
    src/xdr/xdrfile_trr.h
    src/xdr/xdrfile_xtc.c
    src/xdr/xdrfile_xtc.h
    src/commonNeighboursAnalysis.F90
    src/messagesModule.F90
    src/actionModule.F90
    src/actions.F90
    src/addAtoms.F90
    src/addCentres.F90
    src/addDummyParticle.F90
    src/addElements.F90
    src/addFile.F90
    src/addMolecules.F90
    src/alignMolecule.F90
    src/changeCoordinates.F90
    src/computeSteinhardt.F90
    src/createNewCell.F90
    src/createSurface.F90
    src/dcdfort_common.F90
    src/dcdfort_reader.F90
    src/dcdfort_writer.F90
    src/defineVariables.F90
    src/deleteAtoms.F90
    src/densityMap1D.F90
    src/densityMap2D.F90
    src/densityMap3D.F90
    src/distances.F90
    src/dumpCoordinates.F90
    src/elementsModule.F90
    src/extractClusters.F90
    src/extractFramesByProperty.F90
    src/extractProperties.F90
    src/filesModule.F90
    src/gpta.F90
    src/help.F90
    src/inertiaTensor.F90
    src/legendre.F90
    src/meanSquareDisplacement.F90
    src/molecularProperties.F90
    src/mergeFiles.F90
    src/neighboursModule.F90
    src/parseCommandLine.F90
    src/plumed.F90
    src/property.F90
    src/radialPairDistributionFunction.F90
    src/random_module.F90
    src/ranking_module.F90
    src/readCoordinatesCIF.F90
    src/readCoordinatesMisc.F90
    src/readCoordinatesModule.F90
    src/readWriteCoordinatesDCD.F90
    src/readWriteCoordinatesGULP.F90
    src/readWriteCoordinatesGromacs.F90
    src/readWriteCoordinatesLAMMPS.F90
    src/readWriteCoordinatesPDB.F90
    src/readWriteCoordinatesPDBx.F90
    src/replaceMolecules.F90
    src/residenceTime.F90
    src/resizeArray.F90
    src/scattering_factors.F90
    src/selectAtoms.F90
    src/setAtomAttributes.F90
    src/solvationShell.F90
    src/sortModule.F90
    src/spherical_harmonics.F90
    src/steinhardt.F90
    src/stringsModule.F90
    src/superimposeMolecules.F90
    src/symmetryModule.F90
    src/systemComposition.F90
    src/systemModule.F90
    src/testTemplate.F90
    src/tools.F90
    src/topology.F90
    src/variablesModule.F90
    src/wigner3j.F90
    src/writeCoordinates.F90
    src/writeTopologyPSF.F90
    src/writeTopologyTemplateXML.F90
    src/xray.F90
    src/gsas/genhkl.for
    src/gsas/sglatc.for
    src/gsas/sglcen.for
    src/gsas/sglpak.for
    src/gsas/sgmtml.for
    src/gsas/sgoprn.for
    src/gsas/sgrmat.for
    src/gsas/sgroupnp.for
    src/gsas/sgtrfc.for
)

# Add the executable
add_executable(gpta ${SOURCES})

# Change the output name based on the build type
set_target_properties(gpta PROPERTIES
    DEBUG_POSTFIX "_dbg"
)

# Set the default build type if none is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# Set the Fortran compiler 
if(MPI)
    target_compile_definitions(gpta PRIVATE GPTA_MPI=1)

    # Find MPI package
    find_package(MPI REQUIRED)
    
    # Set include directories
    include_directories(${MPI_INCLUDE_PATH})
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DGPTA_MPI")

    # Specify MPI compiler wrappers
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

if(OMP)
    target_compile_definitions(gpta PRIVATE GPTA_OMP=1)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DGPTA_OMP -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
endif()

if(PLUMED)
    target_compile_definitions(gpta PRIVATE GPTA_PLUMED=1)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DGPTA_PLUMED")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lplumed")
endif()

if(XDR)
    target_compile_definitions(gpta PRIVATE GPTA_XDR=1)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DGPTA_XDR")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -g")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -funroll-all-loops")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wpedantic -fbounds-check -fbacktrace")

message(STATUS "Flags  | " ${CMAKE_Fortran_FLAGS})
#message(STATUS "Linker | " ${CMAKE_EXE_LINKER_FLAGS})
#message(STATUS "Release| " ${CMAKE_Fortran_FLAGS_RELEASE})
#message(STATUS "Debug  | " ${CMAKE_Fortran_FLAGS_DEBUG})

