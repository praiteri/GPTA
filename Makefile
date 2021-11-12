# Linux-specific options
VERSION=4.0beta

EXEDIR = bin
EXE = gpta.x
MPIEXE = gpta_mpi.x
HOST := $(shell hostname)

# Fortran source code
SRCF90  := $(wildcard src/*.F90)
SRCFOR  := $(wildcard src/gsas/*.for)
SRCCXX  := $(wildcard src/*.cpp)
SRCC    := $(wildcard src/xdr/*.c)

OBJF90  := $(patsubst src/%,%,$(SRCF90:.F90=.o)) 
OBJFOR  := $(patsubst src/gsas/%,%,$(SRCFOR:.for=.o))
OBJCXX  := $(patsubst src/%,%,$(SRCCXX:.cpp=.o))
OBJC    := $(patsubst src/xdr/%,%,$(SRCC:.c=.o))
OBJECTS := $(OBJF90) $(OBJX90) $(OBJFOR) $(OBJC)

FFLAGS_for = -std=legacy
LDFLAGS  = -lstdc++
FFLAGS = -ffree-line-length-none -ffpe-trap=zero,overflow,invalid,underflow -finit-real=zero -finit-integer=0


# Gromacs binary coordinates files
xdr=y
ifeq ($(xdr),y)
  CPPFLAGS += -DGPTA_XDR
endif

# OMP parallelisation
omp=y
ifeq ($(omp),y)
  CPPFLAGS += -DGPTA_OMP
  FFLAGS += -fopenmp
endif

# PLUMED interface
plumed=no
ifeq ($(plumed),y)
  $(info Compiling PLUMED Interface)
  CPPFLAGS += -DGPTA_PLUMED
  LDFLAGS += -lplumed
endif

# OpenMM interface
openMM=no
ifeq ($(openmm),y)
  CPPFLAGS += -DGPTA_OPENMM
  $(info Compiling openMM Interface)
  OBJECTS += $(OBJCXX)
  INCLUDE_DIR_OPENMM := -I$(OPENMM_INCLUDE_PATH)
  ifeq ($(cuda),y)
    CPPFLAGS += -DGPTA_CUDA
    ifeq ($(HOST),$(filter $(HOST),topaz-1 topaz-2))
      LDFLAGS += /pawsey/centos7.6/devel/binary/cuda/10.1/lib64/stubs/libcuda.so
    endif
    LDFLAGS += -lOpenMMAmoebaCUDA -lOpenMMCUDA
  endif
  LDFLAGS += -lOpenMMAmoeba -lOpenMM -lOpenMMAmoebaReference
endif

# Debugging flags
dbg=no
ifeq ($(dbg),y)
  FFLAGS += -g -fbounds-check -fbacktrace -Wall
  CFLAGS += -g -funroll-all-loops -O3
else
  FFLAGS += -g -funroll-all-loops -O3
  CFLAGS += -g -funroll-all-loops -O3
endif

# Parallel version
ifeq ($(MAKECMDGOALS),mpi)
  export MPI=yes
endif
ifeq ($(MPI),yes)
  FC=mpif90
  CC=mpicc
  CXX=mpic++
  CPPFLAGS += -DGPTA_MPI
  OBJECTS := $(addprefix obj_mpi/,$(OBJECTS))
else
  FC  := gfortran
  CC  := gcc
  CXX := g++
  OBJECTS := $(addprefix obj/,$(OBJECTS))
endif

help:
	@echo '----- build type -------------------------------------------'
	@echo '     serial  :  serial build'
	@echo '        mpi  :  parallel build'
	@echo '----- build options ----------------------------------------'
	@echo '    omp=y/n  :  openMP options (default y)'
	@echo '    xdr=y/n  :  GROMACS trajectory files (default y)'
	@echo '    dbg=y/n  :  debugging options (default n)'
	@echo ' plumed=y/n  :  link to PLUMED library (default n)'
	@echo ' openmm=y/n  :  link to openMM library (default n)'
	@echo '   cuda=y/n  :  enable GPU use for openMM'
	@echo '----- build options ----------------------------------------'
	@echo '      clean  :  clean serial build'
	@echo '  clean-mpi  :  clean parallel build'
	@echo '  clean-all  :  clean all builds'
	@echo '  distclean  :  clean everything'
	@echo '------------------------------------------------------------'
  
# Start up
start:
	@if [ ! -d bin ]; then mkdir bin; fi
	@if [ ! -d obj ]; then mkdir obj; fi
	@rm -f src/*mod

start_mpi:
	@if [ ! -d bin ]; then mkdir bin; fi
	@if [ ! -d obj_mpi ]; then mkdir obj_mpi; fi
	@rm -f src/*mod

# Executables
serial: start
	$(MAKE) 	F90FLAGS="$(FFLAGS) -J./obj" 	CXXFLAGS="$(CFLAGS) $(INCLUDE_DIR_OPENMM)" 	$(EXE)
	mv $(EXE) $(EXEDIR)/

mpi: start_mpi
	$(MAKE) 	F90FLAGS="$(FFLAGS) -J./obj_mpi" 	CXXFLAGS="$(CFLAGS) $(INCLUDE_DIR_OPENMM)" 	$(MPIEXE)
	mv $(MPIEXE) $(EXEDIR)/

clean-all:
	rm -fr obj obj_mpi

clean:
	rm -fr obj

clean-mpi:
	rm -fr obj_mpi

distclean:
	rm -fr obj obj_mpi bin


#
# Executable generation
#
$(EXE): $(OBJECTS)
	$(FC) $(OBJECTS) $(FFLAGS) $(LDFLAGS) -o $(EXE)

$(MPIEXE): $(OBJECTS)
	$(FC) $(OBJECTS) $(FFLAGS) $(LDFLAGS) -o $(MPIEXE)

#
# Serial dependencies
#
obj/actionModule.o: ./src/actionModule.F90 obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/actionModule.F90 -c -o $@
obj/actions.o: ./src/actions.F90 obj/actionModule.o  obj/addAtoms.o  obj/amoeba.o  obj/setAtomAttributes.o  obj/createSurface.o  obj/deleteAtoms.o  obj/densityMap2D.o  obj/densityMap3D.o  obj/densityMap1D.o  obj/writeCoordinates.o  obj/extractClusters.o  obj/extractFramesByProperty.o  obj/extractProperties.o  obj/meanSquareDisplacement.o  obj/messagesModule.o  obj/changeCoordinates.o  obj/molecularProperties.o  obj/topology.o  obj/neighboursModule.o  obj/openMM.o  obj/plumed.o  obj/radialPairDistributionFunction.o  obj/replaceMolecules.o  obj/residenceTime.o  obj/solvationShell.o  obj/stringsModule.o  obj/systemModule.o  obj/testTemplate.o  obj/variablesModule.o  obj/xray.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/actions.F90 -c -o $@
obj/addAtoms.o: ./src/addAtoms.F90 obj/addCentres.o  obj/addDummyParticle.o  obj/addElements.o  obj/addMolecules.o  obj/actionModule.o  obj/distances.o  obj/messagesModule.o  obj/neighboursModule.o  obj/stringsModule.o  obj/systemModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addAtoms.F90 -c -o $@
obj/addCentres.o: ./src/addCentres.F90 obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addCentres.F90 -c -o $@
obj/addDummyParticle.o: ./src/addDummyParticle.F90 obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addDummyParticle.F90 -c -o $@
obj/addElements.o: ./src/addElements.F90 obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/random_module.o  obj/readCoordinatesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addElements.F90 -c -o $@
obj/addMolecules.o: ./src/addMolecules.F90 obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/random_module.o  obj/readCoordinatesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addMolecules.F90 -c -o $@
obj/align_molecules.o: ./src/align_molecules.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/align_molecules.F90 -c -o $@
obj/amoeba.o: ./src/amoeba.F90 obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/openMM.o  obj/property.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/amoeba.F90 -c -o $@
obj/changeCoordinates.o: ./src/changeCoordinates.F90 obj/align_molecules.o  obj/distances.o  obj/elementsModule.o  obj/messagesModule.o  obj/neighboursModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/changeCoordinates.F90 -c -o $@
obj/createSurface.o: ./src/createSurface.F90 obj/ranking_module.o  obj/filesModule.o  obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/createSurface.F90 -c -o $@
obj/dcdfort_common.o: ./src/dcdfort_common.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/dcdfort_common.F90 -c -o $@
obj/dcdfort_reader.o: ./src/dcdfort_reader.F90 obj/dcdfort_common.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/dcdfort_reader.F90 -c -o $@
obj/dcdfort_writer.o: ./src/dcdfort_writer.F90 obj/dcdfort_common.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/dcdfort_writer.F90 -c -o $@
obj/defineVariables.o: ./src/defineVariables.F90 obj/elementsModule.o  obj/messagesModule.o  obj/neighboursModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/defineVariables.F90 -c -o $@
obj/deleteAtoms.o: ./src/deleteAtoms.F90 obj/messagesModule.o  obj/neighboursModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/deleteAtoms.F90 -c -o $@
obj/densityMap1D.o: ./src/densityMap1D.F90 obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/property.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/densityMap1D.F90 -c -o $@
obj/densityMap2D.o: ./src/densityMap2D.F90 obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/property.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/densityMap2D.F90 -c -o $@
obj/densityMap3D.o: ./src/densityMap3D.F90 obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/densityMap3D.F90 -c -o $@
obj/distances.o: ./src/distances.F90 obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/distances.F90 -c -o $@
obj/dumpCoordinates.o: ./src/dumpCoordinates.F90 obj/elementsModule.o  obj/messagesModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/dumpCoordinates.F90 -c -o $@
obj/elementsModule.o: ./src/elementsModule.F90 obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/elementsModule.F90 -c -o $@
obj/extractClusters.o: ./src/extractClusters.F90 obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/extractClusters.F90 -c -o $@
obj/extractFramesByProperty.o: ./src/extractFramesByProperty.F90 obj/filesModule.o  obj/messagesModule.o  obj/property.o  obj/resizeArray.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/extractFramesByProperty.F90 -c -o $@
obj/extractProperties.o: ./src/extractProperties.F90 obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/property.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/extractProperties.F90 -c -o $@
obj/filesModule.o: ./src/filesModule.F90 obj/messagesModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/filesModule.F90 -c -o $@
obj/gpta.o: ./src/gpta.F90 obj/actionModule.o  obj/distances.o  obj/elementsModule.o  obj/messagesModule.o  obj/neighboursModule.o  obj/property.o  obj/random_module.o  obj/readCoordinatesModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/gpta.F90 -c -o $@
obj/help.o: ./src/help.F90 obj/addAtoms.o  obj/amoeba.o  obj/setAtomAttributes.o  obj/createSurface.o  obj/deleteAtoms.o  obj/densityMap2D.o  obj/densityMap3D.o  obj/densityMap1D.o  obj/writeCoordinates.o  obj/extractClusters.o  obj/extractProperties.o  obj/meanSquareDisplacement.o  obj/messagesModule.o  obj/changeCoordinates.o  obj/molecularProperties.o  obj/topology.o  obj/openMM.o  obj/plumed.o  obj/radialPairDistributionFunction.o  obj/replaceMolecules.o  obj/residenceTime.o  obj/solvationShell.o  obj/testTemplate.o  obj/variablesModule.o  obj/xray.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/help.F90 -c -o $@
obj/inertiaTensor.o: ./src/inertiaTensor.F90 obj/elementsModule.o  obj/systemModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/inertiaTensor.F90 -c -o $@
obj/meanSquareDisplacement.o: ./src/meanSquareDisplacement.F90 obj/filesModule.o  obj/messagesModule.o  obj/property.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/meanSquareDisplacement.F90 -c -o $@
obj/messagesModule.o: ./src/messagesModule.F90 obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/messagesModule.F90 -c -o $@
obj/molecularProperties.o: ./src/molecularProperties.F90 obj/actionModule.o  obj/distances.o  obj/elementsModule.o  obj/filesModule.o  obj/messagesModule.o  obj/property.o  obj/resizeArray.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/molecularProperties.F90 -c -o $@
obj/neighboursModule.o: ./src/neighboursModule.F90 obj/ranking_module.o  obj/actionModule.o  obj/distances.o  obj/messagesModule.o  obj/resizeArray.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/neighboursModule.F90 -c -o $@
obj/openMM.o: ./src/openMM.F90 obj/elementsModule.o  obj/messagesModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/openMM.F90 -c -o $@
obj/parseCommandLine.o: ./src/parseCommandLine.F90 obj/actionModule.o  obj/help.o  obj/messagesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/parseCommandLine.F90 -c -o $@
obj/plumed.o: ./src/plumed.F90 obj/elementsModule.o  obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/plumed.F90 -c -o $@
obj/property.o: ./src/property.F90 obj/resizeArray.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/property.F90 -c -o $@
obj/radialPairDistributionFunction.o: ./src/radialPairDistributionFunction.F90 obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/radialPairDistributionFunction.F90 -c -o $@
obj/random_module.o: ./src/random_module.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/random_module.F90 -c -o $@
obj/ranking_module.o: ./src/ranking_module.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/ranking_module.F90 -c -o $@
obj/readCoordinatesCIF.o: ./src/readCoordinatesCIF.F90 obj/messagesModule.o  obj/symmetryModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readCoordinatesCIF.F90 -c -o $@
obj/readCoordinatesMisc.o: ./src/readCoordinatesMisc.F90 obj/elementsModule.o  obj/stringsModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readCoordinatesMisc.F90 -c -o $@
obj/readCoordinatesModule.o: ./src/readCoordinatesModule.F90 obj/dcdfort_common.o  obj/dcdfort_reader.o  obj/actionModule.o  obj/readCoordinatesCIF.o  obj/distances.o  obj/elementsModule.o  obj/filesModule.o  obj/readWriteCoordinatesGULP.o  obj/messagesModule.o  obj/neighboursModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o  obj/xdr.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readCoordinatesModule.F90 -c -o $@
obj/readWriteCoordinatesDCD.o: ./src/readWriteCoordinatesDCD.F90 obj/messagesModule.o  obj/systemModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesDCD.F90 -c -o $@
obj/readWriteCoordinatesGULP.o: ./src/readWriteCoordinatesGULP.F90 obj/messagesModule.o  obj/resizeArray.o  obj/stringsModule.o  obj/symmetryModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesGULP.F90 -c -o $@
obj/readWriteCoordinatesGromacs.o: ./src/readWriteCoordinatesGromacs.F90 obj/elementsModule.o  obj/messagesModule.o  obj/stringsModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesGromacs.F90 -c -o $@
obj/readWriteCoordinatesLAMMPS.o: ./src/readWriteCoordinatesLAMMPS.F90 obj/filesModule.o  obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesLAMMPS.F90 -c -o $@
obj/readWriteCoordinatesPDB.o: ./src/readWriteCoordinatesPDB.F90 obj/elementsModule.o  obj/messagesModule.o  obj/systemModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesPDB.F90 -c -o $@
obj/replaceMolecules.o: ./src/replaceMolecules.F90 obj/align_molecules.o  obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/neighboursModule.o  obj/openMM.o  obj/property.o  obj/readCoordinatesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/replaceMolecules.F90 -c -o $@
obj/residenceTime.o: ./src/residenceTime.F90 obj/filesModule.o  obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/residenceTime.F90 -c -o $@
obj/resizeArray.o: ./src/resizeArray.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/resizeArray.F90 -c -o $@
obj/scattering_factors.o: ./src/scattering_factors.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/scattering_factors.F90 -c -o $@
obj/selectAtoms.o: ./src/selectAtoms.F90 obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/selectAtoms.F90 -c -o $@
obj/setAtomAttributes.o: ./src/setAtomAttributes.F90 obj/elementsModule.o  obj/messagesModule.o  obj/openMM.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/setAtomAttributes.F90 -c -o $@
obj/solvationShell.o: ./src/solvationShell.F90 obj/align_molecules.o  obj/distances.o  obj/elementsModule.o  obj/filesModule.o  obj/messagesModule.o  obj/readCoordinatesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/solvationShell.F90 -c -o $@
obj/sortModule.o: ./src/sortModule.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/sortModule.F90 -c -o $@
obj/stringsModule.o: ./src/stringsModule.F90 obj/messagesModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/stringsModule.F90 -c -o $@
obj/symmetryModule.o: ./src/symmetryModule.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/symmetryModule.F90 -c -o $@
obj/systemComposition.o: ./src/systemComposition.F90 obj/messagesModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/systemComposition.F90 -c -o $@
obj/systemModule.o: ./src/systemModule.F90 obj/messagesModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/systemModule.F90 -c -o $@
obj/testTemplate.o: ./src/testTemplate.F90 obj/ranking_module.o  obj/distances.o  obj/filesModule.o  obj/messagesModule.o  obj/property.o  obj/resizeArray.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/testTemplate.F90 -c -o $@
obj/tools.o: ./src/tools.F90 obj/ranking_module.o  obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/tools.F90 -c -o $@
obj/topology.o: ./src/topology.F90 obj/sortModule.o  obj/distances.o  obj/elementsModule.o  obj/filesModule.o  obj/messagesModule.o  obj/neighboursModule.o  obj/resizeArray.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/topology.F90 -c -o $@
obj/variablesModule.o: ./src/variablesModule.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/variablesModule.F90 -c -o $@
obj/writeCoordinates.o: ./src/writeCoordinates.F90 obj/filesModule.o  obj/readWriteCoordinatesGULP.o  obj/readWriteCoordinatesLAMMPS.o  obj/messagesModule.o  obj/neighboursModule.o  obj/openMM.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o  obj/xdr.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/writeCoordinates.F90 -c -o $@
obj/writeTopologyPSF.o: ./src/writeTopologyPSF.F90 obj/elementsModule.o  obj/systemModule.o  obj/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/writeTopologyPSF.F90 -c -o $@
obj/xdr.o: ./src/xdr.F90 obj/messagesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/xdr.F90 -c -o $@
obj/xray.o: ./src/xray.F90 obj/ranking_module.o  obj/elementsModule.o  obj/filesModule.o  obj/messagesModule.o  obj/stringsModule.o  obj/systemModule.o  obj/variablesModule.o  obj/scattering_factors.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/xray.F90 -c -o $@
obj/interfaceOpenMM.o: ./src/interfaceOpenMM.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/interfaceOpenMM.cpp -c -I. -o $@
obj/genhkl.o: ./src/gsas/genhkl.for
	$(FC) -O3 src/gsas/genhkl.for $(FFLAGS_for) -c -I. -o $@
obj/sglatc.o: ./src/gsas/sglatc.for
	$(FC) -O3 src/gsas/sglatc.for $(FFLAGS_for) -c -I. -o $@
obj/sglcen.o: ./src/gsas/sglcen.for
	$(FC) -O3 src/gsas/sglcen.for $(FFLAGS_for) -c -I. -o $@
obj/sglpak.o: ./src/gsas/sglpak.for
	$(FC) -O3 src/gsas/sglpak.for $(FFLAGS_for) -c -I. -o $@
obj/sgmtml.o: ./src/gsas/sgmtml.for
	$(FC) -O3 src/gsas/sgmtml.for $(FFLAGS_for) -c -I. -o $@
obj/sgoprn.o: ./src/gsas/sgoprn.for
	$(FC) -O3 src/gsas/sgoprn.for $(FFLAGS_for) -c -I. -o $@
obj/sgrmat.o: ./src/gsas/sgrmat.for
	$(FC) -O3 src/gsas/sgrmat.for $(FFLAGS_for) -c -I. -o $@
obj/sgroupnp.o: ./src/gsas/sgroupnp.for
	$(FC) -O3 src/gsas/sgroupnp.for $(FFLAGS_for) -c -I. -o $@
obj/sgtrfc.o: ./src/gsas/sgtrfc.for
	$(FC) -O3 src/gsas/sgtrfc.for $(FFLAGS_for) -c -I. -o $@
obj/xdrfile.o: ./src/xdr/xdrfile.c
	$(CC) $(CFLAGS) src/xdr/xdrfile.c -c -I. -o $@
obj/xdrfile_trr.o: ./src/xdr/xdrfile_trr.c
	$(CC) $(CFLAGS) src/xdr/xdrfile_trr.c -c -I. -o $@
obj/xdrfile_xtc.o: ./src/xdr/xdrfile_xtc.c
	$(CC) $(CFLAGS) src/xdr/xdrfile_xtc.c -c -I. -o $@
obj_mpi/actionModule.o: ./src/actionModule.F90 obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/actionModule.F90 -c -o $@
obj_mpi/actions.o: ./src/actions.F90 obj_mpi/actionModule.o  obj_mpi/addAtoms.o  obj_mpi/amoeba.o  obj_mpi/setAtomAttributes.o  obj_mpi/createSurface.o  obj_mpi/deleteAtoms.o  obj_mpi/densityMap2D.o  obj_mpi/densityMap3D.o  obj_mpi/densityMap1D.o  obj_mpi/writeCoordinates.o  obj_mpi/extractClusters.o  obj_mpi/extractFramesByProperty.o  obj_mpi/extractProperties.o  obj_mpi/meanSquareDisplacement.o  obj_mpi/messagesModule.o  obj_mpi/changeCoordinates.o  obj_mpi/molecularProperties.o  obj_mpi/topology.o  obj_mpi/neighboursModule.o  obj_mpi/openMM.o  obj_mpi/plumed.o  obj_mpi/radialPairDistributionFunction.o  obj_mpi/replaceMolecules.o  obj_mpi/residenceTime.o  obj_mpi/solvationShell.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/testTemplate.o  obj_mpi/variablesModule.o  obj_mpi/xray.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/actions.F90 -c -o $@
obj_mpi/addAtoms.o: ./src/addAtoms.F90 obj_mpi/addCentres.o  obj_mpi/addDummyParticle.o  obj_mpi/addElements.o  obj_mpi/addMolecules.o  obj_mpi/actionModule.o  obj_mpi/distances.o  obj_mpi/messagesModule.o  obj_mpi/neighboursModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addAtoms.F90 -c -o $@
obj_mpi/addCentres.o: ./src/addCentres.F90 obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addCentres.F90 -c -o $@
obj_mpi/addDummyParticle.o: ./src/addDummyParticle.F90 obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addDummyParticle.F90 -c -o $@
obj_mpi/addElements.o: ./src/addElements.F90 obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/random_module.o  obj_mpi/readCoordinatesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addElements.F90 -c -o $@
obj_mpi/addMolecules.o: ./src/addMolecules.F90 obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/random_module.o  obj_mpi/readCoordinatesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/addMolecules.F90 -c -o $@
obj_mpi/align_molecules.o: ./src/align_molecules.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/align_molecules.F90 -c -o $@
obj_mpi/amoeba.o: ./src/amoeba.F90 obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/openMM.o  obj_mpi/property.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/amoeba.F90 -c -o $@
obj_mpi/changeCoordinates.o: ./src/changeCoordinates.F90 obj_mpi/align_molecules.o  obj_mpi/distances.o  obj_mpi/elementsModule.o  obj_mpi/messagesModule.o  obj_mpi/neighboursModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/changeCoordinates.F90 -c -o $@
obj_mpi/createSurface.o: ./src/createSurface.F90 obj_mpi/ranking_module.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/createSurface.F90 -c -o $@
obj_mpi/dcdfort_common.o: ./src/dcdfort_common.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/dcdfort_common.F90 -c -o $@
obj_mpi/dcdfort_reader.o: ./src/dcdfort_reader.F90 obj_mpi/dcdfort_common.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/dcdfort_reader.F90 -c -o $@
obj_mpi/dcdfort_writer.o: ./src/dcdfort_writer.F90 obj_mpi/dcdfort_common.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/dcdfort_writer.F90 -c -o $@
obj_mpi/defineVariables.o: ./src/defineVariables.F90 obj_mpi/elementsModule.o  obj_mpi/messagesModule.o  obj_mpi/neighboursModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/defineVariables.F90 -c -o $@
obj_mpi/deleteAtoms.o: ./src/deleteAtoms.F90 obj_mpi/messagesModule.o  obj_mpi/neighboursModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/deleteAtoms.F90 -c -o $@
obj_mpi/densityMap1D.o: ./src/densityMap1D.F90 obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/property.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/densityMap1D.F90 -c -o $@
obj_mpi/densityMap2D.o: ./src/densityMap2D.F90 obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/property.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/densityMap2D.F90 -c -o $@
obj_mpi/densityMap3D.o: ./src/densityMap3D.F90 obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/densityMap3D.F90 -c -o $@
obj_mpi/distances.o: ./src/distances.F90 obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/distances.F90 -c -o $@
obj_mpi/dumpCoordinates.o: ./src/dumpCoordinates.F90 obj_mpi/elementsModule.o  obj_mpi/messagesModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/dumpCoordinates.F90 -c -o $@
obj_mpi/elementsModule.o: ./src/elementsModule.F90 obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/elementsModule.F90 -c -o $@
obj_mpi/extractClusters.o: ./src/extractClusters.F90 obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/extractClusters.F90 -c -o $@
obj_mpi/extractFramesByProperty.o: ./src/extractFramesByProperty.F90 obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/property.o  obj_mpi/resizeArray.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/extractFramesByProperty.F90 -c -o $@
obj_mpi/extractProperties.o: ./src/extractProperties.F90 obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/property.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/extractProperties.F90 -c -o $@
obj_mpi/filesModule.o: ./src/filesModule.F90 obj_mpi/messagesModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/filesModule.F90 -c -o $@
obj_mpi/gpta.o: ./src/gpta.F90 obj_mpi/actionModule.o  obj_mpi/distances.o  obj_mpi/elementsModule.o  obj_mpi/messagesModule.o  obj_mpi/neighboursModule.o  obj_mpi/property.o  obj_mpi/random_module.o  obj_mpi/readCoordinatesModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/gpta.F90 -c -o $@
obj_mpi/help.o: ./src/help.F90 obj_mpi/addAtoms.o  obj_mpi/amoeba.o  obj_mpi/setAtomAttributes.o  obj_mpi/createSurface.o  obj_mpi/deleteAtoms.o  obj_mpi/densityMap2D.o  obj_mpi/densityMap3D.o  obj_mpi/densityMap1D.o  obj_mpi/writeCoordinates.o  obj_mpi/extractClusters.o  obj_mpi/extractProperties.o  obj_mpi/meanSquareDisplacement.o  obj_mpi/messagesModule.o  obj_mpi/changeCoordinates.o  obj_mpi/molecularProperties.o  obj_mpi/topology.o  obj_mpi/openMM.o  obj_mpi/plumed.o  obj_mpi/radialPairDistributionFunction.o  obj_mpi/replaceMolecules.o  obj_mpi/residenceTime.o  obj_mpi/solvationShell.o  obj_mpi/testTemplate.o  obj_mpi/variablesModule.o  obj_mpi/xray.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/help.F90 -c -o $@
obj_mpi/inertiaTensor.o: ./src/inertiaTensor.F90 obj_mpi/elementsModule.o  obj_mpi/systemModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/inertiaTensor.F90 -c -o $@
obj_mpi/meanSquareDisplacement.o: ./src/meanSquareDisplacement.F90 obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/property.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/meanSquareDisplacement.F90 -c -o $@
obj_mpi/messagesModule.o: ./src/messagesModule.F90 obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/messagesModule.F90 -c -o $@
obj_mpi/molecularProperties.o: ./src/molecularProperties.F90 obj_mpi/actionModule.o  obj_mpi/distances.o  obj_mpi/elementsModule.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/property.o  obj_mpi/resizeArray.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/molecularProperties.F90 -c -o $@
obj_mpi/neighboursModule.o: ./src/neighboursModule.F90 obj_mpi/ranking_module.o  obj_mpi/actionModule.o  obj_mpi/distances.o  obj_mpi/messagesModule.o  obj_mpi/resizeArray.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/neighboursModule.F90 -c -o $@
obj_mpi/openMM.o: ./src/openMM.F90 obj_mpi/elementsModule.o  obj_mpi/messagesModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/openMM.F90 -c -o $@
obj_mpi/parseCommandLine.o: ./src/parseCommandLine.F90 obj_mpi/actionModule.o  obj_mpi/help.o  obj_mpi/messagesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/parseCommandLine.F90 -c -o $@
obj_mpi/plumed.o: ./src/plumed.F90 obj_mpi/elementsModule.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/plumed.F90 -c -o $@
obj_mpi/property.o: ./src/property.F90 obj_mpi/resizeArray.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/property.F90 -c -o $@
obj_mpi/radialPairDistributionFunction.o: ./src/radialPairDistributionFunction.F90 obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/radialPairDistributionFunction.F90 -c -o $@
obj_mpi/random_module.o: ./src/random_module.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/random_module.F90 -c -o $@
obj_mpi/ranking_module.o: ./src/ranking_module.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/ranking_module.F90 -c -o $@
obj_mpi/readCoordinatesCIF.o: ./src/readCoordinatesCIF.F90 obj_mpi/messagesModule.o  obj_mpi/symmetryModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readCoordinatesCIF.F90 -c -o $@
obj_mpi/readCoordinatesMisc.o: ./src/readCoordinatesMisc.F90 obj_mpi/elementsModule.o  obj_mpi/stringsModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readCoordinatesMisc.F90 -c -o $@
obj_mpi/readCoordinatesModule.o: ./src/readCoordinatesModule.F90 obj_mpi/dcdfort_common.o  obj_mpi/dcdfort_reader.o  obj_mpi/actionModule.o  obj_mpi/readCoordinatesCIF.o  obj_mpi/distances.o  obj_mpi/elementsModule.o  obj_mpi/filesModule.o  obj_mpi/readWriteCoordinatesGULP.o  obj_mpi/messagesModule.o  obj_mpi/neighboursModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o  obj_mpi/xdr.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readCoordinatesModule.F90 -c -o $@
obj_mpi/readWriteCoordinatesDCD.o: ./src/readWriteCoordinatesDCD.F90 obj_mpi/messagesModule.o  obj_mpi/systemModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesDCD.F90 -c -o $@
obj_mpi/readWriteCoordinatesGULP.o: ./src/readWriteCoordinatesGULP.F90 obj_mpi/messagesModule.o  obj_mpi/resizeArray.o  obj_mpi/stringsModule.o  obj_mpi/symmetryModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesGULP.F90 -c -o $@
obj_mpi/readWriteCoordinatesGromacs.o: ./src/readWriteCoordinatesGromacs.F90 obj_mpi/elementsModule.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesGromacs.F90 -c -o $@
obj_mpi/readWriteCoordinatesLAMMPS.o: ./src/readWriteCoordinatesLAMMPS.F90 obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesLAMMPS.F90 -c -o $@
obj_mpi/readWriteCoordinatesPDB.o: ./src/readWriteCoordinatesPDB.F90 obj_mpi/elementsModule.o  obj_mpi/messagesModule.o  obj_mpi/systemModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/readWriteCoordinatesPDB.F90 -c -o $@
obj_mpi/replaceMolecules.o: ./src/replaceMolecules.F90 obj_mpi/align_molecules.o  obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/neighboursModule.o  obj_mpi/openMM.o  obj_mpi/property.o  obj_mpi/readCoordinatesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/replaceMolecules.F90 -c -o $@
obj_mpi/residenceTime.o: ./src/residenceTime.F90 obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/residenceTime.F90 -c -o $@
obj_mpi/resizeArray.o: ./src/resizeArray.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/resizeArray.F90 -c -o $@
obj_mpi/scattering_factors.o: ./src/scattering_factors.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/scattering_factors.F90 -c -o $@
obj_mpi/selectAtoms.o: ./src/selectAtoms.F90 obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/selectAtoms.F90 -c -o $@
obj_mpi/setAtomAttributes.o: ./src/setAtomAttributes.F90 obj_mpi/elementsModule.o  obj_mpi/messagesModule.o  obj_mpi/openMM.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/setAtomAttributes.F90 -c -o $@
obj_mpi/solvationShell.o: ./src/solvationShell.F90 obj_mpi/align_molecules.o  obj_mpi/distances.o  obj_mpi/elementsModule.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/readCoordinatesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/solvationShell.F90 -c -o $@
obj_mpi/sortModule.o: ./src/sortModule.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/sortModule.F90 -c -o $@
obj_mpi/stringsModule.o: ./src/stringsModule.F90 obj_mpi/messagesModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/stringsModule.F90 -c -o $@
obj_mpi/symmetryModule.o: ./src/symmetryModule.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/symmetryModule.F90 -c -o $@
obj_mpi/systemComposition.o: ./src/systemComposition.F90 obj_mpi/messagesModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/systemComposition.F90 -c -o $@
obj_mpi/systemModule.o: ./src/systemModule.F90 obj_mpi/messagesModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/systemModule.F90 -c -o $@
obj_mpi/testTemplate.o: ./src/testTemplate.F90 obj_mpi/ranking_module.o  obj_mpi/distances.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/property.o  obj_mpi/resizeArray.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/testTemplate.F90 -c -o $@
obj_mpi/tools.o: ./src/tools.F90 obj_mpi/ranking_module.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/tools.F90 -c -o $@
obj_mpi/topology.o: ./src/topology.F90 obj_mpi/sortModule.o  obj_mpi/distances.o  obj_mpi/elementsModule.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/neighboursModule.o  obj_mpi/resizeArray.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/topology.F90 -c -o $@
obj_mpi/variablesModule.o: ./src/variablesModule.F90
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/variablesModule.F90 -c -o $@
obj_mpi/writeCoordinates.o: ./src/writeCoordinates.F90 obj_mpi/filesModule.o  obj_mpi/readWriteCoordinatesGULP.o  obj_mpi/readWriteCoordinatesLAMMPS.o  obj_mpi/messagesModule.o  obj_mpi/neighboursModule.o  obj_mpi/openMM.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o  obj_mpi/xdr.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/writeCoordinates.F90 -c -o $@
obj_mpi/writeTopologyPSF.o: ./src/writeTopologyPSF.F90 obj_mpi/elementsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/writeTopologyPSF.F90 -c -o $@
obj_mpi/xdr.o: ./src/xdr.F90 obj_mpi/messagesModule.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/xdr.F90 -c -o $@
obj_mpi/xray.o: ./src/xray.F90 obj_mpi/ranking_module.o  obj_mpi/elementsModule.o  obj_mpi/filesModule.o  obj_mpi/messagesModule.o  obj_mpi/stringsModule.o  obj_mpi/systemModule.o  obj_mpi/variablesModule.o  obj_mpi/scattering_factors.o 
	$(FC) $(CPPFLAGS) $(F90FLAGS) src/xray.F90 -c -o $@
obj_mpi/interfaceOpenMM.o: ./src/interfaceOpenMM.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/interfaceOpenMM.cpp -c -I. -o $@
obj_mpi/genhkl.o: ./src/gsas/genhkl.for
	$(FC) -O3 src/gsas/genhkl.for $(FFLAGS_for) -c -I. -o $@
obj_mpi/sglatc.o: ./src/gsas/sglatc.for
	$(FC) -O3 src/gsas/sglatc.for $(FFLAGS_for) -c -I. -o $@
obj_mpi/sglcen.o: ./src/gsas/sglcen.for
	$(FC) -O3 src/gsas/sglcen.for $(FFLAGS_for) -c -I. -o $@
obj_mpi/sglpak.o: ./src/gsas/sglpak.for
	$(FC) -O3 src/gsas/sglpak.for $(FFLAGS_for) -c -I. -o $@
obj_mpi/sgmtml.o: ./src/gsas/sgmtml.for
	$(FC) -O3 src/gsas/sgmtml.for $(FFLAGS_for) -c -I. -o $@
obj_mpi/sgoprn.o: ./src/gsas/sgoprn.for
	$(FC) -O3 src/gsas/sgoprn.for $(FFLAGS_for) -c -I. -o $@
obj_mpi/sgrmat.o: ./src/gsas/sgrmat.for
	$(FC) -O3 src/gsas/sgrmat.for $(FFLAGS_for) -c -I. -o $@
obj_mpi/sgroupnp.o: ./src/gsas/sgroupnp.for
	$(FC) -O3 src/gsas/sgroupnp.for $(FFLAGS_for) -c -I. -o $@
obj_mpi/sgtrfc.o: ./src/gsas/sgtrfc.for
	$(FC) -O3 src/gsas/sgtrfc.for $(FFLAGS_for) -c -I. -o $@
obj_mpi/xdrfile.o: ./src/xdr/xdrfile.c
	$(CC) $(CFLAGS) src/xdr/xdrfile.c -c -I. -o $@
obj_mpi/xdrfile_trr.o: ./src/xdr/xdrfile_trr.c
	$(CC) $(CFLAGS) src/xdr/xdrfile_trr.c -c -I. -o $@
obj_mpi/xdrfile_xtc.o: ./src/xdr/xdrfile_xtc.c
	$(CC) $(CFLAGS) src/xdr/xdrfile_xtc.c -c -I. -o $@
