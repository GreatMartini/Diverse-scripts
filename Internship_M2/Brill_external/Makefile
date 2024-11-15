ifneq (,$(wildcard ./ExternalInitialData.cpp))
	NOT_FOUND = 1
endif
ifneq (,$(wildcard ./ExternalInitialData.h))
	NOT_FOUND = 1
endif
define n


endef
ifndef NOT_FOUND
$(error ExternalInitialData.cpp or ExternalInitialData.h not found. $n \Any of the providen external initial data classes can be used. $n \Rename the proper files, set path to external libraries and compile again.)
endif
SHELL = /bin/sh
HDF5_HOME=/software/lamod/bin/hdf5/1.10
SAMRAI_HOME=/software/lamod/bin/SAMRAI/4.1
MPI_HOME=/software/lamod/bin/openmpi/1.10
LORENE_HOME=/path/to/Lorene
#Uncomment next line for external ReprimAnd EOS (https://zenodo.org/record/3785075  Requires installation)
#EXTERNAL_EOS_HOME=/paht/of/RePrimAnd
ifeq ($(CXX),g++)
  PERFFLAGS=-fno-math-errno -fno-signaling-nans -fomit-frame-pointer -O3
endif
ifeq ($(CXX),gcc)
  PERFFLAGS=-fno-math-errno -fno-signaling-nans -fomit-frame-pointer -O3
endif
ifeq ($(CXX),icc)
  PERFFLAGS=-O3 -xHost -ipo -fno-math-errno -fno-signaling-nans -fp-model precise -fp-model source
endif
ifeq ($(CXX),icpc)
  PERFFLAGS=-O3 -xHost -ipo -fno-math-errno -fno-signaling-nans -fp-model precise -fp-model source
endif
ifdef EXTERNAL_EOS_HOME
  PERFFLAGS += -DEXTERNAL_EOS
  EXTERNAL_INCLUDE += -I$(EXTERNAL_EOS_HOME)/include
  EXTERNAL_LIB += -L$(EXTERNAL_EOS_HOME)/lib/x86_64-linux-gnu -lRePrimAnd
endif

CFLAGS= -std=c++14 -I. $(EXTERNAL_INCLUDE) -I$(MPI_HOME)/include -I$(HDF5_HOME)/include -I$(SAMRAI_HOME)/include -I$(LORENE_HOME)/C++/Include -I$(LORENE_HOME)/Export/C++/Include
LIBS= $(EXTERNAL_LIB) -L$(LORENE_HOME)/Lib -L$(MPI_HOME)/lib -L$(HDF5_HOME)/lib -L$(SAMRAI_HOME)/lib -lSAMRAI_algs -lSAMRAI_appu -lSAMRAI_geom -lSAMRAI_hier -lSAMRAI_math -lSAMRAI_mesh -lSAMRAI_pdat -lSAMRAI_solv -lSAMRAI_tbox -lSAMRAI_xfer -lmpi -lmpi_cxx -ldl -lpthread -lutil -lhdf5 -lgsl -lgslcblas -lX11 -lfftw3 -lstdc++ -lm

DEPS =Commons.h ExternalInitialData.h TimeInterpolator.h LagrangianPolynomicRefine.h IntegrateDataWriter.h PointDataWriter.h RefineSchedule.h RefineAlgorithm.h StandardRefineTransactionFactory.h RefineClasses.h RefineCopyTransaction.h RefineTimeTransaction.h TimeRefinementIntegrator.h SlicerDataWriter.h SphereDataWriter.h Problem.h Functions.h MainRestartData.h

OBJ =Commons.o ExternalInitialData.o TimeInterpolator.o LagrangianPolynomicRefine.o IntegrateDataWriter.o PointDataWriter.o RefineSchedule.o RefineAlgorithm.o StandardRefineTransactionFactory.o RefineClasses.o RefineCopyTransaction.o RefineTimeTransaction.o TimeRefinementIntegrator.o SlicerDataWriter.o SphereDataWriter.o SAMRAIConnector.o Problem.o Functions.o MainRestartData.o

BrillWaves: $(OBJ)
	$(CXX) $(PERFFLAGS) -o $@ $^ $(CFLAGS) $(LIBS)

ExternalInitialData.o: ExternalInitialData.cpp ExternalInitialData.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

TimeInterpolator.o: TimeInterpolator.C TimeInterpolator.h TimeInterpolateOperator.h Commons.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

LagrangianPolynomicRefine.o: LagrangianPolynomicRefine.C LagrangianPolynomicRefine.h RefineTimeTransaction.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

IntegrateDataWriter.o: IntegrateDataWriter.C IntegrateDataWriter.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

PointDataWriter.o: PointDataWriter.C PointDataWriter.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

RefineSchedule.o: RefineSchedule.C RefineSchedule.h RefineClasses.h TimeInterpolateOperator.h RefineTimeTransaction.h RefineCopyTransaction.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

RefineAlgorithm.o: RefineAlgorithm.C RefineAlgorithm.h TimeInterpolateOperator.h RefineClasses.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

StandardRefineTransactionFactory.o: StandardRefineTransactionFactory.C StandardRefineTransactionFactory.h RefineClasses.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

RefineClasses.o: RefineClasses.C RefineClasses.h TimeInterpolateOperator.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

RefineCopyTransaction.o: RefineCopyTransaction.C RefineCopyTransaction.h RefineClasses.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

RefineTimeTransaction.o: RefineTimeTransaction.C RefineTimeTransaction.h RefineClasses.h TimeInterpolateOperator.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

TimeRefinementIntegrator.o: TimeRefinementIntegrator.C TimeRefinementIntegrator.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

SlicerDataWriter.o: SlicerDataWriter.C SlicerDataWriter.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

SphereDataWriter.o: SphereDataWriter.C SphereDataWriter.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

SAMRAIConnector.o: SAMRAIConnector.cpp Problem.h MainRestartData.h TimeRefinementIntegrator.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

Problem.o: Problem.cpp Problem.h Functions.h RefineClasses.h Commons.h TimeInterpolateOperator.h TimeRefinementIntegrator.h RefineSchedule.h RefineAlgorithm.h StandardRefineTransactionFactory.h RefineTimeTransaction.h MainRestartData.h ExternalInitialData.h TimeInterpolator.h LagrangianPolynomicRefine.h IntegrateDataWriter.h PointDataWriter.h SlicerDataWriter.h SphereDataWriter.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

Functions.o: Functions.cpp Functions.h Commons.h 
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

Commons.o: Commons.cpp Commons.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

MainRestartData.o: MainRestartData.cpp MainRestartData.h
	$(CXX) $(PERFFLAGS) -c -o $@ $< $(CFLAGS)

clean:
	/bin/rm *.o BrillWaves
