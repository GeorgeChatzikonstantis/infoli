#!/bin/bash
CSHELL=sh

# --- MACROS
# define program name
EXE= infoli.x
EXEC= stable_version.x

# define C source files
SRC= infoli.c
SRCC= stable_version.c

# define C header files
HDR= infoli.h

# define C object files
OBJ= infoli.o
OBJC= stable_version.o

# common directives
COMP_FLAGS= -w -c -O3 -vec-report6 -opt-subscript-in-range -opt-report-phase ipo_inl -lm
LINK_FLAGS= -O3 -lm

# mic directives
MIC= -mmic

# version directives
OMP_FLAGS= -openmp
HYBRID_FLAGS= -mt_mpi -openmp

# subdirectories
MPI_DIR= MPI/
OMP_DIR= OpenMP/
HYBRID_DIR= Hybrid_New/
VER_DIR= VersionVerification/

# --- TARGETS
all_xeon: mpi_xeon omp_xeon hybrid_xeon mpi_xeon_old

all_phi: mpi_phi omp_phi hybrid_phi

mpi_xeon: ${MPI_DIR}${SRC} ${MPI_DIR}${HDR}
	@echo #
	@echo "Building For Pure MPI Version on Xeon"
	mpiicc ${MPI_DIR}${SRC} -o ${MPI_DIR}${OBJ} ${COMP_FLAGS}
	mpiicc ${MPI_DIR}${OBJ} -o ${MPI_DIR}${EXE} ${LINK_FLAGS}

mpi_phi: ${MPI_DIR}${SRC} ${MPI_DIR}${HDR}
	@echo #
	@echo "Building For Pure MPI Version on Phi"
	mpiicc ${MPI_DIR}${SRC} -o ${MPI_DIR}${OBJ} ${MIC} ${COMP_FLAGS}
	mpiicc ${MPI_DIR}${OBJ} -o ${MPI_DIR}${EXE} ${MIC} ${LINK_FLAGS}

omp_xeon: ${OMP_DIR}${SRC} ${OMP_DIR}${HDR}
	@echo #
	@echo "Building For Pure OpenMP Version on Xeon"
	icc ${OMP_DIR}${SRC} -o ${OMP_DIR}${OBJ} ${OMP_FLAGS} ${COMP_FLAGS}
	icc ${OMP_DIR}${OBJ} -o ${OMP_DIR}${EXE} ${OMP_FLAGS} ${LINK_FLAGS}

omp_phi: ${OMP_DIR}${SRC} ${OMP_DIR}${HDR}
	@echo #
	@echo "Building For Pure OpenMP Version on Phi"
	icc ${OMP_DIR}${SRC} -o ${OMP_DIR}${OBJ} ${MIC} ${OMP_FLAGS} ${COMP_FLAGS}
	icc ${OMP_DIR}${OBJ} -o ${OMP_DIR}${EXE} ${MIC} ${OMP_FLAGS} ${LINK_FLAGS}

hybrid_xeon: ${HYBRID_DIR}${SRC} ${HYBRID_DIR}${HDR}
	@echo #
	@echo "Building For Hybrid Version on Xeon"
	mpiicc ${HYBRID_DIR}${SRC} -o ${HYBRID_DIR}${OBJ} ${HYBRID_FLAGS} ${COMP_FLAGS}
	mpiicc ${HYBRID_DIR}${OBJ} -o ${HYBRID_DIR}${EXE} ${HYBRID_FLAGS} ${LINK_FLAGS}

hybrid_phi: ${HYBRID_DIR}${SRC} ${HYBRID_DIR}${HDR}
	@echo #
	@echo "Building For Hybrid Version on Phi"
	mpiicc ${HYBRID_DIR}${SRC} -o ${HYBRID_DIR}${OBJ} ${MIC} ${HYBRID_FLAGS} ${COMP_FLAGS}
	mpiicc ${HYBRID_DIR}${OBJ} -o ${HYBRID_DIR}${EXE} ${MIC} ${HYBRID_FLAGS} ${LINK_FLAGS}

test_mpi: ${VER_DIR}${SRCC} ${VER_DIR}${HDR}
	@echo #
	@echo "Building For Stable Version on Xeon"
	icc ${VER_DIR}${SRCC} -o ${VER_DIR}${OBJC} ${COMP_FLAGS}
	icc ${VER_DIR}${OBJC} -o ${VER_DIR}${EXEC} ${LINK_FLAGS}
	./VersionVerification/check.sh mpi_xeon

test_omp: ${VER_DIR}${SRCC} ${VER_DIR}${HDR}
	@echo #
	@echo "Building For Stable Version on Xeon"
	icc ${VER_DIR}${SRCC} -o ${VER_DIR}${OBJC} ${COMP_FLAGS}
	icc ${VER_DIR}${OBJC} -o ${VER_DIR}${EXEC} ${LINK_FLAGS}
	./VersionVerification/check.sh omp_xeon

test_hybrid: ${VER_DIR}${SRCC} ${VER_DIR}${HDR}
	@echo #
	@echo "Building For Stable Version on Xeon"
	icc ${VER_DIR}${SRCC} -o ${VER_DIR}${OBJC} ${COMP_FLAGS}
	icc ${VER_DIR}${OBJC} -o ${VER_DIR}${EXEC} ${LINK_FLAGS}
	./VersionVerification/check.sh hybrid_xeon
clean:
	@echo #
	@echo "-- Cleaning Src --"
	rm -f ${MPI_DIR}*.o ${MPI_DIR}*.x ${OMP_DIR}*.o ${OMP_DIR}*.x ${HYBRID_DIR}*.o ${HYBRID_DIR}*.x
	rm -f ${MPI_DIR}*.optrpt ${OMP_DIR}*.optrpt ${HYBRID_DIR}*.optrpt
