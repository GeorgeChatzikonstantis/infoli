#!/bin/bash
CSHELL=sh

# --- MACROS
# define program name
EXE= infoli.x

# define C source files
SRC= infoli.c

# define C header files
HDR= infoli.h

# define C object files
OBJ= infoli.o

# common directives
COMP_FLAGS= -w -c -O3 -vec-report6 -opt-subscript-in-range -opt-report-phase ipo_inl -lm
LINK_FLAGS= -O3 -lm

# mic directives
MIC= -mmic

# define mic power measurement files
MICPWRDIR=../lib/
MICPWRH= ${MICPWRDIR}/power/mic_power.h
MICPWRSO= mic_power
IMFSO= imf
INTLCSO= intlc
IRNGSO= irng
SVMLSO= svml

# version directives
OMP_FLAGS= -openmp
HYBRID_FLAGS= -mt_mpi -openmp

# subdirectories
OMP_DIR= openmp/
HYBRID_DIR= hybrid/
REG_DIR= regression_test/

# bin directory
BIN_DIR= ../bin/

# --- TARGETS
all: hybrid_xeon hybrid_phi omp_xeon omp_phi

omp_xeon: ${OMP_DIR}${SRC} ${OMP_DIR}${HDR}
	@echo #
	@echo "Building For Pure OpenMP Version on Xeon"
	icc ${OMP_DIR}${SRC} -o ${OMP_DIR}${OBJ} ${OMP_FLAGS} ${COMP_FLAGS}
	icc ${OMP_DIR}${OBJ} -o ${BIN_DIR}${EXE} ${OMP_FLAGS} ${LINK_FLAGS}

omp_phi: ${OMP_DIR}${SRC} ${OMP_DIR}${HDR}
	@echo #
	@echo "Building For Pure OpenMP Version on Phi"
	icc ${OMP_DIR}${SRC} -o ${OMP_DIR}${OBJ} ${MIC} ${OMP_FLAGS} ${COMP_FLAGS} -I${MICPWRDIR}
	icc ${OMP_DIR}${OBJ} -o ${OMP_DIR}${EXE} ${MIC} ${OMP_FLAGS} ${LINK_FLAGS} -L${MICPWRDIR} \
	-l${MICPWRSO} -l${IMFSO} -l${INTLCSO} -l${IRNGSO} -l${SVMLSO}

hybrid_xeon: ${HYBRID_DIR}${SRC} ${HYBRID_DIR}${HDR}
	@echo #
	@echo "Building For Hybrid Version on Xeon"
	mpiicc ${HYBRID_DIR}${SRC} -o ${HYBRID_DIR}${OBJ} ${HYBRID_FLAGS} ${COMP_FLAGS}
	mpiicc ${HYBRID_DIR}${OBJ} -o ${BIN_DIR}${EXE} ${HYBRID_FLAGS} ${LINK_FLAGS}

hybrid_phi: ${HYBRID_DIR}${SRC} ${HYBRID_DIR}${HDR}
	@echo #
	@echo "Building For Hybrid Version on Phi"
	mpiicc ${HYBRID_DIR}${SRC} -o ${HYBRID_DIR}${OBJ} ${MIC} ${HYBRID_FLAGS} ${COMP_FLAGS}
	mpiicc ${HYBRID_DIR}${OBJ} -o ${BIN_DIR}${EXE} ${MIC} ${HYBRID_FLAGS} ${LINK_FLAGS}

clean:
	@echo #
	@echo "-- Cleaning Src --"
	rm -f ${OMP_DIR}*.o ${OMP_DIR}*.x ${HYBRID_DIR}*.o ${HYBRID_DIR}*.x
	rm -f ${OMP_DIR}*.optrpt ${HYBRID_DIR}*.optrpt
