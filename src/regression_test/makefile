#!/bin/bash
CSHELL=sh

# --- MACROS
# define program name
EXE= infoli.x
EXE2= checker.x

# define C source files
SRC= stable_version.c
SRC2= regression_checker.c

# define C header files
HDR= infoli.h

# define C object files
OBJ= infoli.o
OBJ2=regression_checker.o

# common directives
COMP_FLAGS= -w -c -O3 -lm
LINK_FLAGS= -O3 -lm

# regression check directive
REG_TEST_FLAGS= -DG_CAL_FROM_FILE

# --- TARGETS

default: ${SRC} ${HDR}
	@echo #
	@echo "Building Regression Test Executable and Checker"
	icc ${SRC} -o ${OBJ} ${REG_TEST_FLAGS} ${COMP_FLAGS}
	icc ${OBJ} -o ${EXE} ${REG_TEST_FLAGS} ${LINK_FLAGS}
	icc ${SRC2} -o ${OBJ2} ${COMP_FLAGS}
	icc ${OBJ2} -o ${EXE2} ${LINK_FLAGS}

clean:
	@echo #
	@echo "-- Cleaning Regression Test Dir --"
	rm -f ${EXE} ${OBJ} ${EXE2} ${OBJ2}
