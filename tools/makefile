CSHELL=sh

# --- MACROS
#  # define program name
CONN_EXE= conn_generator.x
CURR_EXE= curr_generator.x

# define C source files
CONN_SRC= conn_generator.c
CURR_SRC= curr_generator.c

# common directives
FLAGS= -w -O3 -lm

# subdirectories
CONN_DIR= connectivity/
GAUSS3D_DIR= 3d_gaussian/
GAUSS3DCOMP_DIR= 3d_gaussian_compressed/
GAUSS2D_DIR= 2d_gaussian/
NEIGH2D_DIR= 2d_neighbors/
PCT_DIR= percentage/
COUNT_DIR= 3d_synapse_count/

CURR_DIR= iapp/
RANDOM_DIR= random/
SINGLE_DIR= single_global_spike/
LOBES_DIR= two_lobes/
WEIBULL_DIR= weibull/

# --- TARGETS
all:	all_conn all_curr

all_conn:	gauss3d gauss3dcomp gauss2d neigh2d pct count

all_curr:	random single lobes weibull

gauss3d: ${CONN_DIR}${GAUSS3D_DIR}${CONN_SRC}
	icc ${CONN_DIR}${GAUSS3D_DIR}${CONN_SRC} -o ${CONN_DIR}${GAUSS3D_DIR}${CONN_EXE} ${FLAGS}

gauss3dcomp: ${CONN_DIR}${GAUSS3DCOMP_DIR}${CONN_SRC}
	icc ${CONN_DIR}${GAUSS3DCOMP_DIR}${CONN_SRC} -o ${CONN_DIR}${GAUSS3DCOMP_DIR}${CONN_EXE} ${FLAGS}

gauss2d: ${CONN_DIR}${GAUSS2D_DIR}${CONN_SRC}
	icc ${CONN_DIR}${GAUSS2D_DIR}${CONN_SRC} -o ${CONN_DIR}${GAUSS2D_DIR}${CONN_EXE} ${FLAGS}

neigh2d: ${CONN_DIR}${NEIGH2D_DIR}${CONN_SRC}
	icc ${CONN_DIR}${NEIGH2D_DIR}${CONN_SRC} -o ${CONN_DIR}${NEIGH2D_DIR}${CONN_EXE} ${FLAGS}

pct: ${CONN_DIR}${PCT_DIR}${CONN_SRC}
	icc ${CONN_DIR}${PCT_DIR}${CONN_SRC} -o ${CONN_DIR}${PCT_DIR}${CONN_EXE} ${FLAGS}

count:	${CONN_DIR}${COUNT_DIR}${CONN_SRC}
	icc ${CONN_DIR}${COUNT_DIR}${CONN_SRC} -o ${CONN_DIR}${COUNT_DIR}${CONN_EXE} ${FLAGS}

random: ${CURR_DIR}${RANDOM_DIR}${CURR_SRC}
	icc ${CURR_DIR}${RANDOM_DIR}${CURR_SRC} -o ${CURR_DIR}${RANDOM_DIR}${CURR_EXE} ${FLAGS}

single: ${CURR_DIR}${SINGLE_DIR}${CURR_SRC}
	icc ${CURR_DIR}${SINGLE_DIR}${CURR_SRC} -o ${CURR_DIR}${SINGLE_DIR}${CURR_EXE} ${FLAGS}

lobes: ${CURR_DIR}${LOBES_DIR}${CURR_SRC}
	icc ${CURR_DIR}${LOBES_DIR}${CURR_SRC} -o ${CURR_DIR}${LOBES_DIR}${CURR_EXE} ${FLAGS}

weibull: ${CURR_DIR}${WEIBULL_DIR}${CURR_SRC}
	icc ${CURR_DIR}${WEIBULL_DIR}${CURR_SRC} -o ${CURR_DIR}${WEIBULL_DIR}${CURR_EXE} ${FLAGS}

clean:
	@echo #
	@echo "-- Cleaning Tools --"
	rm -f ${CONN_DIR}${GAUSS3D_DIR}*.x ${CONN_DIR}${GAUSS2D_DIR}*.x ${CONN_DIR}${NEIGH2D_DIR}*.x ${CONN_DIR}${PCT_DIR}*.x
	rm -f ${CONN_DIR}${GAUSS3DCOMP_DIR}*.x ${CONN_DIR}${COUNT_DIR}*.x
	rm -f ${CURR_DIR}${RANDOM_DIR}*.x ${CURR_DIR}${SINGLE_DIR}*.x ${CURR_DIR}${LOBES_DIR}*.x ${CURR_DIR}${WEIBULL_DIR}*.x
