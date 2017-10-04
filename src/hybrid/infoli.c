/*
 *
 * Copyright (c) 2012, Neurasmus B.V., The Netherlands,
 * web: www.neurasmus.com email: voltage@neurasmus.com
 *
 * Any use or reproduction in whole or in parts is prohibited
 * without the written consent of the copyright owner.
 *
 * All Rights Reserved.
 *
 *
 * Author: Sebastian Isaza
 * Created: 19-01-2012
 * Modified: 07-08-2012
 *
 * Description: Top source file of the Inferior Olive model, originally written
 * in Matlab by Jornt De Gruijl. It contains the implementation of all functions.
 * The main function allocates the necessary memory, initializes the system
 * state and runs the model calculations.
 *
 */

/*
 * we assume that dim1 refers to the number of ROWS
 * and dim2 refers to COLUMNS (cell network size).
 */


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include "infoli.h"

int core_id, cores, cellCount, core_offset;
int IO_NETWORK_DIM1, IO_NETWORK_DIM2, IO_NETWORK_SIZE;
float CONN_PROBABILITY;
struct timeval tic, toc, intime;

void subtract_and_add (struct timeval *Z,struct timeval*x,struct timeval *y) {

	struct timeval result;
      
       /* Compute the time remaining to wait.
          tv_usec is certainly positive. */
	result.tv_sec = y->tv_sec - x->tv_sec;
	result.tv_usec = y->tv_usec - x->tv_usec;
	   
	Z->tv_sec += result.tv_sec;
	Z->tv_usec += result.tv_usec;
	    
	return;
}

char fpeek(FILE *stream)
{
	char c;

	c = (char)fgetc(stream);
	ungetc(((int)c), stream);

	return c;
}

int main(int argc, char *argv[]){
	
	/* declaration of variables
	 * necessary for program flow
	 */

	int i=0, j, k=0, line_count, targetCore, global_cell_id, print_granularity=1;
	char c;
	char *inFileName;
	char outFileName[100];
	FILE *pInFile, *coreF;
	char conFile[200],core_file[100];
	FILE *pConFile,*pOutFile;
	mod_prec* iAppArray= NULL;
	int simStep = 0, inputFromFile = 0, initSteps,total_simulation_steps;
	float simTime = 0;
	int seedvar;
	char tempbuf[100];		
	mod_prec iApp, voltage;
	mod_prec *gcal;
	int *init_steps;
	long double seconds;
	int simulation_array_ID,target_cell, local_cell, requested_neighbour;
	int package_index, coded_package_index, decoded_package_index;

	cellCompParams cellParamsPtr;

	/* declaration of variables
	 * necessary for ion channel
	 * computations
	 */

	mod_prec q_inf, tau_q, dq_dt;
	mod_prec alpha_r, beta_r, r_inf, tau_r, dr_dt;
	mod_prec alpha_s, beta_s, s_inf, tau_s, ds_dt;
	mod_prec dCa_dt, f, V, I_c_storage;
	mod_prec I_sd, I_CaH_temp, I_K_Ca, I_ld, I_h, dVd_dt;

	mod_prec k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt;
	mod_prec m_inf, h_inf, tau_h, dh_dt;
	mod_prec n_inf, p_inf, tau_n, tau_p, dn_dt, dp_dt;
	mod_prec alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s;
	mod_prec I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt;

	mod_prec m_inf_a, h_inf_a, tau_h_a, dh_dt_a;
	mod_prec alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a;
	mod_prec I_Na_a, I_la, I_sa, I_K_a, dVa_dt;

	/* end of declarations
	 */ 

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &core_id);
	MPI_Comm_size(MPI_COMM_WORLD, &cores);
	MPI_Request* s_request = (MPI_Request*) malloc(cores*sizeof(MPI_Request));      //array of request handles for MPI_Isend, one for each core in the network
	MPI_Request* r_request = (MPI_Request*) malloc(cores*sizeof(MPI_Request));      //array of request handles for MPI_Irecv, one for each core in the network

	/* Process command line arguments
	 * Case argc = 3 then a one-pulse input is stimulated.
	 * Otherwise we receive input from a specified file in argv[3] and simulation runs accordingly
	 * in the case of inputFromFile, we will also need a buffer to hold the stimulus
	 * for each cell on each step
  	 */

	char host_name[100];
	gethostname(host_name, 100);
	printf("Rank %d: Running on %s.\n", core_id, host_name);
	
	if(argc == 4) {
		inputFromFile = 0;
	}
	else if(argc == 5) {
		inputFromFile = 1;
		inFileName = argv[2];
		pInFile = fopen(inFileName,"r");
		if(pInFile==NULL) {
			printf("Error: Couldn't open %s\n", inFileName);
			exit(EXIT_FAILURE);
		}
		iAppArray= (mod_prec*) _mm_malloc(cellCount*(sizeof(mod_prec)), 64);
	}	
	else {
		printf("Error: Too many arguments.\nUsage: ./InferiorOlive <#nrns> <density ([0,1])> <simtime>\nor ./InferiorOlive <#nrns> <density ([0,1])> <simtime> <stimulus_file>\n");
		exit(EXIT_FAILURE);
	}

	/* outFileName is file for output
	 */
        sprintf(outFileName,"InferiorOlive_Output%d.txt", core_id);
        
        IO_NETWORK_SIZE = atoi(argv[1]);
        CONN_PROBABILITY = atof(argv[2]);
        simTime = atof(argv[3]);

        /* compute how many grid cells 
	 * are assigned to each core
	*/

        cellCount= IO_NETWORK_SIZE/cores;

	/* compute offset each core applies
	 * to its assigned (local) cells
	 * in order to calculate a cell id
	 * relevant to the entire (global) network
	 */

	core_offset = core_id*cellCount;

	/* PRINTING is true in case we want to write the
	 * output to a specified file (outFileName)
	 */

	if (PRINTING) {
		pOutFile = fopen(outFileName,"w+");
		if(pOutFile==NULL){
			printf("Error: Couldn't create %s\n", outFileName);
			exit(EXIT_FAILURE);
		}
		print_granularity=100;
		#ifdef G_CAL_FROM_FILE
			print_granularity=1;
		#endif
	}

	/* Channel declaration and mem allocation 
	 * 
	 */

	int* cellID = (int*) _mm_malloc(cellCount*sizeof(int), 64);
	mod_prec* V_dend = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Hcurrent_q = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Calcium_r = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_s = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* I_CaH = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Ca2Plus = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* iAppIn = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* I_c = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* V_soma = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* g_CaL = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Sodium_m = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Sodium_h = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Calcium_k = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Calcium_l = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_n = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_p = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_x_s = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* V_axon = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Sodium_m_a = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Sodium_h_a = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
	mod_prec* Potassium_x_a = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);

	if (!(Potassium_x_a)) {
		printf("Error: Couldn't allocate memory for channels\n");
		exit(EXIT_FAILURE);
	}

	/* Initialise network with appropriate values.
	*/

	for (i=0;i<cellCount;i++) {
		cellID[i] = i;
		//Initial dendritic parameters
		V_dend[i] = -60;
		Calcium_r[i] = 0.0112788;// High-threshold calcium
		Potassium_s[i] = 0.0049291;// Calcium-dependent potassium
		Hcurrent_q[i] = 0.0337836;// H current
		Ca2Plus[i] = 3.7152;// Calcium concentration
		I_CaH[i] = 0.5;// High-threshold calcium current
		//Initial somatic parameters
		g_CaL[i] = 0.68; //default arbitrary value but it should be randomized per cell
		V_soma[i] = -60;
		Sodium_m[i] = 1.0127807;// Sodium (artificial)
		Sodium_h[i] = 0.3596066;
		Potassium_n[i] = 0.2369847;// Potassium (delayed rectifier)
		Potassium_p[i] = 0.2369847;
		Potassium_x_s[i] = 0.1;// Potassium (voltage-dependent)
		Calcium_k[i] = 0.7423159;// Low-threshold calcium
		Calcium_l[i] = 0.0321349;
		// Initial axonal parameters
		V_axon[i] = -60;
		//sisaza: Sodium_m_a doesn't have a state, therefore this assignment doesn'thave any effect
		Sodium_m_a[i] = 0.003596066;// Sodium (thalamocortical)
		Sodium_h_a[i] = 0.9;
		Potassium_x_a[i] = 0.2369847;// Potassium (transient)
	}

	#ifdef G_CAL_FROM_FILE
		read_g_CaL_from_file(g_CaL);
	#endif

	/* cellCompParams contains
	 * connection information for the NW
	 */


	//1D-array storing how many connections each neuron has
	cellParamsPtr.total_amount_of_neighbours = (int*) _mm_malloc(cellCount*sizeof(int), 64);
	for (i=0;i<cellCount;i++)
		cellParamsPtr.total_amount_of_neighbours[i] = 0;	//needed bugfix initialization
	//2D-array storing the ids of each connection
	cellParamsPtr.neighId = (int**) _mm_malloc(cellCount*sizeof(int*), 64);
	//2D-array storing the conductances of each connection
	cellParamsPtr.neighConductances = (mod_prec**) _mm_malloc(cellCount*sizeof(mod_prec*), 64);

	if(cellParamsPtr.neighConductances==NULL){
		printf("Error: Couldn't malloc for cellParamsPtr\n");
		exit(EXIT_FAILURE);
	}

	mod_prec cond_value= 0.04;	//default value for conductance in synapses (microSiemens)

	/* we shall generate a connectivity map
	 * based on the probability given by the user
	 * Note that this map handles the core's cells
	 * connections to the entire network
	 */

	int sender_cell, receiver_cell;
	int* conn_gen_buffer=(int*) calloc(IO_NETWORK_SIZE, sizeof(int));	//temp buffer for storing bonds
	float rndm;

	/* we generate rng checks against the probability
	 * supplied by the user to calculate connection counts
 	 * We store the results of each rng check in a temp buffer
 	 * After we are done, we allocate memory for dendritic data structs
 	 * and fill it with neighbour (bonds) ids
 	 */

	srand ( time(NULL) + core_id );
	short* cellsNeeded = (short*) calloc(IO_NETWORK_SIZE, sizeof(short));	//buffer marking cells in other cores necessary to this core

	for (receiver_cell=0; receiver_cell<cellCount; receiver_cell++) {
		global_cell_id = receiver_cell + core_offset;
		for (sender_cell=0;sender_cell<IO_NETWORK_SIZE;sender_cell++) {

			if (sender_cell == global_cell_id)
				;		//no self-feeding connections allowed
			else {
				rndm = ((float) rand()) / ((float) RAND_MAX);	//generate rng and compare to probability
				if (rndm <= CONN_PROBABILITY) {
					cellParamsPtr.total_amount_of_neighbours[receiver_cell]++;  //increase neighbour count
					conn_gen_buffer[sender_cell]++;	//mark that we formed a bond with this cell
					if ((sender_cell/cellCount)!=core_id) //if this cell does not belong to core
						cellsNeeded[sender_cell]  = 1; //mark it
				}
			}
		}

		//allocate enough space now that we know how many neighbours this receiving cell has
		cellParamsPtr.neighConductances[receiver_cell] = (mod_prec*) _mm_malloc(cellParamsPtr.total_amount_of_neighbours[receiver_cell]*sizeof(mod_prec), 64);
		cellParamsPtr.neighId[receiver_cell] = (int*) _mm_malloc(cellParamsPtr.total_amount_of_neighbours[receiver_cell]*sizeof(int), 64);

		//run through the temporary buffer and fill the data structs with the bonds' info
		i=0;	//this temporary variable will help fill the data structs
		for (sender_cell=0;sender_cell<IO_NETWORK_SIZE;sender_cell++) {
			if (conn_gen_buffer[sender_cell]==0)
				;	//skip this cell, it is not a bond
			else {
				cellParamsPtr.neighConductances[receiver_cell][i]=cond_value;
				cellParamsPtr.neighId[receiver_cell][i]=sender_cell;
				i++;
			}
			if (i>cellParamsPtr.total_amount_of_neighbours[receiver_cell])
				break;
		}

		//reset the buffer for the next receiver cell
		memset(conn_gen_buffer, 0, IO_NETWORK_SIZE*sizeof(int));
	}

	printf("%d: ", core_id);
/*	for (receiver_cell=0; receiver_cell<cellCount; receiver_cell++) {
		printf("| %d | ", receiver_cell);
		for (i=0; i<cellParamsPtr.total_amount_of_neighbours[receiver_cell]; i++)
			printf("%d ", cellParamsPtr.neighId[receiver_cell][i]);
	}
*/	printf("\n");

	/* connections mapping
	 * for the core's cells
	 * to the entire network
 	 * is now complete
 	 */

	/* cores will exchange knowledge
	 * of each sub-network's connections
	 * in order to sync the sub-network's
	 * connectivity needs
	 *
 	 * scatter the cellsNeeded matrix
	 * so that other cores know which cells are needed
	 */

	short** cellsToSend = (short**) malloc(cores*sizeof(short*));
	for (i=0; i<cores; i++) {
		cellsToSend[i] = (short*) malloc(cellCount*sizeof(short));
		MPI_Scatter(cellsNeeded, cellCount, MPI_SHORT, cellsToSend[i], cellCount, MPI_SHORT, i, MPI_COMM_WORLD);
	}

	/* process received matrix so that the number of cells
	 * necessary for exchanging per core is known
	 */

	int* packagesToSend = (int*) calloc(cores, sizeof(int));
	for (i=0; i<cores; i++) {
		for (j=0; j<cellCount; j++) {
			if (cellsToSend[i][j] > 0)
				packagesToSend[i]++;
		}
	}
	packagesToSend[core_id]=0;	//no need for packages to same-core
	int* packagesToReceive = (int*) calloc(cores, sizeof(int));
	for (i=0; i<cores; i++) {
		MPI_Scatter(packagesToSend, 1, MPI_INT, &packagesToReceive[i], 1, MPI_INT, i, MPI_COMM_WORLD);
	}

	/* notify the other cores the (global) indexes
	 * of each cells that they will be receiving
	 */
	int** packagesIndexToSend = (int**) malloc(cores*sizeof(int*));
	for (i=0; i<cores; i++) {
		packagesIndexToSend[i] = (int*) calloc(packagesToSend[i], sizeof(int));
		k=0;
		for (j=0; j<cellCount; j++) {
			global_cell_id = j + core_offset;
			if (cellsToSend[i][j] > 0) {
				packagesIndexToSend[i][k] = global_cell_id;
				k++;
			}
			if (k>(packagesToSend-1))
				break;
		}
	}
	int** packagesIndexToReceive = (int**) malloc(cores*sizeof(int*));
	for (i=0; i<cores; i++) {
		packagesIndexToReceive[i] = (int*) calloc(packagesToReceive[i], sizeof(int));
		if (i!=core_id) {
			MPI_Isend(packagesIndexToSend[i], packagesToSend[i], MPI_INT, i, 0, MPI_COMM_WORLD, &s_request[i]);
			MPI_Irecv(packagesIndexToReceive[i], packagesToReceive[i], MPI_INT, i, 0, MPI_COMM_WORLD, &r_request[i]);
		} else {
			r_request[i]= MPI_REQUEST_NULL;
		}
	}

	/* create buffers that will send and hold
	 * the respective voltages in each sim step
	 */

	mod_prec** packagesVoltToSend = (mod_prec**) _mm_malloc(cores*sizeof(mod_prec*), 64);
	for (i=0; i<cores; i++)
		packagesVoltToSend[i] = (mod_prec*) _mm_malloc(packagesToSend[i]*sizeof(mod_prec), 64);
	mod_prec** packagesVoltToReceive = (mod_prec**) _mm_malloc(cores*sizeof(mod_prec*), 64);
	for (i=0; i<cores; i++) {
		if (i==core_id)
			packagesVoltToReceive[i] = (mod_prec*) _mm_malloc(cellCount*sizeof(mod_prec), 64);
		else
			packagesVoltToReceive[i] = (mod_prec*) _mm_malloc(packagesToReceive[i]*sizeof(mod_prec), 64);
	}

	/* sync asynchronous calls to ensure
	 * core syncing
	 */
	syncing(r_request);
	cleanup_requests(r_request);

/*	printf("%d: ", core_id);
        for (i=0; i<cores; i++) {
                printf("| %d <-  | ", i);
                for (j=0; j<packagesToReceive[i]; j++)
                        printf("%d ", packagesIndexToReceive[i][j]);
        }
        printf("\n");
*/

	/* free up memory structs not necessary
	 * anymore (simulation body only needs
	 * the packages data structs)
	 */
	for (i=0;i<cores;i++)
		free(cellsToSend[i]);
	free(cellsToSend);
	free(cellsNeeded);

	/* create a secondary struct
	 * which allows the core's cells
	 * to know the location (index) of their
	 * connections' data in each of the
	 * received packages from other cores.
	 * This index matching allows
	 * skipping parses of the received 
	 * packages, accelerating post-communication
	 * processing but increasing memory footprint
	 * The index will be "codified" in a way that allows
	 * its value to indicate both  which package to look
	 * for (i.e. which core it came from) and the
	 * exact position within the package (i.e. index)
	 */

	int** packagesIndexMatcher = (int**) _mm_malloc(cellCount*sizeof(int*), 64);
	int packageNumber;	//which package to match index with, corresponds to which core sent the package
	int packageIndex;	//which index of the package we are currently examinating
	int neighbour;		//which connection of the receiver cell we are trying to match
	for (receiver_cell=0; receiver_cell<cellCount; receiver_cell++) { 
		packagesIndexMatcher[receiver_cell] = (int*) _mm_malloc(cellParamsPtr.total_amount_of_neighbours[receiver_cell]*sizeof(int), 64);
		packageNumber=0;
		packageIndex=0;
		neighbour=0;

		while ((packageNumber<cores)&&(neighbour<cellParamsPtr.total_amount_of_neighbours[receiver_cell])) {

			packagesIndexMatcher[receiver_cell][neighbour]=-1;

			if (packageNumber==core_id) {
				packageNumber++;
				if (!(packageNumber<cores))
					break;
				packageIndex=0;
			}
			
			if ((cellParamsPtr.neighId[receiver_cell][neighbour]/cellCount)==core_id) {
				local_cell=cellParamsPtr.neighId[receiver_cell][neighbour];
				packagesIndexMatcher[receiver_cell][neighbour]=local_cell;
				neighbour++;
			} else {

				if (cellParamsPtr.neighId[receiver_cell][neighbour]==packagesIndexToReceive[packageNumber][packageIndex]) {
					packagesIndexMatcher[receiver_cell][neighbour]=packageIndex+packageNumber*cellCount;
					neighbour++;
				}

				packageIndex++;
				if (packageIndex>(packagesToReceive[packageNumber]-1)) {
					packageNumber++;
					packageIndex=0;
				}
			}

		}

		while (neighbour<cellParamsPtr.total_amount_of_neighbours[receiver_cell]) {

			packagesIndexMatcher[receiver_cell][neighbour]=-1;
			if ((cellParamsPtr.neighId[receiver_cell][neighbour]/cellCount)==core_id) {
				local_cell=cellParamsPtr.neighId[receiver_cell][neighbour];
				packagesIndexMatcher[receiver_cell][neighbour]=local_cell;
			}
			neighbour++;
		}

	}


	/* cores are now up-to-date with
	 * each other's communicational needs.
	 * Communication for main simulation
	 * body is properly set up
	 */

	//Initialize output file IF enabled
	#ifndef G_CAL_FROM_FILE
		if (PRINTING) {
			sprintf(tempbuf, "Execution Time for Simulation in ms:              \n");
				fputs(tempbuf, pOutFile);
			sprintf(tempbuf, "#simSteps Input(Iapp) Output(V_axon)");
			fputs(tempbuf, pOutFile);
			for (i=0;i<cellCount;i++) {
				sprintf(tempbuf, "%d ", cellID[i]);
				fputs(tempbuf, pOutFile);
			}
			sprintf(tempbuf, "\n");
			fputs(tempbuf, pOutFile);
		}
	#endif

	//Initialize g_CaL
	seedvar = time(NULL);		//seedvar will be time- and core-dependent now!
	srand(seedvar);

	if (RAND_INIT)
		for(i=0;i<cellCount;i++)
			g_CaL[i] = 0.6+(0.2*(rand()%100)/100);

	//random initialization process
	if (RAND_INIT) {
		for(i=0;i<cellCount;i++) {
			initSteps = rand()%(int)ceil(100/DELTA);
			initSteps = initSteps | 0x00000001;//make initialization steps odd

			for(j=0;j<initSteps;j++){
				//Arrange input
				iAppIn[i] = 0;//No stimulus
				
//				CompDend(&cellParamsPtr[i], 1);  PLACEHOLDER RANDOMIZATION, NEEDS RE-CODING
//				CompSoma(&cellParamsPtr[i]);
//				CompAxon(&cellParamsPtr[i]);
			}

		}
	}

	/* start of the simulation
	 * In case we want to read the stimulus from file inputFromFile = true
	 */
	gettimeofday(&tic, NULL);

	/* 	WARNING, recent changes have made inputFromFile currently unusable-
 	 *  	related code will be commented out until fixed and commited
 	 */

/*	if(inputFromFile){
		simStep = 0;


		while(ReadFileLine(pInFile, iAppArray)) {
			
			simulation_array_ID = simStep%2;
			if (PRINTING) {
				sprintf(tempbuf, "%d %.2f ", (simStep+1), simStep*0.05); // start @ 1 because skipping initial values
				fputs(tempbuf, pOutFile);
			}
			
		#pragma omp parallel for shared(cellParamsPtr, cellPtr) private(iAppArray, target_cell, tempbuf, i, requested_neighbour) firstprivate(simulation_array_ID)	
			for (target_cell=0;target_cell<cellCount;target_cell++) {
				
				for (i=0; i<cellParamsPtr[target_cell].total_amount_of_neighbours; i++){
                                        requested_neighbour = cellParamsPtr[target_cell].neighId[i];
                                        cellParamsPtr[target_cell].neighVdend[i] = cellPtr[simulation_array_ID][requested_neighbour].dend.V_dend;
				}

				cellParamsPtr[target_cell].iAppIn = iAppArray[target_cell];
				cellParamsPtr[target_cell].prevCellState = &cellPtr[simulation_array_ID][target_cell];
				cellParamsPtr[target_cell].newCellState = &cellPtr[simulation_array_ID^1][target_cell];

				CompDend(&cellParamsPtr[target_cell], 0);
				CompSoma(&cellParamsPtr[target_cell]);
				CompAxon(&cellParamsPtr[target_cell]);

				if (PRINTING) {
					sprintf(tempbuf, "%.16f ", cellPtr[(simulation_array_ID)^1][target_cell].axon.V_axon);
					fputs(tempbuf, pOutFile);
				}

			}
			if (PRINTING) {
				sprintf(tempbuf, "\n");
				fputs(tempbuf,pOutFile);
			}

			simStep++;
		}
		
	}
	else {	
*/		total_simulation_steps = (int)(simTime/DELTA);

		for(simStep=0;simStep<total_simulation_steps;simStep++) {
								
			if ((simStep>=20000)&&(simStep<20500-1))
				iApp = 6; 
			else
		   		iApp = 0;
			
			#ifndef G_CAL_FROM_FILE
				if (PRINTING) {
					sprintf(tempbuf, " %d %.2f ", simStep+1, iApp);
					fputs(tempbuf, pOutFile);
				}
			#endif


	/* 	We shall perform core-communication to exchange voltage packages
	 * 	We first fill the voltage-buffer needed to be sent over
	 */ 

			for (targetCore=0; targetCore<cores; targetCore++) {
				if (targetCore==core_id) {
					//fill the receiving Voltage buffer with our updated local cells
					memcpy(packagesVoltToReceive[core_id], V_dend, cellCount*sizeof(mod_prec));
					r_request[targetCore]=MPI_REQUEST_NULL;
				} else {
					for (package_index=0; package_index<packagesToSend[targetCore]; package_index++) {
						requested_neighbour = packagesIndexToSend[targetCore][package_index]-core_offset;
						packagesVoltToSend[targetCore][package_index] = V_dend[requested_neighbour];
					}

	/*	We then call MPI functions for sending the voltage-buffer and receiving
	 *	our respective voltage-buffer
	 */
					if (packagesToSend[targetCore]>0)
						MPI_Isend(packagesVoltToSend[targetCore], packagesToSend[targetCore], MPI_FLOAT, targetCore, 0, MPI_COMM_WORLD, &s_request[targetCore]);
					if (packagesToReceive[targetCore]>0)
						MPI_Irecv(packagesVoltToReceive[targetCore], packagesToReceive[targetCore], MPI_FLOAT, targetCore, 0, MPI_COMM_WORLD, &r_request[targetCore]);
				}
			}
			syncing(r_request);
			cleanup_requests(r_request);

	/*	First Loop takes care of Gap Junction Functions
 	*/			

		#pragma omp parallel for shared (cellParamsPtr, packagesIndexMatcher, packagesVoltToReceive, iApp, iAppIn, V_dend, I_c, \
		pOutFile) private(target_cell, i, requested_neighbour, targetCore, f, V, coded_package_index, decoded_package_index, \
		voltage, I_c_storage) firstprivate(simulation_array_ID)
			for (target_cell=0;target_cell<cellCount;target_cell++) {

				if ((simStep==0)&&(target_cell==0))
					printf("Rank %d: Executing with %d threads.\n", core_id, omp_get_num_threads());

				//Feeding of Input Current
				iAppIn[target_cell] = iApp;

	/*	Gathering of the data concerning Vdend of neighours, then
 	*	computing and storing the incoming Ic from neighbours
 	*/

				I_c_storage = 0;
				__assume_aligned(cellParamsPtr.neighConductances[target_cell], 64);
				__assume_aligned(V_dend, 64);
				__assume_aligned(packagesIndexMatcher[target_cell], 64);
			#pragma ivdep
				for (i=0; i<cellParamsPtr.total_amount_of_neighbours[target_cell]; i++){
					//get the codified index of the cell's neighbour
					coded_package_index = packagesIndexMatcher[target_cell][i];
					//decode the package (core) you need to read
					targetCore = coded_package_index/cellCount;
					//decode the index of the package-to-read
					decoded_package_index = coded_package_index%cellCount;
					//search the neighbour's voltage in the receiving voltage buffer
					voltage = packagesVoltToReceive[targetCore][decoded_package_index];
					//compute neighbour's incoming current
					V = V_dend[target_cell] - voltage;
					f = 0.8f * expf(-1*powf(V, 2)/100) + 0.2f;// SCHWEIGHOFER 2004 VERSION
					I_c_storage += cellParamsPtr.neighConductances[target_cell][i] * f * V;
                                }	
				I_c[target_cell] = I_c_storage;
			}

	/* Second Loop computes the new state (voltages, channels etc.)
 	*/ 

			__assume_aligned(V_dend, 64);
			__assume_aligned(Hcurrent_q, 64);
			__assume_aligned(Calcium_r, 64);
			__assume_aligned(Potassium_s, 64);
			__assume_aligned(I_CaH, 64);
			__assume_aligned(Ca2Plus, 64);
			__assume_aligned(I_c, 64);
			__assume_aligned(iAppIn, 64);

			__assume_aligned(V_soma, 64);
			__assume_aligned(g_CaL, 64);
			__assume_aligned(Sodium_m, 64);
			__assume_aligned(Sodium_h, 64);
			__assume_aligned(Calcium_k, 64);
			__assume_aligned(Calcium_l, 64);
			__assume_aligned(Potassium_n, 64);
			__assume_aligned(Potassium_p, 64);
			__assume_aligned(Potassium_x_s, 64);

			__assume_aligned(V_axon, 64);
			__assume_aligned(Sodium_m_a, 64);
			__assume_aligned(Sodium_h_a, 64);
			__assume_aligned(Potassium_x_a, 64);
		#pragma omp parallel for simd shared (cellParamsPtr, V_dend, Hcurrent_q, Calcium_r, Potassium_s, I_CaH, Ca2Plus,\
		iApp, iAppIn, I_c, V_soma, g_CaL, Sodium_m, Sodium_h, Calcium_k, Calcium_l, Potassium_n, Potassium_p, Potassium_x_s,\
		V_axon, Sodium_m_a, Sodium_h_a, Potassium_x_a, pOutFile) private(target_cell, tempbuf, q_inf, tau_q, dq_dt,\
		alpha_r, beta_r, r_inf, tau_r, dr_dt, alpha_s, beta_s, s_inf, tau_s, ds_dt, dCa_dt, I_sd, I_CaH_temp, I_K_Ca, I_ld,\
		I_h, dVd_dt, k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt, m_inf, h_inf, tau_h, dh_dt, n_inf, p_inf, tau_n, tau_p, dn_dt,\
		dp_dt, alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s, I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt,\
		m_inf_a, h_inf_a, tau_h_a, dh_dt_a, alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a, I_Na_a, I_la, I_sa, I_K_a,\
		dVa_dt)	firstprivate(simulation_array_ID)
			for (target_cell=0;target_cell<cellCount;target_cell++) {

//				~DENDRITIC COMPUTATIONS~

//				Dend H Current Calcs
				q_inf = 1 /(1 + expf((V_dend[target_cell] + 80) / 4));
			        tau_q = 1 /(expf(-0.086f * V_dend[target_cell] - 14.6f) + expf(0.070f * V_dend[target_cell] - 1.87f));
			        dq_dt = (q_inf - Hcurrent_q[target_cell]) / tau_q;
				Hcurrent_q[target_cell] = DELTA * dq_dt + Hcurrent_q[target_cell];

//				Dend Ca Current Calcs
				alpha_r = 1.7f / (1 + expf( -(V_dend[target_cell] - 5) / 13.9f));
		       		beta_r = 0.02f * (V_dend[target_cell] + 8.5f) / (expf((V_dend[target_cell] + 8.5f) / 5) - 1);
	  		     	r_inf = alpha_r / (alpha_r + beta_r);
	       	 		tau_r = 5 / (alpha_r + beta_r);
	       		 	dr_dt = (r_inf - Calcium_r[target_cell]) / tau_r;
	    			Calcium_r[target_cell] = DELTA * dr_dt + Calcium_r[target_cell];

//				Dend K Current Calcs
				alpha_s = min((0.00002f * Ca2Plus[target_cell]), 0.01f);
			        beta_s = 0.015f;
			        s_inf = alpha_s / (alpha_s + beta_s);
			        tau_s = 1 / (alpha_s + beta_s);
			        ds_dt = (s_inf - Potassium_s[target_cell]) / tau_s;
			        Potassium_s[target_cell] = DELTA * ds_dt + Potassium_s[target_cell];

//				Dend Cal Current Calcs
				dCa_dt = -3 * I_CaH[target_cell] - 0.075f * Ca2Plus[target_cell];
				Ca2Plus[target_cell] = DELTA * dCa_dt + Ca2Plus[target_cell];

//				Dendritic Voltage And Current Calcs
				I_sd = (G_INT / (1 - P1)) * (V_dend[target_cell] - V_soma[target_cell]);
				I_CaH_temp =G_CAH * powf(Calcium_r[target_cell], 2) * (V_dend[target_cell] - V_CA);
				I_K_Ca =G_K_CA * Potassium_s[target_cell] * (V_dend[target_cell] - V_K);
				I_ld =G_LD * (V_dend[target_cell] - V_L);
				I_h=G_H * Hcurrent_q[target_cell] * (V_dend[target_cell] - V_H);
				dVd_dt = (-(I_CaH_temp + I_sd+ I_ld + I_K_Ca + I_c[target_cell] + I_h) + iAppIn[target_cell]) / C_M;
				I_CaH[target_cell] = I_CaH_temp;


//				~SOMATIC COMPUTATIONS~

//				Soma Calcium Calcs
				k_inf = (1 / (1 + expf(-1 * (V_soma[target_cell] + 61) / 4.2f)));
			        l_inf = (1 / (1 + expf(( V_soma[target_cell] + 85.5f) / 8.5f)));
			        tau_k = 1;
			        tau_l = ((20 * expf((V_soma[target_cell] + 160) / 30) / (1 + expf((V_soma[target_cell] + 84) / 7.3f))) +35);
			        dk_dt = (k_inf - Calcium_k[target_cell]) / tau_k;
			        dl_dt = (l_inf - Calcium_l[target_cell]) / tau_l;
			        Calcium_k[target_cell] = DELTA * dk_dt + Calcium_k[target_cell];
			        Calcium_l[target_cell] = DELTA * dl_dt + Calcium_l[target_cell];

//				Soma Sodium Calcs
				m_inf = 1 / (1 + (expf((-30 - V_soma[target_cell])/ 5.5f)));
				h_inf = 1 / (1 + (expf((-70 - V_soma[target_cell])/-5.8f)));
				tau_h = 3 * expf((-40 - V_soma[target_cell])/33);
				dh_dt = (h_inf - Sodium_h[target_cell])/tau_h;
				Sodium_m[target_cell] = m_inf;
				Sodium_h[target_cell] = Sodium_h[target_cell] + DELTA * dh_dt;

//				Soma Potassium Calcs
				n_inf = 1 / (1 + expf( ( -3 - V_soma[target_cell]) /10));
			        p_inf = 1 / (1 + expf( (-51 - V_soma[target_cell]) / -12));
			        tau_n = 5 + (47 * expf( -(-50 - V_soma[target_cell]) /900));
			        tau_p = tau_n;
			        dn_dt = (n_inf - Potassium_n[target_cell]) / tau_n;
			        dp_dt = (p_inf - Potassium_p[target_cell]) / tau_p;
			        Potassium_n[target_cell] = DELTA * dn_dt + Potassium_n[target_cell];
				Potassium_p[target_cell] = DELTA * dp_dt + Potassium_p[target_cell];

//				Soma Potassium X Calcs
				alpha_x_s = 0.13f * (V_soma[target_cell] + 25) / (1 - expf(-(V_soma[target_cell] + 25) / 10));
			        beta_x_s= 1.69f * expf(-0.0125f * (V_soma[target_cell] + 35));
				x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
				tau_x_s = 1 / (alpha_x_s + beta_x_s);
			        dx_dt_s = (x_inf_s - Potassium_x_s[target_cell]) / tau_x_s;
			        Potassium_x_s[target_cell] = 0.05f * dx_dt_s + Potassium_x_s[target_cell];

//				Somatic Voltage And Current Calcs
				I_ds= (G_INT / P1) * (V_soma[target_cell] - V_dend[target_cell]);
				I_CaL = g_CaL[target_cell] * powf(Calcium_k[target_cell], 3) * Calcium_l[target_cell] * (V_soma[target_cell] - V_CA);
				I_Na_s= G_NA_S * powf(Sodium_m[target_cell], 3) * Sodium_h[target_cell] * (V_soma[target_cell] - V_NA);
				I_ls= G_LS * (V_soma[target_cell] - V_L);
				I_Kdr_s = G_KDR_S * powf(Potassium_n[target_cell], 4) * (V_soma[target_cell] - V_K);
				I_K_s = G_K_S * powf(Potassium_x_s[target_cell], 4) * (V_soma[target_cell] - V_K);
				I_as= (G_INT / (1 - P2)) * (V_soma[target_cell] - V_axon[target_cell]);
				dVs_dt = (-(I_CaL + I_ds+ I_as + I_Na_s + I_ls + I_Kdr_s + I_K_s)) / C_M;


//				~AXONAL COMPUTATIONS~

//				Axon Sodium Calcs
				m_inf_a = 1 / (1 + (expf((-30 - V_axon[target_cell])/ 5.5f)));
				h_inf_a = 1 / (1 + (expf((-60 - V_axon[target_cell])/-5.8f)));
				tau_h_a = 1.5f * expf((-40 - V_axon[target_cell])/33);
				dh_dt_a = (h_inf_a - Sodium_h_a[target_cell])/tau_h_a;
				Sodium_m_a[target_cell] = m_inf_a;
				Sodium_h_a[target_cell] = Sodium_h_a[target_cell] + DELTA * dh_dt_a;

//				Axon Potassium Calcs
				alpha_x_a = 0.13f * (V_axon[target_cell] + 25) / (1 - expf(-(V_axon[target_cell] + 25) / 10));
			        beta_x_a= 1.69f * expf(-0.0125f * (V_axon[target_cell] + 35));
			        x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a);
			        tau_x_a = 1 / (alpha_x_a + beta_x_a);
			        dx_dt_a = (x_inf_a - Potassium_x_a[target_cell]) / tau_x_a;
			        Potassium_x_a[target_cell] = 0.05f * dx_dt_a + Potassium_x_a[target_cell];

//				Axonal Voltage And Current Calcs
				I_Na_a= G_NA_A * powf(Sodium_m_a[target_cell], 3) * Sodium_h_a[target_cell] * (V_axon[target_cell] - V_NA);
				I_la= G_LA* (V_axon[target_cell] - V_L);
				I_sa= (G_INT / P2) * (V_axon[target_cell] - V_soma[target_cell]);
				I_K_a = G_K_A * powf(Potassium_x_a[target_cell], 4) * (V_axon[target_cell] - V_K);
				dVa_dt = (-(I_K_a + I_sa + I_la + I_Na_a)) / C_M;


//				~NEW VOLTAGES~

				V_dend[target_cell] = DELTA * dVd_dt + V_dend[target_cell];
				V_soma[target_cell] = DELTA * dVs_dt + V_soma[target_cell];
				V_axon[target_cell] = DELTA * dVa_dt + V_axon[target_cell];

//				~END OF COMPUTATIONS~
 
			}


			if (PRINTING&&((simStep%print_granularity)==0)) {

			#pragma omp parallel for simd shared (V_axon, pOutFile) private(target_cell, tempbuf)
				for (target_cell=0;target_cell<cellCount;target_cell++) {
					#ifndef G_CAL_FROM_FILE
						sprintf(tempbuf, "%d : %.8f ", target_cell+1, V_axon[target_cell]);
					#else
						sprintf(tempbuf, "%.8f ", V_axon[target_cell]);
					#endif
					fputs(tempbuf, pOutFile);
				}

				sprintf(tempbuf, "\n");
				fputs(tempbuf, pOutFile);
			}
		
		}

//	}

	gettimeofday(&toc, NULL);
		
	/* SIM END
	 * Free  memory and close files
	 */

	if (PRINTING) {
		fclose (pOutFile);
		chmod(outFileName,0x01B6);
	}

	if(inputFromFile)
		fclose (pInFile);

	/* Execution Time for the Sim
	 */

	printf("Execution Time for Simulation: %0.2f ms.\n", ((toc.tv_sec*1000.0 + ((float)(toc.tv_usec)/1000.0)) - \
		(tic.tv_sec*1000.0 + ((float)(tic.tv_usec)/1000.0))));
	if (PRINTING) {
		pOutFile = fopen (outFileName, "rb+");
		sprintf(tempbuf, "Execution Time for Simulation in ms: %0.2f\n", \
			((toc.tv_sec*1000.0 + ((float)(toc.tv_usec)/1000.0)) - \
			(tic.tv_sec*1000.0 + ((float)(tic.tv_usec)/1000.0))));
		#ifndef G_CAL_FROM_FILE
			fputs(tempbuf, pOutFile);
		#endif
		fclose (pOutFile);
	}

	MPI_Finalize();
	return 0;
}

int ReadFileLine(FILE *pInFile, mod_prec *iAppArray){

	char c= fpeek(pInFile);
	if (c==EOF)
		return 0;
	
	char *strNumber;
	int bufSize = cellCount*20;	//more than enough but will do for nao
	char* buffer = (char*) _mm_malloc(bufSize*sizeof(char), 64);

	int floats_ignored = 0;
	int floats_needed = cellCount;
	int useful_element_found = 0;
	//int i will stand for elements already processed, hence it starts at 0 even after the first strtok
	int i = 0;
	
	//Get one line
	if(fgets(buffer, bufSize, pInFile)){
		//Convert the ASCII string of one element to a double precision floating point value
		strNumber = strtok(buffer," ");
		i = 0;

		//printf("Line:\n");
		while ((strNumber != NULL) && (i<IO_NETWORK_SIZE)){
			if (i>=floats_ignored)
				useful_element_found = 1;

			if (i>=floats_ignored+floats_needed)
				useful_element_found = 0;

			if (useful_element_found)
				iAppArray[i-floats_ignored] = atof(strNumber);	//atof() should change if using integers or fixed point
			//printf ("(%s) %0.2f ", strNumber, iAppArray[i]);
			strNumber = strtok(NULL, " ");
			i++; 
		}
		//printf("i: %d\n", i);
		if(i<IO_NETWORK_SIZE){
			//BUG: if only one element is missing but the line ends in a space, the error is not detected
			printf("Error: Input line doesn't have enough elements, only %d\n", i);
			exit(EXIT_FAILURE);
		}
		free(buffer);
		return 1;//success
	}else{
		if(!feof(pInFile)){
			printf("Error: Reading from input file didn't finish successfully\n");
			exit(EXIT_FAILURE);
		}
		free(buffer);
		return 0;//end of file
	}
}

void read_g_CaL_from_file(mod_prec* g_cal) {

	int i;

	FILE* fd = fopen("gcal_file.txt","r");
	for (i=0;i<cellCount;i++) 
		fscanf(fd, "%f ", &g_cal[i]);
	fclose(fd);

	return;

}

inline mod_prec min(mod_prec a, mod_prec b){

	return (a < b) ? a : b;
}

void syncing(MPI_Request* request) {

/* a function that first makes sure all requests passed on to it are completed
 * and then calls an MPI_Barrier to make sure all processes are in sync. We expect
 * one request from each core employed by the app
 */

        int i, done = 0, flag, no_of_reqs = cores;

        while (!done) {
                done = 1;
                for (i=0;i<no_of_reqs;i++) {
                        if (request[i] == MPI_REQUEST_NULL)
                                ;
                        else {
                                MPI_Test(&(request[i]), &flag, MPI_STATUS_IGNORE);
                                if (!flag)
                                        done = 0;
                        }
                }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        return;

}

void cleanup_requests(MPI_Request* request) {

/* a function that makes sure the resources used up by the mpi communication phase
 * are cleaned up properly
 */

        int i, no_of_reqs = cores;

        for (i=0; i< no_of_reqs; i++) {
                if (request[i] == MPI_REQUEST_NULL)
                        ;
                else
                        MPI_Request_free(&request[i]);
        }
        return;

}
