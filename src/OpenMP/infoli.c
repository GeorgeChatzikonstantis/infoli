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
#include "infoli.h"
#include <omp.h>

int core_id, cores, cellCount;
int IO_NETWORK_DIM1, IO_NETWORK_DIM2, IO_NETWORK_SIZE;
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

	int i=0, j, k=0, line_count,l, p, q, x, y, targetCore;
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
	int simulation_array_ID,target_cell, requested_neighbour;

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


	/* outFileName is file for output
	 * conFile is input of the cell's connections
	 */
	sprintf(outFileName,"InferiorOlive_Output.txt");		
	sprintf(conFile,"cellConnections.txt");

	IO_NETWORK_SIZE = atoi(argv[1]); 
	
	/* compute how many grid cells 
	 * are assigned to each core
	 */

        cellCount= IO_NETWORK_SIZE;

	/* Process command line arguments
	 * Case argc = 3 then a one-pulse input is stimulated.
	 * Otherwise we receive input from a specified file in argv[3] and simulation runs accordingly
	 * in the case of inputFromFile, we will also need a buffer to hold the stimulus
	 * for each cell on each step
  	 */
	
	if(argc == 2) {
		inputFromFile = 0;
	}
	else if(argc == 3) {
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
		printf("Error: Too many arguments.\nUsage: ./InferiorOlive <dim1> <dim2> <Iapp_input_file> or ./InferiorOlive <dim1> <dim2>\n");
		exit(EXIT_FAILURE);
	}

	/* PRINTING is true in case we want to write the
	 * output to a specified file (outFileName)
	 */

	if (PRINTING) {
		pOutFile = fopen(outFileName,"w+");
		if(pOutFile==NULL){
			printf("Error: Couldn't create %s\n", outFileName);
			exit(EXIT_FAILURE);
		}
		sprintf(tempbuf, "#simSteps Time(ms) Input(Iapp) Output(V_axon)");
		fputs(tempbuf, pOutFile);
	}

	pConFile = fopen(conFile, "r");
        if(pConFile==NULL){
                printf("Error: Couldn't create %s\n", conFile);
                exit(EXIT_FAILURE);
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

	if (G_CAL_FROM_FILE)
		read_g_CaL_from_file(g_CaL);

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

	mod_prec cond_value= 0;

	if (COMPRESSED_MAP==0) { //non-sparse matrix format for conn. matrix
		int line_counter;

	/* initial passing of the connections file
 	 * necessary to calculate connection counts
 	 */

		for (line_counter=0;line_counter<IO_NETWORK_SIZE;line_counter++) {
			for (i=0; i<IO_NETWORK_SIZE; i++) {

				fscanf(pConFile, "%f ", &cond_value);
				if (cond_value < 0.0000001)
					;		//this connection is considered not existing if conductance value is too small
				else
					cellParamsPtr.total_amount_of_neighbours[i]++;

			}
		}

		for (i=0; i<IO_NETWORK_SIZE; i++) {
			cellParamsPtr.neighConductances[i] = (mod_prec*) _mm_malloc(cellParamsPtr.total_amount_of_neighbours[i]*sizeof(mod_prec), 64);
			cellParamsPtr.neighId[i] = (int*) _mm_malloc(cellParamsPtr.total_amount_of_neighbours[i]*sizeof(int), 64);
			//we reset the amount of neighbours to re-use it in the next segment
			cellParamsPtr.total_amount_of_neighbours[i] = 0;
		}

	/* now the file is re-processed to plan out the
 	 * mapping of the connections
 	 */

		rewind(pConFile);
		for (line_counter=0;line_counter<IO_NETWORK_SIZE;line_counter++) {
			for (i=0; i<IO_NETWORK_SIZE; i++) {

				fscanf(pConFile, "%f ", &cond_value);
				if (cond_value < 0.0000001)
					;
				else {
					cellParamsPtr.neighId[i][cellParamsPtr.total_amount_of_neighbours[i]] = line_counter; //which cell sends this voltage to us
					cellParamsPtr.neighConductances[i][cellParamsPtr.total_amount_of_neighbours[i]] = cond_value;	//what conductance we use to calculate its impact
					cellParamsPtr.total_amount_of_neighbours[i]++;
				}
			}
		}

	/* connections mapping for non-sparse matrix format
 	 * is now complete
 	 */

	}
	else {	//sparse-matrix format for conn. matrix
		char c = fpeek(pConFile);
		int send_cell, rec_cell;

		while (c!=EOF) {

			c = fpeek(pConFile);
			if (c==EOF)
				break;
			fscanf(pConFile, "%d %d %f\n", &send_cell, &rec_cell, &cond_value);
			if (cond_value < 0.0000001)
				;
			else
				cellParamsPtr.total_amount_of_neighbours[rec_cell]++;
		}
		for (i=0; i<IO_NETWORK_SIZE; i++) {
			cellParamsPtr.neighConductances[i] = (mod_prec*) _mm_malloc(cellParamsPtr.total_amount_of_neighbours[i]*sizeof(mod_prec), 64);
			cellParamsPtr.neighId[i] = (int*) _mm_malloc(cellParamsPtr.total_amount_of_neighbours[i]*sizeof(int), 64);
			cellParamsPtr.total_amount_of_neighbours[i] = 0;
		}

		rewind(pConFile);
		c = 'a';

		while (c!=EOF) {

			c = fpeek(pConFile);
			if (c==EOF)
				break;
			fscanf(pConFile, "%d %d %f\n", &send_cell, &rec_cell, &cond_value);
			if (cond_value < 0.0000001)                     //this connection is considered not existing if conductance value is too small
				;
			else {
				cellParamsPtr.neighId[rec_cell][cellParamsPtr.total_amount_of_neighbours[rec_cell]] = send_cell;
				cellParamsPtr.neighConductances[rec_cell][cellParamsPtr.total_amount_of_neighbours[rec_cell]] = cond_value;
				cellParamsPtr.total_amount_of_neighbours[rec_cell]++;
			}
		}

	/* connections mapping for sparse matrix format
	 * is now complete
	 */

printf("%d %d\n", cellParamsPtr.total_amount_of_neighbours[0], cellParamsPtr.total_amount_of_neighbours[1]);
	}
	if (PRINTING) {
		for (i=0;i<cellCount;i++) {
			sprintf(tempbuf, "%d ", cellID[i]);
			fputs(tempbuf, pOutFile);
		}
		sprintf(tempbuf, "\n");
		fputs(tempbuf, pOutFile);
	}
	
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
*/		simTime = SIMTIME; 
		total_simulation_steps = (int)(simTime/DELTA);

		for(simStep=0;simStep<total_simulation_steps;simStep++) {
								
			if ((simStep>=20000)&&(simStep<20500-1))
				iApp = 6; 
			else
		   		iApp = 0;
			
			if (PRINTING) {
				sprintf(tempbuf, " %d %.2f ", simStep+1, iApp);
				fputs(tempbuf, pOutFile);
			}

	/*	First Loop takes care of Gap Junction Functions
 	*/			

		#pragma omp parallel for shared (cellParamsPtr, iAppIn, V_dend, I_c, pOutFile) private(iApp, target_cell, i, \
		requested_neighbour, f, V, voltage, I_c_storage) firstprivate(simulation_array_ID)
			for (target_cell=0;target_cell<cellCount;target_cell++) {

				//Feeding of Input Current
				iAppIn[target_cell] = iApp;

	/*	Gathering of the data concerning Vdend of neighours, then
 	*	computing and storing the incoming Ic from neighbours
 	*/

				I_c_storage = 0;
				__assume_aligned(cellParamsPtr.neighId[target_cell], 64);
				__assume_aligned(cellParamsPtr.neighConductances[target_cell], 64);
				__assume_aligned(V_dend, 64);
			#pragma ivdep
				for (i=0; i<cellParamsPtr.total_amount_of_neighbours[target_cell]; i++){
					requested_neighbour = cellParamsPtr.neighId[target_cell][i];
					voltage = V_dend[requested_neighbour];
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
		iAppIn, I_c, V_soma, g_CaL, Sodium_m, Sodium_h, Calcium_k, Calcium_l, Potassium_n, Potassium_p, Potassium_x_s,\
		V_axon, Sodium_m_a, Sodium_h_a, Potassium_x_a, pOutFile) private(iApp, target_cell, tempbuf, q_inf, tau_q, dq_dt,\
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

/*				This I/O part is commented out because it ruins vectorization
 
				if (PRINTING) {
					sprintf(tempbuf, "%d : %.8f ", target_cell+1, V_axon[target_cell]);
					fputs(tempbuf, pOutFile);
				}

*/			}

			if (PRINTING) {
				sprintf(tempbuf, "\n");
				fputs(tempbuf, pOutFile);
			}
		
		}

//	}
		
	/* Free  memory and close files
	 */

	if (PRINTING) {
		fclose (pOutFile);
		chmod(outFileName,0x01B6);
	}

	if(inputFromFile)
		fclose (pInFile);

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
