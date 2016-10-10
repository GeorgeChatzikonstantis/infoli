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
#include "infoli_v3.h"
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
	cellState **cellPtr;
	cellCompParams *cellParamsPtr;
	int seedvar;
	char tempbuf[100];		
	mod_prec iApp, voltage;
	mod_prec *gcal;
	int *init_steps;
	long double seconds;
	int simulation_array_ID,target_cell, requested_neighbour;

	/* declaration of variables
	 * necessary for ion channel
	 * computations
	 */

	mod_prec q_inf, tau_q, dq_dt;
	mod_prec alpha_r, beta_r, r_inf, tau_r, dr_dt;
	mod_prec alpha_s, beta_s, s_inf, tau_s, ds_dt;
	mod_prec dCa_dt, f, V, I_c;
	mod_prec I_sd, I_CaH, I_K_Ca, I_ld, I_h, dVd_dt;

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

	
	/* CellPtr is a 2-D array of size 2*CellCount 
	 * and containd cell states used in every simulation 
	 * step accordingly
	 */	
	cellPtr = _mm_malloc(2*sizeof(cellState *), 64);
	if(cellPtr==NULL){
		printf("Error: Couldn't malloc for cellPtr\n");
		exit(EXIT_FAILURE);
	}
	
	cellPtr[0] = _mm_malloc(cellCount*sizeof(cellState), 64);
	cellPtr[1] = _mm_malloc(cellCount*sizeof(cellState), 64);
	if ((!cellPtr[0])||(!cellPtr[1])) {
		printf("Error: Couldn't malloc the array for cellStates\n");
		exit(EXIT_FAILURE);
	}

	/* cellCompParams struct is used to update a cell state.
	 * It contains the current flowing to this specific cell
	 * voltages from communicating cells , and pointers to the previous
	 * cell state and to the next cell state
	 * We allocate cellCount cellCompParams for all the core's cells
	 * Pointer to previous and next states are pointers to CellPtr[0] and CellPtr[1]
	 * elements.
	 */
	cellParamsPtr = _mm_malloc(cellCount*sizeof(cellCompParams), 64);
	if(cellParamsPtr==NULL){
		printf("Error: Couldn't malloc for cellParamsPtr\n");
		exit(EXIT_FAILURE);
	}
	for (i=0;i<cellCount;i++)					//initial amount of neighbours for each cell is 0 (bug fix in case the cell stays isolated)
		cellParamsPtr[i].total_amount_of_neighbours = 0;

	int line_counter;
	mod_prec cond_value= 0;

	if (COMPRESSED_MAP==0) { //non-sparse matrix format for conn. matrix
		int line_counter;

	/* initial passing of the connections file
 	 * necessary to calculate neighbours' counts
 	 */

		for (line_counter=0;line_counter<IO_NETWORK_SIZE;line_counter++) {
			for (i=0; i<IO_NETWORK_SIZE; i++) {

				fscanf(pConFile, "%f ", &cond_value);
				if (cond_value < 0.0000001)
					;		//this connection is considered not existing if conductance value is too small
				else
					cellParamsPtr[i].total_amount_of_neighbours++;

			}
		}

		for (i=0; i<IO_NETWORK_SIZE; i++) {
			cellParamsPtr[i].neighConductances = (mod_prec*) _mm_malloc(cellParamsPtr[i].total_amount_of_neighbours*sizeof(mod_prec), 64);
			cellParamsPtr[i].neighId = (int*) _mm_malloc(cellParamsPtr[i].total_amount_of_neighbours*sizeof(int), 64);
			//we reset the amount of neighbours to re-use it in the next segment
			cellParamsPtr[i].total_amount_of_neighbours = 0;
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
					cellParamsPtr[i].neighId[cellParamsPtr[i].total_amount_of_neighbours] = line_counter;//which cell sends this voltage to us
					cellParamsPtr[i].neighConductances[cellParamsPtr[i].total_amount_of_neighbours] = cond_value;	//what conductance we use to calculate its impact
					cellParamsPtr[i].total_amount_of_neighbours++;
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
				cellParamsPtr[rec_cell].total_amount_of_neighbours++;
		}

		for (i=0; i<IO_NETWORK_SIZE; i++) {
			cellParamsPtr[i].neighConductances = (mod_prec*) _mm_malloc(cellParamsPtr[i].total_amount_of_neighbours*sizeof(mod_prec), 64);
			cellParamsPtr[i].neighId = (int*) _mm_malloc(cellParamsPtr[i].total_amount_of_neighbours*sizeof(int), 64);
			cellParamsPtr[i].total_amount_of_neighbours = 0;
		}

		rewind(pConFile);
		c= 'a';

		while (c!=EOF) {

			c = fpeek(pConFile);
			if (c==EOF)
				break;
			fscanf(pConFile, "%d %d %f\n", &send_cell, &rec_cell, &cond_value);
			if (cond_value < 0.0000001)                     //this connection is considered not existing if conductance value is too small
				;
			else {
				cellParamsPtr[rec_cell].neighId[cellParamsPtr[rec_cell].total_amount_of_neighbours] = send_cell;
				cellParamsPtr[rec_cell].neighConductances[cellParamsPtr[rec_cell].total_amount_of_neighbours] = cond_value;
				cellParamsPtr[rec_cell].total_amount_of_neighbours++;
			}
		}

	/* connections mapping for sparse matrix format
	 * is now complete
	 */

	}

printf("%d %d\n", cellParamsPtr[0].total_amount_of_neighbours, cellParamsPtr[1].total_amount_of_neighbours);
		
	/* Initialise cellPtr[0] with appropriate values.
	 */

	initState(cellPtr[0]);
	
	if (PRINTING) {
		for (i=0;i<cellCount;i++) {
			sprintf(tempbuf, "%d ", cellPtr[0][i].cellID);
			fputs(tempbuf, pOutFile);
		}
		sprintf(tempbuf, "\n");
		fputs(tempbuf, pOutFile);
	}
	
	//Initialize g_CaL
	seedvar = time(NULL);		//seedvar will be time- and core-dependent now!
	srand(seedvar);

	for(i=0;i<cellCount;i++){
		cellPtr[1][i].soma.g_CaL = cellPtr[0][i].soma.g_CaL;
		if (RAND_INIT) {
			cellPtr[0][i].soma.g_CaL = 0.6+(0.2*(rand()%100)/100);
			cellPtr[1][i].soma.g_CaL = cellPtr[0][i].soma.g_CaL;
		}
	}

	//random initialization process
	if (RAND_INIT) {
		for(i=0;i<cellCount;i++) {
			initSteps = rand()%(int)ceil(100/DELTA);
			initSteps = initSteps | 0x00000001;//make it odd, so that the final state is in prevCellState

			for(j=0;j<initSteps;j++){
				//Arrange inputs
				cellParamsPtr[i].iAppIn = 0;//No stimulus
				cellParamsPtr[i].prevCellState = &cellPtr[j%2][i];
				cellParamsPtr[i].newCellState = &cellPtr[(j%2)^1][i];

//				CompDend(&cellParamsPtr[i], 1);
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
								
			simulation_array_ID = simStep%2;
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

		#pragma omp parallel for shared (cellParamsPtr, cellPtr, pOutFile) private(iApp, target_cell, i, requested_neighbour, \
		f, V, voltage, I_c) firstprivate(simulation_array_ID)
			for (target_cell=0;target_cell<cellCount;target_cell++) {

//		Feeding of Input Current and Swapping of New and Old NW States
				cellParamsPtr[target_cell].iAppIn = iApp;
				cellParamsPtr[target_cell].prevCellState = &cellPtr[simulation_array_ID][target_cell];
				cellParamsPtr[target_cell].newCellState = &cellPtr[simulation_array_ID^1][target_cell];

	/*	Gathering of the data concerning Vdend of neighours, then
 	*	computing and storing the incoming Ic from neighbours
 	*/

				I_c = 0;
				__assume_aligned(cellParamsPtr, 64);
				__assume_aligned(cellParamsPtr[target_cell].neighId, 64);
				__assume_aligned(cellParamsPtr[target_cell].neighConductances, 64);
				__assume_aligned(cellPtr, 64);
				__assume_aligned(cellPtr[0], 64);
				__assume_aligned(cellPtr[1], 64);
			#pragma ivdep
				for (i=0; i<cellParamsPtr[target_cell].total_amount_of_neighbours; i++){
					requested_neighbour = cellParamsPtr[target_cell].neighId[i];
					voltage = cellPtr[simulation_array_ID][requested_neighbour].dend.V_dend;
					V = cellParamsPtr[target_cell].prevCellState->dend.V_dend - voltage;
					f = 0.8f * expf(-1*powf(V, 2)/100) + 0.2f;// SCHWEIGHOFER 2004 VERSION
					I_c += cellParamsPtr[target_cell].neighConductances[i] * f * V;
                                }	
				cellParamsPtr[target_cell].I_c = I_c;
			}

	/* Second Loop computes the new state (voltages, channels etc.)
 	*/ 

			__assume_aligned(cellParamsPtr, 64);
			__assume_aligned(cellParamsPtr[target_cell].neighId, 64);
			__assume_aligned(cellParamsPtr[target_cell].neighConductances, 64);
			__assume_aligned(cellPtr, 64);
			__assume_aligned(cellPtr[0], 64);
			__assume_aligned(cellPtr[1], 64);
		#pragma omp parallel for simd shared (cellParamsPtr, cellPtr, pOutFile) private(iApp, target_cell, tempbuf, q_inf, \
		tau_q, dq_dt, alpha_r, beta_r, r_inf, tau_r, dr_dt, alpha_s, beta_s, s_inf, tau_s, ds_dt, dCa_dt, I_sd, I_CaH, I_K_Ca, \
		I_ld, I_h, dVd_dt, k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt, m_inf, h_inf, tau_h, dh_dt, n_inf, p_inf, tau_n, tau_p, \
		dn_dt, dp_dt, alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s, I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt, \
		m_inf_a, h_inf_a, tau_h_a, dh_dt_a, alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a, I_Na_a, I_la, I_sa, I_K_a, dVa_dt) \
		firstprivate(simulation_array_ID)
			for (target_cell=0;target_cell<cellCount;target_cell++) {

//				~DENDRITIC COMPUTATIONS~

//				Dend H Current Calcs
				q_inf = 1 /(1 + expf((cellParamsPtr[target_cell].prevCellState->dend.V_dend + 80) / 4));
			        tau_q = 1 /(expf(-0.086f * cellParamsPtr[target_cell].prevCellState->dend.V_dend - 14.6f) + expf(0.070f * cellParamsPtr[target_cell].prevCellState->dend.V_dend - 1.87f));
			        dq_dt = (q_inf - cellParamsPtr[target_cell].prevCellState->dend.Hcurrent_q) / tau_q;
			        cellParamsPtr[target_cell].newCellState->dend.Hcurrent_q = DELTA * dq_dt + cellParamsPtr[target_cell].prevCellState->dend.Hcurrent_q;

//				Dend Ca Current Calcs
				alpha_r = 1.7f / (1 + expf( -(cellParamsPtr[target_cell].prevCellState->dend.V_dend - 5) / 13.9f));
		       		beta_r = 0.02f * (cellParamsPtr[target_cell].prevCellState->dend.V_dend + 8.5f) / (expf((cellParamsPtr[target_cell].prevCellState->dend.V_dend + 8.5f) / 5) - 1);
	  		     	r_inf = alpha_r / (alpha_r + beta_r);
	       	 		tau_r = 5 / (alpha_r + beta_r);
	       		 	dr_dt = (r_inf - cellParamsPtr[target_cell].prevCellState->dend.Calcium_r) / tau_r;
	    			cellParamsPtr[target_cell].newCellState->dend.Calcium_r = DELTA * dr_dt + cellParamsPtr[target_cell].prevCellState->dend.Calcium_r;

//				Dend K Current Calcs
				alpha_s = min((0.00002f * cellParamsPtr[target_cell].prevCellState->dend.Ca2Plus), 0.01f);
			        beta_s = 0.015f;
			        s_inf = alpha_s / (alpha_s + beta_s);
			        tau_s = 1 / (alpha_s + beta_s);
			        ds_dt = (s_inf - cellParamsPtr[target_cell].prevCellState->dend.Potassium_s) / tau_s;
			        cellParamsPtr[target_cell].newCellState->dend.Potassium_s = DELTA * ds_dt + cellParamsPtr[target_cell].prevCellState->dend.Potassium_s;

//				Dend Cal Current Calcs
				dCa_dt = -3 * cellParamsPtr[target_cell].prevCellState->dend.I_CaH - 0.075f * cellParamsPtr[target_cell].prevCellState->dend.Ca2Plus;
				cellParamsPtr[target_cell].newCellState->dend.Ca2Plus = DELTA * dCa_dt + cellParamsPtr[target_cell].prevCellState->dend.Ca2Plus;

//				Dendritic Voltage And Current Cacls
				I_sd = (G_INT / (1 - P1)) * (cellParamsPtr[target_cell].prevCellState->dend.V_dend - cellParamsPtr[target_cell].prevCellState->soma.V_soma);
				I_CaH =G_CAH * powf(cellParamsPtr[target_cell].newCellState->dend.Calcium_r, 2) * (cellParamsPtr[target_cell].prevCellState->dend.V_dend - V_CA);
				I_K_Ca =G_K_CA * cellParamsPtr[target_cell].newCellState->dend.Potassium_s * (cellParamsPtr[target_cell].prevCellState->dend.V_dend - V_K);
				I_ld =G_LD * (cellParamsPtr[target_cell].prevCellState->dend.V_dend - V_L);
				I_h=G_H * cellParamsPtr[target_cell].newCellState->dend.Hcurrent_q * (cellParamsPtr[target_cell].prevCellState->dend.V_dend - V_H);
				dVd_dt = (-(I_CaH + I_sd+ I_ld + I_K_Ca + cellParamsPtr[target_cell].I_c + I_h) + cellParamsPtr[target_cell].iAppIn) / C_M;
				cellParamsPtr[target_cell].newCellState->dend.V_dend = DELTA * dVd_dt + cellParamsPtr[target_cell].prevCellState->dend.V_dend;
				cellParamsPtr[target_cell].newCellState->dend.I_CaH = I_CaH;


//				~SOMATIC COMPUTATIONS~

//				Soma Calcium Calcs
				k_inf = (1 / (1 + expf(-1 * (cellParamsPtr[target_cell].prevCellState->soma.V_soma + 61) / 4.2f)));
			        l_inf = (1 / (1 + expf(( cellParamsPtr[target_cell].prevCellState->soma.V_soma + 85.5f) / 8.5f)));
			        tau_k = 1;
			        tau_l = ((20 * expf((cellParamsPtr[target_cell].prevCellState->soma.V_soma + 160) / 30) / (1 + expf((cellParamsPtr[target_cell].prevCellState->soma.V_soma + 84) / 7.3f))) +35);
			        dk_dt = (k_inf - cellParamsPtr[target_cell].prevCellState->soma.Calcium_k) / tau_k;
			        dl_dt = (l_inf - cellParamsPtr[target_cell].prevCellState->soma.Calcium_l) / tau_l;
			        cellParamsPtr[target_cell].newCellState->soma.Calcium_k = DELTA * dk_dt + cellParamsPtr[target_cell].prevCellState->soma.Calcium_k;
			        cellParamsPtr[target_cell].newCellState->soma.Calcium_l = DELTA * dl_dt + cellParamsPtr[target_cell].prevCellState->soma.Calcium_l;

//				Soma Sodium Calcs
				m_inf = 1 / (1 + (expf((-30 - cellParamsPtr[target_cell].prevCellState->soma.V_soma)/ 5.5f)));
				h_inf = 1 / (1 + (expf((-70 - cellParamsPtr[target_cell].prevCellState->soma.V_soma)/-5.8f)));
				tau_h = 3 * expf((-40 - cellParamsPtr[target_cell].prevCellState->soma.V_soma)/33);
				dh_dt = (h_inf - cellParamsPtr[target_cell].prevCellState->soma.Sodium_h)/tau_h;
				cellParamsPtr[target_cell].newCellState->soma.Sodium_m = m_inf;
				cellParamsPtr[target_cell].newCellState->soma.Sodium_h = cellParamsPtr[target_cell].prevCellState->soma.Sodium_h + DELTA * dh_dt;

//				Soma Potassium Calcs
				n_inf = 1 / (1 + expf( ( -3 - cellParamsPtr[target_cell].prevCellState->soma.V_soma) /10));
			        p_inf = 1 / (1 + expf( (-51 - cellParamsPtr[target_cell].prevCellState->soma.V_soma) / -12));
			        tau_n = 5 + (47 * expf( -(-50 - cellParamsPtr[target_cell].prevCellState->soma.V_soma) /900));
			        tau_p = tau_n;
			        dn_dt = (n_inf - cellParamsPtr[target_cell].prevCellState->soma.Potassium_n) / tau_n;
			        dp_dt = (p_inf - cellParamsPtr[target_cell].prevCellState->soma.Potassium_p) / tau_p;
			        cellParamsPtr[target_cell].newCellState->soma.Potassium_n = DELTA * dn_dt + cellParamsPtr[target_cell].prevCellState->soma.Potassium_n;
				cellParamsPtr[target_cell].newCellState->soma.Potassium_p = DELTA * dp_dt + cellParamsPtr[target_cell].prevCellState->soma.Potassium_p;

//				Soma Potassium X Calcs
				alpha_x_s = 0.13f * (cellParamsPtr[target_cell].prevCellState->soma.V_soma + 25) / (1 - expf(-(cellParamsPtr[target_cell].prevCellState->soma.V_soma + 25) / 10));
			        beta_x_s= 1.69f * expf(-0.0125f * (cellParamsPtr[target_cell].prevCellState->soma.V_soma + 35));
				x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
				tau_x_s = 1 / (alpha_x_s + beta_x_s);
			        dx_dt_s = (x_inf_s - cellParamsPtr[target_cell].prevCellState->soma.Potassium_x_s) / tau_x_s;
			        cellParamsPtr[target_cell].newCellState->soma.Potassium_x_s = 0.05f * dx_dt_s + cellParamsPtr[target_cell].prevCellState->soma.Potassium_x_s;

//				Somatic Voltage And Current Calcs
				I_ds= (G_INT / P1) * (cellParamsPtr[target_cell].prevCellState->soma.V_soma - cellParamsPtr[target_cell].prevCellState->dend.V_dend);
				I_CaL = cellParamsPtr[target_cell].prevCellState->soma.g_CaL * powf(cellParamsPtr[target_cell].newCellState->soma.Calcium_k, 3) * cellParamsPtr[target_cell].newCellState->soma.Calcium_l * (cellParamsPtr[target_cell].prevCellState->soma.V_soma - V_CA);
				I_Na_s= G_NA_S * powf(cellParamsPtr[target_cell].newCellState->soma.Sodium_m, 3) * cellParamsPtr[target_cell].newCellState->soma.Sodium_h * (cellParamsPtr[target_cell].prevCellState->soma.V_soma - V_NA);
				I_ls= G_LS * (cellParamsPtr[target_cell].prevCellState->soma.V_soma - V_L);
				I_Kdr_s = G_KDR_S * powf(cellParamsPtr[target_cell].newCellState->soma.Potassium_n, 4) * (cellParamsPtr[target_cell].prevCellState->soma.V_soma - V_K);
				I_K_s = G_K_S * powf(cellParamsPtr[target_cell].newCellState->soma.Potassium_x_s, 4) * (cellParamsPtr[target_cell].prevCellState->soma.V_soma - V_K);
				I_as= (G_INT / (1 - P2)) * (cellParamsPtr[target_cell].prevCellState->soma.V_soma - cellParamsPtr[target_cell].prevCellState->axon.V_axon);
				dVs_dt = (-(I_CaL + I_ds+ I_as + I_Na_s + I_ls + I_Kdr_s + I_K_s)) / C_M;
				cellParamsPtr[target_cell].newCellState->soma.V_soma = DELTA * dVs_dt + cellParamsPtr[target_cell].prevCellState->soma.V_soma;


//				~AXONAL COMPUTATIONS~

//				Axon Sodium Calcs
				m_inf_a = 1 / (1 + (expf((-30 - cellParamsPtr[target_cell].prevCellState->axon.V_axon)/ 5.5f)));
				h_inf_a = 1 / (1 + (expf((-60 - cellParamsPtr[target_cell].prevCellState->axon.V_axon)/-5.8f)));
				tau_h_a = 1.5f * expf((-40 - cellParamsPtr[target_cell].prevCellState->axon.V_axon)/33);
				dh_dt_a = (h_inf_a - cellParamsPtr[target_cell].prevCellState->axon.Sodium_h_a)/tau_h_a;
				cellParamsPtr[target_cell].newCellState->axon.Sodium_m_a = m_inf_a;
				cellParamsPtr[target_cell].newCellState->axon.Sodium_h_a = cellParamsPtr[target_cell].prevCellState->axon.Sodium_h_a + DELTA * dh_dt_a;

//				Axon Potassium Calcs
				alpha_x_a = 0.13f * (cellParamsPtr[target_cell].prevCellState->axon.V_axon + 25) / (1 - expf(-(cellParamsPtr[target_cell].prevCellState->axon.V_axon + 25) / 10));
			        beta_x_a= 1.69f * expf(-0.0125f * (cellParamsPtr[target_cell].prevCellState->axon.V_axon + 35));
			        x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a);
			        tau_x_a = 1 / (alpha_x_a + beta_x_a);
			        dx_dt_a = (x_inf_a - cellParamsPtr[target_cell].prevCellState->axon.Potassium_x_a) / tau_x_a;
			        cellParamsPtr[target_cell].newCellState->axon.Potassium_x_a = 0.05f * dx_dt_a + cellParamsPtr[target_cell].prevCellState->axon.Potassium_x_a;

//				Axonal Voltage And Current Calcs
				I_Na_a= G_NA_A * powf(cellParamsPtr[target_cell].newCellState->axon.Sodium_m_a, 3) * cellParamsPtr[target_cell].newCellState->axon.Sodium_h_a * (cellParamsPtr[target_cell].prevCellState->axon.V_axon - V_NA);
				I_la= G_LA* (cellParamsPtr[target_cell].prevCellState->axon.V_axon - V_L);
				I_sa= (G_INT / P2) * (cellParamsPtr[target_cell].prevCellState->axon.V_axon - cellParamsPtr[target_cell].prevCellState->soma.V_soma);
				I_K_a = G_K_A * powf(cellParamsPtr[target_cell].newCellState->axon.Potassium_x_a, 4) * (cellParamsPtr[target_cell].prevCellState->axon.V_axon - V_K);
				dVa_dt = (-(I_K_a + I_sa + I_la + I_Na_a)) / C_M;
				cellParamsPtr[target_cell].newCellState->axon.V_axon = DELTA * dVa_dt + cellParamsPtr[target_cell].prevCellState->axon.V_axon;


//				~END OF COMPUTATIONS~

/*				This I/O part is commented out because it ruins vectorization
 
				if (PRINTING) {
					sprintf(tempbuf, "%d : %.8f ", target_cell+1, cellPtr[(simulation_array_ID)^1][target_cell].axon.V_axon);
					fputs(tempbuf, pOutFile);
				}

*/			}

			if (PRINTING) {
				sprintf(tempbuf, "\n");
				fputs(tempbuf, pOutFile);
			}
		
		}

		simStep = total_simulation_steps;//so that simStep in the end has the exact value of how many steps we had in this sim, regardless of input method (useful to know which cellPtr has what)
//	}
		
	/* Free  memory and close files
	 */

	if (PRINTSTATE) {
		char resultFileName[50];
		sprintf(resultFileName, "lastStateDump.txt");
		printState(cellPtr[simStep%2], resultFileName);		//simStep%2 here should refer to the cellPtr which has the last state of the network that we calculated
	}

//	_mm_free(cellPtr[0]);
//	_mm_free(cellPtr[1]);
//	_mm_free(cellPtr);
//	_mm_free(cellParamsPtr);
	
	if (PRINTING) {
		fclose (pOutFile);
		chmod(outFileName,0x01B6);
	}

	if(inputFromFile)
		fclose (pInFile);

	return 0;
}


//Initialization Function, important ! --------------------------

void initState(cellState *cellPtr){

	int i, j;
	cellState initState;
	//Initial dendritic parameters
	initState.dend.V_dend = -60;
	initState.dend.Calcium_r = 0.0112788;// High-threshold calcium
	initState.dend.Potassium_s = 0.0049291;// Calcium-dependent potassium
	initState.dend.Hcurrent_q = 0.0337836;// H current
	initState.dend.Ca2Plus = 3.7152;// Calcium concentration
	initState.dend.I_CaH = 0.5;// High-threshold calcium current
	//Initial somatic parameters
	initState.soma.g_CaL = 0.68; //default arbitrary value but it should be randomized per cell
	initState.soma.V_soma = -60;
	initState.soma.Sodium_m = 1.0127807;// Sodium (artificial)
	initState.soma.Sodium_h = 0.3596066;
	initState.soma.Potassium_n = 0.2369847;// Potassium (delayed rectifier)
	initState.soma.Potassium_p = 0.2369847;
	initState.soma.Potassium_x_s = 0.1;// Potassium (voltage-dependent)
	initState.soma.Calcium_k = 0.7423159;// Low-threshold calcium
	initState.soma.Calcium_l = 0.0321349;
	// Initial axonal parameters
	initState.axon.V_axon = -60;
	//sisaza: Sodium_m_a doesn't have a state, therefore this assignment doesn'thave any effect
	initState.axon.Sodium_m_a = 0.003596066;// Sodium (thalamocortical)
	initState.axon.Sodium_h_a = 0.9;
	initState.axon.Potassium_x_a = 0.2369847;// Potassium (transient)

//	initState.cell_x=0;
//	initState.cell_y=0; //Initial irrelevant value of compartment's ID and coords
//	initState.cellID= core_id * cellCount;	//Compute the cellID of the first cell in this core
	initState.cellID= 0;	//Compute the cellID of the first cell in this core

	//Copy init state to all cell states and calculate their coords
	for(j=0;j<cellCount;j++){
//		initState.cell_x= initState.cellID / IO_NETWORK_DIM2;		//we assume that dim1 refers to the number of ROWS and dim2 refers to COLUMNS !!!!
//		initState.cell_y= initState.cellID % IO_NETWORK_DIM2;	
		memcpy(&cellPtr[j], &initState, sizeof(cellState));
		initState.cellID ++;				//next cell's ID is increased by 1
	}

	if (G_CAL_FROM_FILE)
		read_g_CaL_from_file(cellPtr);
	
	return;
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

void read_g_CaL_from_file(cellState* cellPtr) {

	int i;

	FILE* fd = fopen("gcal_file.txt","r");
	for (i=0;i<cellCount;i++) 
		fscanf(fd, "%f ", &cellPtr[i].soma.g_CaL);
	fclose(fd);

	return;

}

void printState(cellState* cellPtr, char* filename) {

	FILE* fd = fopen(filename, "w");
	mod_prec* s = (mod_prec*) malloc(16*sizeof(mod_prec));
	int i, j;

	for (i=0; i<cellCount; i++) {

		s[0] = cellPtr[i].soma.V_soma;
		s[1] = cellPtr[i].soma.Sodium_m;
		s[2] = cellPtr[i].soma.Potassium_n;
		s[3] = cellPtr[i].soma.Potassium_x_s;
		s[4] = cellPtr[i].soma.Calcium_k;
		s[5] = cellPtr[i].soma.Calcium_l;
		s[6] = cellPtr[i].dend.V_dend;
		s[7] = cellPtr[i].dend.Calcium_r;
		s[8] = cellPtr[i].dend.Potassium_s;
		s[9] = cellPtr[i].dend.Hcurrent_q;
		s[10] = cellPtr[i].dend.Ca2Plus;
		s[11] = cellPtr[i].dend.I_CaH;
		s[12] = cellPtr[i].axon.V_axon;
		s[13] = cellPtr[i].axon.Sodium_m_a;
		s[14] = cellPtr[i].axon.Sodium_h_a;
		s[15] = cellPtr[i].axon.Potassium_x_a;

		for (j=0; j<16; j++)
			fprintf(fd, "%.8f ", s[j]);
		fprintf(fd, "\n");

	}

	fclose(fd);
	return;

}

inline mod_prec min(mod_prec a, mod_prec b){

	return (a < b) ? a : b;
}
