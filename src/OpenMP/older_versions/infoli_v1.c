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
#include "infoli_v1.h"
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
	
	/* RCCE_init everytime we want 
     * to use the library 
     * RCCC_ue returns id of the cor
     * RCCC_num_ues is the number of the cores used at this time
     * dimensions of the cell grid network ( IO network )
	 * 
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
	//packet* pointer;	
	mod_prec iApp, voltage;
	mod_prec *gcal;
	int *init_steps;
	long double seconds;
	int simulation_array_ID,target_cell, requested_neighbour;
	communication_node *communicating_cell,*next_in_list,*communication_list_head;
	//communicating_packet *packet_in = malloc(sizeof(communicating_packet));
	//communicating_packet *packet_out = malloc(sizeof(communicating_packet));
	sending_node *sending_list_head = NULL;
	receiving_node *receiving_list_head = NULL;

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
		iAppArray= (mod_prec*) malloc(cellCount*(sizeof(mod_prec)));
	}	
	else {
		printf("Error: Too many arguments.\nUsage: ./InferiorOlive <dim1> <dim2> <Iapp_input_file> or ./InferiorOlive <dim1> <dim2>\n");
		exit(EXIT_FAILURE);
	}

//	printf("OK 1\n");
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
	cellPtr = malloc(2*sizeof(cellState *));
	if(cellPtr==NULL){
		printf("Error: Couldn't malloc for cellPtr\n");
		exit(EXIT_FAILURE);
	}
	
	cellPtr[0] = malloc(cellCount*sizeof(cellState));
	cellPtr[1] = malloc(cellCount*sizeof(cellState));
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
	cellParamsPtr = malloc(cellCount*sizeof(cellCompParams));
	if(cellParamsPtr==NULL){
		printf("Error: Couldn't malloc for cellParamsPtr\n");
		exit(EXIT_FAILURE);
	}
	for (i=0;i<cellCount;i++)					//initial amount of neighbours for each cell is 0 (bug fix in case the cell stays isolated)
		cellParamsPtr[i].total_amount_of_neighbours = 0;

	int line_counter;
	mod_prec cond_value= 0;

	if (COMPRESSED_MAP==0) {
		int line_counter;
		for (line_counter=0;line_counter<IO_NETWORK_SIZE;line_counter++) {
			for (i=0; i<IO_NETWORK_SIZE; i++) {

				fscanf(pConFile, "%f ", &cond_value);
				if (cond_value < 0.0000001)			//this connection is considered not existing if conductance value is too small
					;
				else {
				//part of the code handling RECEIVING and noting which of my cells needs input from which other cells, from ANY core

					if (cellParamsPtr[i].total_amount_of_neighbours == 0) {                  //if this is the first neighbour, initialize buffers
						cellParamsPtr[i].neighVdend = NULL;
						cellParamsPtr[i].neighConductances = NULL;
						cellParamsPtr[i].neighId = NULL;
					}

					cellParamsPtr[i].neighId = allocate_space_int(cellParamsPtr[i].neighId, cellParamsPtr[i].total_amount_of_neighbours);
					cellParamsPtr[i].neighId[cellParamsPtr[i].total_amount_of_neighbours] = line_counter;//which cell sends this voltage to us (GLOBAL ID)
					cellParamsPtr[i].neighConductances = allocate_space_mod(cellParamsPtr[i].neighConductances, cellParamsPtr[i].total_amount_of_neighbours);
					cellParamsPtr[i].neighConductances[cellParamsPtr[i].total_amount_of_neighbours] = cond_value;	//what conductance we use to calculate its impact
					//allocate space for storing this voltage
					cellParamsPtr[i].neighVdend = allocate_space_mod(cellParamsPtr[i].neighVdend, cellParamsPtr[i].total_amount_of_neighbours);
					cellParamsPtr[i].total_amount_of_neighbours++;
				}
			}
		}
	}
	else {
		char c = fpeek(pConFile);
		int send_cell, rec_cell;
		while (c!=EOF) {
			c = fpeek(pConFile);
			if (c==EOF)
				break;
			fscanf(pConFile, "%d %d %f\n", &send_cell, &rec_cell, &cond_value);
			if (cond_value < 0.0000001)                     //this connection is considered not existing if conductance value is too small
				;
			else {

				if (cellParamsPtr[rec_cell].total_amount_of_neighbours == 0) {                  //if this is the first neighbour, initialize buffers
					cellParamsPtr[rec_cell].neighVdend = NULL;
					cellParamsPtr[rec_cell].neighConductances = NULL;
					cellParamsPtr[rec_cell].neighId = NULL;
				}
				cellParamsPtr[rec_cell].neighId = allocate_space_int(cellParamsPtr[rec_cell].neighId, cellParamsPtr[rec_cell].total_amount_of_neighbours);
				cellParamsPtr[rec_cell].neighId[cellParamsPtr[rec_cell].total_amount_of_neighbours] = send_cell;
				cellParamsPtr[rec_cell].neighConductances = allocate_space_mod(cellParamsPtr[rec_cell].neighConductances, cellParamsPtr[rec_cell].total_amount_of_neighbours);
				cellParamsPtr[rec_cell].neighConductances[cellParamsPtr[rec_cell].total_amount_of_neighbours] = cond_value;
				cellParamsPtr[rec_cell].neighVdend = allocate_space_mod(cellParamsPtr[rec_cell].neighVdend, cellParamsPtr[rec_cell].total_amount_of_neighbours);
				cellParamsPtr[rec_cell].total_amount_of_neighbours++;

			}
		}
	}
		
	/* Initialise cellPtr[0] with appropriate values.
	 */

printf("%d %d\n", cellParamsPtr[0].total_amount_of_neighbours, cellParamsPtr[1].total_amount_of_neighbours);

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

				CompDend(&cellParamsPtr[i], 1);
				CompSoma(&cellParamsPtr[i]);
				CompAxon(&cellParamsPtr[i]);
			}

		}
	}
	

	/* start of the simulation
	 * In case we want to read the stimulus from file inputFromFile = true
	 */

	if(inputFromFile){
		simStep = 0;

		/* Read full lines until end of file. 
		 * Every iteration (line) is one simulation step.
		 */
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
		simTime = SIMTIME; 
		total_simulation_steps = (int)(simTime/DELTA);
				
		for(simStep=0;simStep<total_simulation_steps;simStep++) {
								
			simulation_array_ID = simStep%2;
			if ((simStep>=20000)&&(simStep<20500-1))
			//if ((simStep>=20000)&&((simStep%20000)>=0)&&((simStep%20000)<500))		//re-occuring spikes after first sec
				iApp = 6; 
			else
		   		iApp = 0;
			
			if (PRINTING) {
				sprintf(tempbuf, " %d %.2f ", simStep+1, iApp);
				fputs(tempbuf, pOutFile);
			}
			
		#pragma omp parallel for shared(cellParamsPtr, cellPtr, pOutFile) private(iApp, target_cell, tempbuf, i, requested_neighbour) firstprivate(simulation_array_ID)
			for (target_cell=0;target_cell<cellCount;target_cell++) {
				/* we simulate a hardcoded input pulse here 
				 * that differs from step to step 
				 */

			#pragma ivdep
				for (i=0; i<cellParamsPtr[target_cell].total_amount_of_neighbours; i++){
                                        requested_neighbour = cellParamsPtr[target_cell].neighId[i];
					voltage = cellPtr[simulation_array_ID][requested_neighbour].dend.V_dend;
                                        cellParamsPtr[target_cell].neighVdend[i] = voltage;
//					cellParamsPtr[target_cell].neighVdend[i] = cellPtr[simulation_array_ID][cellParamsPtr[target_cell].neighId[i]].dend.V_dend;
                                }

				cellParamsPtr[target_cell].iAppIn = iApp;			
				cellParamsPtr[target_cell].prevCellState = &cellPtr[simulation_array_ID][target_cell];
				cellParamsPtr[target_cell].newCellState = &cellPtr[simulation_array_ID^1][target_cell];
	
				CompDend(&cellParamsPtr[target_cell], 0);
				CompSoma(&cellParamsPtr[target_cell]);
				CompAxon(&cellParamsPtr[target_cell]);
/*				
				if (PRINTING) {
					sprintf(tempbuf, "%d : %.8f ", target_cell+1, cellPtr[(simulation_array_ID)^1][target_cell].axon.V_axon);
					fputs(tempbuf, pOutFile);
				}

	*/		}
			
			if (PRINTING) {
				sprintf(tempbuf, "\n");
				fputs(tempbuf, pOutFile);
			}
		
		}

		simStep = total_simulation_steps;//so that simStep in the end has the exact value of how many steps we had in this sim, regardless of input method (useful to know which cellPtr has what)
	}
		
	/* Free  memory and close files
	 */

	if (PRINTSTATE) {
		char resultFileName[50];
		sprintf(resultFileName, "lastStateDump.txt");
		printState(cellPtr[simStep%2], resultFileName);		//simStep%2 here should refer to the cellPtr which has the last state of the network that we calculated
	}

//	free(cellPtr[0]);
//	free(cellPtr[1]);
//	free(cellPtr);
//	free(cellParamsPtr);
	/* reminder, maybe free communication list */
	
	
	if (PRINTING) {
		fclose (pOutFile);
		chmod(outFileName,0x01B6);
	}

	if(inputFromFile)
		fclose (pInFile);

	return 0;
}




//DENDRITIC COMPUTATIONAL PART ------------------

void CompDend(cellCompParams *cellParamsPtr, int randomness){

	struct channelParams chPrms;
	struct dendCurrVoltPrms chComps;

	//printf("Dendrite ");

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->dend.V_dend;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->dend.Hcurrent_q;
	chPrms.newComp1 = &cellParamsPtr->newCellState->dend.Hcurrent_q;
	//Compute
	DendHCurr(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->dend.V_dend;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->dend.Calcium_r;
	chPrms.newComp1 = &cellParamsPtr->newCellState->dend.Calcium_r;
	//Compute
	DendCaCurr(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->dend.Potassium_s;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->dend.Ca2Plus;
	chPrms.newComp1 = &cellParamsPtr->newCellState->dend.Potassium_s;
	//Compute
	DendKCurr(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->dend.Ca2Plus;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->dend.I_CaH;
	chPrms.newComp1 = &cellParamsPtr->newCellState->dend.Ca2Plus;
	//Compute
	DendCal(&chPrms);

	/* ANDREAS change here the last parameter of IcNeighbors, we no longer use number_of_neighbours 
	 * instead we have an index in struct cellParameters which tells us where the cell's neighbouring
	 * voltages end in the array
	 */

	//in random initialization mode, cells run in closed circuit mode so neighboring current equals zero
	if (randomness==1)
		chComps.iC = 0;
	else
		chComps.iC = IcNeighbors(cellParamsPtr->neighVdend, cellParamsPtr->neighConductances, cellParamsPtr->prevCellState->dend.V_dend, cellParamsPtr->total_amount_of_neighbours);
	
	chComps.iApp = &cellParamsPtr->iAppIn;
	chComps.vDend = &cellParamsPtr->prevCellState->dend.V_dend;
	chComps.newVDend = &cellParamsPtr->newCellState->dend.V_dend;
	chComps.vSoma = &cellParamsPtr->prevCellState->soma.V_soma;
	chComps.q = &cellParamsPtr->newCellState->dend.Hcurrent_q;
	chComps.r = &cellParamsPtr->newCellState->dend.Calcium_r;
	chComps.s = &cellParamsPtr->newCellState->dend.Potassium_s;
	chComps.newI_CaH = &cellParamsPtr->newCellState->dend.I_CaH;
	DendCurrVolt(&chComps);

	return;
}

void DendHCurr(struct channelParams *chPrms){

	mod_prec q_inf, tau_q, dq_dt, q_local;

	//Get inputs
	mod_prec prevV_dend = *chPrms->v;
	mod_prec prevHcurrent_q = *chPrms->prevComp1;

	// Update dendritic H current component
	q_inf = 1 /(1 + exp((prevV_dend + 80) / 4));
	tau_q = 1 /(exp(-0.086 * prevV_dend - 14.6) + exp(0.070 * prevV_dend - 1.87));
	dq_dt = (q_inf - prevHcurrent_q) / tau_q;
	q_local = DELTA * dq_dt + prevHcurrent_q;
	//Put result
	*chPrms->newComp1 = q_local;

	return;
}

void DendCaCurr(struct channelParams *chPrms){

	mod_prec alpha_r, beta_r, r_inf, tau_r, dr_dt, r_local;

	//Get inputs
	mod_prec prevV_dend = *chPrms->v;
	mod_prec prevCalcium_r = *chPrms->prevComp1;

	// Update dendritic high-threshold Ca current component
	alpha_r = 1.7 / (1 + exp( -(prevV_dend - 5) / 13.9));
	beta_r = 0.02 * (prevV_dend + 8.5) / (exp((prevV_dend + 8.5) / 5) - 1);
	r_inf = alpha_r / (alpha_r + beta_r);
	tau_r = 5 / (alpha_r + beta_r);
	dr_dt = (r_inf - prevCalcium_r) / tau_r;
	r_local = DELTA * dr_dt + prevCalcium_r;
	//Put result
	*chPrms->newComp1 = r_local;

	return;
}

void DendKCurr(struct channelParams *chPrms){

	mod_prec alpha_s, beta_s, s_inf, tau_s, ds_dt, s_local;

	//Get inputs
	mod_prec prevPotassium_s = *chPrms->prevComp1;
	mod_prec prevCa2Plus = *chPrms->prevComp2;

	// Update dendritic Ca-dependent K current component
	alpha_s = min((0.00002*prevCa2Plus), 0.01);
	beta_s = 0.015;
	s_inf = alpha_s / (alpha_s + beta_s);
	tau_s = 1 / (alpha_s + beta_s);
	ds_dt = (s_inf - prevPotassium_s) / tau_s;
	s_local = DELTA * ds_dt + prevPotassium_s;
	//Put result
	*chPrms->newComp1 = s_local;

	return;
}

//Consider merging DendCal into DendKCurr since DendCal's output doesn't go to DendCurrVolt but to DendKCurr
void DendCal(struct channelParams *chPrms){

	mod_prec dCa_dt, Ca2Plus_local;

	//Get inputs
	mod_prec prevCa2Plus = *chPrms->prevComp1;
	mod_prec prevI_CaH = *chPrms->prevComp2;

	// update Calcium concentration
	dCa_dt = -3 * prevI_CaH - 0.075 * prevCa2Plus;
	Ca2Plus_local = DELTA * dCa_dt + prevCa2Plus;
	//Put result
	*chPrms->newComp1 = Ca2Plus_local;//This state value is read in DendKCurr

	return;
}

void DendCurrVolt(struct dendCurrVoltPrms *chComps){

	//Loca variables
	mod_prec I_sd, I_CaH, I_K_Ca, I_ld, I_h, dVd_dt;

	//Get inputs
	mod_prec I_c = chComps->iC;
	mod_prec I_app = *chComps->iApp;
	mod_prec prevV_dend = *chComps->vDend;
	mod_prec prevV_soma = *chComps->vSoma;
	mod_prec q = *chComps->q;
	mod_prec r = *chComps->r;
	mod_prec s = *chComps->s;

	// DENDRITIC CURRENTS

	// Soma-dendrite interaction current I_sd
	I_sd = (G_INT / (1 - P1)) * (prevV_dend - prevV_soma);
	// Inward high-threshold Ca current I_CaH
	I_CaH=G_CAH * r * r * (prevV_dend - V_CA);
	// Outward Ca-dependent K current I_K_Ca
	I_K_Ca =G_K_CA * s * (prevV_dend - V_K);
	// Leakage current I_ld
	I_ld =G_LD * (prevV_dend - V_L);
	// Inward anomalous rectifier I_h
	I_h=G_H * q * (prevV_dend - V_H);

	dVd_dt = (-(I_CaH + I_sd+ I_ld + I_K_Ca + I_c + I_h) + I_app) / C_M;

	//Put result (update V_dend)
	*chComps->newVDend = DELTA * dVd_dt + prevV_dend;
	*chComps->newI_CaH = I_CaH;//This is a state value read in DendCal
	return;
}

mod_prec IcNeighbors(mod_prec *neighVdend, mod_prec *neighConductances, mod_prec prevV_dend, int neighbors){

	int i;
	mod_prec f, V, Cond, I_c=0;
	for(i=0;i<neighbors;i++){
		V = prevV_dend - neighVdend[i];
		f = 0.8 * exp(-1*pow(V, 2)/100) + 0.2;// SCHWEIGHOFER 2004 VERSION
		Cond = neighConductances[i];
		I_c = I_c + (Cond * f * V);
//		V = 1+1;
//		f = 3+3;
	}

	return I_c;
}

//SOMATIC COMPUTATIONAL PART -----------------

void CompSoma(cellCompParams *cellParamsPtr){

	struct channelParams chPrms;
	struct somaCurrVoltPrms chComps;

	// update somatic components
	// SCHWEIGHOFER:

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->soma.V_soma;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->soma.Calcium_k;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->soma.Calcium_l;
	chPrms.newComp1 = &cellParamsPtr->newCellState->soma.Calcium_k;
	chPrms.newComp2 = &cellParamsPtr->newCellState->soma.Calcium_l;
	//Compute
	SomaCalcium(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->soma.V_soma;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->soma.Sodium_m;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->soma.Sodium_h;
	chPrms.newComp1 = &cellParamsPtr->newCellState->soma.Sodium_m;
	chPrms.newComp2 = &cellParamsPtr->newCellState->soma.Sodium_h;
	//Compute
	SomaSodium(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->soma.V_soma;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->soma.Potassium_n;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->soma.Potassium_p;
	chPrms.newComp1 = &cellParamsPtr->newCellState->soma.Potassium_n;
	chPrms.newComp2 = &cellParamsPtr->newCellState->soma.Potassium_p;
	//Compute
	SomaPotassium(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->soma.V_soma;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->soma.Potassium_x_s;
	chPrms.newComp1 = &cellParamsPtr->newCellState->soma.Potassium_x_s;
	//Compute
	SomaPotassiumX(&chPrms);

	chComps.g_CaL = &cellParamsPtr->prevCellState->soma.g_CaL;
	chComps.vDend = &cellParamsPtr->prevCellState->dend.V_dend;
	chComps.vSoma = &cellParamsPtr->prevCellState->soma.V_soma;
	chComps.newVSoma = &cellParamsPtr->newCellState->soma.V_soma;
	chComps.vAxon = &cellParamsPtr->prevCellState->axon.V_axon;
	chComps.k = &cellParamsPtr->newCellState->soma.Calcium_k;
	chComps.l = &cellParamsPtr->newCellState->soma.Calcium_l;
	chComps.m = &cellParamsPtr->newCellState->soma.Sodium_m;
	chComps.h = &cellParamsPtr->newCellState->soma.Sodium_h;
	chComps.n = &cellParamsPtr->newCellState->soma.Potassium_n;
	chComps.x_s = &cellParamsPtr->newCellState->soma.Potassium_x_s;
	SomaCurrVolt(&chComps);

	return;
}

void SomaCalcium(struct channelParams *chPrms){

	mod_prec k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt, k_local, l_local;

	//Get inputs
	mod_prec prevV_soma = *chPrms->v;
	mod_prec prevCalcium_k = *chPrms->prevComp1;
	mod_prec prevCalcium_l = *chPrms->prevComp2;

	k_inf = (1 / (1 + exp(-1 * (prevV_soma + 61) / 4.2)));
	l_inf = (1 / (1 + exp(( prevV_soma + 85.5) / 8.5)));
	tau_k = 1;
	tau_l = ((20 * exp((prevV_soma + 160) / 30) / (1 + exp((prevV_soma + 84) / 7.3))) +35);
	dk_dt = (k_inf - prevCalcium_k) / tau_k;
	dl_dt = (l_inf - prevCalcium_l) / tau_l;
	k_local = DELTA * dk_dt + prevCalcium_k;
	l_local = DELTA * dl_dt + prevCalcium_l;
	//Put result
	*chPrms->newComp1= k_local;
	*chPrms->newComp2= l_local;

	return;
}

void SomaSodium(struct channelParams *chPrms){

	mod_prec m_inf, h_inf, tau_h, dh_dt, m_local, h_local;

	//Get inputs
	mod_prec prevV_soma = *chPrms->v;
	//mod_prec prevSodium_m = *chPrms->prevComp1;
	mod_prec prevSodium_h = *chPrms->prevComp2;

	// RAT THALAMOCORTICAL SODIUM:
	m_inf = 1 / (1 + (exp((-30 - prevV_soma)/ 5.5)));
	h_inf = 1 / (1 + (exp((-70 - prevV_soma)/-5.8)));
	tau_h = 3 * exp((-40 - prevV_soma)/33);
	dh_dt = (h_inf - prevSodium_h)/tau_h;
	m_local = m_inf;
	h_local = prevSodium_h + DELTA * dh_dt;
	//Put result
	*chPrms->newComp1 = m_local;
	*chPrms->newComp2 = h_local;

	return;
}

void SomaPotassium(struct channelParams *chPrms){

	mod_prec n_inf, p_inf, tau_n, tau_p, dn_dt, dp_dt, n_local, p_local;

	//Get inputs
	mod_prec prevV_soma = *chPrms->v;
	mod_prec prevPotassium_n = *chPrms->prevComp1;
	mod_prec prevPotassium_p = *chPrms->prevComp2;

	// NEOCORTICAL
	n_inf = 1 / (1 + exp( ( -3 - prevV_soma) /10));
	p_inf = 1 / (1 + exp( (-51 - prevV_soma) / -12));
	tau_n = 5 + (47 * exp( -(-50 - prevV_soma) /900));
	tau_p = tau_n;
	dn_dt = (n_inf - prevPotassium_n) / tau_n;
	dp_dt = (p_inf - prevPotassium_p) / tau_p;
	n_local = DELTA * dn_dt + prevPotassium_n;
	p_local = DELTA * dp_dt + prevPotassium_p;
	//Put result
	*chPrms->newComp1 = n_local;
	*chPrms->newComp2 = p_local;

	return;
}

void SomaPotassiumX(struct channelParams *chPrms){

	mod_prec alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s, x_s_local;

	//Get inputs
	mod_prec prevV_soma = *chPrms->v;
	mod_prec prevPotassium_x_s = *chPrms->prevComp1;

	// Voltage-dependent (fast) potassium
	alpha_x_s = 0.13 * (prevV_soma + 25) / (1 - exp(-(prevV_soma + 25) / 10));
	beta_x_s= 1.69 * exp(-0.0125 * (prevV_soma + 35));
	x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
	tau_x_s = 1 / (alpha_x_s + beta_x_s);
	dx_dt_s = (x_inf_s - prevPotassium_x_s) / tau_x_s;
	x_s_local = 0.05 * dx_dt_s + prevPotassium_x_s;
	//Put result
	*chPrms->newComp1 = x_s_local;

	return;
}

void SomaCurrVolt(struct somaCurrVoltPrms *chComps){

	//Local variables
	mod_prec I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt;

	//Get inputs
	mod_prec g_CaL = *chComps->g_CaL;
	mod_prec prevV_dend = *chComps->vDend;
	mod_prec prevV_soma = *chComps->vSoma;
	mod_prec prevV_axon = *chComps->vAxon;
	mod_prec k = *chComps->k;
	mod_prec l = *chComps->l;
	mod_prec m = *chComps->m;
	mod_prec h = *chComps->h;
	mod_prec n = *chComps->n;
	mod_prec x_s = *chComps->x_s;

	// SOMATIC CURRENTS

	// Dendrite-soma interaction current I_ds
	I_ds= (G_INT / P1) * (prevV_soma - prevV_dend);
	// Inward low-threshold Ca current I_CaL
	I_CaL = g_CaL * k * k * k * l * (prevV_soma - V_CA); //k^3
	// Inward Na current I_Na_s
	I_Na_s= G_NA_S * m * m * m * h * (prevV_soma - V_NA);
	// Leakage current I_ls
	I_ls= G_LS * (prevV_soma - V_L);
	// Outward delayed potassium current I_Kdr
	I_Kdr_s = G_KDR_S * n * n * n * n * (prevV_soma - V_K); // SCHWEIGHOFER
	// I_K_s
	I_K_s = G_K_S * pow(x_s, 4) * (prevV_soma - V_K);
	// Axon-soma interaction current I_as
	I_as= (G_INT / (1 - P2)) * (prevV_soma - prevV_axon);

	dVs_dt = (-(I_CaL + I_ds+ I_as + I_Na_s + I_ls + I_Kdr_s + I_K_s)) / C_M;
	*chComps->newVSoma = DELTA * dVs_dt + prevV_soma;

	return;
}

//AXONAL COMPUTATIONAL PART -------------------

void CompAxon(cellCompParams *cellParamsPtr){

	struct channelParams chPrms;
	struct axonCurrVoltPrms chComps;

	// update somatic components
	// SCHWEIGHOFER:

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->axon.V_axon;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->axon.Sodium_h_a;
	chPrms.newComp1 = &cellParamsPtr->newCellState->axon.Sodium_h_a;
	chPrms.newComp2 = &cellParamsPtr->newCellState->axon.Sodium_m_a;
	//Compute
	AxonSodium(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->axon.V_axon;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->axon.Potassium_x_a;
	chPrms.newComp1 = &cellParamsPtr->newCellState->axon.Potassium_x_a;
	//Compute
	AxonPotassium(&chPrms);

	//Get inputs
	chComps.vSoma = &cellParamsPtr->prevCellState->soma.V_soma;
	chComps.vAxon = &cellParamsPtr->prevCellState->axon.V_axon;
	chComps.newVAxon = &cellParamsPtr->newCellState->axon.V_axon;
	chComps.m_a = &cellParamsPtr->newCellState->axon.Sodium_m_a;
	chComps.h_a = &cellParamsPtr->newCellState->axon.Sodium_h_a;
	chComps.x_a = &cellParamsPtr->newCellState->axon.Potassium_x_a;
	AxonCurrVolt(&chComps);

	return;
}

void AxonSodium(struct channelParams *chPrms){

	mod_prec m_inf_a, h_inf_a, tau_h_a, dh_dt_a, m_a_local, h_a_local;

	//Get inputs
	mod_prec prevV_axon = *chPrms->v;
	mod_prec prevSodium_h_a = *chPrms->prevComp1;

	// Update axonal Na components
	// NOTE: current has shortened inactivation to account for high
	// firing frequencies in axon hillock
	m_inf_a = 1 / (1 + (exp((-30 - prevV_axon)/ 5.5)));
	h_inf_a = 1 / (1 + (exp((-60 - prevV_axon)/-5.8)));
	tau_h_a = 1.5 * exp((-40 - prevV_axon)/33);
	dh_dt_a = (h_inf_a - prevSodium_h_a)/tau_h_a;
	m_a_local = m_inf_a;
	h_a_local = prevSodium_h_a + DELTA * dh_dt_a;
	//Put result
	*chPrms->newComp1 = h_a_local;
	*chPrms->newComp2 = m_a_local;

	return;
}

void AxonPotassium(struct channelParams *chPrms){

	mod_prec alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a, x_a_local;

	//Get inputs
	mod_prec prevV_axon = *chPrms->v;
	mod_prec prevPotassium_x_a = *chPrms->prevComp1;

	// D'ANGELO 2001 -- Voltage-dependent potassium
	alpha_x_a = 0.13 * (prevV_axon + 25) / (1 - exp(-(prevV_axon + 25) / 10));
	beta_x_a= 1.69 * exp(-0.0125 * (prevV_axon + 35));
	x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a);
	tau_x_a = 1 / (alpha_x_a + beta_x_a);
	dx_dt_a = (x_inf_a - prevPotassium_x_a) / tau_x_a;
	x_a_local = 0.05 * dx_dt_a + prevPotassium_x_a;
	//Put result
	*chPrms->newComp1 = x_a_local;

	return;
}

void AxonCurrVolt(struct axonCurrVoltPrms *chComps){

	//Local variable
	mod_prec I_Na_a, I_la, I_sa, I_K_a, dVa_dt;

	//Get inputs
	mod_prec prevV_soma = *chComps->vSoma;
	mod_prec prevV_axon = *chComps->vAxon;
	mod_prec m_a = *chComps->m_a;
	mod_prec h_a = *chComps->h_a;
	mod_prec x_a = *chComps->x_a;

	// AXONAL CURRENTS
	// Sodium
	I_Na_a= G_NA_A* m_a * m_a * m_a * h_a * (prevV_axon - V_NA);
	// Leak
	I_la= G_LA* (prevV_axon - V_L);
	// Soma-axon interaction current I_sa
	I_sa= (G_INT / P2) * (prevV_axon - prevV_soma);
	// Potassium (transient)
	I_K_a = G_K_A * pow(x_a, 4) * (prevV_axon - V_K);
	dVa_dt = (-(I_K_a + I_sa + I_la + I_Na_a)) / C_M;
	*chComps->newVAxon = DELTA * dVa_dt + prevV_axon;

	return;
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

//	gettimeofday(&tic, NULL);
	char c= fpeek(pInFile);
	if (c==EOF)
		return 0;
	
	char *strNumber;
//	int bufSize = cores*cellCount*20;	//more than enough but will do for nao
	int bufSize = cellCount*20;	//more than enough but will do for nao
	char* buffer = (char*) malloc(bufSize*sizeof(char));

//	int floats_ignored = core_id*cellCount;
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
//		gettimeofday(&toc, NULL);
//		subtract_and_add(&intime, &tic, &toc);
		return 1;//success
	}else{
		if(!feof(pInFile)){
			printf("Error: Reading from input file didn't finish successfully\n");
			exit(EXIT_FAILURE);
		}
		free(buffer);
//		gettimeofday(&toc, NULL);
  //              subtract_and_add(&intime, &tic, &toc);
		return 0;//end of file
	}
}

void read_g_CaL_from_file(cellState* cellPtr) {

	int i;

	FILE* fd = fopen("gcal_file.txt","r");
	for (i=0;i<cellCount;i++) 
		fscanf(fd, "%lf ", &cellPtr[i].soma.g_CaL);
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
			fprintf(fd, "%.8lf ", s[j]);
		fprintf(fd, "\n");

	}

	fclose(fd);
	return;

}

inline mod_prec min(mod_prec a, mod_prec b){

	return (a < b) ? a : b;
}

//sending_node* Make_Core_Communication_List_newest_format (char *filename, cellCompParams *core_cells) {

/*
 * Fair Warning: This function needs heavy testing and is currently being designed as a way to use conductivity matrixes, which seem to be mandatory now, instead of the old
 * communication format. A conductivity matrix basically describes that cell x in rows sends to cell y in columns via a conductivity value z = conductivity_matrix[x][y].
 * The matrix does not have to be symmetrical. Actually the only limit assumed here is that conductivity_matrix[i][i]==0 (no self-feeding, that would make no sense).
 * This function builds a sending list and the necessary receiving buffers of the cellCompParams structure. The new sending list should behave much better since now
 * all cells that need to be sent from one core to the other are BUNDLED TOGETHER, so that we can send all of them in one go (and not sending
 * duplicates of the same cell to the same core either). Read below for functions completing communication based on the structures designed here.
*/
/*	FILE *input_file = fopen(filename,"r");
	sending_node *temporary = NULL, *list_head = NULL;
	int i, k, my_cell, my_sending_cell, core_to_send, line_counter, my_cell_lower_bound = cellCount*core_id, my_cell_upper_bound = my_cell_lower_bound + cellCount;
	int* my_neighbour_count = (int*) calloc(cellCount, sizeof(int));
	mod_prec cond_value;

	//this list has one node for every core employed, detailing which cells must be sent to each core. The node that corresponds to this core will be kept blank
	list_head = (sending_node*) malloc(1*sizeof(sending_node));
	list_head->target_core = 0;
	list_head->next = NULL;
	list_head->total_cells_to_send = 0;
	list_head->cells_to_send = NULL;
	list_head->voltages_to_send = NULL;
	temporary = list_head;

	for (i=1;i<cores;i++) {
		temporary->next = (sending_node*) malloc(1*sizeof(sending_node));
		temporary = temporary->next;
		temporary->target_core = i;
		temporary->next = NULL;
            	temporary->total_cells_to_send = 0;
               	temporary->cells_to_send = NULL;
		temporary->voltages_to_send = NULL;
	}

	//it is important to keep in mind that the matrix is read as if cell in ROW is sending information to cell in COLUMN!

	for (line_counter=0;line_counter<IO_NETWORK_SIZE;line_counter++) {

		temporary = list_head;				//before examining every (supposedly sending) line, we reset temporary to the beginning of our sending list
		for (i=0; i<IO_NETWORK_SIZE; i++) {

			fscanf(input_file, "%lf ", &cond_value);
			if (cond_value==0)						//this connection is considered not existing if conductance = 0
				;
			else {

			//part of the code handling RECEIVING and noting which of my cells needs input from which other cells, from ANY core

				if ((i>=my_cell_lower_bound)&&(i<my_cell_upper_bound)) {				//these are this core's cells' COLUMNS (incoming)
					my_cell = i%cellCount;
					if (my_neighbour_count[my_cell]==0) {                  //if this is the first neighbour, initialize buffers
						core_cells[my_cell].neighVdend = NULL;
						core_cells[my_cell].neighConductances = NULL;
						core_cells[my_cell].neighId = NULL;
					}

					core_cells[my_cell].neighId = allocate_space_int(core_cells[my_cell].neighId, my_neighbour_count[my_cell]);
					core_cells[my_cell].neighId[my_neighbour_count[my_cell]] = line_counter;		//which cell sends this voltage to us (GLOBAL ID)

					core_cells[my_cell].neighConductances = allocate_space_mod(core_cells[my_cell].neighConductances, my_neighbour_count[my_cell]);
					core_cells[my_cell].neighConductances[my_neighbour_count[my_cell]] = cond_value;	//what conductance we use to calculate its impact

					//allocate space for storing this voltage
					core_cells[my_cell].neighVdend = allocate_space_mod(core_cells[my_cell].neighVdend, my_neighbour_count[my_cell]);

					my_neighbour_count[my_cell]++;						//how many neighbours this cell has so far (from ANY core)
					core_cells[my_cell].total_amount_of_neighbours = my_neighbour_count[my_cell];
				}

			//part of the code handling SENDING and bundling together my cells I need to send each core
			
			//NOT NEEDED
				if ((line_counter>=my_cell_lower_bound)&&(line_counter<my_cell_upper_bound)) {		//these are this core's cells' ROWS (outgoing)
					core_to_send = i/cellCount;							//which core I need to send to
					while ((temporary!=NULL)&&(temporary->target_core!=core_to_send))		//search for the right list node
						temporary = temporary->next;

					if (core_to_send==core_id)							//Obviously I do not send to my own core
						;
					else if (temporary!=NULL) {
						k = temporary->total_cells_to_send;
						my_sending_cell = line_counter%cellCount;

						if ((k>0)&&(temporary->cells_to_send[k-1]==my_sending_cell))		//I have already marked this cell to send core_to_send
							;
						else {
							//add this cell to the core_to_send's cells_to_send-array ->WITH LOCAL ID<-
							temporary->cells_to_send = allocate_space_int(temporary->cells_to_send, temporary->total_cells_to_send);
							temporary->cells_to_send[k] = my_sending_cell;
							temporary->voltages_to_send = allocate_space_mod(temporary->voltages_to_send, temporary->total_cells_to_send);
							temporary->total_cells_to_send++;
						}
					}
				}

			//end of the code handling proper list and buffer creation concerning communication

			}
		}

	}

//	for (i=0; i<cellCount; i++)
//		printf("%d) [%d] = (%d-%lf,%d-%lf,%d-%lf,%d-%lf)\n", core_id, i, core_cells[i].neighId[0], core_cells[i].neighConductances[0], core_cells[i].neighId[1], core_cells[i].neighConductances[1], core_cells[i].neighId[2], core_cells[i].neighConductances[2], core_cells[i].neighId[3], core_cells[i].neighConductances[3]);

	fclose(input_file);
	return list_head;

}*/


//NOT NEEDED

//receiving_node* reckon_phase(sending_node *sending_list_head, cellCompParams *core_cells) {

/* In this function, we will make the opposite of the sending_list: We will create a list detailing what
 * information we will receive from each other core. Since cell mapping and connections do not change,
 * we will create this list once in the beginning so that during simulation, only voltages are to be exchanged
 * This list will be created via some information exchange between cores (although it could be built during
 * the conductivity matrix parsing, this seems simpler)
 */
/*	receiving_node *r_temporary = NULL, *receiving_list_head = NULL; 
	sending_node *s_temporary = NULL;
        int i, k, targetCore;
	MPI_Request s_request;	//irrelevant
	MPI_Request* request1 = (MPI_Request*) malloc(cores*sizeof(MPI_Request));
	MPI_Request* request2 = (MPI_Request*) malloc(cores*sizeof(MPI_Request));

	receiving_list_head = (receiving_node*) malloc(1*sizeof(receiving_node));
	receiving_list_head->target_core = 0;
	receiving_list_head->next = NULL;
	receiving_list_head->total_cells_to_receive = 0;
	receiving_list_head->cells_to_receive = NULL;
	receiving_list_head->voltages_to_receive = NULL;
	r_temporary = receiving_list_head;

	for (i=1;i<cores;i++) {
		r_temporary->next = (receiving_node*) malloc(1*sizeof(receiving_node));
		r_temporary = r_temporary->next;
		r_temporary->target_core = i;
		r_temporary->next = NULL;
		r_temporary->total_cells_to_receive = 0;
		r_temporary->cells_to_receive = NULL;
		r_temporary->voltages_to_receive = NULL;
	}

	s_temporary = sending_list_head;
	r_temporary = receiving_list_head;

	//phase 1: inform cores how many cells need to be exchanged

	k = 0;
	while ((r_temporary!=NULL)&&(s_temporary!=NULL)) {
		targetCore = r_temporary->target_core;
		if (targetCore == core_id)
			request1[k] = MPI_REQUEST_NULL;
		else {
			MPI_Isend(&(s_temporary->total_cells_to_send), 1, MPI_INT, targetCore, 0, MPI_COMM_WORLD, &s_request);
			MPI_Irecv(&(r_temporary->total_cells_to_receive), 1, MPI_INT, targetCore, 0, MPI_COMM_WORLD, &request1[k]);
		}
		r_temporary = r_temporary->next;
		s_temporary = s_temporary->next;
		k++;
	}

	syncing(request1);		//phase 1 complete, waiting for sync

	// phase 2: allocate buffers to hold the incoming cells and exchange cell ids to know who sends what
	// WARNING: the receiving buffer here gets filled with ids as the sending core sends them - they are NOT global ids
	// and thus will NEED "translation"

	s_temporary = sending_list_head;
        r_temporary = receiving_list_head;

	k = 0;
	while ((r_temporary!=NULL)&&(s_temporary!=NULL)) {
                targetCore = r_temporary->target_core;
                if (targetCore == core_id)
                        request2[k] = MPI_REQUEST_NULL;
                else {

			if (s_temporary->total_cells_to_send != 0)
				MPI_Isend(&(s_temporary->cells_to_send[0]), s_temporary->total_cells_to_send, MPI_INT, targetCore, 0, MPI_COMM_WORLD, &s_request);

			if (r_temporary->total_cells_to_receive != 0) {
				r_temporary->cells_to_receive = (int*) malloc(r_temporary->total_cells_to_receive * sizeof(int));
				r_temporary->voltages_to_receive = (mod_prec*) malloc(r_temporary->total_cells_to_receive * sizeof(mod_prec));
				MPI_Irecv(&(r_temporary->cells_to_receive[0]), r_temporary->total_cells_to_receive, MPI_INT, targetCore, 0, MPI_COMM_WORLD, &request2[k]);
			} else {
				request2[k] = MPI_REQUEST_NULL;
			}
		}
                r_temporary = r_temporary->next;
                s_temporary = s_temporary->next;
		k++;
        }

	syncing(request2);		//phase 2 complete, waiting for sync

	return receiving_list_head;

}*/

void perform_communication_step(cellCompParams* params, cellState* cells) {

/* The function where the magic happens: using the structures created during initialization, exchange necessary dendritic voltage
 * between all cores. This time, we use non-blocking functions and voltages are bundled together per core-target (and this time,
 * no duplicates when a core needs a voltage for more than one of its own cells).
 */

//	sending_node* s_temp;
//	receiving_node* r_temp;
	int* processing_pointer = (int*) calloc(cellCount, sizeof(int));
	int i, j, k=0, cell_id, translated_cell_id, requested_neighbour, my_requested_cell;

	//phase 3: we will now process what we received from all cores in one passing of our structs altogether
	//i check for every core, which of the cells are needed by which of my cells
//	r_temp = receiving_list_head;
//	while (r_temp!=NULL) {
//		if (r_temp->total_cells_to_receive == 0)
//			;
//		else
//			for (i=0; i<r_temp->total_cells_to_receive; i++) {		//we will run through the list we received
//				cell_id = r_temp->cells_to_receive[i];
//				translated_cell_id = cell_id + cellCount*(r_temp->target_core);			//global id of the cell under examination

				for (j=0; j<cellCount; j++) {							//we shall examine each one of our cells to see which one needs it
					if (processing_pointer[j]>=params[j].total_amount_of_neighbours)
						;
					else {
						requested_neighbour = params[j].neighId[processing_pointer[j]];		//id of the next neighbour this cell needs
						//while (((requested_neighbour/cellCount)==core_id)&&(requested_neighbour!=-1)) {			//if the requested neighbour belongs to this core
						while ((requested_neighbour!=-1)) {			//if the requested neighbour belongs to this core
							//this is what i have to do from now on
						//	my_requested_cell = requested_neighbour%cellCount;		//which of my cells is needed
							params[j].neighVdend[processing_pointer[j]] = cells[requested_neighbour].dend.V_dend;	//the voltage of the cell requested

							processing_pointer[j]++;
							if (processing_pointer[j]<params[j].total_amount_of_neighbours)
								requested_neighbour = params[j].neighId[processing_pointer[j]];	//request next neighbour global id
							else
								requested_neighbour = -1;						//if we are done then make sure we do not request any more cells
						}
				//		if (requested_neighbour==translated_cell_id) {						//we found a cell we need
				//			params[j].neighVdend[processing_pointer[j]] = r_temp->voltages_to_receive[i];
				//			processing_pointer[j]++;
				//		}
					}
				}
//			}

//		r_temp = r_temp->next;
//	}

	//phase 4: check whether all cells have filled up their neighVdend buffers (some of them may still need info from our own cells)
	
/*	for (j=0; j<cellCount; j++) {
		if (processing_pointer[j]>=params[j].total_amount_of_neighbours)
			;
		else {
			requested_neighbour = params[j].neighId[processing_pointer[j]];
			while (((requested_neighbour/cellCount)==core_id)&&(requested_neighbour!=-1)) {

				my_requested_cell = requested_neighbour%cellCount;
				params[j].neighVdend[processing_pointer[j]] = cells[my_requested_cell].dend.V_dend;

				processing_pointer[j]++;
				if (processing_pointer[j]<params[j].total_amount_of_neighbours)
					requested_neighbour = params[j].neighId[processing_pointer[j]];
				else
					requested_neighbour = -1;
			}
		}
	}*/
	
	//all tasks have been completed and processed so we can move out of the communication function
	return;

}

int* allocate_space_int(int* pointer, int existing_slots) {

	int new_total_slots = existing_slots + 1;
	int* new_pointer = (int*) realloc(pointer, new_total_slots*sizeof(int));
	return new_pointer;
}

mod_prec* allocate_space_mod(mod_prec* pointer, int existing_slots) {

	int new_total_slots = existing_slots + 1;
        mod_prec* new_pointer = (mod_prec*) realloc(pointer, new_total_slots*sizeof(mod_prec));
	return new_pointer;
}

//void syncing(MPI_Request* request) {

/* a function that first makes sure all requests passed on to it are completed
 * and then calls an MPI_Barrier to make sure all processes are in sync. We expect
 * one request from each core employed by the app
 */

/*	int i, done = 0, flag, no_of_reqs = cores;

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

}*/

