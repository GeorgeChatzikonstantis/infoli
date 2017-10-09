/*
 *
 * Copyright (c) 2012, Neurasmus B.V., The Netherlands,
 * web: www.neurasmus.com email: info@neurasmus.com
 *
 * Any use reproduction in whole or in parts is prohibited
 * without the written consent of the copyright owner.
 *
 * All Rights Reserved.
 *
 *
 * Author: Sebastian Isaza
 * Created: 10-04-2012
 * Modified: 06-06-2012
 *
 * Description : Top header file of the Inferior Olive model. It contains the
 * constant model conductances, the data structures that hold the cell state and
 * the function prototypes.
 *
 */

//#include <mic_power.h>

#ifndef MAIN_H_
#define MAIN_H_
/*** MACROS ***/
#define RAND_INIT 0 // make it zero to facilitate debugging , 0 for debugging / 1 for random states

//IO network size is IO_NETWORK_DIM1*IO_NETWORK_DIM2
#define IAPP_MAX_CHARS 6 //	2 integer, the dot, 2 decimals and the delimiter 
#define PRINTING 1	 //	flag enabling or disabling output , axon's voltage is the output at every step

int core_id, cores, cellCount;
int IO_NETWORK_DIM1, IO_NETWORK_DIM2, IO_NETWORK_SIZE;
float CONN_PROBABILITY;
struct timeval tic, toc;

/*** TYPEDEFS AND STRUCTS***/
//typedef double mod_prec;
typedef float mod_prec;                 //BE VERY CAREFUL TO CHECK ALL DAMNED SCANFS TO BE SURE YOU SCAN FOR SINGLE-POINT ACCURACY, KNOWN ISSUE WITH COND VALUES) AND MPI_TYPES

// Cell properties, biological properties for the cells
mod_prec DELTA;    // 0.05 milli sec = 50 micro sec
//Conductance for neighbors' coupling
mod_prec CONDUCTANCE;     
// Capacitance
mod_prec C_M;
// Somatic conductances (mS/cm2)
mod_prec G_NA_S;        // Na gate conductance (=90 in Schweighofer code, 70 in paper) 120 too little
mod_prec G_KDR_S;       // K delayed rectifier gate conductance (alternative value: 18)
mod_prec G_K_S; // Voltage-dependent (fast) potassium
mod_prec G_LS;  // Leak conductance (0.015)
// Dendritic conductances (mS/cm2)
mod_prec G_K_CA;        // Potassium gate conductance (35)
mod_prec G_CAH; // High-threshold Ca gate conductance (4.5)
mod_prec G_LD;  // Dendrite leak conductance (0.015)
mod_prec G_H;   // H current gate conductance (1.5) (0.15 in SCHWEIGHOFER 2004)
// Axon hillock conductances (mS/cm2)
mod_prec G_NA_A;// Na gate conductance (according to literature: 100 to 200 times as big as somatic conductance)
mod_prec G_NA_R;        // Na (resurgent) gate conductance
mod_prec G_K_A; // K voltage-dependent
mod_prec G_LA;  // Leak conductance
// Cell morphology
mod_prec P1;    // Cell surface ratio soma/dendrite (0.2)
mod_prec P2;    // Cell surface ratio axon(hillock)/soma (0.1)
mod_prec G_INT; // Cell internal conductance (0.13)
// Reversal potentials
mod_prec V_NA;  // Na reversal potential (55)
mod_prec V_K;   // K reversal potential
mod_prec V_CA;  // Ca reversal potential (120)
mod_prec V_H;   // H current reversal potential
mod_prec V_L;   // leak current

/* struct that encapsulates the network's
* stable parameters (connection conductances
* and Ids)
*/

typedef __attribute__((aligned(64))) struct cellCompParams{
	int *total_amount_of_neighbours;
	mod_prec **neighConductances;
	int **neighId;
}cellCompParams;

/*** FUNCTION PROTOTYPES ***/

int ReadFileLine(FILE *, mod_prec *);
void read_g_CaL_from_file(mod_prec *);
inline mod_prec min(mod_prec a, mod_prec b);

#endif /* MAIN_H_ */
