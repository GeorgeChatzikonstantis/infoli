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
#define G_CAL_FROM_FILE 0		//	flag enabling or disabling user defined initial settings of gcal


// Cell properties , biological properties for the cells ( irrevelant for now )
#define DELTA 0.05f // o.05 milli sec = 50 micro sec
//Conductance for neighbors' coupling
#define CONDUCTANCE 0.04f 
// Capacitance
#define C_M 1
// Somatic conductances (mS/cm2)
#define G_NA_S 150      // Na gate conductance (=90 in Schweighofer code, 70 in paper) 120 too little
#define G_KDR_S 9    // K delayed rectifier gate conductance (alternative value: 18)
#define G_K_S 5      // Voltage-dependent (fast) potassium
#define G_LS 0.016f  // Leak conductance (0.015)
// Dendritic conductances (mS/cm2)
#define G_K_CA 35       // Potassium gate conductance (35)
#define G_CAH 4.5f     // High-threshold Ca gate conductance (4.5)
#define G_LD 0.016f   // Dendrite leak conductance (0.015)
#define G_H 0.125f    // H current gate conductance (1.5) (0.15 in SCHWEIGHOFER 2004)
// Axon hillock conductances (mS/cm2)
#define G_NA_A 240      // Na gate conductance (according to literature: 100 to 200 times as big as somatic conductance)
#define G_NA_R 0      // Na (resurgent) gate conductance
#define G_K_A 20      // K voltage-dependent
#define G_LA 0.016f  // Leak conductance
// Cell morphology
#define P1 0.25f        // Cell surface ratio soma/dendrite (0.2)
#define P2 0.15f      // Cell surface ratio axon(hillock)/soma (0.1)
#define G_INT 0.13f       // Cell internal conductance (0.13)
// Reversal potentials
#define V_NA 55       // Na reversal potential (55)
#define V_K -75       // K reversal potential
#define V_CA 120       // Ca reversal potential (120)
#define V_H -43       // H current reversal potential
#define V_L 10       // leak current


/*** TYPEDEFS AND STRUCTS***/
//typedef double mod_prec;
typedef float mod_prec;			//BE VERY CAREFUL TO CHECK ALL DAMNED SCANFS TO BE SURE YOU SCAN FOR SINGLE-POINT ACCURACY, KNOWN ISSUE WITH COND VALUES) AND MPI_TYPES

/* struct that encaptulates the network's 
*  collection of channels and voltages
*/


/* struct that encapsulates the network's
* stable parameters (connection conductances
* and Ids)
*/

typedef __attribute__((align(64))) struct cellCompParams{
	int *total_amount_of_neighbours;
	mod_prec **neighConductances;
	int **neighId;
}cellCompParams;

/*** FUNCTION PROTOTYPES ***/

int ReadFileLine(FILE *, mod_prec *);
void read_g_CaL_from_file(mod_prec *);
inline mod_prec min(mod_prec a, mod_prec b);

#endif /* MAIN_H_ */
