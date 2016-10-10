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


#ifndef MAIN_H_
#define MAIN_H_
/*** MACROS ***/
#define RAND_INIT 0 // make it zero to facilitate debugging , 0 for debugging / 1 for random states
#define SIMTIME 10  // in ms, for when no input file is provided , time in msec when you do not have an input file
//IO network size is IO_NETWORK_DIM1*IO_NETWORK_DIM2
#define IAPP_MAX_CHARS 6 //	2 integer, the dot, 2 decimals and the delimiter 
#define PRINTING 1	 //	flag enabling or disabling output , axon's voltage is the output at every step
#define PRINTSTATE 0	//	flag enabling dumping of the final states of the neurons (search for it in results folder)
#define G_CAL_FROM_FILE 0		//	flag enabling or disabling user defined initial settings of gcal

// Cell properties , biological properties for the cells ( irrevelant for now )
#define DELTA 0.05 // o.05 milli sec = 50 micro sec
//Conductance for neighbors' coupling
#define CONDUCTANCE 0.04 
// Capacitance
#define C_M 1
// Somatic conductances (mS/cm2)
#define G_NA_S 150      // Na gate conductance (=90 in Schweighofer code, 70 in paper) 120 too little
#define G_KDR_S 9.0    // K delayed rectifier gate conductance (alternative value: 18)
#define G_K_S 5      // Voltage-dependent (fast) potassium
#define G_LS 0.016  // Leak conductance (0.015)
// Dendritic conductances (mS/cm2)
#define G_K_CA 35       // Potassium gate conductance (35)
#define G_CAH 4.5     // High-threshold Ca gate conductance (4.5)
#define G_LD 0.016   // Dendrite leak conductance (0.015)
#define G_H 0.125    // H current gate conductance (1.5) (0.15 in SCHWEIGHOFER 2004)
// Axon hillock conductances (mS/cm2)
#define G_NA_A 240      // Na gate conductance (according to literature: 100 to 200 times as big as somatic conductance)
#define G_NA_R 0      // Na (resurgent) gate conductance
#define G_K_A 20      // K voltage-dependent
#define G_LA 0.016  // Leak conductance
// Cell morphology
#define P1 0.25        // Cell surface ratio soma/dendrite (0.2)
#define P2 0.15      // Cell surface ratio axon(hillock)/soma (0.1)
#define G_INT 0.13       // Cell internal conductance (0.13)
// Reversal potentials
#define V_NA 55       // Na reversal potential (55)
#define V_K -75       // K reversal potential
#define V_CA 120       // Ca reversal potential (120)
#define V_H -43       // H current reversal potential
#define V_L 10       // leak current


/*** TYPEDEFS AND STRUCTS***/
typedef double mod_prec;	//BE VERY CAREFUL TO CHECK ALL DAMNED SCANFS TO BE SURE YOU SCAN FOR DOUBLE-POINT ACCURACY, KNOWN ISSUE WITH COND VALUES) AND MPI_TYPES

/* compartment structs 
*
*/

typedef struct dend{
	mod_prec V_dend;
	mod_prec Hcurrent_q;
	mod_prec Calcium_r;
	mod_prec Potassium_s;
	mod_prec I_CaH;
	mod_prec Ca2Plus;
} dend;

typedef struct soma{
	mod_prec g_CaL;
	mod_prec V_soma;
	mod_prec Sodium_m;
	mod_prec Sodium_h;
	mod_prec Calcium_k;
	mod_prec Calcium_l;
	mod_prec Potassium_n;
	mod_prec Potassium_p;
	mod_prec Potassium_x_s;
} soma;


typedef struct axon{
	mod_prec V_axon;
	mod_prec Sodium_m_a;
	mod_prec Sodium_h_a;
	mod_prec Potassium_x_a;
} axon;

/* struct that represents a cell
 *
 */
typedef struct cellState{
	int cellID;
//	int cell_x;
//	int cell_y;
	struct dend dend;
	struct soma soma;
	struct axon axon;
} cellState; 

typedef struct cellCompParams{
	mod_prec iAppIn;
	mod_prec *neighVdend;		
	mod_prec *neighConductances;
	int *neighId;
	cellState *prevCellState;
	cellState *newCellState;
	int index_of_neighVdend;
	int total_amount_of_neighbours;
}cellCompParams;

/* might be a little irrelevant for the time 
 * passing parameters in a nice way
 */

typedef struct channelParams{
	mod_prec *v;
	mod_prec *prevComp1, *prevComp2;
	mod_prec *newComp1, *newComp2;
} channelParams;

typedef struct dendCurrVoltPrms{
	mod_prec *iApp;
	mod_prec iC;
	mod_prec *vDend;
	mod_prec *vSoma;
	mod_prec *q, *r, *s;
	mod_prec *newVDend;
	mod_prec *newI_CaH;
} dendCurrVoltPrms;

typedef struct somaCurrVoltPrms{
	mod_prec *g_CaL;
	mod_prec *vSoma;
	mod_prec *vDend;
	mod_prec *vAxon;
	mod_prec *k, *l, *m, *h, *n, *x_s;
	mod_prec *newVSoma;
} somaCurrVoltPrms;

typedef struct axonCurrVoltPrms{
	mod_prec *vSoma;
	mod_prec *vAxon;
	mod_prec *m_a, *h_a, *x_a;
	mod_prec *newVAxon;
} axonCurrVoltPrms;

/*** FUNCTION PROTOTYPES ***/


void perform_communication_step(cellCompParams *, cellState *);

void CompDend(cellCompParams *, int);
void DendHCurr(struct channelParams *);
void DendCaCurr(struct channelParams *);
void DendKCurr(struct channelParams *);
void DendCal(struct channelParams *);
void DendCurrVolt(struct dendCurrVoltPrms *);
mod_prec IcNeighbors(mod_prec *, mod_prec *, mod_prec, int);

void CompSoma(cellCompParams *);
void SomaCalcium(struct channelParams *);
void SomaSodium(struct channelParams *);
void SomaPotassium(struct channelParams *);
void SomaPotassiumX(struct channelParams *);
void SomaCurrVolt(struct somaCurrVoltPrms *);

void CompAxon(cellCompParams *);
void AxonSodium(channelParams *);
void AxonPotassium(channelParams *);
void AxonCurrVolt(axonCurrVoltPrms *);

void initState(cellState *);
int ReadFileLine(FILE *, mod_prec *);
void read_g_CaL_from_file(cellState *);
void printState(cellState *, char *);
inline mod_prec min(mod_prec a, mod_prec b);
int* allocate_space_int(int *, int);
mod_prec* allocate_space_mod(mod_prec *, int);

#endif /* MAIN_H_ */
