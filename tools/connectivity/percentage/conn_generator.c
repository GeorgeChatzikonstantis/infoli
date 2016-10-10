#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DEF_GAP_JUNCTION 0.04

int main (int argc,char **argv) {

	int neurons, i, j;
	float rndm, percentage;
	char filename[30];
	sprintf(filename,"cellConnections.txt");
	FILE *out_file;
	float **global_connection_table;

	if (argc!=3) {
		printf("Error: Usage ./conn_generator.x #neurons density_pct%\n");
		return 0;
	}

	neurons = atoi(argv[1]);
	percentage = atof(argv[2]);
	if (percentage>1) {
		printf("Warning: Connectivity Density assumed at 100%\n");
		percentage = 1;
	}
	out_file = fopen(filename, "w"); 
	global_connection_table = (float**) malloc(sizeof(float*)*neurons);
	for (i=0; i<neurons; i++) {
		global_connection_table[i] = (float*) calloc(neurons, sizeof(float));
	}
	srand ( time(NULL) );
	
	for (i=0; i<neurons; i++) {
		for (j=0; j<neurons; j++) {
			if (i==j)
				;
			else {
				rndm = ((float) rand()) / ((float) RAND_MAX);
				if (rndm <= percentage)
					global_connection_table[i][j] = DEF_GAP_JUNCTION;
			}
			fprintf(out_file, "%0.2f ", global_connection_table[i][j]);
		}
		fprintf(out_file, "\n");
	}
	fclose(out_file);
	chmod(out_file,0x01B6);
	
	return 0;

}
