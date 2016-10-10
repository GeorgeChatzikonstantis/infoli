#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

int main(int argc, char *argv[]){

	if (argc!=5) {
		printf("Error!\nUsage: ./create_stim network_X network_Y lamda total_time_in_ms.\n");
		return 0;
	}
	int network_x = atoi(argv[1]);			//network dimension 1
	int network_y = atoi(argv[2]);			//network dimension 2
	float lamda = atof(argv[3]);			//how many spikes per simulated second of brain time we expect to have
	int total_steps = (int)((atof(argv[4]))*20);	//experiment duration in steps, delta time is 0.05 so each ms is 20 steps

	char currFileName[100];
	sprintf(currFileName,"currentInput.txt");
	FILE* fd= fopen(currFileName, "w+");
	
	int** stim_matrix = (int**) malloc((sizeof(int*))*network_x);		//this matrix shows whether a stimulus is active or not, set at zero initially
	int i, j, t= 0;
	for (i=0;i<network_x;i++)
		stim_matrix[i] = (int*) calloc(network_y, (sizeof(int)));	//calloc ensures all cells begin with zero stimulus

/*
 * we are going to assign a chance to each cell, based on lamda, on each step to create a stimulus. whenever that is done,
 * we will set the stim_matrix to 500 which signifies 500 steps of active stimulus - the stim_matrix will be reduced by 1 each step if >0
*/

	float chance = lamda / (20*1000);			//lamda spikes per sec, meaning lamda spikes per 20k steps
	srand(time(NULL));
	float random_number;

	while (t<total_steps) {

		for (i=0;i<network_x;i++)
			for (j=0;j<network_y;j++) {

				if (stim_matrix[i][j]==0) {
					random_number = ((float) rand()) / RAND_MAX;
					if (random_number<chance) {
						fprintf(fd, "6.000 ");
						stim_matrix[i][j]=499;
					} else {
						fprintf(fd, "0.000 ");
					}
				} else {
					fprintf(fd, "6.000 ");
					stim_matrix[i][j]--;
				}

			}
		fprintf(fd, "\n");
		t++;
	}

	fclose(fd);
	return 0;

}
