#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

float MTTF;

int weibull_spike_fired() {

	float chance;				//will use this float to calculate chances for a spike to fire
	float random_num;			//random number in [0,1] that compares against chance to decide whether spike is fired
		
	chance = 1 - expf(-(0.05/MTTF));	//we check every step, therefore, dt=1step=0.05 ms
	random_num = ((float) rand()) / RAND_MAX;

	if (random_num < chance)
		return 1;
	else return 0;
	
}

int main(int argc, char *argv[]){

	if (argc!=6) {
		printf("Error!\nUsage: ./stim_generator.x network_X network_Y MTTF_in_ms spike_width experiment_duration_in_ms.\n");
		return 0;
	}

	int network_x = atoi(argv[1]);			//network dimension 1
	int network_y = atoi(argv[2]);			//network dimension 2
	MTTF = atof(argv[3]);				//Mean Time To Fire in ms
	float spike_duration = atof(argv[4]);		//spike duration in ms
	int spike_steps = (int)(spike_duration*20);	//spike duration in steps, delta is assumed at 0.05ms
	float exp_duration = atof(argv[5]);		//experiment duration in ms
	int total_steps = (int)(exp_duration*20);	//experiment duration in steps
	int i, j, t= 0;
	
	char currFileName[100];
	sprintf(currFileName,"currentInput.txt");
	FILE* fd= fopen(currFileName, "w+");

	char fileName2[100];
	sprintf(fileName2,"stimulusReport.txt");
	FILE* fd2= fopen(fileName2, "w+");

	int** stim_matrix = (int**) malloc((sizeof(int*))*network_x);		//this matrix shows whether a stimulus is generated somewhere
	for (i=0;i<network_x;i++)
		stim_matrix[i] = (int*) calloc(network_y, (sizeof(int)));

	int** counter = (int**) malloc((sizeof(int*))*network_x);
	for (i=0;i<network_x;i++)
		counter[i] = (int*) calloc(network_y, (sizeof(int)));
	
	int line_message_size = 6 * network_x * network_y + 2;			//calculation of each line's byte size
	char* line_message_buffer = (char*) malloc(line_message_size*sizeof(char));

	srand(time(NULL));

	while (t < total_steps) {

		strcpy(line_message_buffer, "");	//clean your output buffer
		for (i=0;i<network_x;i++)
			for (j=0;j<network_y;j++) {

				if (stim_matrix[i][j] > 0) {		//ongoing stimulus for this cell
					strcat(line_message_buffer, "6.000 ");
					stim_matrix[i][j]--;		//count down the stim duration
				} else {
					strcat(line_message_buffer, "0.000 ");
					if (weibull_spike_fired()) {		//one check per step for spike
						stim_matrix[i][j]=spike_steps;	//stimulus activation
						counter[i][j]++;
					}
				}
			
			}
		strcat(line_message_buffer, "\n");
		fprintf(fd, "%s", line_message_buffer);
		t++;

	}

	fprintf(fd2, "Spikes fired (MTTF:%.2fms, Duration:%.2fms)\n\n", MTTF, exp_duration);
	for (i=0;i<network_x;i++) {
		for (j=0;j<network_y;j++)
			fprintf(fd2, "%d\t", counter[i][j]);
		fprintf(fd2, "\n");
	}

	fclose(fd);
	fclose(fd2);
	return 0;
}
