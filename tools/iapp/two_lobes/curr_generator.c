#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_POISSON_VALUES 10	//always set this less than 40, things might get nasty in EXTREMELY rare cases where you request more than 40 spikes per sec

int factorial(int n) {

 	int retval = 1, i;
	if (n<0)
		printf("Warning: Negative factorial calculation demand detected.\n");
 	for (i = n; i > 1; --i)
 		retval *= i;
 	return retval;
}

void detect_lobe_spike_creation(int x, int y, int** stim_matrix, float frequency, float phase, int t) {

	if (frequency==0) {
		stim_matrix[x][y] = 0;
		return;
	}
	int T = (int)(20000.0 / frequency);			//find how many steps transpire between each spike stimulation
	int offset_steps = (int)((phase/360.0)*T);
	if ((t>=offset_steps)&&(((t-offset_steps)%T)==0)) {
		stim_matrix[x][y] = 500;
	} else
		stim_matrix[x][y] = 0;
	return;

}

void detect_noise_spike_creation(int x, int y, int*** spikes, int** spikes_created, int** stim_matrix, int t) {

	int sim_sec = t/20000;
	int window_spikes_left = spikes[x][y][sim_sec] - spikes_created[x][y];
	if ((spikes[x][y][sim_sec]==0)||(window_spikes_left==0)) {		//no spikes (left) for this observation window in this cell!
		stim_matrix[x][y] = 0;
		return;
	}

/* each second is broken down in frames according to how many spikes are generated within this second. Then we assign an increasing chance while we are
   in this frame to stimulate a spike, which ends when we reach near the frame's end in which case we force the spike to happen. No check needed for now
   concerning the frame's width relevant to the window, we can have 10 frames (spikes) max within 20,000 steps, hence at a frame width of at least 2,000
   steps, more than enough for a stimulus to be created here. */

	int spike_frame_beginning = (sim_sec + ( ((float)(spikes_created[x][y])) / spikes[x][y][sim_sec] ))*20000;
	int spike_frame_ending = (sim_sec + ( ((float)(spikes_created[x][y]+1)) / spikes[x][y][sim_sec] ))*20000 -1;

//	if (t==spike_frame_beginning)
//		printf("(%d,%d):%d-%d\n", x, y, spike_frame_beginning, spike_frame_ending);

	if (t<spike_frame_beginning) {						//we are not in the next spike frame yet!
		stim_matrix[x][y] = 0;
		return;
	}
	if (t>= (spike_frame_ending-501)) {
		stim_matrix[x][y] = 500;					//we need to manually set the stimulus before we run out of time for this frame!
		spikes_created[x][y]++;
		return;
	}
	
	float random_number = ((float) rand()) / RAND_MAX;
	float chance_to_get_spike = ((float)(t - spike_frame_beginning)) / ((float)(spike_frame_ending - 500 - spike_frame_beginning));
	if (random_number < chance_to_get_spike) {
		stim_matrix[x][y] = 500;					//we generated a spike within this spike frame!
		spikes_created[x][y]++;
	} else {
		stim_matrix[x][y] = 0;					//we generated a spike within this spike frame!
	}
	return;
}

int main (int argc, char* argv[]) {

	if (argc!=11) {
		printf("Error!\nUsage: ./input_generate network_X network_Y lobe1_radius lobe2_radius\nlobes_distance noise_lamda f1 f2 relative_phase_degrees duration_s.\n");
		return 0;
	}
	int dim_1 = atoi(argv[1]);
	int dim_2 = atoi(argv[2]);
	int r1 = atoi(argv[3]);
	int r2 = atoi(argv[4]);
	int d = atoi(argv[5]);
	float lamda = atof(argv[6]);
	float f1 = atof(argv[7]);
	float f2 = atof(argv[8]);
	float phase = atof(argv[9]);
	if ((phase < 0)||(phase > 360)) {
		printf("Error!\nPhase degrees needs to be within the [0, 360] interval.\n");
		return 0;
	}
	int experiment_secs = atoi(argv[10]);
	int total_steps = experiment_secs*20000;

	char currFileName[100];
	sprintf(currFileName,"currentInput.txt");
	FILE* fd= fopen(currFileName, "w+");

	int i, j, k, t = 0, distance_x, distance1_y, distance2_y;
	float eucl_distance;
	int c1_y = ((dim_2 - d) / 2) -1;
	int c2_y = c1_y + d + 1;
	int c_x = dim_1 /2;
	
	char** grid = (char**) malloc(dim_1*sizeof(char*));		//create a grid indicating which cell has lobe 1-lobe 2-noise properties
	for (i=0; i<dim_1; i++) {
		grid[i] = (char*) malloc(dim_2*sizeof(char));

		distance_x = abs(i-c_x);
		for (j=0;j<dim_2;j++) {

			grid[i][j] = 'n';			//noise cells
			distance1_y = abs(j-c1_y);
			distance2_y = abs(j-c2_y);

			eucl_distance = sqrt(pow(distance_x,2)+pow(distance1_y,2));
			if (eucl_distance <= (float) r1)
				grid[i][j] = '1';		//lobe 1 cells
			eucl_distance = sqrt(pow(distance_x,2)+pow(distance2_y,2));
			if (eucl_distance <= (float)  r2) {
				if (grid[i][j] == '1') {
					printf("Error!\nExperiment parameters lead to lobe cells overlapping!\n");
					return 0;
				}
				grid[i][j] = '2';		//lobe 2 cells
			}
		}
	}

	int** stim_matrix = (int**) malloc((sizeof(int*))*dim_1);		//this matrix shows whether a stimulus is active or not, set at zero initially
	for (i=0;i<dim_1;i++)
		stim_matrix[i] = (int*) calloc(dim_2, (sizeof(int)));		//calloc ensures all cells begin with zero stimulus

/* based on our lamda, we need a series of values denoting the poisson distribution values for probability denoting how possible
   it is to have UP TO x spikes in each time window (window default length: 1 sec). We shall name these values as poisson_values. Due to infoli's
   behaviour, its lamda is expected to be between 0.5 and 2, hence we will evaluate a set amount of values. Mind that these values are CUMULATIVE
   possibilities, thus creating a "bucket" system which eases the next step.*/

	float* poisson_values = malloc(MAX_POISSON_VALUES*sizeof(float));
	float total_poisson_value = 0.0;
	for (i=0;i<MAX_POISSON_VALUES;i++) {
		total_poisson_value += (pow(lamda,i))/((expf(lamda))*(factorial(i)));
		poisson_values[i] = total_poisson_value;
	}

/* based on these values, we can now calculate, for each noise cell, how many spikes it shall generate in every window of 1 sec brain time.
   We do this by generating a random number and start comparing it with the first probability bucket, in increasing index fashion. If 
   rnd <= poisson[i], then we got i spikes for this cell's window. If not, then we carry on our comparisons, up till we reach max i. At this point
   we no longer compare and just assign max i spikes for this window -highly unlikely-. Ask a mathematician whether this is a valid method.. */

	int*** spikes = (int***) malloc(dim_1*sizeof(int**));
	srand(time(NULL));
	float random_number;

	for (i=0;i<dim_1;i++) {
		spikes[i] = (int**) malloc(dim_2*sizeof(int*));
		for (j=0;j<dim_2;j++) {
			spikes[i][j] = (int*) malloc(experiment_secs*sizeof(int));
			if (grid[i][j]=='n')						//we only generate these values for noise cells!
				for (t=0;t<experiment_secs;t++) {

					random_number = ((float) rand()) / RAND_MAX;
					k=0;
					while ((random_number > poisson_values[k])&&(k < MAX_POISSON_VALUES))
						k++;
					spikes[i][j][t] = k;

				}
		}
	}
	free(poisson_values);

/*	for (t=0;t<experiment_secs;t++) {
		for (i=0;i<dim_1;i++) {
			for (j=0;j<dim_2;j++)
				printf("%d ", spikes[i][j][t]);
			printf("\n");
		}
		printf("\n");
	} */

	int** spikes_created = (int**) malloc(dim_1*sizeof(int*));		//matrix that records how many spikes have been raised by every cell for this sec
	for (i=0;i<dim_1;i++)
		spikes_created[i] = (int*) calloc(dim_2, sizeof(int));

	int chars_per_line = (dim_1 * dim_2 * 6) + 3;
	char* line_message_buffer = (char*) malloc(chars_per_line*sizeof(char));	//buffer holding this step's currents, needs 4 digits+1 dot+1 space per cell
											//and a special |\n series of char at the end of the line and \0 at str's end

	t=0;
	while (t<total_steps) {

		if ((t%20000)==0)				//we need to clear the spikes_created matrix in the beginning of each sec
			for (i=0;i<dim_1;i++)
				for (j=0;j<dim_2;j++)
					spikes_created[i][j] = 0;

		strcpy(line_message_buffer, "");
		for (i=0;i<dim_1;i++)
			for (j=0;j<dim_2;j++) {

				if (stim_matrix[i][j]!=0) {
					strcat(line_message_buffer, "6.000 ");
					stim_matrix[i][j]--;
				} else {
					strcat(line_message_buffer, "0.000 ");
					switch (grid[i][j]) {
						case 'n': {
							detect_noise_spike_creation(i, j, spikes, spikes_created, stim_matrix, t);
							break;
						}
						case '1': {
							detect_lobe_spike_creation(i, j, stim_matrix, f1, 0, t);
							break;
						}
						case '2': {
							detect_lobe_spike_creation(i, j, stim_matrix, f2, phase, t);
							break;
						}
						default: {
							printf("This shouldn't be happening...\n");
							break;
						}
					}
				}

			}

		strcat(line_message_buffer, "\n");
		fprintf(fd, "%s", line_message_buffer);
		t++;
	}

	fclose(fd);
	return 1;

}
