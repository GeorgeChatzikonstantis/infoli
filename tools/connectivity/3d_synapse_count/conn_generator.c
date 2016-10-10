#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

void Gaussian_probability_function(double*, int, int, int);
int calculate_synapse_count(double*, int, int, int);
int coin_flip(double);
void clear_buffer(float*, int);

int main(int argc, char *argv[]){

  	FILE  *fd_out;
	char outFileName[100];
	int i, j, k, x, y, z, connections_counter=0;
	int N, target_synapse_count, rng, count, sparse=0;
	float avg_counter=0.0;
		
	if (argc > 4){
		x  = atoi(argv[1]);
		y  = atoi(argv[2]);
		z  = atoi(argv[3]);
		N = x * y* z;
		target_synapse_count = atoi(argv[4]);
	}
	if (argc > 5)
		sparse = atoi(argv[5]);
	if ((argc<5)||(argc>6)){
	  	printf("Error: Wrong arguments:\nUsage: ./conn_generator.x x y z average_synapse_count <sparse_matrix>\n");
		exit(1);
	}
	
	if (target_synapse_count > 0.95*N) {
		printf("Error:\nTarget Synapse Count must be less than 95%% of NW size\n");
                exit(1);
	}

	sprintf(outFileName,"cellConnections.txt");

	fd_out = fopen(outFileName,"w+");
	if(fd_out==NULL){
	   	printf("Error: Couldn't create %s\n", outFileName);
		exit(1);
	}

	int* synapse_count = (int*) calloc(N, sizeof(int));
	double* synapse_count_probabilities = (double*) calloc((N+1), sizeof(double));
	Gaussian_probability_function(synapse_count_probabilities, (N+1), target_synapse_count, 1);

	srand(time(NULL));
	for (i=0; i<N; i++)
		synapse_count[i] = calculate_synapse_count(synapse_count_probabilities, (N+1), 1, (i%2));
	free(synapse_count_probabilities);

	float *connections = (float *) calloc(N, sizeof(float));
	for (i=0; i<N; i++) {
		count = 0;
		while (count < synapse_count[i]) {
			rng = (int) (((float) rand() / (float) RAND_MAX) * N);
			if ((rng<N)&&(rng!=i)&&(connections[rng]==0)) {
				connections[rng] = 0.04;
				count++;
			}
		}

		if (sparse==0) {
			for (j=0; j<N; j++) {
				fprintf(fd_out,"%0.2lf ", connections[j]);
				if (connections[j]>0)
					connections_counter++;
			}
			fprintf(fd_out,"\n");
		} else {
			for (j=0; j<N; j++)
				if (connections[j]>0) {
					fprintf(fd_out,"%d %d %0.2lf\n", i, j, connections[j]);
					connections_counter++;
				}
		}

		clear_buffer(connections, N);
	}

	avg_counter = ((float)connections_counter) / ((float)N);
	printf("Average Synapse Count Per Neuron: %0.2f\n", avg_counter);
	fclose(fd_out);

	return 0;
}

void Gaussian_probability_function(double* probs, int size, int mean, int sigma) {

	double xd, yd;
	int i;

	xd = 1 / (sigma*sqrt(2*M_PI));
        yd = -(1 / (2 * pow(sigma, 2)));

	for (i=0; i< size; i++)
		probs[i] =  xd * exp(pow((i-mean),2)*yd);

	return;
}

int calculate_synapse_count(double* probs, int size, int sigma, int reverse) {

	int i=0, count=0;
	double max_chance = 1/(sigma*sqrt(2*M_PI));
	double rng = ((double) rand() / (double) RAND_MAX) * max_chance;
	
	if (reverse==0) {
		for (i=0; i<size; i++)
			if (rng<probs[i])
				count = i;
	} else {
		for (i=(size-1); i>=0; i--)
			if (rng<probs[i])
				count = i;
	}
		
	return count;

}

int coin_flip(double chance) {

	double rng = (double) rand() / (double) RAND_MAX;
	if (rng < chance)
		return 1;
	else
		return 0;

}

void clear_buffer(float* buffer, int size) {

	int i;
	for (i=0; i<size; i++)
		buffer[i]=0;
	
}
