#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

/* gaussian_distribution() function
* 
* 
* 
* 
* 
*/


void gaussian_distribution(int, int, int, int, int, int,  double **, double *);

int main(int argc, char *argv[]){

  	FILE  *fd_out;
	int N;
	char outFileName[100];
	double mean, variance, average_connections_ratio=0.0;
	float numb_of_cells;

	int numb_neurons, i, j, k, x, y, z, counter;
	int max_dist_square; 
	double xd, yd;
	
	
	if (argc == 6){
		x  = atoi(argv[1]);
		y  = atoi(argv[2]);
		z  = atoi(argv[3]);
		mean = atoi(argv[4]);
		variance = atoi(argv[5]);
	}
	else {
	  	printf("Error: Not enough arguments: \nUsage: ./CubeGaussian <x> <y> <z> <mean value> <variance>");
		exit(1);
	}

//	N = cbrt(numb_of_cells);

	sprintf(outFileName,"cellConnections.txt");

	fd_out = fopen(outFileName,"w+");
	if(fd_out==NULL){
	   	printf("Error: Couldn't create %s\n", outFileName);
		exit(1);
	}

//	numb_neurons = N*N*N;
	numb_neurons = x * y* z;
	
	//step = 1 for the distances
//	max_dist_square = 3*pow((N-1),2);
	max_dist_square = pow((x-1),2) + pow((y-1),2) + pow((z-1),2);
	double *probabilities = (double*)malloc((max_dist_square+1) * sizeof(double));

	 xd = 1 / (variance*sqrt(2*M_PI));
         yd = -(1 / (2 * pow(variance, 2)));
	
	//Compute the probabilities for every distance
	for (i=1; i<= max_dist_square; i++){
		probabilities[i] =  2 * xd * exp(pow(sqrt(i)-mean,2)*yd);
	}

//	for (i=1; i<= max_dist_square; i++){
//		printf("%lf\n", probabilities[i]) ;
//	}
	
	double **connections = (double **)malloc(numb_neurons * sizeof(double *));
	for (i=0; i<numb_neurons; i++)
	  	connections[i] = (double *)malloc(numb_neurons * sizeof(double));
	float *connections_counter = (float*)calloc(numb_neurons, sizeof(float));	
	
	for (i=0; i<numb_neurons; i++){
	  	for (j=0; j<numb_neurons; j++)
		  	connections[i][j] = 0.0;
	}

	//u_random = fopen("/dev/urandom","r");

	//Compute the connections for every neuron
	srand(time(NULL));
	for (i=0; i<x; i++)
		for(j=0; j<y; j++)
		  	for(k=0; k<z; k++)
				gaussian_distribution(i, j, k, x, y, z, connections, probabilities);
	
	//Write the connectivity array to the output file
	for (i=0; i<numb_neurons; i++){
	  	for (j=0; j<numb_neurons; j++){
		  	fprintf(fd_out,"%lf", connections[i][j]);
			fprintf(fd_out," ");
			if (connections[i][j]>0)
				connections_counter[i]++;
		}
		fprintf(fd_out,"\n");
	}
	
	fclose(fd_out);

	if (numb_neurons<=1)
		printf("Network Too Small: %d cells.\n", numb_neurons);
	else {
		for (i=0; i<numb_neurons; i++) {
			connections_counter[i] = connections_counter[i] / ((float)(numb_neurons-1));
			average_connections_ratio += connections_counter[i];
		}
		average_connections_ratio = (100 * average_connections_ratio) / ((float)numb_neurons);
		printf("%0.3f\n", average_connections_ratio);
	}

	for (i=0; i<numb_neurons; i++){
	  	free(connections[i]);
	}

	free(connections);
	free(probabilities);

	return 0;

}


void gaussian_distribution(int ii, int jj, int kk,  int xx, int yy, int zz, double **connections, double *probs)
{	

  	int i, j, k, cell_to_check, cell_evaluating, N_square, distance_without_root; 
	double random_number, probability;
	
	N_square = xx*yy;
	cell_evaluating = kk * N_square + jj * xx + ii;
	
	//Compute the connections for every neuron
	for (i=0; i<xx; i++)
		for(j=0; j<yy; j++)
		  	for(k=0; k<zz; k++){
							
			  	cell_to_check = k * N_square + j * xx + i; 
				distance_without_root = pow(i-ii,2) + pow(j-jj,2) + pow(k-kk,2);
				probability = probs[distance_without_root];
				
//				connections[cell_evaluating][cell_to_check] = 2 * probability * 0.04;
				random_number = (double) rand() / (double) RAND_MAX;

				if (random_number <= probability)
					connections[cell_evaluating][cell_to_check] = 0.04;

			}

	return;
}
