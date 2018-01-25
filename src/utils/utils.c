#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

int cellCount, IO_NETWORK_DIM1, IO_NETWORK_DIM2, IO_NETWORK_SIZE;

typedef float mod_prec;

void print_usage() {

	printf("Usage for infoli simulator:\n");
	printf("/path/to/bin/infoli.x -n <integer> -p <[0,1]> -t <float>\n");
	printf("-i <filename> -c <filename> -o <filename>\n");
	printf("Argument summary:\n");
	printf("-n denotes number of neurons to be simulated. Default value: 1000.\n");
	printf("-p denotes network density in [0,1] format. Default value: 0.5 (=50%%).\n");
	printf("-t denotes simulation time in milliseconds. Default value: 5000.\n");
	printf("-i denotes a filename for input stimulus. By default, stimulus is hardcoded.\n");
	printf("-c denotes a filename for configuration file. By default, configuration is read from default.conf.\n");
	printf("-o denotes a filename for simulation output. Default value: InferiorOlive_Output.txt.\n");

}

void print_connections(int* conn_array, int cell_id, int connections_no) {

	printf("Cell %d has %d connections:", cell_id, connections_no);
	int i;
	for (i=0; i<connections_no;i++) {
		if ((i%10)==0)
			printf("\n");
		printf("%d ", conn_array[i]);
	}
	printf("\n");

	return;

}

char fpeek(FILE *stream)
{
	char c;

	c = (char)fgetc(stream);
	ungetc(((int)c), stream);

	return c;
}

void removeSubstring(char *s, const char *toremove) {

	s = strstr(s,toremove);
	memmove(s,s+strlen(toremove),1+strlen(s+strlen(toremove)));

}

void stopAtSubstring(char *s, const char *toremove) {

	char* temp = strstr(s,toremove);
	if (temp!=NULL) {
		temp[0] = '\0';
		temp = NULL;
	}

}

int ReadFileLine(FILE *pInFile, mod_prec *iAppArray){

	char c= fpeek(pInFile);
	if (c==EOF)
		return 0;
	
	char *strNumber;
	int bufSize = cellCount*20;	//more than enough but will do for nao
	char* buffer = (char*) _mm_malloc(bufSize*sizeof(char), 64);

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
		return 1;//success
	}else{
		if(!feof(pInFile)){
			printf("Error: Reading from input file didn't finish successfully\n");
			exit(EXIT_FAILURE);
		}
		free(buffer);
		return 0;//end of file
	}
}

void read_g_CaL_from_file(mod_prec* g_cal) {

	int i;

	FILE* fd = fopen("gcal_file.txt","r");
	for (i=0;i<cellCount;i++) 
		fscanf(fd, "%f ", &g_cal[i]);
	fclose(fd);

	return;

}

void read_parameters(char* paramsFileName, mod_prec* values) {

	int i, linesToRead = 42;
	char* valueName= (char*) malloc(100*sizeof(char));

	FILE* fd = fopen(paramsFileName,"r");
	if(fd==NULL) {
		printf("Error: Couldn't open %s\n", paramsFileName);
		exit(EXIT_FAILURE);
	}

	for (i=0; i<linesToRead; i++) {
		if ((fscanf(fd, "%s = %f\n", valueName, &values[i]))!=2) {
			printf("Configuration Error: Improper Parameter %s in file %s, line %d.\n", valueName, paramsFileName, i+1);
			printf("For a proper Configuration file format, please see default.conf in bin folder.\n");
			exit(EXIT_FAILURE);
		}
	}

	fclose(fd);
	free(valueName);

	return;
}

int conn_marking_uniform(int* conn_buffer, short* external_cell, int global_cell_id, int local_nw_size, int nw_size, double probability, int core_id) {

	int total_amount_connections=0, sender_cell;
	float rndm;

	for (sender_cell=0;sender_cell<nw_size;sender_cell++) {
		if (sender_cell == global_cell_id)
			;               //no self-feeding connections allowed
		else {
			rndm = ((float) rand()) / ((float) RAND_MAX);   //generate rng and compare to probability
			if (rndm <= probability) {
				total_amount_connections++;  //increase neighbour count
				conn_buffer[sender_cell]++; //mark that we formed a bond with this cell
				if ((sender_cell/local_nw_size)!=core_id) //if this cell does not belong to core
					external_cell[sender_cell]  = 1; //mark it
			}
		}

	}
	return total_amount_connections;
}

int generate_gaussian_number(int lower_bound, int higher_bound) {

	float x = 0, rndm;
	int range = higher_bound - lower_bound;
	int i, reps=20;

	for (i=0;i<reps;i++) {
		rndm = ((float) rand()) / ((float) RAND_MAX);
		x += rndm;
	}
	x /= reps; //generate a normally-distributed number in [0,1]
	x *= range; //generate a normally-distribute number in [0,range]
	x += lower_bound; //generate a normally-distributed number in [lower_bound, higher_bound]
	
	return (round(x));

}

int conn_marking_gaussian_1D(int* conn_buffer, short* external_cell, int global_cell_id, int local_nw_size, int nw_size, double probability, int core_id) {

	int total_amount_connections=0, sender_cell;
	float mean=0, deviation=1.0, scale=367.647;
	float gaussian_probability;
	float rndm;
	

	float average_connections=nw_size*probability;
	if (probability==0)
		total_amount_connections=0;
	else
		total_amount_connections=generate_gaussian_number(round(average_connections*0.9), round(average_connections*1.1));

	float x = 1 / (deviation*sqrt(2*M_PI));
	float y = -(1 / (2 * pow(deviation, 2)));

	float distance_1D=0.0;
	int connections_found=0;
	while (connections_found<total_amount_connections) {
		for (sender_cell=0; sender_cell<nw_size; sender_cell++) {
			if (connections_found>=total_amount_connections)
				break;
			if (sender_cell == global_cell_id)
				;		//no self-feeding connections allowed
			else if (conn_buffer[sender_cell])
				;		//cell already marked
			else {
				distance_1D = ((float)(abs(sender_cell-global_cell_id)))/scale;
				gaussian_probability = x * exp(pow(distance_1D-mean,2)*y); //if ((core_id==0)&&((sender_cell%100)==0)) printf("%.2f:%.3f%% ", distance_1D, 100*gaussian_probability);
				rndm = ((float) rand()) / ((float) RAND_MAX);   //generate rng and compare to probability
				if (rndm <= gaussian_probability) {
					conn_buffer[sender_cell]++; //mark that we formed a bond with this cell
					connections_found++;  //increase amount of marks
					if ((sender_cell/local_nw_size)!=core_id) //if this cell does not belong to core
						external_cell[sender_cell]  = 1; //mark it
				}
			}
		}

	}
	return connections_found;

}

int conn_marking_gaussian_3D(int* conn_buffer, short* external_cell, int global_cell_id, int local_nw_size, int nw_size, double probability, int core_id) {

	int total_amount_connections=0, sender_cell;
	float mean=0, deviation=1.0, scale=20; //scale 7.15 makes it so that neuron 1000 units away has 10% chance
	float gaussian_probability;
	float rndm;

	int base_3D= ((int)(pow(nw_size,0.33333)))+1;
	int base_3D_squared=base_3D*base_3D;

	int my_x=global_cell_id/base_3D_squared;
	int my_y=(global_cell_id%base_3D_squared)/base_3D;
	int my_z=(global_cell_id%base_3D_squared)%base_3D;

	int sender_x, sender_y, sender_z;

	float average_connections=nw_size*probability;
	if (probability==0)
		total_amount_connections=0;
	else
		total_amount_connections=generate_gaussian_number(round(average_connections*0.9), round(average_connections*1.1));

        float const1 = 2.5 / (deviation*sqrt(2*M_PI));
        float const2 = -(1 / (2 * pow(deviation, 2)));

        float distance_3D=0.0;
        int connections_found=0;
        while (connections_found<total_amount_connections) {
                for (sender_cell=0; sender_cell<nw_size; sender_cell++) {
                        if (connections_found>=total_amount_connections)
                                break;
                        if (sender_cell == global_cell_id)
                                ;               //no self-feeding connections allowed
                        else if (conn_buffer[sender_cell])
                                ;               //cell already marked
                        else {
				sender_x=sender_cell/base_3D_squared;
				sender_y=(sender_cell%base_3D_squared)/base_3D;
				sender_z=(sender_cell%base_3D_squared)%base_3D;

                                distance_3D = (sqrt(pow(my_x-sender_x,2)+pow(my_y-sender_y,2)+pow(my_z-sender_z,2)))/scale;
                                gaussian_probability = const1 * exp(pow(distance_3D-mean,2)*const2);
		//	if (((sender_cell%1)==0)&&(global_cell_id==0)) printf("x:%d y:%d z:%d->%.2f, %.3f%%\n", sender_x, sender_y, sender_z, distance_3D, 100*gaussian_probability);
                                rndm = ((float) rand()) / ((float) RAND_MAX);   //generate rng and compare to probability
                                if (rndm <= gaussian_probability) {
                                        conn_buffer[sender_cell]++; //mark that we formed a bond with this cell
                                        connections_found++;  //increase amount of marks
                                        if ((sender_cell/local_nw_size)!=core_id) //if this cell does not belong to core
                                                external_cell[sender_cell]  = 1; //mark it
                                }
                        }
                }

        }
        return connections_found;

}

int conn_marking_nearest_1D(int* conn_buffer, short* external_cell, int global_cell_id, int local_nw_size, int nw_size, double probability, int core_id) {

	int total_amount_connections=0, sender_cell, connections_found=0, step;
	float average_connections=nw_size*probability;

	if (probability==0)
		total_amount_connections=0;
	else
		total_amount_connections=generate_gaussian_number(round(average_connections*0.9), round(average_connections*1.1));

	step=1;
	while ((connections_found < total_amount_connections)&&(step<nw_size)) {

		sender_cell=global_cell_id+step;
		if (sender_cell<nw_size) {
			conn_buffer[sender_cell]++; //mark that we formed a bond with this cell
			connections_found++;  //increase amount of marks
			if ((sender_cell/local_nw_size)!=core_id) //if this cell does not belong to core
				external_cell[sender_cell]  = 1; //mark it
		}

		sender_cell=global_cell_id-step;
		if (sender_cell>=0) {
			conn_buffer[sender_cell]++; //mark that we formed a bond with this cell
			connections_found++;  //increase amount of marks
			if ((sender_cell/local_nw_size)!=core_id) //if this cell does not belong to core
				external_cell[sender_cell]  = 1; //mark it
		}

		step++;
	}

	return connections_found;
}

inline mod_prec min(mod_prec a, mod_prec b){

	return (a < b) ? a : b;
}
