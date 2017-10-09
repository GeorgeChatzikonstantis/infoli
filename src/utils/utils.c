#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

int core_id, cores, cellCount;
int IO_NETWORK_DIM1, IO_NETWORK_DIM2, IO_NETWORK_SIZE;

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

char fpeek(FILE *stream)
{
	char c;

	c = (char)fgetc(stream);
	ungetc(((int)c), stream);

	return c;
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

inline mod_prec min(mod_prec a, mod_prec b){

	return (a < b) ? a : b;
}
