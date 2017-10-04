#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

int main(int argc, char *argv[]){

	FILE *pFileBase, *pFileTest;
	char fileBaseName[100];
	char fileTestName[100];
	float error, relative_error, max_error=0.0;
	float stored_v1, stored_v2;
	int index, line, stored_line, stored_index;

	if (argc!=2) {
		printf("Checker Error: Incorrect arguments.\nUsage: ./checker <nw_size>\n");
		exit(EXIT_FAILURE);
	}
	int nw_size = atoi(argv[1]);
	float* matrixBase = (float*) calloc(nw_size, sizeof(float));
	float* matrixTest = (float*) calloc(nw_size, sizeof(float));

	sprintf(fileBaseName,"InferiorOlive_Output_Baseline.txt");
	sprintf(fileTestName,"InferiorOlive_Output.txt");

	pFileBase = fopen(fileBaseName,"r");
	pFileTest = fopen(fileTestName,"r");

	line=1;
	while ((ReadFileLine(pFileBase, &matrixBase[0], nw_size))&&(ReadFileLine(pFileTest, &matrixTest[0], nw_size))) {
		for (index=0; index<nw_size; index++) {
			error = fabsf(matrixBase[index]-matrixTest[index]);
			relative_error = fabsf(error/matrixBase[index]);
			if (relative_error > max_error) {
				max_error = relative_error;
				stored_v1 = matrixBase[index];
				stored_v2 = matrixTest[index];
				stored_line = line;
				stored_index = index;
			}
		}
		line++;
	}
	printf("Maximum error was found %0.4f%%.\n", max_error);
	if (max_error>0) {
		printf("Maximum discrepancy was found at line %d, neuron %d:\n", stored_line, stored_index);
		printf("Correct value: %0.8f - Runtime value: %0.8f.", stored_v1, stored_v2);
	}

	fclose(pFileBase);
	fclose(pFileTest);

	return 1;
}

char fpeek(FILE *stream) {

	char c;

	c = (char)fgetc(stream);
	ungetc(((int)c), stream);

	return c;
}

int ReadFileLine(FILE *pFile, float *matrix, int size){

	char c= fpeek(pFile);
	if (c==EOF)
		return 0;
	char *strNumber;
	int bufSize = size*30;
	char* buffer = (char*) calloc(bufSize, sizeof(char));
	int i = 0;

	//Get one line
	if(fgets(buffer, bufSize, pFile)){
	//Convert the ASCII string of one element to a floating point value
		strNumber = strtok(buffer," ");
		i = 0;

		while ((strNumber != NULL) && (i<size)){

			//store number
                        matrix[i] = atof(strNumber);
			//get next
			strNumber = strtok(NULL, " ");
			i++;
		}

		if(i<size){
		//BUG: if only one element is missing but the line ends in a space, the error is not detected
			printf("Error: Input line doesn't have enough elements, only %d\n", i);
			exit(EXIT_FAILURE);
		}

		free(buffer);
		return 1;//success
	}else{

		if(!feof(pFile)){
			printf("Error: Reading from input file didn't finish successfully\n");
			exit(EXIT_FAILURE);
		}

		free(buffer);
		return 0;//end of file
	}
}
