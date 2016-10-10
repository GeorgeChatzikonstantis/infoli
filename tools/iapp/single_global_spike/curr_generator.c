#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char* argv[]) {

	if (argc!=4) {
		printf("Error!\nUsage: ./stim_generator.x network_X network_Y duration_ms.\n");
		return 0;
	}
	int dim_1 = atoi(argv[1]);
	int dim_2 = atoi(argv[2]);
	int experiment_msecs = atoi(argv[3]);
	int total_steps = experiment_msecs*20;

	char currFileName[100];
	sprintf(currFileName,"currentInput.txt");
	FILE* fd= fopen(currFileName, "w+");

	int i, j, t=0;
	int line_message_size = 6 * dim_1 * dim_2 + 2;
	char* line_message_buffer = (char*) malloc(line_message_size*sizeof(char));

	while (t<total_steps) {

		strcpy(line_message_buffer, "");
		for (i=0;i<dim_1;i++)
			for (j=0;j<dim_2;j++)

				if ((t>20000)&&(t<20500))
					strcat(line_message_buffer, "6.000 ");
				else
					strcat(line_message_buffer, "0.000 ");


		strcat(line_message_buffer, "\n");
		fprintf(fd, "%s", line_message_buffer);
		t++;
	}

	fclose(fd);
	return 1;

}
