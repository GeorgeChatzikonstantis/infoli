#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

struct connection {
	int cell_id;
	struct connection *next;
};
typedef struct connection node;

void Evaluate_connections (int, int, FILE*, FILE*, int, int, double, double, node***, double, double*);

int main (int argc,char **argv) {

	int rows, columns, i, j, max_distance, seed;
	char filename[30];
	sprintf(filename,"cellConnections.txt");
	FILE *out_file, *u_random;
	double x, y, middle_value, variance, *evaluated_possibilities;
	node ***global_connection_table;

	if (argc!=5) {
		printf("Error: Usage ./random_connectivity rows columns middle_value variance\n");
		return 0;
	}

	rows = atoi(argv[1]);
	columns = atoi(argv[2]);
	middle_value = atoi(argv[3]);
	variance = atoi(argv[4]);

	out_file = fopen(filename, "w"); 
	u_random = fopen("/dev/urandom","r");

	global_connection_table = malloc(sizeof(node**)*rows);
	for (i=0;i<rows;i++) {
		global_connection_table[i] = malloc(sizeof(node*)*columns);
		for (j=0;j<columns;j++) 
				global_connection_table[i][j] = NULL;
	}

	max_distance = pow(rows-1,2) + pow(columns-1,2);

	evaluated_possibilities = calloc(max_distance+1, sizeof(double));


	x = 1 / (variance*sqrt(2*M_PI));
	y = -(1 / (2 * pow(variance, 2)));

	for (i=1;i<max_distance;i++)
		evaluated_possibilities[i] = x * exp(pow(sqrt(i)-middle_value,2)*y);
	
	/* Calculate connections 
	 * for each cell 
	 */

	fread(&seed, sizeof(int), 1, u_random);
        srand(seed);

	for (i=0;i<rows;i++)
		for (j=0;j<columns;j++)
			Evaluate_connections(i, j, out_file, u_random, rows, columns, x, y, global_connection_table, middle_value, evaluated_possibilities);
	
	return 0;

}


/* Change file format from now on 
 * file will be formatted as follows : number of line corresponds to number of cell,for example line 1 is for cell 0 etc
 * first number is a number of cell then the next is 1/2/3 (cell_line will send to/receive from/BOTH from cell described here)
 */

void Evaluate_connections (int cell_x, int cell_y, FILE* out_file,FILE* u_random,int rows, int columns, double x, double y,node ***global_connection_table, double middle_value, double *evaluated_possibilities) {

	int i, j, seed, cell_to_check, cell_evaluating, distance_without_root;
	double possibility, random_number;
	char cell_name[10];
	node *temporary, *next;

	for (i=0;i<rows;i++) 
		for (j=0;j<columns;j++) {

			cell_to_check = i*columns + j;
			cell_evaluating = cell_x*columns + cell_y;

			if (cell_to_check > cell_evaluating) {


				distance_without_root = pow(i-cell_x,2)+pow(j-cell_y,2);
				possibility = evaluated_possibilities[distance_without_root];
					
				random_number = (double) rand() / (double) RAND_MAX;

				if (random_number <= possibility) {

					next = malloc( sizeof(node));
					next->cell_id = cell_to_check;
					temporary = global_connection_table[cell_x][cell_y];
					next->next = temporary;
					global_connection_table[cell_x][cell_y] = next;

					sprintf(cell_name,"%d ",cell_to_check);
					fwrite(cell_name, sizeof(char), strlen(cell_name), out_file);
					fwrite("B ", sizeof(char), 2, out_file);
					continue;
				}

				continue;	
			}

			temporary = global_connection_table[i][j];
			while (temporary) {
				if (temporary->cell_id == cell_evaluating) {
					sprintf(cell_name,"%d ",cell_to_check);
					fwrite(cell_name, sizeof(char), strlen(cell_name), out_file);
					fwrite("B ", sizeof(char), 2, out_file);
					break;
				}
				temporary = temporary->next;
			}
		}

	fwrite("|\n", sizeof(char), 2, out_file);
	

}
