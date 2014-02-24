/*row_wise.c - Divide en bloques de rayas por filas para A y en bloqus de rayas por filas para B*/
/*El número de procesadores debe ser divisible por el orden de las matrices n*/
#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define MAX 10000

void Read_matrix (char *prompt, float *matrix, int n);
void Print_matrix (char *prompt, float *matrix, int n);

main(int argc, char * argv[]) {
	srand(time(NULL));
	int 	my_rank, size, stage, temp;
	int 	n, local_n, i, j, k, source, dest, q, ind;
	float	*matrix_A;
	float	*matrix_B;
	float	*matrix_C;
	float   *local_A;
	float   *local_B;
	float	*local_C;
	double	start, end;
	MPI_Datatype	column;
	MPI_Status	status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (my_rank == 0) {
		printf("\t\t****************************************************\n");
		printf("\t\t*   Descomposición en Bloques de Rayas por Filas   *\n");
		printf("\t\t****************************************************\n\n");
	}

	if (my_rank == 0) {
		matrix_A = (float *)malloc(MAX*MAX*sizeof(float));
		matrix_B = (float *)malloc(MAX*MAX*sizeof(float));
		matrix_C = (float *)malloc(MAX*MAX*sizeof(float));
		if (argc == 2) {
			sscanf(argv[1], "%d", &n);
		} else if (my_rank == 0) {
			printf("¿ Cuál es el orden de las matrices ? \n");
			scanf("%d", &n);
		}
		local_n = n / size;
		/* Se lee la matriz A */
		Read_matrix ("Ingrese A :", matrix_A, n);
		Print_matrix ("Se leyó A :", matrix_A, n);

		/* Se lee la matriz B */
		Read_matrix ("Ingrese B :", matrix_B, n);
		Print_matrix ("Se leyó B :", matrix_B, n);

	}

	MPI_Bcast(&local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	local_A = (float *) malloc(MAX*MAX*sizeof(float));
	local_B = (float *) malloc(MAX*MAX*sizeof(float));
	local_C = (float *) malloc(MAX*MAX*sizeof(float));

/******************************************************************************/
/* Distribuir los bloques de filas de A y los bloques de filas de B entre 
   los procesadores del comunicador global */

	// Enviar los bloques de filas de A y B a todos los procesadores
	MPI_Scatter(matrix_A, local_n*n, MPI_FLOAT, local_A, local_n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Scatter(matrix_B, local_n*n, MPI_FLOAT, local_B, local_n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
/*****************************************************************************/
/* Algoritmo de Descomposición por Bloques de Rayas por filas para la matriz A 
   y por columnas para la matriz B */
	
	source = (my_rank - 1 + size) % size;
	dest = (my_rank + 1) % size;
	q = n / local_n;
	start = MPI_Wtime();
	for (stage = 0; stage < q; stage++) {
		ind = (my_rank - stage + size) % size;
		for (j = 0; j < local_n; j++) {
			for (i = 0; i < n; i++) {
				for (k = 0; k < local_n; k++) {
					local_C[i + j*n] += local_A[local_n*ind + k + j*n] * local_B[i + k*n];
				}
			}
		}
		MPI_Sendrecv_replace(local_B, local_n*n, MPI_FLOAT, dest, 0, source, 0, MPI_COMM_WORLD, &status);
	}
	end = MPI_Wtime();

/*****************************************************************************/

	// Recolectar los bloques local_C de cada procesador en matrix_C en le procesador 0
	MPI_Gather(local_C, local_n*n, MPI_FLOAT, matrix_C, local_n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	if (my_rank == 0) {
		Print_matrix ("El producto es C : ", matrix_C, n);
	}

	if (my_rank == 0)
		printf("n : %d\nDesc. por filas : %f segundos.\n", n, end - start);
	
	MPI_Finalize();

} /* main */

void Read_matrix (char *prompt, float *matrix, int n) {
	int i, j;
	printf("%s\n", prompt);
	for(i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				matrix[i*n + j] = 1 + rand()%5;
			}
	}
} /* Read_matrix */

void Print_matrix (char *prompt, float *matrix, int n)  {
	int i, j;
	printf("%s\n", prompt);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%4.2f ", matrix[i*n + j]);
		}
		printf("\n");
	}
	printf("\n");
} /* Print_matrix */
