#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include <random>

#define EPS 0.0001

int* CreateAndFillVector(int n) {
	int i;
	int* vector = (int*)malloc(n * sizeof(int));
	for (i = 0; i < n; i++) {
		vector[i] = 0 + rand() % 100; //от нуля до 100
	}
	return vector;
}

int AverageValueOfVector(int* vector, int n) {
	int i;
	int ave = 0;
	int *temp = vector;
	for (i = 0; i < n; i++) {
		ave = ave + vector[i];
		}
	ave = ave / n;
	return ave;
}

void DeleteVector(int* vector) {
	free(vector);
}

void PrintVector(int *vector, int n) {
	int i;
	for (i = 0; i < n; i++) {
		printf_s("%5d ", vector[i]);
		printf_s("\n");
	}
	printf_s("\n");
}

using namespace std;
int main(int argc, char* argv[]) {
	double time1, time2, delta_time_1, delta_time_2;
	MPI_Status status;
	FILE *f = NULL;
	int i, n,k, i1, i2, AveSum = 0, Ave = 0;
	int *vector = NULL;
	int ProcNum, ProcRank;
	int dataSize, bufferSize;
	int deltaSize;
	int *vector1 = NULL;
	if (argc >= 2) {
		n = atoi(argv[1]);
	}
	else {
		printf_s("Error with argv: argc!=2\n");
	}

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	n = 300;
	vector1 = (int *)malloc(n * sizeof(int));
	if (ProcRank == 0) {
		vector = CreateAndFillVector(n);
		vector1 = vector;
		printf_s("n = %d", n);
		PrintVector(vector, n);
		time1 = MPI_Wtime();
	}
	
	dataSize = (int)(n / (ProcNum)); 
	deltaSize = n - dataSize*(ProcNum);

	bufferSize = dataSize;

	if (ProcRank < (deltaSize + 1) && ProcRank != 0 ){
		bufferSize = dataSize + 1;
	}
	if (ProcRank == 0) {
		int *temp_start_vector = vector;

		for (i = 1; i < ProcNum; i++) {
			if (i < deltaSize + 1) {
				MPI_Send(temp_start_vector, bufferSize + 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				temp_start_vector = temp_start_vector + bufferSize + 1;
			}
			else {
				MPI_Send(temp_start_vector, bufferSize, MPI_INT, i, 0, MPI_COMM_WORLD);
				temp_start_vector = temp_start_vector + bufferSize;
			}
		}
	}

	k = n / ProcNum;
	i1 = k *   ProcRank;
	i2 = k * (ProcRank + 1);
	if (ProcRank == ProcNum - 1) i2 = n;
	for (int i = i1; i < i2; i++)
		AveSum = AveSum + vector[i];

	if (ProcRank == 0) {
		Ave = AveSum;
		for (int i = 1; i < ProcNum; i++) {
			MPI_Recv(&AveSum, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			Ave = Ave + AveSum;
		}
	}
	else // Все процессы отсылают свои частичные суммы
		MPI_Send(&AveSum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

		Ave = Ave / n;
		time2 = MPI_Wtime();
		delta_time_1 = time2 - time1;
		printf_s("MPI:\n");
		printf_s("	Time = %f\n", delta_time_1);
		printf_s("	Threads = %d\n", ProcNum);
		printf_s("	AverageValue = %d\n\n", Ave);
		if (delta_time_1 > EPS) {
			fopen_s(&f, "../../log/parallel.txt", "a");
			fprintf_s(f, "%f %d %d\n", delta_time_1, n, ProcNum);
			fflush(f);
			fclose(f);
		}
		time1 = MPI_Wtime();
		printf_s("ProcRank = %d", ProcRank);
		Ave = AverageValueOfVector(vector1, n);
		time2 = MPI_Wtime();
		delta_time_2 = time2 - time1;
		printf_s("Consistent implementation:\n");
		printf_s("	Time = %f\n", delta_time_2);
		printf_s("	AverageValue = %d\n", Ave);
		if (delta_time_2>EPS) {
			fopen_s(&f, "../../log/consistent.txt", "a");
			fprintf(f, "%f %d\n", delta_time_2, n);
			fflush(f);
			fclose(f);
		}
		printf_s("\nOptimal choice: %d", Ave);
		if (delta_time_1 < delta_time_2)
		{
			printf_s("parallel implementation.\n");
		}
		else
		{
			printf_s("consistent implementation.\n");
		}
		DeleteVector(vector1);
		DeleteVector(vector);
	MPI_Finalize();
	return 0;
}