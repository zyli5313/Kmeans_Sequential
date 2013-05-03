/*
 ============================================================================
 Name        : openMPI-Kmeans.c
 Author      : Guanyu Wang, Zeyuan Li
 Version     :
 Copyright   :
 Description : Parallel K-means method using openMPI
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "seq-Kmeans.h"
#define MY_MAXITER 1000

float** kmeans_read(char *fname, int nline, int ndim) {
	float **data;
	char *token;
	int i = 1, num = 0;

	//printf("nline:%d\tndim:%d\n", *nline, ndim);

	FILE * fp;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	char cwd[1024];
	getcwd(cwd, sizeof(cwd));
	printf("Current working dir: %s\n", cwd);

	fp = fopen(fname, "r");
	if (fp == NULL )
		exit(EXIT_FAILURE);

	data = (float **) malloc(nline * sizeof(float *));
	data[0] = (float *) malloc(nline * ndim * sizeof(float));
	for (i = 1; i < nline; i++)
		data[i] = data[i - 1] + ndim;

	while ((read = getline(&line, &len, fp)) != -1) {
		int j = 0;
		token = strtok(line, " ,");
		while (token != NULL ) {
			data[num][j++] = atof(token);
//				printf("%s\t", token);
			token = strtok(NULL, " ,");
		}
		num++;
//			printf("%d\n", num);
//			printf("Retrieved line of length %zu :\n", read);
//			printf("%s", line);
	}

	return data;
}

int kmeans_write(char *outputfilename, int numberofLocalData,
		int numberofTotalData, int numberofClusters, int numberofCoordinates,
		float **clusters, int *localMemebership, int ranktooutput) {
	FILE * fp;
	int i, j;

	printf(
			"Start writing results (cluster centroid and membership) of K=%d cluster centers to file \"%s\"\n",
			numberofClusters, outputfilename);

//	struct stat st = {0};
//
//	if (stat("/some/directory", &st) == -1) {
//		mkdir("/some/directory", 0777);
//	}

	fp = fopen(outputfilename, "w");
	if (fp == NULL ) {
		printf("Error: cannot access the outputfile: %s\n", outputfilename);
		exit(EXIT_FAILURE);
	}
	else {
		for (i = 0; i < numberofClusters; i++) {
			char str[32];
			fprintf(fp, "%d ", i);
			//				MPI_File_write(mpif, str, strlen(str), MPI_CHAR, &mpistatus);
			for (j = 0; j < numberofCoordinates; j++) {
				fprintf(fp, "%f ", clusters[i][j]);

			}
			//				MPI_File_write(mpif, "\n", 1, MPI_CHAR, &mpistatus);
			fprintf(fp, "\n");
		}
	}
	printf("Finish writing clusters\n");

	//wait for the other memebership from other processors*/
	int outputCount = 0;
	//print out local membership
	for (j = 0; j < numberofLocalData; j++)
		fprintf(fp, "%d %d\n", outputCount++, localMemebership[j]);
	fclose(fp);

	return 0;
}

//for find the nearest neighbor in the given set;
int find_NN(float *datapoint, float ** neighborset, int numberofNeighber,
		int numberofCoordinates) {
	int i, j;
	int nearest_neighbor = -1;
	float distance, mindist;
	mindist = FLT_MAX;
	for (i = 0; i < numberofNeighber; i++) {
		distance = 0.0;
		for (j = 0; j < numberofCoordinates; j++) {
			distance += (datapoint[j] - neighborset[i][j])
					* (datapoint[j] - neighborset[i][j]);
			if (distance < mindist) {
				nearest_neighbor = i;
			}
		}
	}
	return nearest_neighbor;
}

int kmeans(float **data, int numberofClusters, int numberofCoordinates,
		int numberofTotalData, float stopthreshold, int *membership,
		float **clusters) {
	float **updatedClusters;
	int *updatedClusterSize;

	int i, j;
	//initialization
	//malloc space for pointers
	updatedClusterSize = (int *) calloc(numberofClusters, sizeof(int));

	updatedClusters = (float **) malloc(numberofClusters * sizeof(float*));
	updatedClusters[0] = (float *) calloc(
			numberofClusters * numberofCoordinates, sizeof(float));

	if (!updatedClusterSize || !updatedClusters || !updatedClusters[0]) {
		printf("Error: Cannot calloc space for the new cluster variables");
		exit(1);
	}

	//reset memeber ship
	membership[0] = -1;
	for (i = 1; i < numberofTotalData; i++) {
		updatedClusters[i] = updatedClusters[i - 1] + numberofCoordinates;
		membership[i] = -1;
	}

	//get the total data number

	float delta;
	delta = FLT_MAX;
	int index, differences;
	int iterations;
	iterations = 0;
	while (delta > stopthreshold && iterations < MY_MAXITER) {
		iterations++;
//		may use the Wtime to record computing time
//		double time = MPI_Wtime();

		delta = 0.0;
		for (i = 0; i < numberofTotalData; i++) {
			index = find_NN(data[i], clusters, numberofClusters,
					numberofCoordinates);
			if (index < 0) {
				printf("Error: mistake in finding nearest cluster.");
				exit(1);
			}
			if (index != membership[i])
				differences++;
			membership[i] = index;
			updatedClusterSize[index]++;
			for (j = 0; j < numberofCoordinates; j++) {
				updatedClusters[index][j] += data[i][j];
			}
		}

		//compute the new cluster center
		for (i = 0; i < numberofClusters; i++) {
			for (j = 0; j < numberofCoordinates; j++) {
				if (updatedClusterSize[i] > 0) {
					clusters[i][j] /= updatedClusterSize[i];
				}
				updatedClusters[i][j] = 0;
			}
			updatedClusterSize[i] = 0;
		}
		delta = differences / (double) numberofTotalData;
	}

	free(updatedClusters);
	free(updatedClusters[0]);
	free(updatedClusterSize);
	return 1;
}
