/*
 * runKMeans.c
 *
 *  Created on: May 1, 2013
 *      Author: lzy
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include "seq-Kmeans.h"


void printUsage(char *name) {
	fprintf(stderr, "Usage: %s -p dataType -i inFile -o outFile -k #OfClusters -t threshold -l #ofLinesInInputFile -d #dimension\n", name);
}

int main(int argc, char **argv) {
	extern char *optarg;
	extern int optind, optopt;
	int c, ncluster = 4, nline, totalLine, ndim, i = 0, j = 0;
	int type;
	char *inFile, *outFile;
	float thres = 0.01;

	float **data;	// input data in this process
	float **centroid;	// all cluster centroids
	int *label;	// for each data point, find its new class label

	clock_t stime, etime, 	// whole system time
		stimeCluster, etimeCluster;	// maximum cluster time among all processes

	if(argc != 15) {
		printUsage(argv[0]);
		exit(EXIT_FAILURE);
	}

	while ((c = getopt(argc, argv, "p:i:o:k:t:l:d:")) != EOF) {
		switch (c) {
		case 'p':
			type = atoi(optarg);
			break;
		case 'i':
			inFile = optarg;
			break;
		case 'o':
			outFile = optarg;
			break;
		case 'k':
			ncluster = atoi(optarg);
			break;
		case 't':
			thres = atof(optarg);
			break;
		case 'l':
			nline = atof(optarg);
			break;
		case 'd':
			ndim = atof(optarg);
			break;
		default:
			printUsage(argv[0]);
			exit(EXIT_FAILURE);
		}
	}

	stime = time(NULL);

	printf("init done\n");

	// read input data (each process has a portion of all data)
	totalLine = nline;
	data = kmeans_read(inFile, nline, ndim);
	printf("read data done.\n");

	// initialize cluster centers
	centroid = (float **) malloc(ncluster * sizeof(float *));	// pointer to each line
	centroid[0] = (float *) malloc(ncluster * ndim * sizeof(float));
	for(i = 1; i < ncluster; i++)
		centroid[i] = centroid[i-1] + ndim;

	for(i = 0; i < ncluster; i++)
		for(j = 0; j < ndim; j++)
			centroid[i][j] = data[i][j];
	printf("init cluster center done\n");

	// do kmeans calculation
	stimeCluster = time(NULL);
	label = (int *) malloc(nline * sizeof(int));

	//printf("rank:%d prior kmeans\n", rank);
	kmeans(type, data, ncluster, ndim, nline, thres, label, centroid);
	printf("kmeans done\n");
	etimeCluster = time(NULL);

	kmeans_write(outFile, nline, totalLine, ncluster, ndim , centroid, label, 0);

	free(label);
	free(centroid[0]);
	free(centroid);
	etime = time(NULL);

	// performance report
	printf("System time: %f secs\tClustering time: %f secs\n", (double)(etime-stime), (double)(etimeCluster-stimeCluster));

	return EXIT_SUCCESS;
}
