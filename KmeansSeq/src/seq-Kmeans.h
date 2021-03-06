/*
 * seq-Kmeans.h
 *
 *  Created on: May 1, 2013
 *      Author: guanyuw
 */

#ifndef SEQ_KMEANS_H_
#define NORMDATA 0
#define DNADATA 1
#define DNATYPENUM 4
#define SEQ_KMEANS_H_


float** kmeans_read(char*, int, int);
int     kmeans(int, float**, int, int, int, float, int*, float**);
int     kmeans_write(char*, int, int, int, int, float**, int*, int);


#endif /* SEQ_KMEANS_H_ */
