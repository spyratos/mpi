#ifndef SEQ_CAMEDOIDS
#define SEQ_CAMEDOIDS

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <sys/time.h>
#include <float.h>


void init_clusters(int nclusters, int nelements,int centroids[]);
double euclid (int n,double** data,int index1,double* index2);
int find_nearest_cluster(int numClusters,int ncoords,int object,double **objects,double **cluster);
double seq_camedoids(double **objects,int numObjs,int numCoords,int numClusters,int nCands,int *membership,double **clusters);  

#endif


