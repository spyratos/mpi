#ifndef FILE_IO
#define FILE_IO

#include <stdio.h>
#include <stdlib.h>
#include <string.h>     
#include <sys/types.h>  
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>    
#include <errno.h>
#include <assert.h>

#define MAX_CHAR_PER_LINE 128

#define NONE 0
#define FIRST 1
#define LAST 2
#define BOTH 3 

int * clusters_read(char *filename,int numclusters);
double** file_read(char *filename,int  *numObjs,int  *numCoords, int lines_to_skip,int attr_to_skip);
void print_results(char *filename,double **objects, int *clusterid,double **clusters, int numObjs, int numCoords, int numclusters);      

#endif

