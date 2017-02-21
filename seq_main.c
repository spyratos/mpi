#include "seq_camedoids.h"
#include "file_io.h"



static void usage(char *argv0) {
    char *help =
        "Usage: %s [switches] -i filename -n num_clusters -k num_cands\n"
        "       -f filename      : file containing data to be clustered\n"
        "       -c num_clusters  : number of clusters (K must > 1)\n"
		"       -i filename      : file containing the initial centroids \n"
        "       -l lines_to_skip : lines to be ignored from the beggining of the file\n"
        "       -a attr_to_skip  : attributes to be ignored (options 1:first, 2:last, 3: first & last)\n"
		"       -k num_cands     : number of candidates\n"
        "       -h               : print this help information\n";
    fprintf(stderr, help, argv0);
    exit(-1);
}



int main(int argc, char **argv) {
          
    
    int     i,j,opt;
	extern char   *optarg;
    int     numCoords, numObjs,numclusters=2;
	int numCands=100;
    double **objects;    /* [numObjs][numCoords] data objects */   
    double error;
    int *clusters;  	/* [numClusters] cluster medoids */
	double **medoids;
	int *clusterid;      /* [numObjs] membership */
    double time2=0, time3;
	char *input_file,*clusters_file;
    unsigned long start,end;
	struct timeval  tv1,tv2;	
    int lines_to_skip=0,attr_to_skip=0,clustersfromfile=0;
    
	while ( (opt=getopt(argc,argv,"f:c:i:l:a:t:k:h"))!= EOF) {
        switch (opt) {
            case 'f': input_file=optarg;
                      break;
            case 'c': numclusters = atoi(optarg);
                      break;
		    case 'i': clustersfromfile=1;
			          clusters_file=optarg;
					  break;
            case 'l': lines_to_skip = atoi(optarg);
                      break;
            case 'a': attr_to_skip = atoi(optarg);
			          break;
			case 'k': numCands=atoi(optarg);
			          break;
            case 'h': 
            default: usage(argv[0]);
                      break;
        }
    }    
	    
	if(numclusters<=1)
    {
		printf("Too few clusters\n");
	    return;
	}
	
	/* read data points from file ------------------------------------------*/
    objects=file_read(input_file, &numObjs, &numCoords,lines_to_skip,attr_to_skip); 
    	
    if(objects==NULL)
         return;

	printf("::Objects loaded::\n");
    
	
	/* initialize cluster medoids ------------------------------------------*/
	clusterid = (int*)malloc(numObjs * sizeof(int));
	
	medoids=(double**)malloc(numclusters * sizeof(double*));
	for(i=0;i<numclusters;i++)
	        medoids[i]=(double*)malloc(numCoords * sizeof(double));
	
	if(clustersfromfile==0)
	{
		clusters = (int*)malloc(numclusters * sizeof(int));
		
		gettimeofday(&tv1, NULL);
	 
		init_clusters(numclusters,numObjs,clusters); 
   
		gettimeofday(&tv2, NULL);
		
		start = (unsigned long)(tv1.tv_usec + tv1.tv_sec * 1000000);
		end = (unsigned long)(tv2.tv_usec + tv2.tv_sec * 1000000);

		time2=(double)((end - start) / 1000000.0);
	}
	else
	{ 
		clusters=clusters_read(clusters_file,numclusters); 		
	}
	
    if(clusters==NULL)
          return;

	for(i=0;i<numclusters;i++)
		for(j=0;j<numCoords;j++)
	        medoids[i][j]=objects[clusters[i]][j];
    	  
     		
	printf("::Clusters initialized::\n");

	/* clustering  ------------------------------------------*/
    gettimeofday(&tv1, NULL);
    
	error=seq_camedoids(objects,numObjs,numCoords,numclusters,numCands,clusterid,medoids);
    
	gettimeofday(&tv2, NULL);
	
	start = (unsigned long)(tv1.tv_usec + tv1.tv_sec * 1000000);
	end = (unsigned long)(tv2.tv_usec + tv2.tv_sec * 1000000);

	time3=(double)((end - start) / 1000000.0);
    
	print_results(input_file,objects,clusterid,medoids,numObjs,numCoords,numclusters);

	printf("::Clustering done::\n\n\n");

    /* output performance and clean ------------------------------------------*/
	printf("--- Clustering info ---\n");
	printf("File: %s\n",input_file);
	printf("Initial clusters: ");
	if(clustersfromfile==1)
		printf("inserted from file\n");
	else
		printf("random\n");
	printf("Skipped lines: %d\n",lines_to_skip);
	printf("Ignored attribute: ");
	switch (attr_to_skip) {
            case NONE: printf("none\n");
                      break;
            case FIRST: printf("first\n");
                      break;
            case LAST: printf("last\n");
                      break;
            case BOTH: printf("first and last\n");
                      break;
            default: printf("\n");
        }
    printf("Objects: %d\n",numObjs);
    printf("Attributes: %d\n",numCoords);
    printf("Clusters: %d\n",numclusters);  
    printf("Number of candidates: %d\n\n",numCands);	
	
 
	printf("--- Results: ---\n");
	printf("Error: %f\n",error);
	printf("Time for clusters' initialization: %lf\n", time2);
	printf("Time for clustering: %lf\n", time3);
	printf("Total time: %lf\n",time2+time3);
	
	
  
    free(objects);
	free(clusters);
	free(clusterid);

   
    return(0);
}

