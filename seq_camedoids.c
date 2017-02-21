#include "seq_camedoids.h"



/* Euclid distance between two multi-dimensional points  */
double euclid (int n,/* no. objects */
			   double** data,/*[numObjs][numCoords] objects */
			   int index1,/*position of the first element*/
			   double* index2/*attributes's vector of the second element*/)
{
  double result = 0.;
  int i;

  for (i = 0; i < n; i++)
    { 
         double term = data[index1][i] - index2[i];

         result += term*term;

    }
 
  result=sqrt(result);

  return result;
}


/* Initialize medoids with random elements of the dataset */
void init_clusters(int nclusters,/* no. clusters */
				   int nelements,/* no. objects */
				   int *centroids /*[numClusters] initial medoids position */)
{

  int i,j,k;

  srand(time(NULL));	   
  
  for(i=0;i<nclusters;i++)
		centroids[i]=rand()%(nelements-1);

}



/* find the cluster id that has min distance to object */
int find_nearest_cluster(int numClusters, /* no. clusters */
						 int ncoords, /* no. coordinates */
						 int object, /* position of the specifi object */
						 double **objects /*[numObjs][numCoords] */,
						 double **cluster /*[numClusters][numCoords] medoids*/)       
{
    int   index, i;
    double dist, min_dist;

    index    = 0;
    min_dist = euclid(ncoords,objects,object,cluster[0]); 

    for (i=1; i<numClusters; i++) {
        
	    dist=euclid(ncoords,objects,object,cluster[i]); 
   
		if (dist < min_dist)
		{ 
            min_dist = dist;
            index    = i;
        }

    }
	
    return(index);
}



double seq_camedoids(double **objects /* in: [numObjs][numCoords]*/,
					int numObjs /* no. objects */,
					int numCoords /* no. features */,
					int numClusters /* no. clusters */,
					int nCands /* no. candidates */,
					int *membership /* out: [numObjs] */,
					double **clusters /* out: [numClusters][numCoords] medoids (position of medoids into the dataset) */)          
{
    int    i,j,k,rep, index, loop=0, delta;          
    int pos;
	int     *clusterSize; /* [numClusters]: no. objects assigned in each new cluster */
	
    double **clusters_means;  /* [numClusters][numCoords]  cluster's means */
	
	double error=0,*errors,d;
	
    int **medoid_cand;   /* [numClusters][nCands]  medoid candidates for each cluster */
	double **cand_sum_dis; /* [numClusters][nCands]  accumulative distance of each medoid candidate to all other elements of its cluster*/
	double **cand_dis;   /* [numClusters][nCands]  distance of each medoid candidate to the mean of its cluster*/
	double *max_cand;    /* [numClusters] candidate with the max distance of every other candidate in the same cluster*/
	int *cand;           /* [numClusters] counter of candidates that have already been found for every cluster*/
	int *max_cand_pos;   /* [numClusters] posistion of the candidate with the max distance of every other candidate in the same cluster*/
    double ***final_cand; /* [numClusters][nCands][numCoords] the final medoids' candidates for every cluster*/
	
	cand=(int*)calloc(numClusters , sizeof(int));

    max_cand_pos=(int*)calloc(numClusters , sizeof(int));
  
    medoid_cand=(int**)malloc(numClusters * sizeof(int*));
    for(j=0;j<numClusters;j++)
          medoid_cand[j]=(int*)calloc(nCands , sizeof(int));

	final_cand=(double***)malloc(numClusters * sizeof(double**));
    for(j=0;j<numClusters;j++)
    {
        final_cand[j]=(double**)malloc(nCands * sizeof(double*));
		for(i=0;i<nCands;i++)
			final_cand[j][i]=(double*)malloc(numCoords * sizeof(double));
    }

	cand_dis=(double**)malloc(numClusters * sizeof(double*));
    for(j=0;j<numClusters;j++)
        cand_dis[j]=(double*)calloc(nCands , sizeof(double));

	
    max_cand=(double*)calloc(numClusters , sizeof(double));
	
	cand_sum_dis=(double**)malloc(numClusters * sizeof(double*));
	for(i=0;i<numClusters;i++)
    {
	   cand_sum_dis[i]=(double*)calloc(nCands , sizeof(double));

	}
	

	errors=(double*)malloc(numClusters * sizeof(double));

	
	clusters_means    = (double**) malloc(numClusters * sizeof(double*));
	for (i=0; i<numClusters; i++)
        clusters_means[i] = (double*) calloc(numCoords , sizeof(double));
	

    clusterSize = (int*) calloc(numClusters , sizeof(int));

    

    for (i=0; i<numObjs; i++) membership[i] = -1;
    

	do {
        
			delta=0;//counter gia tin loopa
            error=0;
			
           
            for (i=0; i<numObjs; i++)
			{
                /* find the array index of nestest cluster center */   
                index = find_nearest_cluster(numClusters, numCoords,i,objects, clusters);
     
	            /* if membership changes, increase delta by 1 */
                if (membership[i] != index) delta += 1.0;
       
	            /* assign the membership to object i */
                membership[i] = index;
                
				/* update new cluster center : sum of objects located within */
				clusterSize[index]++;
			
				for (j=0; j<numCoords; j++)
                    clusters_means[index][j] += objects[i][j];
					    
            }
               
			/* average the sum and replace old cluster center with newClusters */
			for (i=0; i<numClusters; i++) 
			{

				if(clusterSize[i]>1)
				{
					for(k=0;k<numCoords; k++)//gia na orisoume candidates
						clusters_means[i][k]/=clusterSize[i];
                }
			}				 

			/* find nCands candidates for every cluster  */	
			for (i = 0; i < numObjs; i++)
			{
 
				d=0;
    	  
		        /*calculate the distance of the current element of its cluster mean*/
				for (j = 0; j < numCoords; j++)
				{ 
					double term = objects[i][j] - clusters_means[membership[i]][j];

					d += term*term;
				}
				
				d=sqrt(d);
				
                 				
				if(cand[membership[i]]<nCands) /*if candidates of the cluster < nCands */
				{/*add the object to the cluster's candidate */
						medoid_cand[membership[i]][cand[membership[i]]]=i;
						cand_dis[membership[i]][cand[membership[i]]]=d;
						if(d>max_cand[membership[i]])
						{
							max_cand[membership[i]]=d;
							max_cand_pos[membership[i]]=cand[membership[i]];
						}
						cand[membership[i]]++;
				}
				else if(d<max_cand[membership[i]]) /*or distance of the current element smaller than the biggest distance of cluster's candidates so far*/
				{/*add the object to the cluster's candidate */
						medoid_cand[membership[i]][max_cand_pos[membership[i]]]=i;
						cand_dis[membership[i]][max_cand_pos[membership[i]]]=d;
		 
						max_cand[membership[i]]=cand_dis[membership[i]][0];
						max_cand_pos[membership[i]]=0;
						for(j=1;j<nCands;j++)
							if(cand_dis[membership[i]][j]>max_cand[membership[i]])
							{
								max_cand[membership[i]]=cand_dis[membership[i]][j];  
								max_cand_pos[membership[i]]=j;
							}	
				}
					
				
			}
			
			for(i=0;i<numClusters;i++)
				for (j = 0; j <nCands; j++)
			        for(k=0;k<numCoords;k++)
						final_cand[i][j][k]=objects[medoid_cand[i][j]][k];
  
	        /* calculate the accumulative distance of every candidate to the other objects of its cluster*/
		   	for (i = 0; i < numObjs; i++)
			{
				for (j = 0; j <nCands; j++)
				{ 
					cand_sum_dis[membership[i]][j]+=euclid(numCoords,objects,i,final_cand[membership[i]][j]); 			
				}
			}

		    /* set the candidate with the smallest accumulative distance as the medoid of its cluster */
			for(i=0;i<numClusters;i++)
			{ 
				
				errors[i]=cand_sum_dis[i][0];
                pos=0;
				
				
				for(j=1;j<nCands;j++)
				{
					if(cand_sum_dis[i][j]<errors[i])
					{
						errors[i]=cand_sum_dis[i][j];
						pos=j;
					}	

				}
				
				for(j=0;j<numCoords;j++)
						clusters[i][j]=final_cand[i][pos][j];
						
			}   
               
			 
			/* initialize for the next loop*/
            for(i=0;i<numClusters;i++)
			{
				for(j=0;j<numCoords;j++)
                    clusters_means[i][j]=0;					
            }
					
			for(j=0;j<numClusters;j++)
			{
				cand[j]=0;
				max_cand[j]=0;
				for(k=0;k<nCands;k++)
					cand_dis[j][k]=0;
		
			}
              		
            for(i=0;i<numClusters;i++)
			{
				clusterSize[i]=0;
                for(j=0;j<nCands;j++)
                    cand_sum_dis[i][j]=0;					
            }
			
	} while (delta!=0 && loop++ < 500); /*until no change to the membership of loops>500 */
	
	for(i=0;i<numClusters;i++)
		error+=errors[i];
    
	for(j=0;j<numClusters;j++)
    {
		for(i=0;i<nCands;i++)
           free(final_cand[j][i]);
	    free(final_cand[j]);
    }   
    free(final_cand);
	for(i=0;i<numClusters;i++)
           free(clusters_means[i]);
	free(clusters_means);
	free(clusterSize);
	for(i=0;i<numClusters;i++)
           free(medoid_cand[i]);
	free(medoid_cand);
	for(i=0;i<numClusters;i++)
           free(cand_sum_dis[i]);
	free(cand_sum_dis);
	free(errors);
	free(max_cand_pos);
	for(i=0;i<numClusters;i++)
           free(cand_dis[i]);
	free(cand_dis);
	free(max_cand);
	free(cand);
	
    return error;
}
