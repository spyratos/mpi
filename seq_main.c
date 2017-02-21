#include "file_io.h"
#include "seq_camedoids.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define FLAG 0

static void usage(char *argv0) {
  char *help =
      "Usage: %s [switches] -i filename -n num_clusters -k num_cands\n"
      "       -f filename      : file containing data to be clustered\n"
      "       -c num_clusters  : number of clusters (K must > 1)\n"
      "       -i filename      : file containing the initial centroids \n"
      "       -l lines_to_skip : lines to be ignored from the beggining of the "
      "file\n"
      "       -a attr_to_skip  : attributes to be ignored (options 1:first, "
      "2:last, 3: first & last)\n"
      "       -k num_cands     : number of candidates\n"
      "       -h               : print this help information\n";
  fprintf(stderr, help, argv0);
  exit(-1);
}

/******************************************************************/
void p() { printf("\n"); }
void print_array_1d_double(double *A, int height, int width, int rank) {
  int count = 0;
  for (int i = 0; i < height * width; i++) {
    printf("R:%d - %f\t", rank, A[i]);
    count++;
    if (count == width) {
      printf("\n");
      count = 0;
    }
  }
  p();
}

void print_array_1d_int(int *A, int height, int width, int rank) {
  int count = 0;
  for (int i = 0; i < height * width; i++) {
    printf("R:%d - %d\t", rank, A[i]);
    count++;
    if (count == width) {
      printf("\n");
      count = 0;
    }
  }
  p();
}

void print_array(double *A, size_t height, size_t width, int process) {
  int count = 0;
  for (size_t i = 0; i < height; ++i) {
    for (size_t j = 0; j < width; ++j) {

      printf("%f,", A[i * width + j]);
    }
  }
  p();
  p();
}

/*************************************************************************/

int main(int argc, char **argv) {
/***seq_camedoids***/
double error = 0, *errors;

int **medoid_cand;

double **cand_sum_dis;

double **cand_dis;
double *max_cand;
int *cand;
int *max_cand_pos;
/* [numclusters] posistion of the candidate with the max
          distance of every other candidate in the same cluster*/
double ***final_cand;
/* [numclusters][numCands][numCoords] the final medoids'
          candidates for every cluster*/

  // MY VAR initial
  int extraObjs, size, my_rank, sliceSize, sliceSizeExt, target, startSlice,
      endSlice;
  int *a;

  int i, j, opt,pos;
  extern char *optarg;
  int numCoords, numObjs, numclusters = 2;
  int numCands = 100;
  double **objects; /* [numObjs][numCoords] data objects */

  int *clusters; /* [numclusters] cluster medoids */
  // double **medoids;
  int *clusterid, *clusterid_extr_array; /* [numObjs] membership */
  double time2 = 0, time3;
  char *input_file, *clusters_file;
  unsigned long start, end;
  struct timeval tv1, tv2;
  int lines_to_skip = 0, attr_to_skip = 0, clustersfromfile = 0;
  // Command line arguments
  while ((opt = getopt(argc, argv, "f:c:i:l:a:t:k:h")) != EOF) {
    switch (opt) {
    case 'f':
      input_file = optarg;
      break;
    case 'c':
      numclusters = atoi(optarg);
      break;
    case 'i':
      clustersfromfile = 1;
      clusters_file = optarg;
      break;
    case 'l':
      lines_to_skip = atoi(optarg);
      break;
    case 'a':
      attr_to_skip = atoi(optarg);
      break;
    case 'k':
      numCands = atoi(optarg);
      break;
    case 'h':
    default:
      usage(argv[0]);
      break;
    }
  }

  if (numclusters <= 1) {
    printf("Too few clusters\n");
    return -1;
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status status;

  // dataset slicing!
  a = malloc(2 * sizeof(int));
  int *clusterSize, *clusteridAll, *CLUSTER_SIZE, delta, *dlt, *displs,
      *recvcount, sumdelta = -2, k;
  int *flag = FLAG, loop;
  double *distance, *distance_extr_array, *distance_SUM;

  /* read data points from file ------------------------------------------*/
  if (my_rank == 0) {

    objects = file_read(input_file, &numObjs, &numCoords, lines_to_skip,
                        attr_to_skip);

    // printf("\n Rank %d , Calculating Slice Size: %d", my_rank, sliceSize);
    if (objects == NULL)
      return -1;

    printf("::Objects loaded::\n");

    a[0] = numObjs;
    a[1] = numCoords;
  }

  MPI_Bcast(a, 2, MPI_INT, 0, MPI_COMM_WORLD);
  numObjs = a[0];
  numCoords = a[1];

  extraObjs = numObjs % size;
  sliceSize = numObjs / size;
  sliceSizeExt = (numObjs / size) + extraObjs;





  double(*extr_array)[numCoords] = malloc(sizeof(double[(sliceSize)][numCoords]));
  double(*medoids)[numCoords] =
      malloc(sizeof(double[(numclusters)][numCoords]));
  clusterid = (int *)malloc(sliceSize * sizeof(int));
  clusterid_extr_array = (int *)malloc(sliceSizeExt * sizeof(int));
  clusteridAll = (int *)malloc(numObjs * sizeof(int));
  clusterSize = (int *)calloc(numclusters, sizeof(int));
  distance = calloc(sliceSize, sizeof(double));
  distance_extr_array = calloc(sliceSizeExt, sizeof(double));
  distance_SUM = calloc(numObjs, sizeof(double));
  displs = calloc(size, sizeof(int));
  recvcount = calloc(size, sizeof(int));
  CLUSTER_SIZE = (int *)calloc(numclusters, sizeof(int));
  /* medoids = (double **)malloc(numclusters * sizeof(double *));
   for (i = 0; i < numclusters; i++)
     medoids[i] = (double *)malloc(numCoords * sizeof(double));*/
  double(*cluster_distance_sum)[numCoords] =
      malloc(sizeof(double[(numclusters)][numCoords]));
  double(*clusters_means)[numCoords] =
      malloc(sizeof(double[(numclusters)][numCoords]));
  double(*CLUSTER_MEANS)[numCoords] =
      malloc(sizeof(double[(numclusters)][numCoords]));
      if (my_rank == 0 ){

       final_cand = (double ***)malloc(numclusters * sizeof(double **));
       for (j = 0; j < numclusters; j++) {
         final_cand[j] = (double **)malloc(numCands * sizeof(double *));
         for (i = 0; i < numCands; i++)
           final_cand[j][i] = (double *)malloc(numCoords * sizeof(double));
       }

       cand = (int *)calloc(numclusters, sizeof(int));

        medoid_cand = (int **)malloc(numclusters * sizeof(int *));
        for (j = 0; j < numclusters; j++){
          medoid_cand[j] = (int *)calloc(numCands, sizeof(int));
         }
          max_cand_pos = (int *)calloc(numclusters, sizeof(int));
            cand = (int *)calloc(numclusters, sizeof(int));

            max_cand_pos = (int *)calloc(numclusters, sizeof(int));

            medoid_cand = (int **)malloc(numclusters * sizeof(int *));
            for (j = 0; j < numclusters; j++){
              medoid_cand[j] = (int *)calloc(numCands, sizeof(int));
          }


            cand_dis = (double **)malloc(numclusters * sizeof(double *));
            for (j = 0; j < numclusters; j++){
              cand_dis[j] = (double *)calloc(numCands, sizeof(double));
          }
            max_cand = (double *)calloc(numclusters, sizeof(double));

            cand_sum_dis = (double **)malloc(numclusters * sizeof(double *));
            for (i = 0; i < numclusters; i++) {
              cand_sum_dis[i] = (double *)calloc(numCands, sizeof(double));
            }

            errors = (double *)malloc(numclusters * sizeof(double));
      }

  if (my_rank == 0) {
    clusters = clusters_read(clusters_file, numclusters);
    if (clusters == NULL)
      return -1;
    for (i = 0; i < numclusters; i++) {
      for (j = 0; j < numCoords; j++) {
        medoids[i][j] = objects[clusters[i]][j];
      }
    }
  }
  double(*rank0_extr_array)[numCoords] =
      malloc(sizeof(double[sliceSizeExt][numCoords]));

  if (my_rank == 0) {

    // RANDOM MEDOIDS - FIRST TIME
    // MASTER PROCESS CODE - SLicing!
    for (target = 0; target < size; target++) {
      if (target == 0) {
        startSlice = target * sliceSize;
        endSlice = startSlice + sliceSizeExt - 1;
      } else if (target < size) {
        startSlice = target * sliceSize + extraObjs;
        endSlice = startSlice + sliceSize - 1;
      }

      for (i = startSlice; i <= endSlice; i++) {
        for (j = 0; j < numCoords; j++) {
          if (target == 0) {

            rank0_extr_array[i - startSlice][j] = objects[i][j];
          } else if (target > 0 && target < size) {
            extr_array[i - startSlice][j] = objects[i][j];
          }
        }
      }

      if (target > 0 && target < size) {
        MPI_Send(extr_array, sliceSize * numCoords, MPI_DOUBLE, target, 10,
                 MPI_COMM_WORLD);
      }
    }

  } else if (my_rank > 0 && my_rank < size) {

    MPI_Recv(extr_array, sliceSize * numCoords, MPI_DOUBLE, target, 10,
             MPI_COMM_WORLD, &status);

    /* [numclusters][numCoords]  cluster's distance sum */
  }

  // COLLECTIVE COMMUNICATION!
  // for (int duo; duo < 2; duo++) {

  MPI_Barrier(MPI_COMM_WORLD);
  time_t t0, t1;
  if (my_rank == 0) {
    gettimeofday(&t0, 0);
  }
  for (int REC = 0; REC < 500; REC++) {
    MPI_Bcast(medoids, numclusters * numCoords, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {

      // RANDOM MEDOIDS - FIRST TIME
      // MASTER PROCESS CODE - SLicing!
      sumdelta = 0;

      cluster_distance_sum =
          seq_camedoids2(sliceSizeExt, numCoords, rank0_extr_array, numclusters,
                         clusterid_extr_array, clusterSize, medoids, &dlt, &flag);
      delta = dlt;


    } else if (my_rank > 0 && my_rank < size) {

      /* [numclusters][numCoords]  cluster's distance sum */

      cluster_distance_sum =
          seq_camedoids2(sliceSize, numCoords, extr_array, numclusters, clusterid,
                         clusterSize, medoids, &dlt, &flag);

      delta = dlt;

    }
    MPI_Allreduce(cluster_distance_sum, CLUSTER_MEANS, numclusters * numCoords,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(clusterSize, CLUSTER_SIZE, numclusters, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);

    if (my_rank == 0) {

      for (int qq = 0; qq < numclusters; qq++) {
        if (CLUSTER_SIZE[qq] > 1) {
          for (k = 0; k < numCoords; k++)
            CLUSTER_MEANS[qq][k] /= CLUSTER_SIZE[qq];
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&delta, &sumdelta, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Bcast(&sumdelta, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(CLUSTER_MEANS, numclusters * numCoords, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

    if (my_rank == 0) {
      distance_extr_array = mean_distance(sliceSizeExt, numCoords, rank0_extr_array,
                                   clusterid_extr_array, CLUSTER_MEANS);

      for (int dcounter = 1; dcounter < size; dcounter++) {
        recvcount[dcounter] = sliceSize;
        displs[dcounter] = dcounter * sliceSize + extraObjs;
      }
      recvcount[0] = sliceSizeExt;
      displs[0] = 0;
    } else if (my_rank > 0 && my_rank < size) {

      distance =
          mean_distance(sliceSize, numCoords, extr_array, clusterid, CLUSTER_MEANS);
    }

    MPI_Bcast(recvcount, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0) {
      MPI_Gatherv(distance_extr_array, recvcount[my_rank], MPI_DOUBLE, distance_SUM,
                  recvcount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      MPI_Gatherv(clusterid_extr_array, recvcount[my_rank], MPI_INT, clusteridAll,
                  recvcount, displs, MPI_INT, 0, MPI_COMM_WORLD);
    } else if (my_rank > 0 && my_rank < size) {
      MPI_Gatherv(distance, recvcount[my_rank], MPI_DOUBLE, distance_SUM,
                  recvcount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(clusterid, recvcount[my_rank], MPI_INT, clusteridAll,
                  recvcount, displs, MPI_INT, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    if (my_rank == 0) {
     error = 0;

                      for (i = 0; i < numObjs; i++) {



                        if (cand[clusteridAll[i]] < numCands) {

                          medoid_cand[clusteridAll[i]][cand[clusteridAll[i]]] = i;
                          cand_dis[clusteridAll[i]][cand[clusteridAll[i]]] = distance_SUM[i];
                          if (distance_SUM[i] > max_cand[clusteridAll[i]]) {

                            max_cand[clusteridAll[i]] = distance_SUM[i];
                            max_cand_pos[clusteridAll[i]] = cand[clusteridAll[i]];
                          }
                          cand[clusteridAll[i]]++;
                        } else if (distance_SUM[i] < max_cand[clusteridAll[i]])
                        // or distance of the current element smaller than the biggest //distance
                        // of
                        //   cluster's candidates so far
                        {

                          medoid_cand[clusteridAll[i]][max_cand_pos[clusteridAll[i]]] = i;
                          cand_dis[clusteridAll[i]][max_cand_pos[clusteridAll[i]]] = distance_SUM[i];

                          max_cand[clusteridAll[i]] = cand_dis[clusteridAll[i]][0];
                          max_cand_pos[clusteridAll[i]] = 0;
                          for (j = 1; j < numCands; j++)

                            if (cand_dis[clusteridAll[i]][j] > max_cand[clusteridAll[i]]) {

                              max_cand[clusteridAll[i]] = cand_dis[clusteridAll[i]][j];
                              max_cand_pos[clusteridAll[i]] = j;
                            }
                        }
                      }

                      /* 3D */
                      for (i = 0; i < numclusters; i++)
                        for (j = 0; j < numCands; j++)
                          for (k = 0; k < numCoords; k++)
                            final_cand[i][j][k] = objects[medoid_cand[i][j]][k];

                      /* calculate the accumulative distance of every candidate to the other
                       * objects of its cluster*/
                      for (i = 0; i < numObjs; i++) {
                        for (j = 0; j < numCands; j++) {

                          cand_sum_dis[clusteridAll[i]][j] +=
                              euclid_dynamic(objects, numCoords, i, final_cand[clusteridAll[i]][j]);
                        }
                      }

                    

                      for (i = 0; i < numclusters; i++) {

                        errors[i] = cand_sum_dis[i][0];
                        pos = 0;

                        for (j = 1; j < numCands; j++) {
                          if (cand_sum_dis[i][j] < errors[i]) {
                            errors[i] = cand_sum_dis[i][j];
                            pos = j;
                          }
                        }

                        for (j = 0; j < numCoords; j++)
                          medoids[i][j] = final_cand[i][pos][j];
                      }


                      for (j = 0; j < numclusters; j++) {
                        cand[j] = 0;
                        max_cand[j] = 0;
                        for (k = 0; k < numCands; k++)
                          cand_dis[j][k] = 0;
                      }

                      for (i = 0; i < numclusters; i++) {

                        for (j = 0; j < numCands; j++)
                          cand_sum_dis[i][j] = 0;
                      }

                      // to ypologizei kathe fora ara if sumdelta == 0 den ypologizei tipota.

                      for (i = 0; i < numclusters; i++) {
                        error += errors[i];

                      }






                    }


    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < numclusters; i++) {
      clusterSize[i] = 0;
      CLUSTER_SIZE[i] = 0;
      for (j = 0; j < numCoords; j++) {
        CLUSTER_MEANS[i][j] = 0;
        cluster_distance_sum[i][j] = 0;
      }
    }

    // MPI_Bcast(&sumdelta, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /*MPI_Bcast(&sumdelta, 1, MPI_INT, 0, MPI_COMM_WORLD);*/
    if (sumdelta == 0) {
      break;
    }
  }


  if (my_rank == 0) {
    print_results(input_file, numCoords, objects, clusteridAll, medoids,
                  numObjs, numclusters);
  }
  /* output performance and clean
     ------------------------------------------*/
  if (my_rank == 0) {
    printf("::Clustering done::\n\n\n");
    printf("--- Clustering info ---\n");
    printf("File: %s\n", input_file);
    printf("Initial clusters: ");
    if (clustersfromfile == 1)
      printf("inserted from file\n");
    else
      printf("random\n");
    printf("Skipped lines: %d\n", lines_to_skip);
    printf("Ignored attribute: ");
    switch (attr_to_skip) {
    case NONE:
      printf("none\n");
      break;
    case FIRST:
      printf("first\n");
      break;
    case LAST:
      printf("last\n");
      break;
    case BOTH:
      printf("first and last\n");
      break;
    default:
      printf("\n");
    }
    printf("Objects: %d\n", numObjs);
    printf("Attributes: %d\n", numCoords);
    printf("Clusters: %d\n", numclusters);
    printf("Number of candidates: %d\n\n", numCands);

    printf("--- Results: ---\n");
    printf("Error: %f\n", error);
    printf("Time for clusters' initialization: %lf\n", time2);
    printf("Time for clustering: %lf\n", time3);
    gettimeofday(&t1, 0);
    long elapsed = (t1 - t0) * 1000000 + t1 - t0;
    printf("Total time: (%ld) %d seconds\n",  elapsed, elapsed / 1000000);

    free(objects);
    free(clusters);
    free(clusterid);
    free(clusterid_extr_array);
    free(extr_array);
    free(rank0_extr_array);
    free(clusteridAll);
    free(CLUSTER_MEANS);
    free(CLUSTER_SIZE);
    free(cluster_distance_sum);
    free(clusters_means);
    free(distance);
    free(distance_extr_array);
    free(distance_SUM);

     for (j = 0; j < numclusters; j++) {
       for (i = 0; i < numCands; i++)
         free(final_cand[j][i]);
       free(final_cand[j]);
     }
     free(final_cand);

for (i = 0; i < numclusters; i++)
       free(medoid_cand[i]);
     free(medoid_cand);
     for (i = 0; i < numclusters; i++)
       free(cand_sum_dis[i]);
     free(cand_sum_dis);
     free(errors);
     free(max_cand_pos);
     for (i = 0; i < numclusters; i++)
       free(cand_dis[i]);
     free(cand_dis);
     free(max_cand);
     free(cand);
  }

  MPI_Finalize();

  return (0);
}
