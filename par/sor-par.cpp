#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#define TOLERANCE 0.00002       /* termination criterion */
#define MAGIC 0.8               /* magic factor */
#define MY_INIT_TAG 0xbad00
#define MY_FINAL_TAG 0xbad01
#define MY_ITER_TAG 0xbad02
#define MY_PRINT_TAG 0x003

int world_rank;
int P;

int on_master(){
  return world_rank == 0;
}

static int
even(int i)
{
    return !(i & 1);
}


static double
stencil(double **G, int row, int col)
{
    return (G[row - 1][col] + G[row + 1][col] +
            G[row][col - 1] + G[row][col + 1]) / 4.0;
}


static void
alloc_grid(double ***Gptr, int locN, unsigned int N)
{
    double    **G = (double**)malloc(locN * sizeof *G);
    if (G == NULL) {
        printf("malloc failed\n");
        exit(42);
    }

    for (int i = 0; i < locN; i++) {   /* malloc the own range plus one more line */
        /* of overlap on each side */
        G[i] = (double*)malloc(N * sizeof *G[i]);
        if (G[i] == NULL) {
            printf("malloc failed\n");
            exit(42);
        }
    }

    *Gptr = G;
}


static void
init_grid(double **G, int locN, int N)
{
    /* initialize the grid */
    for (int i = 0; i < locN; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0 && !world_rank)
                G[i][j] = 4.56;
            else if (i == (locN - 1) && world_rank == (P - 1))
                G[i][j] = 9.85;
            else if (j == 0)
                G[i][j] = 7.32;
            else if (j == N - 1)
                G[i][j] = 6.88;
            else
                G[i][j] = 0.0;
        }
    }
}


void
print_grid(double **G, int locN, int N)
{
    for (int i = 0; i < locN; i++) {
        for (int j = 0; j < N ; j++) {
            printf("%10.3f ", G[i][j]);
        }
        printf("\n");
    }
}



void __sync_neigh(double **G, unsigned int locN, unsigned int N){
  MPI_Request request, request_2;
  MPI_Status status;
  if(world_rank)
      MPI_Isend(G[1], N, MPI_DOUBLE, world_rank - 1, MY_ITER_TAG, MPI_COMM_WORLD, &request);
  if(world_rank != P - 1)
      MPI_Isend(G[locN - 2], N, MPI_DOUBLE, world_rank + 1, MY_ITER_TAG, MPI_COMM_WORLD, &request_2);
  if(world_rank)
      MPI_Recv(G[0], N, MPI_DOUBLE, world_rank - 1, MY_ITER_TAG, MPI_COMM_WORLD, &status);
  if(world_rank != P - 1)
      MPI_Recv(G[locN - 1], N, MPI_DOUBLE, world_rank + 1, MY_ITER_TAG, MPI_COMM_WORLD, &status);
}

void __sync_matrix(double **G, unsigned int locN, int N){
    MPI_Status status;
  if(!world_rank){
    // Master node
    for(int i = 1; i < P; ++i)
        for(int j = 0; j < locN - 1; ++j)
            MPI_Recv(G[i * (locN - 2) + j + 1], N, MPI_DOUBLE, i, MY_FINAL_TAG, MPI_COMM_WORLD, &status);
  } else {
     for(int j = 0; j < (locN - 1); ++j)
        MPI_Send((double*)G[j + 1], N, MPI_DOUBLE, 0, MY_FINAL_TAG, MPI_COMM_WORLD);
  }
}


double solve(unsigned int N, char print){

  int         ncol, nrow;     /* number of rows and columns */
  double      Gnew;
  double      r;
  unsigned int locN;
  double **G;
  double      omega;
  double start, end;
  double time;
  /* differences btw grid points in iters */
  double      stopdiff;
  double      maxdiff, glob_maxdiff;
  int         iteration; /* counters */
  double      diff;
  /* finally N*N (from argv) array points will be computed */

  /* set up a quadratic grid */
  ncol = nrow = N + 2;
  r = 0.5 * (cos(M_PI / ncol) + cos(M_PI / nrow));
  omega = 2.0 / (1 + sqrt(1 - r * r));
  stopdiff = TOLERANCE / (2.0 - omega);
  omega *= MAGIC;
  locN = N / P + 2; // Need two border lines;
  N += 2;
  if(on_master()){
    printf("Running %d x %d SOR\n", N - 2, N - 2);
    alloc_grid(&G, N, N);
  } else {
    alloc_grid(&G, locN, N);
  }
  init_grid(G, locN, N);

  start = MPI_Wtime();
  /* now do the "real" computation */
  iteration = 0;
  int border = locN - 1;
  int even_border = N/P * (world_rank);
  int big_border = N - 1;
  do {
      maxdiff = 0.0;
      for (int phase = 0; phase < 2; phase++) {
          __sync_neigh(G, locN, N);
          for (int i = 1; i < border ; i++) {
              for (int j = 1 + (even(i + (even_border)) ^ phase); j < big_border; j += 2) {
                  Gnew = stencil(G, i, j);
                  diff = fabs(Gnew - G[i][j]);
                  if (diff > maxdiff) {
                      maxdiff = diff;
                  }
                  G[i][j] = G[i][j] + omega * (Gnew - G[i][j]);

              }
          }
      }
      iteration++;
      //__sync_matrix(G, locN, N);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(&maxdiff, &glob_maxdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if(on_master())
        printf(" ITERATION %d\n", iteration);
  } while (glob_maxdiff > stopdiff);

  if(on_master()){
    end = MPI_Wtime();
    time = end - start;
    printf("SOR took %10.3f seconds\n", time);
    printf("Used %5d iterations, diff is %10.6f, allowed diff is %10.6f\n",
           iteration, glob_maxdiff, stopdiff);
    if (print == 1)
      print_grid(G, N, N);
    return end - start;
  } else
    return -1;

}

int
main(int argc, char *argv[])
{
    int         N;              /* problem size */
    double    **G;              /* the grid */
    double      time;
    char         print = 0;

    // Setting up MPI environment
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    /* set up problem size */
    N = 10;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-print") == 0) {
            print = 1;
        } else {
            if (sscanf(argv[i], "%d", &N) != 1) {
                printf("Positional parameter N must be an int, not '%s'\n",
                        argv[i]);
                exit(42);
            }
        }
    }

    time = solve(N, print);
    MPI_Finalize();
    return 0;
}
