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
alloc_grid(double ***Gptr, int N){
	int coeff = world_rank? P : 1;
    double **G = (double**)malloc((N / coeff + 2)  * sizeof(double*));
    if (G == NULL) {
        printf("malloc failed\n");
        exit(42);
    }

    for (int i = 0; i <= (N / coeff + 2); i++) {   /* malloc the own range plus one more line */
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
init_grid(double **G, int N)
{
	int coeff = world_rank? P : 1;
    /* initialize the grid */
    for (int i = 0; i <= (N / coeff); i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0 && ( world_rank == 0))
                G[i][j] = 4.56;
            else if (i == (N / coeff) && ((world_rank == P - 1) || !world_rank) )
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
print_grid(double **G, int N, int world_rank)
{
    int res;
    int coeff = world_rank? (N/P + 2) : N;
    for (int i = 0; i < (coeff) ; i++) {
        for (int j = 0; j < N ; j++) {
            printf("%10.3f ", G[i][j]);
        }
        printf("\n");
    }
}

void __sync_neigh(double **G, unsigned int N){
  MPI_Request request, request_2;
  MPI_Status status;
  MPI_Barrier(MPI_COMM_WORLD);
  if(world_rank)
      MPI_Isend(G[1], N, MPI_DOUBLE, world_rank - 1, MY_ITER_TAG, MPI_COMM_WORLD, &request);
  if(world_rank != P - 1)
      MPI_Isend(G[N/P], N, MPI_DOUBLE, world_rank + 1, MY_ITER_TAG, MPI_COMM_WORLD, &request_2);
  if(world_rank)
      MPI_Recv(G[0], N, MPI_DOUBLE, world_rank - 1, MY_ITER_TAG, MPI_COMM_WORLD, &status);
  if(world_rank != P - 1)
      MPI_Recv(G[N/P], N, MPI_DOUBLE, world_rank + 1, MY_ITER_TAG, MPI_COMM_WORLD, &status);
  MPI_Barrier(MPI_COMM_WORLD);
}

void __sync_matrix(double **G, unsigned int N){
    MPI_Status status;
  if(!world_rank){
    // Master node
    double new_row[N];
    for(int i = 1; i < P; ++i)
        for(int j = 1; j <= N/P; ++j)
            MPI_Recv(G[i * N/P + j - 1], N, MPI_DOUBLE, i, MY_FINAL_TAG, MPI_COMM_WORLD, &status);
  } else {
     for(int j = 1; j <= N/P; ++j)
        MPI_Send((double*)G[j], N, MPI_DOUBLE, 0, MY_FINAL_TAG, MPI_COMM_WORLD);
  }
}

int main(int argc, char** argv) {
    double **G;
    int         ncol, nrow;     /* number of rows and columns */
    double      Gnew;
    double      r;
    double      omega;
    double      stopdiff;
    double      maxdiff, glob_maxdiff;
    double      diff;
    int         iteration; /* counters */
    double start;
    double end;
    double      time;
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    MPI_Status source;

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    P = world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    unsigned int N = 10;
    if(!world_rank)
    	printf("Running %d x %d SOR\n", N, N);
    ncol = nrow = N;
        r = 0.5 * (cos(M_PI / ncol) + cos(M_PI / nrow));
        omega = 2.0 / (1 + sqrt(1 - r * r));
        stopdiff = TOLERANCE / (2.0 - omega);
        omega *= MAGIC;
    if(!world_rank){
        // Master node
        int msg = 0;
        /* set up a quadratic grid */
        alloc_grid(&G, N);
        init_grid(G, N);

    } else {
        // Slave node
        int msg;
        alloc_grid(&G, N);
        init_grid(G, N);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    __sync_matrix(G, N);
    MPI_Barrier(MPI_COMM_WORLD);

    if(!world_rank){
    	start = MPI_Wtime();
    }


	MPI_Barrier(MPI_COMM_WORLD);
    // Parallel SOR:
    iteration = 0;
    do {
    	int borde_cond = world_rank == 0? 0 : 1;
    	int even_bound = world_rank == 0? 0 : 1;
    	if(world_rank){
        print_grid(G, N, world_rank);
        printf("\n");
      }



      maxdiff = 0.0;
      for (int phase = 0; phase < 2; phase++) {
        //  if(phase == 1){
            __sync_neigh(G, N);
        //  }
          for (int i = 1; i <= (N/P - borde_cond); i++) {
              for (int j = 1 + (even(i + ((P + 1) * N/P)) ^ phase); j < N - 1; j += 2) {
                  Gnew = stencil(G, i, j);
                  diff = fabs(Gnew - G[i][j]);
                  if (diff > maxdiff) {
                      maxdiff = diff;
                  }
                  G[i][j] = G[i][j] + omega * (Gnew - G[i][j]);
              }
          }
      }
    ++iteration;
    MPI_Barrier(MPI_COMM_WORLD);
    __sync_matrix(G, N);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&maxdiff, &glob_maxdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if(!world_rank)
      std::cout<<" ITERATION " << iteration << std::endl;
  } while (glob_maxdiff > stopdiff);


	if(!world_rank){
  	end = MPI_Wtime();
  	time = end - start;
    printf("SOR took %10.3f seconds\n", time);
    printf("Used %5d iterations, diff is %10.6f, allowed diff is %10.6f\n",
    	iteration, maxdiff, stopdiff);
  }
  __sync_matrix(G, N);
	if(!world_rank)
		print_grid(G, N, world_rank);
  // Finalize the MPI environment.
  MPI_Finalize();
  exit(0);
  return 0;
}
