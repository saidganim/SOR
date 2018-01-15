#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#define TOLERANCE 0.00002       /* termination criterion */
#define MAGIC 0.8               /* magic factor */
#define P 2
#define MY_INIT_TAG 0xbad00
#define MY_FINAL_TAG 0xbad01
#define MY_ITER_TAG 0xbad02
#define MY_PRINT_TAG 0x003

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
alloc_grid(double ***Gptr, int N)
{
    double **G = (double**)malloc(N * sizeof(double*));
    if (G == NULL) {
        fprintf(stderr, "malloc failed\n");
        exit(42);
    }

    for (int i = 0; i < N; i++) {   /* malloc the own range plus one more line */
        /* of overlap on each side */
        G[i] = (double*)malloc(N * sizeof *G[i]);
        if (G[i] == NULL) {
            fprintf(stderr, "malloc failed\n");
            exit(42);
        }
    }

    *Gptr = G;
}

static void
init_grid(double **G, int N)
{
    /* initialize the grid */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0)
                G[i][j] = 4.56;
            else if (i == N - 1)
                G[i][j] = 9.85;
            else if (j == 0)
                G[i][j] = 7.32;
            else if (j == N - 1)
                G[i][j] = 6.88;
            else
                G[i][j] = 1.0;
        }
    }
}


void
print_grid(double **G, int N, int world_rank)
{
    int res;
    MPI_Status status;
    if(world_rank)
        MPI_Recv(&res, 1, MPI_INT, 0, MY_PRINT_TAG, MPI_COMM_WORLD, &status);
    printf("BEGINNING OF GRID OF %d\n", world_rank);
    for (int i = 1; i < N - 1; i++) {
        printf("GRID OF %d ", world_rank);
        for (int j = 1; j < N - 1; j++) {
            printf("%10.3f ", G[i][j]);
        }
        printf("\n");
    }
    printf("END OF GRID OF %d\n", world_rank);
    if(!world_rank)
        MPI_Send(&res, 1, MPI_INT, 1, MY_PRINT_TAG, MPI_COMM_WORLD);
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
    struct timeval start;
    struct timeval end;
    double      time;
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    MPI_Status source, status;
    MPI_Request request, request_2;
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    /* set up problem size */
    unsigned int N = 10;
    printf(" CPU FOUND %d\n", world_rank);
    N += 2;    
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
        if(gettimeofday(&start, 0) != 0) {
            fprintf(stderr, "could not do timing\n");
            exit(1);
        }
    } else {
        // Slave node
        int msg;
        alloc_grid(&G, N);
        init_grid(G, N);
    }
     if(world_rank == P - 1){
		memcpy(G[N/P], G[N - 1], N * sizeof(double));
     }
     //print_grid(G, N, world_rank);
    // Real Computation
    iteration = 0;
    do {
        if(world_rank)
            MPI_Isend(G[1], N, MPI_DOUBLE, world_rank - 1, MY_ITER_TAG, MPI_COMM_WORLD, &request);
        if(world_rank != P - 1)
            MPI_Isend(G[N/P - 1], N, MPI_DOUBLE, world_rank + 1, MY_ITER_TAG, MPI_COMM_WORLD, &request_2);
        if(world_rank)
            MPI_Recv(G[0], N, MPI_DOUBLE, world_rank - 1, MY_ITER_TAG, MPI_COMM_WORLD, &status);
        if(world_rank != P - 1)
            MPI_Recv(G[N/P], N, MPI_DOUBLE, world_rank + 1, MY_ITER_TAG, MPI_COMM_WORLD, &status);
        // if(world_rank) MPI_Wait(&request, &status);
        // if(world_rank != P - 1) MPI_Wait(&request_2, &status);
        MPI_Barrier(MPI_COMM_WORLD);
        maxdiff = 0.0;
        for (int phase = 0; phase < 2; phase++) {
            for (int i = 1; i <= N/P - 1; i++) {
                for (int j = 1 + (even(i) ^ phase); j < N - 1; j += 2) {
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
        MPI_Barrier(MPI_COMM_WORLD);
         // if(world_rank == 1)printf("SUCCESS #1\n");
        MPI_Allreduce(&maxdiff, &glob_maxdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         // if(world_rank == 1)printf("SUCCESS #2\n");
        //MPI_Barrier(MPI_COMM_WORLD);
         // if(world_rank == 1)printf("SUCCESS #3\n");
    } while (glob_maxdiff > stopdiff);
    // printf("Running %d x %d SOR on cpu with rank %d out of world=%d\n", N, N, world_rank, world_size);
    
    if(!world_rank){
        // Master node
        double new_row[N];
        for(int i = 1; i < P; ++i){
            double new_row[N];
            for(int j = 1; j < N/P; ++j){
                int n = 0;
                
                // MPI_Probe(1, MY_FINAL_TAG, MPI_COMM_WORLD, &status);
                MPI_Recv((double*)new_row, N, MPI_DOUBLE, i, MY_FINAL_TAG, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &n);
                // 
                memcpy((void*)G[i * N/P + j - 1], (void*)new_row, sizeof(double) * N);    
            }
        }
    } else {
         for(int j = 1; j <= N/P; ++j){
            double a[N];
            for(int k = 0; k < N; ++k){
                a[k] = 1;
            }
            MPI_Send((double*)G[j], N, MPI_DOUBLE, 0, MY_FINAL_TAG, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    print_grid(G, N, world_rank);
    // if(!world_rank)
    // Finalize the MPI environment.
    MPI_Finalize();
    exit(0);
    return 0;
}
