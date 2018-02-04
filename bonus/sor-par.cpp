#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>

#define TOLERANCE 0.00002       /* termination criterion */
#define MAGIC 0.8               /* magic factor */
#define MY_INIT_TAG 0xbad00
#define MY_FINAL_TAG 0xbad01
#define MY_ITER_TAG 0xbad02
#define MY_PRINT_TAG 0x003

#define P_THREADS 2

int world_rank;
int P;
double glob_maxdiff;

pthread_barrier_t thread_fence;

struct thread_args{
  unsigned int N;
  short thread_id;
  double **G;
  unsigned int locN;
  unsigned int thread_quant;
};

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



void __sync_neigh(double **G, unsigned int locN, unsigned int N, unsigned thread_id){
  if(!thread_id)
    return; // It seems like a cheating, but it's better to only 0's thread take care of MPI things :)
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

static inline void barrier_wait(unsigned thread_id){
  while(pthread_barrier_wait(&thread_fence) == EAGAIN){}; // wait
  if(thread_id == 0)
    pthread_barrier_init(&thread_fence, NULL, P_THREADS); // reinitialize barrier
}

void* __solve(void* arg){
    /*
      G** double;
      lockN unsigned
      N unsigned
    */
    struct thread_args& args = *(struct thread_args*)arg;
    double **G = args.G;
    double Gnew;
    unsigned int locN = args.locN, N = args.N;
    unsigned int thread_quant = args.thread_quant,  thread_id = args.thread_id;
    int         ncol, nrow;     /* number of rows and columns */
    double      r;
    double      omega;
    double      stopdiff;
    double      maxdiff;
    int         iteration; /* counters */
    double      diff;
    ncol = nrow = N;
    r = 0.5 * (cos(M_PI / ncol) + cos(M_PI / nrow));
    omega = 2.0 / (1 + sqrt(1 - r * r));
    stopdiff = TOLERANCE / (2.0 - omega);
    omega *= MAGIC;

    iteration = 0;
    int border = locN - 1;
    int even_border = N/P * (P + 1) + N/P/P_THREADS * (P/P_THREADS + 1) + ((world_rank != 0 && thread_id != 0)? 1 : 0); // It's not magic but quite tough
    int big_border = N - 1;
    // do {
        maxdiff = 0.0;
        for (int phase = 0; phase < 2; phase++) {
            __sync_neigh(G, locN, N, thread_id);
            unsigned start_i = (thread_quant * thread_id);
            for (int i = start_i + 1 ; i < (start_i + thread_quant - 3) ; i++) {
                for (int j = 1 + (even(i + (even_border)) ^ phase); j < big_border; j += 2) {
                    // if(thread_id == 0) std::cout<<"J __SOLVE\n";
                    Gnew = stencil(G, i, j);
                    diff = fabs(Gnew - G[i][j]);
                    if (diff > maxdiff) {
                        maxdiff = diff;
                    }
                    G[i][j] = G[i][j] + omega * (Gnew - G[i][j]);

                }
            }
            barrier_wait(thread_id);
        }
        iteration++;
        if(thread_id == 0 ){
          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Allreduce(&maxdiff, &glob_maxdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
        // barrier_wait(thread_id);

        if(on_master() && thread_id == 0)
          printf(" ITERATION %d ::: DIFS (%f, %f)\n", iteration, stopdiff, glob_maxdiff);
    // } while (glob_maxdiff > stopdiff);

    if(on_master() && thread_id == 0){
      printf("Used %5d iterations, diff is %10.6f, allowed diff is %10.6f\n",
             iteration, maxdiff, stopdiff);

    }
    return nullptr;
}

double solve(unsigned N){


  double      Gnew;

  unsigned int locN, thread_quant;
  double **G;
  double start, end;
  double time;
  /* differences btw grid points in iters */
  pthread_t workers[P_THREADS];
  struct thread_args worker_args[P_THREADS];
  pthread_barrier_init(&thread_fence, NULL, P_THREADS);

  /* set up a quadratic grid */
  thread_quant = N / P / P_THREADS + 2;
  locN = N / P + 2; // Need two border lines;
  N += 2;
  if(on_master()){
    start = MPI_Wtime();
    alloc_grid(&G, N, N);
  } else {
    alloc_grid(&G, locN, N);
  }
  init_grid(G, locN, N);

  start = MPI_Wtime();
  /* now do the "real" computation */

  for(int thread_i = 0 ; thread_i < P_THREADS; ++thread_i){
    worker_args[thread_i].N = N;
    worker_args[thread_i].thread_id = thread_i;
    worker_args[thread_i].locN = locN;
    worker_args[thread_i].G = G;
    worker_args[thread_i].thread_quant = thread_quant;
    pthread_create(&workers[thread_i], NULL, &__solve, &worker_args[thread_i]);
  }

  for(int thread_i = 0 ; thread_i < P_THREADS; ++thread_i){
    void **res;
    pthread_join(workers[thread_i], res);
  }

  if(on_master()){
    end = MPI_Wtime();
    return end - start;
  } else
    return -1;

}

int
main(int argc, char *argv[])
{

    int         N;              /* problem size */
    double      time;
    int         print = 0;
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

    time = solve(N);
    if(on_master()){
      // if (print == 1) {
      //   print_grid(G, N, N);
      // }
      std::cout<<" solved ::: " << time <<std::endl;
    }

    MPI_Finalize();
    return 0;
}
