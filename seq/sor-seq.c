/*
 * Successive over relaxation
 * (red-black SOR)
 *
 * Written by Rob van Nieuwpoort, 6-oct-2003
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#define TOLERANCE 0.00002       /* termination criterion */
#define MAGIC 0.8               /* magic factor */

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
    double    **G = malloc(N * sizeof *G);
    if (G == NULL) {
        fprintf(stderr, "malloc failed\n");
        exit(42);
    }

    for (int i = 0; i < N; i++) {   /* malloc the own range plus one more line */
        /* of overlap on each side */
        G[i] = malloc(N * sizeof *G[i]);
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
                G[i][j] = 0.0;
        }
    }
}


void
print_grid(double **G, int N)
{
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            printf("%10.3f ", G[i][j]);
        }
        printf("\n");
    }
}


int
main(int argc, char *argv[])
{
    int         N;              /* problem size */
    int         ncol, nrow;     /* number of rows and columns */
    double      Gnew;
    double      r;
    double      omega;
    /* differences btw grid points in iters */
    double      stopdiff;
    double      maxdiff;
    double      diff;
    double    **G;              /* the grid */
    int         iteration; /* counters */
    struct timeval start;
    struct timeval end;
    double      time;
    int         print = 0;

    /* set up problem size */
    N = 10;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-print") == 0) {
            print = 1;
        } else {
            if (sscanf(argv[i], "%d", &N) != 1) {
                fprintf(stderr, "Positional parameter N must be an int, not '%s'\n",
                        argv[i]);
                exit(42);
            }
        }
    }

    fprintf(stderr, "Running %d x %d SOR\n", N, N);

    N += 2;                     /* add the two border lines */
    /* finally N*N (from argv) array points will be computed */

    /* set up a quadratic grid */
    ncol = nrow = N;
    r = 0.5 * (cos(M_PI / ncol) + cos(M_PI / nrow));
    omega = 2.0 / (1 + sqrt(1 - r * r));
    stopdiff = TOLERANCE / (2.0 - omega);
    omega *= MAGIC;

    alloc_grid(&G, N);
    init_grid(G, N);

    if (gettimeofday(&start, 0) != 0) {
        fprintf(stderr, "could not do timing\n");
        exit(1);
    }

    /* now do the "real" computation */
    iteration = 0;
    do {
        maxdiff = 0.0;
        for (int phase = 0; phase < 2; phase++) {
            for (int i = 1; i < N - 1; i++) {
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
        printf(" ITERATION %d\n", iteration);
    } while (maxdiff > stopdiff);

    if (gettimeofday(&end, 0) != 0) {
        fprintf(stderr, "could not do timing\n");
        exit(1);
    }

    time = (end.tv_sec + (end.tv_usec / 1000000.0)) -
        (start.tv_sec + (start.tv_usec / 1000000.0));

    fprintf(stderr, "SOR took %10.3f seconds\n", time);

    printf("Used %5d iterations, diff is %10.6f, allowed diff is %10.6f\n",
           iteration, maxdiff, stopdiff);

    if (print == 1) {
        print_grid(G, N);
    }

    return 0;
}
