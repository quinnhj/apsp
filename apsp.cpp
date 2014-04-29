#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

const float INF = 100000.0;

void floyd_warshall(int n, int* par, float* dist) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j || dist[i*n + j] == INF) {
                par[i*n + j] = -1;
            } else {
                par[i*n + j] = i;
            }
        }
    }

    float new_dist;
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                new_dist = dist[i*n + k] + dist[k*n + j];
                if (new_dist < dist[i*n + j]) {
                    dist[i*n + j] = new_dist;
                    par[i*n + j] = par[k*n + j];
                }
            }
        }
    }

}


int main( int argc, char **argv )
{    
    //Declaring variables
    double floyd_time, di_time;

    //
    // Handling user input.
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of vertices\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 100 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    /*
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;
    */

    //
    // Setting up the data structures
    //

    // An n*n matrix, where the element in the ith row and jth column
    // represents the weight of the edge from vertex i to j.
    float *graph = (float*) malloc( n * n * sizeof(float));
    
    // Parent and dist arrays for floyd warshall
    float *fw_dist = (float*) malloc( n * n * sizeof(float));
    int *fw_par = (int*) malloc( n * n * sizeof(int));

    // Parent and dist arrays for DI
    float *di_dist = (float*) malloc( n * n * sizeof(float));
    int *di_par = (int*) malloc( n * n * sizeof(int));
    
    // Initialize every edge to have a weight between 0 and 1.
    srand48(5);
    for (int i = 0; i < n * n; i++) {
        //no edge between an edge and itself.
        if (i%n == i/n) {
            graph[i] = INF;
        } else {
            graph[i] = (float)drand48();
        }
    }
    //Print for viewing if small enough
    if (n < 10) {
        printf("\n---Our Graph---\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.4f\t", graph[i*n + j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    //Initialize arrays
    for (int i = 0; i < n*n; i++) {
        fw_dist[i] = graph[i];
        di_dist[i] = graph[i];
        //fw_par[i] = -1;
    }

    //
    // Here we do the fast alg
    //
    di_time = read_timer();
    //Should be call to actual alg once implemented.
    floyd_warshall(n, di_par, di_dist);
    di_time = read_timer() - di_time;


    //Running Floyd Warshall for standard.
    floyd_time = read_timer();
    floyd_warshall(n, fw_par, fw_dist);
    floyd_time = read_timer() - floyd_time;


    //Validating output against Floyd Warshall
    if( find_option( argc, argv, "-no" ) == -1 ) {
        for (int i = 0; i < n*n; i++) {
            if (abs(di_dist[i] - fw_dist[i]) > 0.0001) {
                printf("\n\nFAILURE at %d\n\n", i);
                break;
            }
        }

        //Printing Floyd Warshall Results if Short Enough
        if (n < 10) {
            printf("\n---Floyd Warshall Dist---\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%.4f\t", fw_dist[i*n + j]);
                }
                printf("\n");
            }
            printf("\n");
            printf("\n---Floyd Warshall Par---\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%d\t", fw_par[i*n + j]);
                }
                printf("\n");
            }
            printf("\n");
        }

        //Printing DI Results if Short Enough
        if (n < 10) {
            printf("\n---Demetrescu Italiano Dist---\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%.4f\t", di_dist[i*n + j]);
                }
                printf("\n");
            }
            printf("\n");
            printf("\n---Demetrescu Italiano Par---\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%d\t", di_par[i*n + j]);
                }
                printf("\n");
            }
            printf("\n");
        }

    }

    //Printing out times
    printf("\nTimes:\n");
    printf("Floyd Warshall: %f\n", floyd_time);
    printf("Demetrescu Italiano: %f\n", di_time);

    //
    // Printing summary data
    //
    /*
    if( fsum) 
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);
    */

    //
    // Clearing space
    //
    /*
    if( fsum )
        fclose( fsum );    
    if( fsave )
        fclose( fsave );
    */

   // Freeing memory
   free(graph);
   free(fw_dist);
   free(fw_par);
   free(di_dist);
   free(di_par);
    
    return 0;
}
