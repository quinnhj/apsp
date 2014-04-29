#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
#include <boost/heap/fibonacci_heap.hpp>

using namespace boost::heap;
const float INF = 100000.0;
//Note, reversed comparator, so call "increase" when decreasing key.
fibonacci_heap<vert_pair, compare<vert_comparator> > fib_heap;

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

void di_init(int n, float* dist, float* graph, vert_pair* q_arr, int* p, int* q,
             std::list<int>* L, std::list<int>* R,
             fibonacci_heap<vert_pair, compare<vert_comparator> >::handle_type* handles) {
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dist[i*n + j] = INF;
            p[i*n + j] = -1;
            q[i*n + j] = -1;
            //L and R are already initialized
        }
    }
    
    for (int i = 0; i < n; i++) {
        dist[i*n + i] = 0.0;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (graph[i*n + j] != INF) {
                
                dist[i*n + j] = graph[i*n + j];
                p[i*n + j] = j;
                q[i*n + j] = i;

                vert_pair vert = q_arr[i*n + j];
                vert.u = i;
                vert.v = j;
                vert.dist = graph[i*n + j];
                handles[i*n + j] = fib_heap.push(vert);
            }
        }
    }
}

void di_examine(int n, int u, int v, int w, float* dist, float* graph, vert_pair* q_arr, 
                int* p, int* q, std::list<int>* L, std::list<int>* R,
                fibonacci_heap<vert_pair, compare<vert_comparator> >::handle_type* handles) {
    if (dist[u*n + v] + dist[v*n + w] < dist[u*n + w]) {
        dist[u*n + w] = dist[u*n + v] + dist[v*n + w];
        
        vert_pair vert = q_arr[u*n + w];
        vert.u = u;
        vert.v = w;
        vert.dist = dist[u*n + w];
            
        if (p[u*n + w] == -1) {    
            handles[u*n + w] = fib_heap.push(vert);
        } else {
            fib_heap.increase(handles[u*n + w], vert);
        }

        p[u*n + w] = p[u*n + v];
        q[u*n + w] = q[v*n + w];
    }
}

void di_apsp(int n, float* dist, float* graph, vert_pair* q_arr, 
             int* p, int* q, std::list<int>* L, std::list<int>* R,
             fibonacci_heap<vert_pair, compare<vert_comparator> >::handle_type* handles) {
    
    di_init(n, dist, graph, q_arr, p, q, L, R, handles);
    int u, v;
    //Main loop
    while (!fib_heap.empty()) {
        vert_pair vert = fib_heap.top();
        u = vert.u;
        v = vert.v;
        //Has to happen AFTER setting u and v, or there will be a segfault
        fib_heap.pop();
        L[p[u*n + v]*n + v].push_back(u);
        R[u*n + q[u*n + v]].push_back(v);
        for (std::list<int>::iterator li = L[u*n + q[u*n + v]].begin();
                li != L[u*n + q[u*n + v]].end(); li++) {
            
            di_examine(n, *li, u, v, dist, graph, q_arr, 
                p, q, L, R, handles);
            
        }

        for (std::list<int>::iterator li = R[p[u*n + v]*n + v].begin();
                li != R[p[u*n + v]*n + v].end(); li++) {
            
            di_examine(n, u, v, *li, dist, graph, q_arr, 
                p, q, L, R, handles);
            
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
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 100 );

    //
    // Setting up the data structures
    //

    // An n*n matrix, where the element in the ith row and jth column
    // represents the weight of the edge from vertex i to j.
    float *graph = (float*) malloc( n * n * sizeof(float));
    
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
    
    // Parent and dist arrays for floyd warshall
    float *fw_dist = (float*) malloc( n * n * sizeof(float));
    int *fw_par = (int*) malloc( n * n * sizeof(int));

    // Parent and dist arrays for DI
    float *di_dist = (float*) malloc( n * n * sizeof(float));
    int *di_par = (int*) malloc( n * n * sizeof(int));
   
    // Other arrays for DI
    vert_pair *q_arr = (vert_pair*) malloc (n * n * sizeof(vert_pair));
    int *p = (int*) malloc(n * n * sizeof(int));
    int *q = (int*) malloc(n * n * sizeof(int));
    std::list<int> *L = new std::list<int>[n*n];
    std::list<int> *R = new std::list<int>[n*n];
    fibonacci_heap<vert_pair, compare<vert_comparator> >::handle_type *handles = 
            new fibonacci_heap<vert_pair, compare<vert_comparator> >::handle_type[n*n];

    //Initialize arrays
    for (int i = 0; i < n*n; i++) {
        fw_dist[i] = graph[i];
        di_dist[i] = graph[i];
    }

    //
    // Here we do the fast alg
    //
    
    di_time = read_timer();
    di_apsp(n, di_dist, graph, q_arr, 
            p, q, L, R, handles);
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

    // Freeing memory
    free(graph);
    free(fw_dist);
    free(fw_par);
    free(di_dist);
    free(di_par);
    free(q_arr);
    free(p);
    free(q);
    
    //Not necessary with new?
    //free(L);
    //free(R);
    
    return 0;
}
