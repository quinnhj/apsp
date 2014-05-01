#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
#include <list>
#include <vector>

// Globals for convenience.
const float INF = 100000.0;
float min_edge = 1.0;
int last_bucket = 0;

/*
 * Reference Implementation of Floyd Warshall.
 * It is used for correctness checking, as well as sanity check on performance
 */
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

/*
 * Inserts into our bucket based heap.
 * It doesn't maintain a true heap, but a coarsely grained heap with
 * buckets of size 1/n*n. 
 */
void heap_insert(std::list<int>* bucket_heap, int* b_num, int pair, float val, int n,
        std::list<int>::iterator* iters) {
    
    // Figure out which bucket to insert into
    int idx = (int)(val/min_edge);
    idx = min(idx, n*n);

    // Add to bucket
    bucket_heap[idx].push_back(pair);
    
    // Update bookkeeping so that we can decrement this in constant time.
    b_num[pair] = idx; // array storing bucket numbers
    std::list<int>::iterator iter = bucket_heap[idx].end();
    iter--;
    iters[pair] = iter; // pointer to position in said bucket
}

/*
 * Decrease the value of a key in our heap.
 * We simply delete it from its relevant bucket, then call
 * heap_insert to insert it with a new value.
 */
void heap_decrease(std::list<int>* bucket_heap, int* b_num, int pair, float val, int n,
        std::list<int>::iterator* iters) {
    
    // Deleting 'pair' from its bucket.
    bucket_heap[b_num[pair]].erase(iters[pair]);
    
    // Inserting it into heap with a new value.
    heap_insert(bucket_heap, b_num, pair, val, n, iters);
}

/*
 * Attempt to extract minimum value from the heap.
 * It returns -1 as failure, signifying an empty heap.
 */
int heap_extract(std::list<int>* bucket_heap, int* b_num, int n,
        std::list<int>::iterator* iters) {

    int ret_val;

    // Walk through each bucket in order. The initializing/incrementing
    // with last_bucket ensures that we only ever see each index once.
    for (int i = last_bucket; i < n*n + 1; i++,last_bucket++) {
        if (bucket_heap[i].size() > 0) {
            // If it's a normal bucket
            if (i < n*n) {
                
                // Remove an item from the bucket and return it.
                ret_val = bucket_heap[i].front();
                bucket_heap[i].pop_front();
                return ret_val;

            // Else it's the overflow bucket, which is handled using a fib heap.
            } else {
                // This should be implemented using a fib heap as per the paper.
                // For now enforcing that this scenario can't happen by limiting
                // the minimum edge value.
                return -1;
            }
        }
    }

    // Return -1 if the heap is empty
    return -1;
}

/*
 * Initializer function for the DI APSP algorithm.
 */
void di_init(int n, float* dist, float* graph, int* p, int* q,
             std::vector<int>* L, std::vector<int>* R,
             std::list<int>* bucket_heap, int* b_num,
             std::list<int>::iterator* iters) {
    
    // For every pair of vertices
    for (int i = 0; i < n*n; i++) {
        dist[i] = INF;

        // Because these are always set to positive ints, this is our null
        p[i] = -1;
        q[i] = -1;
        //L and R are already initialized
    }
    
    /* 
     * Uncomment this if you want distance to self to be 0.
     * Otherwise it finds the shortest path out and back to itself. 
     *
    for (int i = 0; i < n; i++) {
        dist[i*n + i] = 0.0;
    }
    */

    // For every edge that exists in the graph.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (graph[i*n + j] != INF) {
                
                dist[i*n + j] = graph[i*n + j];
                p[i*n + j] = j;
                q[i*n + j] = i;

                heap_insert(bucket_heap, b_num, i*n + j, graph[i*n + j], n, iters);
            }
        }
    }
}

/*
 * Function that examines for an opportunity to discover a new shortest path
 * from u to v, using w.
 */
void di_examine(int n, int u, int v, int w, float* dist, float* graph, 
                int* p, int* q, std::vector<int>* L, std::vector<int>* R,
                std::list<int>* bucket_heap, int* b_num,
                std::list<int>::iterator* iters) {
    
    if (dist[u*n + v] + dist[v*n + w] < dist[u*n + w]) {
        dist[u*n + w] = dist[u*n + v] + dist[v*n + w];
        
        // If it's never been inserted, insert, else decrease.
        if (p[u*n + w] == -1) {    
            heap_insert(bucket_heap, b_num, u*n + w, dist[u*n + w], n, iters);
        } else {
            heap_decrease(bucket_heap, b_num, u*n + w, dist[u*n + w], n, iters);
        }

        p[u*n + w] = p[u*n + v];
        q[u*n + w] = q[v*n + w];
    }
}

/*
 * Function that examines across the L and R lists.
 */

void di_lists(int n, int u, int v, float* dist, float* graph,
                int* p, int* q, std::vector<int>* L, std::vector<int>* R,
                std::list<int>* bucket_heap, int* b_num,
                std::list<int>::iterator* iters) {

    int lidx = u*n + q[u*n + v];
    int ridx = p[u*n + v]*n + v;
    std::vector<int>::size_type i = 0;

    for (i = 0; i != L[lidx].size(); i++) {

        di_examine(n, L[lidx][i], u, v, dist, graph,
            p, q, L, R, bucket_heap, b_num, iters);

    }
  
    for (i = 0; i != R[ridx].size(); i++) {

        di_examine(n, u, v, R[ridx][i], dist, graph,
            p, q, L, R, bucket_heap, b_num, iters);

    }
}

/*
 * Primary function for DI APSP.
 */
void di_apsp(int n, float* dist, float* graph,
             int* p, int* q, std::vector<int>* L, std::vector<int>* R,
             std::list<int>* bucket_heap, int* b_num,
             std::list<int>::iterator* iters) {
    
    // Initialize
    di_init(n, dist, graph, p, q, L, R, bucket_heap, b_num, iters);
    int u, v, pair;

    //Main loop
    while (1) {

        // Get min (u,v) pair from the heap.
        // Break the loop if heap was empty.
        pair = heap_extract(bucket_heap, b_num, n, iters);
        if (pair < 0) break;
        u = pair/n;
        v = pair%n;
       
        // Add to relevant lists
        L[p[u*n + v]*n + v].push_back(u);
        R[u*n + q[u*n + v]].push_back(v);

        // Examine through the lists.
        di_lists(n, u, v, dist, graph, 
            p, q, L, R, bucket_heap, b_num, iters);
    }
}


int main( int argc, char **argv )
{    
    //Declaring variables
    double floyd_time, di_time;

    /*
     *Handling user input.
     */
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of vertices\n" );
        printf( "-no turns off all correctness checks\n");
        printf( "-csv outputs in csv format (n, floyd_time, di_time)\n");
        return 0;
    } 
    int n = read_int( argc, argv, "-n", 100 );

    /*
     * Setting up the data structures
     * All data structures are assumed to be row major.
     */

    // An n*n matrix, where the element in the ith row and jth column
    // represents the weight of the edge from vertex i to j.
    float *graph = (float*) malloc( n * n * sizeof(float));
    
    // An n*n matrix, where the element in the ith row and jth column
    // represents the second vertex on the shortest path from u -> v
    int *p = (int*) malloc(n * n * sizeof(int));

    // An n*n matrix, where the element in the ith row and jth column
    // represents the penultimate vertex on the shortest path from u -> v
    int *q = (int*) malloc(n * n * sizeof(int));

    // An n*n matrix, where the element in the ith row and jth column
    // is a list of vertices w for which w -> u -Shortest-Path-> v
    // is known to be a shortest path
    std::vector<int> *L = new std::vector<int>[n*n];
    
    // An n*n matrix, where the element in the ith row and jth column
    // is a list of vertices w for which u -Shortest-Path-> v -> w
    // is known to be a shortest path
    std::vector<int> *R = new std::vector<int>[n*n];

    // An n*n+1 matrix, where each element holds a list of vertex pairs.
    // This is the core data structure for our bucket based heap.
    std::list<int> *bucket_heap = new std::list<int>[n*n + 1];

    // An n*n matrix, where the element in the ith row and jth column
    // is the iterator/pointer to the linked list element holding
    // the vertex pair (i,j) in bucket_heap.
    std::list<int>::iterator *iters = new std::list<int>::iterator[n*n];
    
    // An n*n matrix, where the element in the ith row and jth column
    // is the index of the bucket that holds the vertex pair (i,j) in bucket_heap.
    int *b_num = (int*) malloc(n * n * sizeof(int));
    
    // Parent and dist arrays for floyd warshall
    float *fw_dist = (float*) malloc( n * n * sizeof(float));
    int *fw_par = (int*) malloc( n * n * sizeof(int));

    // Dist array for DI
    float *di_dist = (float*) malloc( n * n * sizeof(float));
 
    /*
     * Initializing values.
     */

    srand48(1337);
    // We set a minimum value to avoid the overflow bucket scenario.
    float min_allowable_value = (1.0/(n*n));
    
    // Initialize every edge to have a weight between min_allowable_value  and 1.
    for (int i = 0; i < n * n; i++) {
        float temp;
        //no edge between an edge and itself.
        if (i%n == i/n) {
            graph[i] = INF;
        } else {
            temp = (float)drand48();
            if (temp < min_allowable_value) temp += min_allowable_value;
            if (temp < min_edge) min_edge = temp;
            graph[i] = temp;
        }
    }
   
    //Initialize dist arrays.
    for (int i = 0; i < n*n; i++) {
        fw_dist[i] = graph[i];
        di_dist[i] = graph[i];
        b_num[i] = 0;
    }

    /*
     * Executing and timing algorithms
     */
    
    // Running Demetrescu Italiano.
    di_time = read_timer();
    di_apsp(n, di_dist, graph, p, q, L, R, bucket_heap, b_num, iters);
    di_time = read_timer() - di_time;

    //Running Floyd Warshall as a standard implementation.
    floyd_time = read_timer();
    floyd_warshall(n, fw_par, fw_dist);
    floyd_time = read_timer() - floyd_time;

    /*
     * Validating DI results against reference Floyd Warshall
     */

    if( find_option( argc, argv, "-no" ) == -1 ) {
        
        //Report a failure if anything varies by more than 0.0001
        for (int i = 0; i < n*n; i++) {
            if ((di_dist[i] - fw_dist[i] > 0.0001) || (fw_dist[i] - di_dist[i] > 0.0001) ) {
                printf("\n\nFAILURE at %d\n\n", i);
                break;
            }
        }

        //Print out matrices if n is small enough.
        if (n < 10) {

            printf("\n---Our Graph---\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%.4f\t", graph[i*n + j]);
                }
                printf("\n");
            }
            printf("\n");
     
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

            printf("\n---Demetrescu Italiano Dist---\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%.4f\t", di_dist[i*n + j]);
                }
                printf("\n");
            }
            printf("\n");
        }

    }

    /*
     * Print out results
     */

    if( find_option( argc, argv, "-csv" ) == -1 ) {
        printf("\nTimes:\n");
        printf("Floyd Warshall: %f\n", floyd_time);
        printf("Demetrescu Italiano: %f\n", di_time);
    } else {
        printf("%d\t%f\t%f\n", n, floyd_time, di_time);
    }

    /*
     * Free up manually allocated memory.
     */

    free(graph);
    free(fw_dist);
    free(fw_par);
    free(di_dist);
    free(p);
    free(q);
    free(b_num);
    
    return 0;
}
