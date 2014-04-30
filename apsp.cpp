#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/binomial_heap.hpp>

#define HEAP fibonacci_heap
#define BUCKET

using namespace boost::heap;
const float INF = 100000.0;
float min_edge = 1.0;
int last_bucket = 0;
bool empty = false;
//Note, reversed comparator, so call "increase" when decreasing key.
HEAP<vert_pair, compare<vert_comparator> > heap;

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

void heap_insert(std::list<int>* bucket_heap, int* b_num, int pair, int val, int n,
        std::list<int>::iterator* iters) {
    int idx = min((int)(val/min_edge), n*n);
    bucket_heap[idx].push_back(pair);
    b_num[pair] = idx;
    if (bucket_heap[idx].size() < 1) {
        printf("SIZE OF ZERO\n");
    }
    std::list<int>::iterator iter = bucket_heap[idx].end();
    iter--;
    iters[pair] = iter;
}


void heap_decrease(std::list<int>* bucket_heap, int* b_num, int pair, int val, int n,
        std::list<int>::iterator* iters) {
    int orig_idx = b_num[pair];
    std::list<int> list = bucket_heap[orig_idx];
    
    //Naive O(N) removal of item from linked list. Slow as hell.
    list.remove(pair);
    //What we should do, but it's not that simple :/
    //list.erase(iters[pair]);
    heap_insert(bucket_heap, b_num, pair, val, n, iters);
}


int heap_extract(std::list<int>* bucket_heap, int* b_num, int n,
        std::list<int>::iterator* iters) {
    int ret_val;
    for (int i = last_bucket; i < n*n + 1; i++,last_bucket++) {
        if (bucket_heap[i].size() > 0) {
            if (i < n*n) {
                ret_val = bucket_heap[i].front();
                bucket_heap[i].pop_front();
                return ret_val;
            } else {
                //This should be implemented using a real heap as per the paper.
                // For now enforcing that this scenario can't happen by limiting
                // the minimum edge value.
                return -1;
            }
        }
    }
    return -1;
}

void di_init(int n, float* dist, float* graph, vert_pair* q_arr, int* p, int* q,
             std::vector<int>* L, std::vector<int>* R,
             std::list<int>* bucket_heap, int* b_num,
             std::list<int>::iterator* iters,
             HEAP<vert_pair, compare<vert_comparator> >::handle_type* handles) {
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dist[i*n + j] = INF;
            p[i*n + j] = -1;
            q[i*n + j] = -1;
            //L and R are already initialized
        }
    }
    
    /* 
     * Uncomment this if you want distance to self to be 0.
     * Otherwise it finds the shortest path out and back to itself. 
     *
    for (int i = 0; i < n; i++) {
        dist[i*n + i] = 0.0;
    }
    */

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
                
                #ifdef BUCKET
                heap_insert(bucket_heap, b_num, i*n + j, graph[i*n + j], n, iters);
                #else
                handles[i*n + j] = heap.push(vert);
                #endif
            }
        }
    }
}

void di_examine(int n, int u, int v, int w, float* dist, float* graph, vert_pair* q_arr, 
                int* p, int* q, std::vector<int>* L, std::vector<int>* R,
                std::list<int>* bucket_heap, int* b_num,
                std::list<int>::iterator* iters,
                HEAP<vert_pair, compare<vert_comparator> >::handle_type* handles) {
    if (dist[u*n + v] + dist[v*n + w] < dist[u*n + w]) {
        dist[u*n + w] = dist[u*n + v] + dist[v*n + w];
        
        vert_pair vert = q_arr[u*n + w];
        vert.u = u;
        vert.v = w;
        vert.dist = dist[u*n + w];
            
        if (p[u*n + w] == -1) {    
            #ifdef BUCKET
            heap_insert(bucket_heap, b_num, u*n + w, dist[u*n + w], n, iters);
            #else
            handles[u*n + w] = heap.push(vert);
            #endif
        } else {
            #ifdef BUCKET
            heap_decrease(bucket_heap, b_num, u*n + w, dist[u*n + w], n, iters);
            #else
            heap.increase(handles[u*n + w], vert);
            #endif
        }

        p[u*n + w] = p[u*n + v];
        q[u*n + w] = q[v*n + w];
    }
}


void di_lists(int n, int u, int v, float* dist, float* graph, vert_pair* q_arr, 
                int* p, int* q, std::vector<int>* L, std::vector<int>* R,
                std::list<int>* bucket_heap, int* b_num,
                std::list<int>::iterator* iters,
                HEAP<vert_pair, compare<vert_comparator> >::handle_type* handles) {

    int lidx = u*n + q[u*n + v];
    int ridx = p[u*n + v]*n + v;
    std::vector<int>::size_type i = 0;

    for (i = 0; i != L[lidx].size(); i++) {

        di_examine(n, L[lidx][i], u, v, dist, graph, q_arr, 
            p, q, L, R, bucket_heap, b_num, iters, handles);

    }
  
    for (i = 0; i != R[ridx].size(); i++) {

        di_examine(n, u, v, R[ridx][i], dist, graph, q_arr, 
            p, q, L, R, bucket_heap, b_num, iters, handles);

    }
}

void di_apsp(int n, float* dist, float* graph, vert_pair* q_arr, 
             int* p, int* q, std::vector<int>* L, std::vector<int>* R,
             std::list<int>* bucket_heap, int* b_num,
             std::list<int>::iterator* iters,
             HEAP<vert_pair, compare<vert_comparator> >::handle_type* handles) {
    
    di_init(n, dist, graph, q_arr, p, q, L, R, bucket_heap, b_num, iters, handles);
    int u, v, pair;
    //Main loop
    #ifdef BUCKET
    while (1) {
    #else
    while (!heap.empty()) {
    #endif

        #ifdef BUCKET
        pair = heap_extract(bucket_heap, b_num, n, iters);
        if (pair < 0) break;
        u = pair/n;
        v = pair%n;
        #else
        vert_pair vert = heap.top();
        u = vert.u;
        v = vert.v;
        //Has to happen AFTER setting u and v, or there will be a segfault
        heap.pop();
        #endif
        
        L[p[u*n + v]*n + v].push_back(u);
        R[u*n + q[u*n + v]].push_back(v);

        di_lists(n, u, v, dist, graph, q_arr, 
            p, q, L, R, bucket_heap, b_num, iters, handles);
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
        printf( "-no turns off all correctness checks\n");
        printf( "-csv outputs in csv format (n, floyd_time, di_time)\n");
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
    srand48(1337);
    float min_allowable_value = (1.0/(n*n));
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
    std::vector<int> *L = new std::vector<int>[n*n];
    std::vector<int> *R = new std::vector<int>[n*n];
    HEAP<vert_pair, compare<vert_comparator> >::handle_type *handles = 
            new HEAP<vert_pair, compare<vert_comparator> >::handle_type[n*n];
    std::list<int> *bucket_heap = new std::list<int>[n*n + 1];
    std::list<int>::iterator *iters = new std::list<int>::iterator[n*n];
    int *b_num = (int*) malloc(n * n * sizeof(int));

    //Initialize arrays
    for (int i = 0; i < n*n; i++) {
        fw_dist[i] = graph[i];
        di_dist[i] = graph[i];
        b_num[i] = 0;
    }

    //
    // Here we do the fast alg
    //
    
    di_time = read_timer();
    di_apsp(n, di_dist, graph, q_arr, 
            p, q, L, R, bucket_heap, b_num, iters, handles);
    di_time = read_timer() - di_time;

    //Running Floyd Warshall for standard.
    floyd_time = read_timer();
    floyd_warshall(n, fw_par, fw_dist);
    floyd_time = read_timer() - floyd_time;

    //Validating output against Floyd Warshall
    if( find_option( argc, argv, "-no" ) == -1 ) {
        for (int i = 0; i < n*n; i++) {
            if ((di_dist[i] - fw_dist[i] > 0.0001) || (fw_dist[i] - di_dist[i] > 0.0001) ) {
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

    if( find_option( argc, argv, "-csv" ) == -1 ) {
        printf("\nTimes:\n");
        printf("Floyd Warshall: %f\n", floyd_time);
        printf("Demetrescu Italiano: %f\n", di_time);
    } else {
        printf("%d\t%f\t%f\n", n, floyd_time, di_time);
    }

    // Freeing memory
    free(graph);
    free(fw_dist);
    free(fw_par);
    free(di_dist);
    free(di_par);
    free(q_arr);
    free(p);
    free(q);
    free(b_num);
    
    return 0;
}
