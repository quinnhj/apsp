#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <omp.h>
#include <list>
#include <vector>
#include <unordered_set>

#define NUM_LOCKS 100

// Globals for convenience.
const float INF = 100000.0;
float min_edge = 1.0;
int last_bucket = 0;
int last_count = 1;
omp_lock_t* heap_locks;
omp_lock_t success_lock;
omp_lock_t delta_lock;
int n_threads = 1;

int* pcount;
int* qcount;
int* Lcount;
int* Rcount;
int* distcount;
int* popped_items;

// Counting variables
int size_sum = 0;
int size_count = 0;

/*
 * Reference Implementation of Floyd Warshall.
 * It is used for correctness checking, as well as sanity check on performance
 */
void floyd_warshall(int n, int* par, float* dist) {
    
    /*
    // Commenting out parent calculations
    // Uncomment if  you want parent matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j || dist[i*n + j] == INF) {
                par[i*n + j] = -1;
            } else {
                par[i*n + j] = i;
            }
        }
    }
    */

    float new_dist;
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                new_dist = dist[i*n + k] + dist[k*n + j];
                if (new_dist < dist[i*n + j]) {
                    dist[i*n + j] = new_dist;
                    
                    // Uncomment this if you want parent matrix
                    //par[i*n + j] = par[k*n + j];
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

    if (idx < last_bucket) {
        printf("INSERTING TOO EARLY, idx = %d\tlast_bucket= %d\n", idx, last_bucket);
    }

    // Add to bucket
    omp_set_lock(heap_locks + (idx % NUM_LOCKS));
    bucket_heap[idx].push_back(pair);
    std::list<int>::iterator iter = bucket_heap[idx].end();
    omp_unset_lock(heap_locks + (idx % NUM_LOCKS));
    
    // Update bookkeeping so that we can decrement this in constant time.
    b_num[pair] = idx; // array storing bucket numbers
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
    omp_set_lock(heap_locks + (b_num[pair] % NUM_LOCKS));
    bucket_heap[b_num[pair]].erase(iters[pair]);
    omp_unset_lock(heap_locks + (b_num[pair] % NUM_LOCKS));
    
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


bool incr_check (int* arr, bool* success, int val, int pos, bool write) {
    
    int newval = write ? val : val + 1;
    bool temp;

    omp_set_lock(&success_lock);
    if (arr[pos] >= newval || !*success) {
        temp = false;
    } else {
        arr[pos] = val;
        temp = true;
    }

    *success = temp;
    omp_unset_lock(&success_lock);

}


bool check_overlap (std::list<int>* bucket_heap, int* b_num, int n,
        std::list<int>::iterator* iters,
        float* dist, int* p, int* q, std::vector<int>* L, std::vector<int>* R,
        int pair, int last_count, bool* success, int* smallest_delta) {

    int u = pair/n;
    int v = pair%n;
    int lidx = u*n + q[u*n + v];
    int ridx = p[u*n + v]*n + v;
    omp_set_lock(&delta_lock);
    int my_smallest_delta = *smallest_delta;
    omp_unset_lock(&delta_lock);
    
    incr_check(distcount, success, last_count, pair, false);
    incr_check(pcount, success, last_count, pair, false);
    incr_check(qcount, success, last_count, pair, false);
    
    incr_check(Lcount, success, last_count, lidx, false);
    incr_check(Lcount, success, last_count, p[pair]*n + v, true);
    incr_check(Rcount, success, last_count, ridx, false);
    incr_check(Rcount, success, last_count, u*n + q[pair], true);
   
    int lsize = (int)L[lidx].size();
    int rsize = (int)R[ridx].size();
    float val;
    
    for (int j = 0; j < lsize; j++) {
        val = dist[pair] + dist[L[lidx][j]*n + u];
        my_smallest_delta = min((int)(val/min_edge), my_smallest_delta);
        incr_check(distcount, success, last_count, L[lidx][j]*n + u, false);
        incr_check(distcount, success, last_count, L[lidx][j]*n + v, true);
        incr_check(pcount, success, last_count, L[lidx][j]*n + v, true);
        incr_check(pcount, success, last_count, L[lidx][j]*n + u, false);
        incr_check(qcount, success, last_count, L[lidx][j]*n + v, true);
    }
   
    for (int j = 0; j < rsize; j++) {
        val = dist[pair] + dist[v*n +  R[ridx][j]];
        my_smallest_delta = min((int)(val/min_edge), my_smallest_delta);
        incr_check(distcount, success, last_count, v*n + R[ridx][j], false);
        incr_check(distcount, success, last_count, u*n + R[ridx][j], true);
        incr_check(pcount, success, last_count, u*n + R[ridx][j], true);
        incr_check(qcount, success, last_count, u*n + R[ridx][j], true);
        incr_check(qcount, success, last_count, v*n + R[ridx][j], false);
    }
    
    omp_set_lock(&delta_lock);
    *smallest_delta = min(my_smallest_delta, *smallest_delta);
    omp_unset_lock(&delta_lock);

} 



/*
 * Attempt to extract as many minimum values from the heap as possible.
 * It returns empty as failure, signifying an empty heap.
 * Otherwise, it returns a vector.
 */
std::vector<int> heap_extract_multiple(std::list<int>* bucket_heap, int* b_num, int n,
        std::list<int>::iterator* iters,
        float* dist, int* p, int* q, std::vector<int>* L, std::vector<int>* R) {

    // Setting up data structures to ensure we don't destroy stuff in parallel.
    std::vector<int> ret_val;
    int smallest_delta = n*n + 1;
    int parallel_smallest_delta = smallest_delta;
    int pair, u, v, set_start_size, num_added, lidx, ridx;
    int rollback_bucket = 0;
    bool success = true;
    bool parallel_success = true;
    last_count += 2;
    
    // Walk through each bucket in order. The initializing/incrementing
    // with last_bucket ensures that we only ever see each index once.
    for (int i = last_bucket; i < n*n + 1; i++, last_bucket++) {
        while (bucket_heap[i].size() > 0) {

            //If it's a normal bucket
            if (i < n*n) {
                int temp_i = i;
                int num_items = 0;
                std::list<int>::iterator listIter;
                
                //Populate popped_items with the next < n_threads items
                while (num_items < n_threads && temp_i < n*n) {
                    //printf("List size: %d\n", (int)bucket_heap[temp_i].size());
                    for(listIter = bucket_heap[temp_i].begin(); 
                            listIter != bucket_heap[temp_i].end() && num_items < n_threads;
                            listIter++, num_items++) {

                        //printf("Adding to popped items, num_items = %d\n", num_items);
                        popped_items[num_items] = *listIter;
                        //printf("Popped item: %d\n", popped_items[num_items]);

                    }
                    temp_i++; // one greater than last bucket
                }

                //Check everything in popped_items in parallel
                #pragma omp parallel for
                for (int a = 0; a < num_items; a++) {
                    //printf("Parallel Check, a = %d\tnum_items = %d\n", a, num_items);
                    check_overlap (bucket_heap, b_num, n, iters, dist, p, q, L, R,
                             popped_items[a], last_count, &parallel_success, &parallel_smallest_delta);
                }
               
                //printf("Finished Parallel Check\n");
                //If everything in this set of num_items was good
                if (parallel_success && temp_i < parallel_smallest_delta) {
                    //printf("Parallel Success!\n");
                    //Add num_items to retval
                    for (int a = 0; a < num_items; a++) {
                        ret_val.push_back(popped_items[a]);
                    }
                    
                    //Pop num_items from lists.
                    temp_i = i;
                    num_items = 0;
                    while (num_items < n_threads && temp_i < n*n) {
                       
                        //Because you have to set INITIAL size, or C++ does stupid shit
                        int size = bucket_heap[temp_i].size();
                        for (int a = 0; a < size && num_items < n_threads; a++) {
                            //printf("Removed item: %d\n", bucket_heap[temp_i].front());
                            bucket_heap[temp_i].pop_front();
                            num_items++;
                        }
                        temp_i++; // one greater than last bucket
                    }
                    
                    //Temp_i is one greater than the last bucket. Over estimate.
                    rollback_bucket = max(temp_i - 3, 0); // Last bucket to have 0 or more items.
                    i = rollback_bucket;
                    last_bucket = i;

                } else {
                    //printf("Parallel Failure...\n");

                    for (int a = i; a < n*n; a++, last_bucket++) {
                        //printf("Innermost loop, a = %d\n", a);
                        while (bucket_heap[a].size() > 0) {
                            if (a < n*n) {
                                /*
                                pair = bucket_heap[a].front();

                                if ((int)ret_val.size() == 0) {
                                    printf("Returning one because empty\n");
                                    ret_val.push_back(pair);
                                    bucket_heap[a].pop_front();
                                } else {
                                    printf("Returning what we had\n");
                                    last_bucket = rollback_bucket;
                                }
                                return ret_val;
                                */
                                    
                                //printf("Success before: %s\n", success ? "True" : "False");
                                pair = bucket_heap[a].front();
                                check_overlap (bucket_heap, b_num, n, iters, dist, p, q, L, R,
                                        pair, last_count, &success, &smallest_delta);
                                //printf("Success after : %s\n", success ? "True" : "False");

                                if (a > smallest_delta || !success) {
                                    
                                    //printf("Serial Failiure. SmallestDelta = %d\n", smallest_delta);
                                    if ((int)ret_val.size() == 0) {
                                        //printf("ret_val.size() == 0\n");
                                        ret_val.push_back(pair);
                                        bucket_heap[a].pop_front();
                                    } else {
                                        //printf("Setting back bucket\n");
                                        last_bucket = rollback_bucket;
                                    }

                                    //last_bucket = 0;
                                    return ret_val;
                                }
                               
                                //printf("Failed, returning: %d\n", pair);
                                ret_val.push_back(pair);
                                bucket_heap[a].pop_front();
                                rollback_bucket = a;
                            
                            }
                        }
                    }
                    
                    //last_bucket = 0;
                    return ret_val;
                }

            // Else it's the overflow bucket, which is handled using a fib heap.
            } else {
                // This should be implemented using a fib heap as per the paper.
                // For now enforcing that this scenario can't happen by limiting
                // the minimum edge value.
                return ret_val;
            }
        }
    }

    // Return -1 if the heap is empty
    return ret_val;
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
    int i = 0;

    int lsize = (int)L[lidx].size();
    //#pragma omp parallel for
    for (i = 0; i < lsize; i++) {

        di_examine(n, L[lidx][i], u, v, dist, graph,
            p, q, L, R, bucket_heap, b_num, iters);

    }
 
    int rsize = (int)R[ridx].size();
    //#pragma omp parallel for
    for (i = 0; i < rsize; i++) {
        
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
    int pair, set_size;
    std::vector<int> working_set;

    //Main loop
    while (1) {

        /*
        // Get min (u,v) pair from the heap.
        // Break the loop if heap was empty.
        pair = heap_extract(bucket_heap, b_num, n, iters);
        if (pair < 0) break;
        u = pair/n;
        v = pair%n;
        */

        // Get as many min (u,v) pairs as we can do in parallel from the heap.
        // Break the loop if heap was empty.
        working_set = heap_extract_multiple(bucket_heap, b_num, n, iters,
                            dist, p, q, L, R);
        set_size = (int) working_set.size();
        if (set_size == 0) {
            break;
        }
    
        size_sum += set_size;
        size_count++;
        
        // For everything in the set, do in parallel.
        // Slow as hell, but if the set returned is valid, it is correct.
        // Way too much overhead at the moment though :(.
        #pragma omp parallel for //schedule(dynamic)
        for (int i = 0; i < set_size; i++) {

            int u = working_set[i]/n;
            int v = working_set[i]%n;

            // Add to relevant lists
            L[p[u*n + v]*n + v].push_back(u);
            R[u*n + q[u*n + v]].push_back(v);

            // Examine through the lists.
            di_lists(n, u, v, dist, graph, 
                p, q, L, R, bucket_heap, b_num, iters);
        
        }
        //#pragma omp barrier

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
        printf( "-t <int> to set number of threads. Default 1\n");
        printf( "-no turns off all correctness checks\n");
        printf( "-csv outputs in csv format (n, floyd_time, di_time)\n");
        return 0;
    } 
    int n = read_int( argc, argv, "-n", 100 );
    n_threads = read_int( argc, argv, "-t", 1);
    omp_set_num_threads(n_threads);
    printf("Running with %d threads.\n", n_threads);
    
    heap_locks = (omp_lock_t*) malloc(NUM_LOCKS * sizeof(omp_lock_t));
    for (int i = 0; i < NUM_LOCKS; i++) {
        omp_init_lock(heap_locks + i);
    }
    omp_init_lock(&success_lock);
    omp_init_lock(&delta_lock);
    
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

    pcount = (int*) malloc(n*n*sizeof(int));
    qcount = (int*) malloc(n*n*sizeof(int));
    Lcount = (int*) malloc(n*n*sizeof(int));
    Rcount = (int*) malloc(n*n*sizeof(int));
    distcount = (int*) malloc(n*n*sizeof(int));
    popped_items = (int*) malloc(n_threads*sizeof(int));

 
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
        pcount[i] = 0;
        qcount[i] = 0;
        Lcount[i] = 0;
        Rcount[i] = 0;
        distcount[i] = 0;
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
        printf("Average size of parallel set: %f\ttotal:%d\n", 
                (float)size_sum/(float)size_count, size_sum);
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
    
    free(pcount);
    free(qcount);
    free(Lcount);
    free(Rcount);
    free(distcount);
    free(popped_items);

    delete[] L;
    delete[] R;
    delete[] bucket_heap;
    delete[] iters;

    for (int i = 0; i < NUM_LOCKS; i++) {
        omp_destroy_lock(heap_locks + i);
    } 
    free(heap_locks);

    return 0;
}
