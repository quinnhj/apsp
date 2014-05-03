#ifndef __APSP_COMMON_H__
#define __APSP_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }
inline float min( float a, float b ) { return a < b ? a : b; }
inline float max( float a, float b ) { return a > b ? a : b; }

typedef struct {
    int u;
    int v;
    float dist;
} vert_pair;

struct vert_comparator {
    bool operator() (const vert_pair & l, const vert_pair & r) const {
        return l.dist > r.dist;
    }
};

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;
const int PER_BIN = 10;

#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timing routines
//
double read_timer( );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
