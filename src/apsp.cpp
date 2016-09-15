#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <limits>

#ifdef Prof
  #include <google/profiler.h>
#endif

using namespace std;

typedef vector< vector<double> > vec2D;
#define MAX_DOUBLE numeric_limits<double>::max()/10

struct Graph
{
    Graph() {}
    void init( int n ){
        W.resize( n, vector<double>(n, MAX_DOUBLE) );
        reset_matx();
    }
    void reset_matx(){
        size_t n=W.size();
        P.clear(); L.clear();
        P.resize( n, vector< vector<double> >( n, vector<double>(n, MAX_DOUBLE) ) );
        L.resize( n, vector< vector<double> >( n, vector<double>(n, MAX_DOUBLE) ) );
        for( size_t s=1; s<n; ++s ){ for( size_t i=1; i<n; ++i ) { L[s][i][i] = 0; } }
    }
    void init_L_P_1_matx();

    double &w( const int &i, const int &o ) { return W[i][o]; }
    double &p( const int &n, const int &i, const int &o ) { return P[n][i][o]; }
    double &l( const int &n, const int &i, const int &o ) { return L[n][i][o]; }

    vec2D W;        // adjacent matrix
    vector< vec2D > P;   // predecessor matrix
    vector< vec2D > L;   // shortest path cost matrix
};

bool read_file( const string& inFileName, Graph &g );
void asap_slow( Graph &g );
void asap_fast( Graph &g );
void ext_sp_slow(int m , Graph &g);
void ext_sp_fast(int m , Graph &g);
void show_stat( int m, Graph &g );
void show_matx(const vec2D &v );
bool is_inf( const double& n ) { return (n >= MAX_DOUBLE); }

int main(int argc, char **argv)
{
#ifdef Prof
    ProfilerStart( "my.prof" );
#endif
    //if( argc != 2 ) { cout << "./asap <in>" << endl; return 1; }

    Graph g;
    //if( read_file( argv[1], g )== false ) { return 1; }
    if( read_file( "../input/gn7.dot", g )== false ) { return 1; }

    asap_slow( g );

    asap_fast( g );

    cout << "\n\n\n\n\n\n\n\n\n" << endl;
#ifdef Prof
    ProfilerStop();
#endif
    return 0;
}

bool read_file( const string& inFileName, Graph &g )
{
    ifstream in( inFileName.c_str() );
    if( !in.is_open() ) { cout << "cannot open file: " << inFileName << endl; return false; }

    int numNodes = atoi( inFileName.substr( inFileName.rfind("n")+1 ).c_str() ) + 1;

    g.init( numNodes );
    string tmpStr;
    bool isStartEdge = false;
    while( in >> tmpStr ){
        if( tmpStr == "{" ) { isStartEdge = true; continue; }
        if( tmpStr == "}" ) { break; }
        if( isStartEdge ){ // process one edge
            int source  = atoi( tmpStr.substr(1).c_str() );
            in >> tmpStr >> tmpStr;
            int sink    = atoi( tmpStr.substr(1).c_str() );
            in >> tmpStr >> tmpStr >> tmpStr;
            int edgeWeight = atoi( tmpStr.substr( 1, tmpStr.rfind( "\"" )-1 ).c_str() );
            //cout << "v" << source << " -> v" << sink << " : " << edgeWeight << endl;
            g.w( source, sink ) = edgeWeight;
        }
    }

    //show_matx( g.W );

    in.close();
    return true;
}

void asap_slow( Graph &g )
{
    size_t n = g.W.size()-1;

    g.reset_matx();
    g.init_L_P_1_matx();

    for( size_t m=2; m<n; ++m ){
        ext_sp_slow( m, g );
    }
    show_stat( n-2, g );
}

void asap_fast( Graph &g )
{
    size_t n = g.W.size()-1;

    g.reset_matx();
    g.init_L_P_1_matx();

    size_t m=2;

    while( m < n-1 ){
        ext_sp_fast( m, g );
        m *= 2;
    }
    show_stat( m/2, g );
}

void ext_sp_slow( int m, Graph &g )
{
    size_t n = g.W.size()-1;

    for( size_t i=1; i<=n; ++i ){
        for( size_t j=1; j<=n; ++j ){
            g.p(m,i,j) = g.p(m-1,i,j);
            for( size_t k=1; k<=n; ++k ){
                double relaxCost = g.l(m-1,i,k) + g.w(k,j);
                if( relaxCost < g.l(m,i,j) ){
                    g.l(m,i,j) = relaxCost;
                    g.p(m,i,j) = k;
                }
            }
        }
    }
    //show_stat( m, g );
}

void ext_sp_fast(int m , Graph &g)
{
    cout << "ext fast " << m << endl;
    size_t n = g.W.size()-1;

    for( size_t i=1; i<=n; ++i ){
        for( size_t j=1; j<=n; ++j ){
            g.p(m,i,j) = g.p(m/2,i,j);
            for( size_t k=1; k<=n; ++k ){
                if( j==k ) { continue; }
                double relaxCost = g.l(m/2,i,k) + g.l(m/2,k,j);
                if( relaxCost < g.l(m,i,j) ){
                    g.l(m,i,j) = relaxCost;
                    g.p(m,i,j) = g.p(m/2,k,j);
                }
            }
        }
    }
    //show_stat( m, g );
}

void show_stat( int m, Graph &g )
{
    cout << "**********************" << endl;
    cout << "--- L(" << m << ") ---" << endl;
    show_matx(g.L[m]);
    cout << "--- P(" << m << ") ---" << endl;
    show_matx(g.P[m]);
    cout << "**********************" << endl;
}

void show_matx( const vec2D &v )
{
    int n = v.size();
    for( int i = 1; i<n; ++i ){
        for( int j = 1; j<n; ++j ){
            if( is_inf( v[i][j] ) ) { cout << "."; }
            else                    { cout << v[i][j]; }
            cout << "\t";
        }
        cout << endl;
    }
    cout << endl;
}


void Graph::init_L_P_1_matx()
{
    size_t n = W.size()-1;
    for( size_t i=1; i<=n; ++i ){
        for( size_t j=1; j<=n; ++j ){
            if( is_inf( w(i,j) ) == false ){
                l(1,i,j) = w(i,j);
                p(1,i,j) = i;
            }
        }
    }
}
