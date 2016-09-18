#ifndef APSP_H
#define APSP_H

#include "retimeDB.h"
#define MAX_APSP 10000.0

struct Cost
{
    Cost( int w_=MAX_APSP, double d_=-MAX_APSP ) { init(w_, d_); }
    void init( int w_, double d_ ){ w=w_; d=d_; }
    Cost operator + ( const Cost & c_ ){ return Cost( w + c_.w, d + c_.d ); }
    bool operator < ( const Cost & c_ ){
        bool isLess1 = w < c_.w;
        bool isLess2 = ( w==c_.w && w!=MAX_APSP && d<c_.d );
        return isLess1 || isLess2;
    }
    Cost &operator = ( const Cost & c_ ){
        if( w==c_.w && d==c_.d ) { return *this; }
        w = c_.w; d = c_.d; return *this;
    }
    friend ostream &operator << ( ostream&, const Cost& );
    int w;
    double d;
};
inline ostream& operator << (ostream& out, const Cost &c){
    out << c.w << "," << -c.d;
    return out;
}

typedef vector< vector<double> > vec2D;
typedef vector< vector<Cost> > vec2D_C; 

#define ForEachEntry( i, j, n ) \
    for( i = 0; i<n; ++i ) for( j=0; j<n; ++j )

struct ApspGraph
{
    ApspGraph( int n, retime::RetimeDB &rt, bool isRetimed ) {
        c.init( MAX_APSP, -MAX_APSP );
        W.resize( n, vector< Cost >(n, c) );
        reset_matx();
        delays.resize(n, -1);
        
        for( size_t i=0; i<n; ++i ){
            retime::Vertex &u = rt.get_v( i );
            for ( size_t j = 0; j<u.outs.size(); ++j ){
                retime::Edge &e = *u.outs[j];
                retime::Vertex &v = *e.t;
                double ew = (isRetimed) ? e.wr : e.w;
                w( u.id, v.id ) = Cost( ew, -u.d );
            }
            delays[ u.id ] = -u.d;
        }
    }
    void reset_matx(){
        size_t n=W.size();
        P.clear(); L.clear();
        P.resize( n, vec2D( n, vector<double>(n, MAX_APSP) ) );
        L.resize( n, vec2D_C( n, vector<Cost>(n, c) ) );
        Cost cost0( 0, 0 );
        for( size_t s=0; s<n; ++s ){ for( size_t i=0; i<n; ++i ) { L[s][i][i] = cost0; } }
    }
    void apsp(){
        size_t n = W.size();
        reset_matx();
        init_L_P_1_matx();
        for( size_t m=2; m<n; ++m ) { ext_sp( m ); }
        size_t i=0,j=0;
        ForEachEntry( i, j, n ){
            l(n-1,i,j) = l(n-1,i,j) + Cost( 0, delays[j] );
        }
    }
    void ext_sp( size_t m ){
        size_t n=W.size(), i=0, j=0;

        ForEachEntry( i, j, n ){
            p(m,i,j) = p(m-1,i,j);
            for( size_t k=0; k<n; ++k ){
                Cost relaxCost = l(m-1,i,k) + w(k,j);
                if( relaxCost < l(m,i,j) ){
                    l(m,i,j) = relaxCost;
                    p(m,i,j) = k;
                }
            }
        }
    }
    void show_matx_pure( const vec2D &v, const string &name )
    {
        int n = v.size();
        cout << " == " << name << " ==" << endl;
        for( int i = 0; i<n; ++i ){
            for( int j = 0; j<n; ++j ){
                if( is_inf( v[i][j] ) )    { cout << "."; }
                else                       { cout << v[i][j]; }
                cout << "\t"; }
            cout << endl; }
        cout << endl;
    }
    void init_L_P_1_matx()
    {
        size_t n = W.size(), i=0, j=0;
        ForEachEntry( i, j, n ){
            if( is_inf( w(i,j) )==false ){
                l(1,i,j) = w(i,j);
                p(1,i,j) = i;
            }
        }
    }
    bool is_inf( const Cost& c ) { return (c.w >= MAX_APSP ); }

    Cost &w( const int &i, const int &o )                   { return W[i][o]; }
    double &p( const int &n, const int &i, const int &o )   { return P[n][i][o]; }
    Cost &l( const int &n, const int &i, const int &o )     { return L[n][i][o]; }
    vector<double>      delays;
    vec2D_C             W;    // adjacent matrix
    vector< vec2D >     P;    // predecessor matrix
    vector< vec2D_C >   L;    // shortest path cost matrix
    Cost c;
};

#endif
