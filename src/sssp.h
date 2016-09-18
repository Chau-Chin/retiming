#ifndef _SSSP_H
#define _SSSP_H

#include "retimeDB.h"
#define MAX_SSSP 1000000.0

using namespace retime;
struct ConstGraph
{
    ConstGraph( RetimeDB &rdb, const double &c_ ) { init( rdb, c_ ); }
    void init( RetimeDB &rdb, const double &c_ );
    void do_bellmanFord_vertexD();
    void relax_one_round();
    bool is_negative_cycle();
 
    Edges ces;
    Vertices cvs;
};

void ConstGraph::init( RetimeDB &rdb, const double &c_ ){
    Edges &es_ = rdb.es;
    Vertices &vs_ = rdb.vs;
    retime::vec2D &D_ = rdb.D, &W_ = rdb.W;
    // vertex
    for( size_t i=0; i<vs_.size(); ++i ){
        Vertex &v = vs_[i];
        cvs.push_back( Vertex( v.id, v.name, MAX_SSSP ) );
    }
    cvs.push_back( Vertex( vs_.size(), "src", 0.0 ) );
    // edge
    for( size_t i=0; i<es_.size(); ++i ){
        Edge &e = es_[i];
        ces.push_back( Edge( e.id, string( e.t->name + "->" + e.s->name ), e.w, &cvs[ e.t->id ], &cvs[ e.s->id ] ) );
    }
    int n = D_.size();
    int m = ces.size();
    for( size_t i=0; i<n; ++i ){
        for( size_t j=0; j<n; ++j ){
            if( D_[i][j] <= c_ ) { continue; }
            Vertex &s = cvs[j], &t = cvs[i];
            ces.push_back( Edge( m, string( s.name + "->" + t.name ), W_[i][j] - 1, &s, &t ) );
            ++m;
            //cout << D_[i][j] << " " << ces.back() << endl;
        }
    }
    Vertex &src = cvs.back();
    for( size_t i=0; i<n; ++i ){
        ces.push_back( Edge( m, string( src.name + "->" + cvs[i].name ), 0, &src, &cvs[i] ) );
        ++m;
    }
    
    for( size_t i=0; i<ces.size(); ++i ){
        Edge &e = ces[i];
        Vertex &s = *e.s, &t = *e.t;
        s.outs.push_back( &e );
        t.ins.push_back( &e );
    }
    //for( size_t i=0; i<ces.size(); ++i ){ cout << ces[i] << endl; }
    //for( size_t i=0; i<cvs.size(); ++i ){ cout << cvs[i] << endl; }
    //cout << c_ << endl; getchar();
}

void ConstGraph::do_bellmanFord_vertexD()
{
    // initialization for d value
    double maxD = MAX_SSSP;
    for( size_t i=0; i<cvs.size()-1; ++i ) { cvs[i].d = maxD; }
    cvs.back().d = 0;
    for( size_t i=1; i<cvs.size(); ++i ) { relax_one_round(); }
}

void ConstGraph::relax_one_round()
{
    for( size_t j=0; j<ces.size(); ++j ){
        Edge &e = ces[j];
        Vertex &u = *e.s, &v = *e.t;
        if( u.d + e.w < v.d ){
            v.d = u.d + e.w;
        }
    }
}

bool ConstGraph::is_negative_cycle()
{
    vector<double> finalD; finalD.clear();
    for( size_t i=0; i<cvs.size(); ++i )
        finalD.push_back( cvs[i].d );
    relax_one_round();
    for( size_t i=0; i<cvs.size(); ++i )
        if( cvs[i].d != finalD[i] )
            return true;
    return false;
}
#endif
