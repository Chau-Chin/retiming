#ifndef _retime_h_
#define _retime_h_

#include "util.h"
namespace retime
{

struct Edge;
struct Vertex;
class RetimeDB;

typedef vector<Edge> 				Edges;
typedef vector<Edge*> 				Eptrs;
typedef vector<Vertex> 				Vertices;
typedef vector<Vertex*> 			Vptrs;
typedef vector< vector<double> >	vec2D;

struct Edge
{
	Edge( int id_, string name_, double w_, Vertex *s_, Vertex *t_ ):
		id( id_ ), name( name_ ), w( w_ ), s( s_ ), t( t_ ) {}
	
	int id;
	string name;
	int w;	// #latches in this edge
	Vertex *s, *t;
    friend ostream &operator << ( ostream&, const Edge& );

  private:
	Edge() {}
};

struct Vertex
{
	Vertex( int id_, string name_, double d_ ):
		id( id_ ), name( name_ ), d( d_ ), r(0) {}

	int 	id;
	string 	name;
	double 	d;	// delay of this vertex
    double  r;  // retime function of this vertex
	Eptrs 	ins, outs;
    friend ostream &operator << ( ostream&, const Vertex& );

  private:
	Vertex() {}
};


class RetimeDB
{
  public:
	RetimeDB(){}

	void    init_from_dot( const string &inName );
	void    compute_W_D();
    void    init_clock_cycle_candidate_C();
    void    compute_opt_clock_cycle_optC();
    bool    is_all_constraints_satisfied_BellmanFord( const double& ); 
    void    compute_retime_function_R();
    double  get_clock_cycle_from_R( vector<double> R_ );

  //private:
    Edge   &get_e( const size_t &i )   { return es[ i ]; }
    Vertex &get_v( const size_t &i ) { return vs[ i ]; }
    double &get_w( const size_t &i, const size_t &j ) { return W[i][j]; }
    double &get_d( const size_t &i, const size_t &j ) { return D[i][j]; }
    vector<double> get_R(){
        int n = vs.size();
        vector<double> retR( n, 0 );
        for(size_t i=0; i<n; ++i) { retR[i] = vs[i].r; }
        return retR;
    }

	Edges 		es;
	Vertices 	vs;
	vec2D		W, D;
    vector<double> C; // candidate cycle time
    double      optC;
};

inline ostream& operator << (ostream& out, const Edge &e){
    out << " edge(" <<e.id<< ") : " << e.name << " : " << e.w
		<< " " << e.s->name << " -> " << e.t->name;
    return out;
}
inline ostream& operator << (ostream& out, const Vertex &v){
    out << " vertex(" <<v.id<< ") : " << v.name << " : " << v.d
		<< " #ins:" << v.ins.size() << " #outs:" << v.outs.size() << " r:" << v.r;
    return out;
}

} // end of namespace retime
#endif
