#include "retimeDB.h"

namespace retime
{

void RetimeDB::init_from_dot( const string &inName )
{
	es.clear(); vs.clear();

    ifstream in( inName.c_str() );
    string tmp, tmp1, tmp2;
	map< string, int > name2Vidmap; name2Vidmap.clear();
	int numVs = 0;

	// read vertices
	while( in >> tmp ){
		if( tmp == "digraph" ) { break; }
		if( tmp != "//" ) { cout << "vertex should begin with //" <<endl; exit( 0 ); }

		in >> tmp1 >> tmp2;
		vs.push_back( Vertex( numVs, tmp1, atof( tmp2.c_str() ) ) );
		++numVs;
	}

	cout << "#vertices : " << numVs << endl;
	for( size_t i=0; i<numVs; ++i ){
		Vertex &v = vs[i];
		name2Vidmap[ v.name ] = v.id;
		//cout << v << endl;
	}

	int numEs = 0;
	// read edges
	while( in >> tmp ){
		if( name2Vidmap.find( tmp ) == name2Vidmap.end() ) { continue; }
		
		in >> tmp1 >> tmp2;
		Vertex *s = &vs[ name2Vidmap[ tmp ] ], *t = &vs[ name2Vidmap[ tmp2 ] ];
		
		in >> tmp >> tmp1 >> tmp2;
		if( tmp2.find("\"")==string::npos ) { cout << "format error in edge description!" << endl; exit(0); }
		tmp = tmp2.substr( tmp2.find("\"")+1, tmp2.rfind("\"")-1 );
		int weight = atoi( tmp.c_str() );
		es.push_back( Edge( numEs, string( s->name + "-" + t->name ), weight, s, t ) );
		++numEs;
	}

	cout << "#edges : " << numEs << endl;
	for( size_t i=0; i<numEs; ++i ){
		Edge &e = es[i];
        e.s->outs.push_back( &e );
        e.t->ins.push_back( &e );
		//cout << e << endl;
	}

	in.close();
}

#define MAX_C 10000.0
//numeric_limits<int>::max()/100000

struct Cost
{
    Cost( int w_=MAX_C, double d_=-MAX_C ) { init(w_, d_); }
    void init( int w_, double d_ ){ w=w_; d=d_; }
    Cost operator + ( const Cost & c_ ){
        return Cost( w + c_.w, d + c_.d );
    }
    bool operator < ( const Cost & c_ ){
        bool isLess1 = w < c_.w;
        bool isLess2 = ( w==c_.w && w!=MAX_C && d<c_.d );
        bool isLess = isLess1 || isLess2;
        //cout << "\t" << *this << " <=? " << c_ << " : " << isLess1 << " " << isLess2 << endl;
        return isLess;
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

typedef vector< vector<Cost> > vec2D_C; 


bool is_inf( const Cost& c ) { return (c.w >= MAX_C ); }
void show_matx( const vec2D_C &v )
{
    int n = v.size();
    for( int i = 0; i<n; ++i ){
        for( int j = 0; j<n; ++j ){
            if( i==0 )              { cout << j-1; }
            else if( j==0 )         { cout << i-1; }
            if( is_inf( v[i][j] ) ) { cout << "."; }
            else                    { cout << v[i][j]; }
            cout << "\t"; }
        cout << endl; }
    cout << endl;
}
void show_matx( const vec2D &v )
{
    int n = v.size();
    for( int i = 0; i<n; ++i ){
        for( int j = 0; j<n; ++j ){
            if( i==0 )                 { cout << j-1; }
            else if( j==0 )            { cout << i-1; }
            if( is_inf( v[i][j] ) )    { cout << "."; }
            else                       { cout << v[i][j]-1; }
            cout << "\t"; }
        cout << endl; }
    cout << endl;
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
struct Graph
{
    Graph( int n ) { init( n ); }
    void init( int n ){
        c.init( MAX_C, -MAX_C );
        W.resize( n, vector< Cost >(n, c) );
        reset_matx();
        delays.resize(n, -1);
    }
    void apsp(){
        size_t n = W.size()-1;
        reset_matx();
        init_L_P_1_matx();
        for( size_t m=2; m<n; ++m ){
            ext_sp( m );
            //show_stat(m);
        }
        for( size_t i=0; i<n; ++i ){
            for( size_t j=0; j<n; ++j ){
                ls(n-1,i,j) = ls(n-1,i,j) + Cost( 0, delays[j] );
            }
        }
        //show_stat(n-1);
    }
    void ext_sp( size_t m ){
        size_t n = W.size()-1;

        for( size_t i=0; i<n; ++i )
            for( size_t j=0; j<n; ++j ){
                ps(m,i,j) = ps(m-1,i,j);
                for( size_t k=0; k<n; ++k ){
                    Cost relaxCost = ls(m-1,i,k) + ws(k,j);
                    //cout << i << " " << j << " " << k 
                    //     << "\t relaxC: " << relaxCost 
                    //     << "\t l'(i,k):" << ls(m-1,i,k) 
                    //     << "\t w(k,j):" << ws(k,j)
                    //     << " < l :" << ls(m,i,j) << endl;
                    //cout << "relax cost for " << i+1 << "->" << j+1 << " " << k+1 << "-iter : " << relaxCost << endl;
                    if( relaxCost < ls(m,i,j) ){
                        ls(m,i,j) = relaxCost;
                        ps(m,i,j) = k+1;
                        //cout << "p[" << i << "][" << j << "] = " << k << endl;
                    }
                }
                //cout << " L:: " << endl;
                //show_matx( L[m] );
                //cout << " P:: " << endl;
                //show_matx( P[m] ); //getchar();
            }
    }
    void show_stat( int m )
    {
        cout << "**********************" << endl;
        //cout << "--- W(" << m << ") ---" << endl;
        //show_matx(W);
        cout << "--- L(" << m << ") ---" << endl;
        show_matx(L[m]);
        cout << "--- P(" << m << ") ---" << endl;
        show_matx(P[m]);
        cout << "**********************" << endl;
    }
    void reset_matx(){
        size_t n=W.size();
        P.clear(); L.clear();
        P.resize( n, vec2D( n, vector<double>(n, MAX_C) ) );
        L.resize( n, vec2D_C( n, vector<Cost>(n, c) ) );
        for( size_t s=1; s<n; ++s ){ for( size_t i=1; i<n; ++i ) { L[s][i][i] = Cost(0,0); } }
    }
    void init_L_P_1_matx()
    {
        size_t n = W.size()-1;
        for( size_t i=0; i<n; ++i ) for( size_t j=0; j<n; ++j )
            if( !is_inf( ws(i,j) ) ){
                ls(1,i,j) = ws(i,j);
                ps(1,i,j) = i+1;
            }
        //show_stat(1);
    }

    Cost &w( const int &i, const int &o )                   { return W[i][o]; }
    Cost &ws( const int &i, const int &o )                  { return W[i+1][o+1]; }
    double &p( const int &n, const int &i, const int &o )   { return P[n][i][o]; }
    double &ps( const int &n, const int &i, const int &o )  { return P[n][i+1][o+1]; }
    Cost &l( const int &n, const int &i, const int &o )     { return L[n][i][o]; }
    Cost &ls( const int &n, const int &i, const int &o )    { return L[n][i+1][o+1]; }
    vector<double>      delays;
    vec2D_C             W;    // adjacent matrix
    vector< vec2D >     P;    // predecessor matrix
    vector< vec2D_C >   L;    // shortest path cost matrix
    Cost c;
};

void RetimeDB::compute_W_D()
{
    size_t n = vs.size(), numEs = es.size();
    W.clear(); W.resize( n, vector<double>(n,-1) );   
    D.clear(); D.resize( n, vector<double>(n,-1) );   
    
    Graph g( n+1 ); // one is for dummy

    // init graph
    for( size_t i=0; i<n; ++i ){
        Vertex &u = get_v( i );
        for ( size_t j = 0; j<u.outs.size(); ++j ){
            Edge &e = *u.outs[j];
            Vertex &v = *e.t;
            g.ws( u.id, v.id ) = Cost( e.w, -u.d );
        }
        g.delays[ u.id ] = -u.d;
    }
    //show_matx( g.W );
    g.apsp();
    
    for( size_t i=0; i<n; ++i ){
        for( size_t j=0; j<n; ++j ){
            get_w(i,j) = g.ls(n-1,i,j).w;
            get_d(i,j) = -g.ls(n-1,i,j).d;
        }
    }
    show_matx_pure( W, "W" );
    show_matx_pure( D, "D" );
}

void RetimeDB::init_clock_cycle_candidate_C()
{
    C.clear();

    map<double,bool> clkMap;
    map<double,bool>::iterator iter;
    int n = D.size();
    for( size_t i=0; i<n; ++i )
        for( size_t j=0; j<n; ++j )
            clkMap[ get_d(i,j) ] = true;

    int cnt = 0;
    for( iter=clkMap.begin(); iter!=clkMap.end(); ++iter )
        C.push_back( iter->first ); 

    int numCandClks = C.size();
    cout << "numCandClks : " << numCandClks << endl;
}

void RetimeDB::compute_opt_clock_cycle_optC()
{
    optC = 0;
    
    for(int i = int(C.size()-1); i>=0; --i){
        if( !is_all_constraints_satisfied_BellmanFord( C[i] ) ){
            optC = C[i+1];
            break;
        }
        cout << "feasible for clk : " << C[i] << endl;
    }
}


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
    vec2D &D_ = rdb.D;
    vec2D &W_ = rdb.W;
    // vertex
    for( size_t i=0; i<vs_.size(); ++i ){
        Vertex &v = vs_[i];
        cvs.push_back( Vertex( v.id, v.name, 100000 ) );
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
    double maxD = 10000;
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
bool RetimeDB::is_all_constraints_satisfied_BellmanFord( const double &c )
{
    ConstGraph cg( *this, c );
    cg.do_bellmanFord_vertexD();
    return !cg.is_negative_cycle();
}

void RetimeDB::compute_retime_function_R()
{
    ConstGraph cg( *this, optC );
    cg.do_bellmanFord_vertexD();
    
    for( size_t i=0; i<vs.size(); ++i )
        vs[i].r = cg.cvs[i].d;
    
    for( size_t i=0; i<vs.size(); ++i )
        cout << vs[i] << endl;
}

double RetimeDB::get_clock_cycle_from_R( vector<double> R_ )
{
    //int m = es.size();
    //vector<double> Wr( m, 0);
    //for( size_t i=0; i<m; ++i ){
    //    Edge &e = es[i];
    //    Vertex &s = *e.s, &t = *e.t;
    //    Wr[i] = e.w + R_[ t.id ] - R_[ s.id ];
    //}
    
}

} // end of namespace retime
