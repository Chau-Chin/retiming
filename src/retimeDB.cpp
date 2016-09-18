#include "retimeDB.h"
#include "apsp.h"
#include "sssp.h"

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

void RetimeDB::compute_W_D(bool isRetimed)
{
    size_t n = vs.size();
    W.clear(); W.resize( n, vector<double>(n,-1) );   
    D.clear(); D.resize( n, vector<double>(n,-1) );   
    
    ApspGraph g( n, *this, isRetimed );

    g.apsp();
    
    for( size_t i=0; i<n; ++i ){
        for( size_t j=0; j<n; ++j ){
            get_w(i,j) = g.l(n-1,i,j).w;
            get_d(i,j) = -g.l(n-1,i,j).d;
        }
    }
    //g.show_matx_pure( W, "W" );
    //g.show_matx_pure( D, "D" );
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

    //int numCandClks = C.size();
    //cout << "numCandClks : " << numCandClks << endl;
}

void RetimeDB::compute_opt_clock_cycle_optC()
{
    optC = 0;
    
    for(int i = int(C.size()-1); i>=0; --i){
        if( !is_all_constraints_satisfied_BellmanFord( C[i] ) ){
            optC = C[i+1];
            break;
        }
        //cout << "feasible for clk : " << C[i] << endl;
    }
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
    
    //for( size_t i=0; i<vs.size(); ++i )
    //    cout << vs[i] << endl;
}

double RetimeDB::get_clock_cycle_from_R()
{
    double clkCycle = 0;
    int m = es.size();
    vector<double> Wr( m, 0);
    for( size_t i=0; i<m; ++i ){
        Edge &e = es[i];
        Vertex &s = *e.s, &t = *e.t;
        e.wr = e.w + t.r - s.r;
        //cout << e << endl;
    }

    compute_W_D( true );
    
    for( size_t i=0; i<vs.size(); ++i ){
        for( size_t j=0; j<vs.size(); ++j ){
            if( get_w(i,j) == 0 ){
                clkCycle = max( clkCycle, get_d(i,j) );
            }
        }
    }
    return clkCycle;
}

} // end of namespace retime
