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
		cout << v << endl;
		name2Vidmap[ v.name ] = v.id;
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
		cout << e << endl;
	}

	in.close();
}
} // end of namespace retime
