#include "util.h"
#include "retimeDB.h"

using namespace std;

int main( int argc, char **argv )
{
	if( argc != 2 ) { cout<<"./rt <in_name>"<<endl; exit( 0 ); }

    retime::RetimeDB rt;
	rt.init_from_dot( argv[1] );

    return 0;
}
