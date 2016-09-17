#include "util.h"
#include "retimeDB.h"

using namespace std;

int main( int argc, char **argv )
{
	if( argc != 2 ) { cout<<"./rt <in_name>"<<endl; exit( 0 ); }

    retime::RetimeDB rt;
	rt.init_from_dot( argv[1] );
    
    rt.compute_W_D();
    //cout << "current cycle before: " << rt.get_clock_cycle_from_R( rt.get_R() ) << endl;
    rt.init_clock_cycle_candidate_C();
    rt.compute_opt_clock_cycle_optC();
    rt.compute_retime_function_R();
    //cout << "cycle time after: " << rt.get_clock_cycle_from_R( rt.get_R() ) << endl;
    return 0;
}
