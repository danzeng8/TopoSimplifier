/**
 * \file   procstatus.cpp
 * \brief  read and check memory limit
 *
 * \author Mario Ruthmaier 
 * \date   2015-05-03
 */

#include "procstatus.h"
#include <string>

u_int ProcStatus::memlimit = numeric_limits<u_int>::max();
u_int ProcStatus::maxusedmem = 0;

ProcStatus::ProcStatus()
{

}

void ProcStatus::setMemLimit( u_int lim )
{
	memlimit = lim;
}

u_int ProcStatus::mem()
{
	unsigned long vsize;
	{
		string ignore;
		ifstream ifs( "/proc/self/stat", std::ios_base::in );
		ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
			>> ignore >> ignore >> vsize; // position 23
	}
	double mb = (double) vsize / (1024 * 1024);
	return ceil( mb );
}

bool ProcStatus::memOK()
{
	u_int mb = mem();
	if( mb > maxusedmem ) maxusedmem = mb;
	if( maxusedmem > memlimit ) {
		cout << "### Memory-Usage too high: " << maxusedmem << " MB\n";
		return false;
	}
	else return true;
}

ProcStatus::~ProcStatus()
{

}
