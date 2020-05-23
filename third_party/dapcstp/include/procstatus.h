/**
 * \file   procstatus.cpp
 * \brief  read and check memory limit
 *
 * \author Mario Ruthmaier 
 * \date   2012-05-03
 */
#ifndef PROCSTATUS_H_
#define PROCSTATUS_H_

#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>

using namespace std;
typedef unsigned int u_int;

/*
 * class for querying current process status
 */
class ProcStatus
{

private:

	static u_int memlimit;

public:

	static u_int maxusedmem;

	ProcStatus();
	virtual ~ProcStatus();

	static void setMemLimit( u_int lim );
	static u_int mem();
	static bool memOK();

};

#endif /* PROCSTATUS_H_ */
