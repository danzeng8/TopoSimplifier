/**
 * \file   sol.h
 * \brief  problem solution class
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#ifndef SOL_H_
#define SOL_H_

#include "inst.h"

class Sol {
public:
	vector<flag_t> arcs, nodes;
	int r = -1;
	// partial solutions are solutions that still require backmapping (may contain anti-parallel arcs)
	bool partial = false;
	weight_t obj = WMAX;
	const Inst& inst;

	Sol(const Inst& inst);
	~Sol();
	Sol(const Sol& S);
	Sol& operator=(const Sol& S);

	// replace solution if given solution is better
	void update(Sol& S);

	// change root to r
	int rootSolution(int r);
	weight_t recomputeObjective();
	bool validate();
};

#endif // SOL_H_
