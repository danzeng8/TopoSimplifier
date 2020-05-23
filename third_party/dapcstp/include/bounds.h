/**
 * \file   bounds.h
 * \brief  procedures for computing lower bounds
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#ifndef BOUNDS_H_
#define BOUNDS_H_

#include "def.h"
#include "ds.h"
#include "inst.h"
#include "sol.h"
#include "options.h"

// dual ascent
template<typename U> weight_t daR(int r, Inst& inst, vector<U>& c, vector<U>& cr, vector<U>& pr, weight_t ub, double eager, Sol* inc, bool heur = false);

#endif // BOUNDS_H_
