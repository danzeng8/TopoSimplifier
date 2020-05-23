/**
 * \file   heur.h
 * \brief  primal heuristics
 *
 * \author Martin Luipersbeck
 * \date   2013-05-21
 */

#ifndef HEUR_H_
#define HEUR_H_

#include "def.h"
#include "inst.h"
#include "sol.h"

void strongprune(Inst& inst, Sol& sol);
Sol primI(int r, Inst& inst, vector<weight_t>& tw);

#endif // HEUR_H_
