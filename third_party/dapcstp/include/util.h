/**
 * \file   util.h
 * \brief  loaders
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */
#ifndef UTIL_H_
#define UTIL_H_

#include "def.h"
#include "inst.h"
#include "sol.h"

Inst load(const char* fn);
Inst loadPCSTP(const char* fn);
Inst loadNWSTP(const char* fn);
Inst loadMWCS(const char* fn);

double getBestKnownBound(const char* fn, const char* boundfile);
void   printBound(Inst& inst, weight_t bound);
void   printBoundPadded(Inst& inst, weight_t bound);
double format(weight_t bound, Inst& inst);

void recoverPartialSol(Sol& sol, Inst& inst1);
int  cntReachable(int r, Sol& sol, Inst& inst);
Sol  genPartialSol(Sol& sol, Inst& inst);
Inst genInst(Sol& sol, Inst& inst1, vector<int>& amap, vector<int>& nmap);

void enlargeStack();

void writeSolution(const char* file, Inst& inst, Sol& sol);

Sol loadSol(const char* fn, Inst& inst);

weight_t oracle(Inst& inst);
void writeInstance(const char* file, Inst& inst);
Sol dmst(Inst& inst, vector<weight_t>& cr);

#endif // UTIL_H_
