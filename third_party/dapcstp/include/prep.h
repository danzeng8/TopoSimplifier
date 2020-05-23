/**
 * \file   prep.h
 * \brief  reduction tests
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#ifndef PREP_H_
#define PREP_H_

#include "def.h"
#include "inst.h"

void costShift(Inst& inst);

int MAcutarc(Inst& inst, vector<flag_t>& ap, vector<int>& lastap);
int MAcutnode(Inst& inst, vector<flag_t> ap, vector<int>& lastap);
int MA(Inst& inst);

int ntd1(Inst& inst);
int ntd2(Inst& inst);

int lc(Inst& inst);

bool nr(Inst& inst);

pair<int,int> bbred(Inst& inst, weight_t lb, weight_t ub, vector<weight_t>& cr, vector<weight_t>& pr);

int APfixing(Inst& inst, vector<flag_t> ap, vector<int>& lastap);

#endif // PREP_H_
