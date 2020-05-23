/**
 * \file   bbnode.cpp
 * \brief  branch-and-bound node
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#include "bbnode.h"

BBNode::BBNode(Inst* inst) {
	v = -1;
	bdir = -1;
	this->inst = inst;
}

BBNode::BBNode(Inst* _inst, int root, vector<int>& fe0) {
	v = root;
	bdir = 1;

	inst = new Inst(*_inst);
	inst->r = root;

	if(root != -1) {
		inst->p[root] = WMAX;
		inst->T[root] = true;
		inst->f1[root] = true;
	}

	for(int ij : fe0) {
		if(inst->fe0[ij]) continue;
		inst->delArc(ij);
		inst->fe0[ij] = true;
	}
}

BBNode::BBNode(const BBNode* b, int var, int bdir) : inst(b->inst), bdir(bdir) {
	lb = b->lb;
	v = var;
	depth = b->depth+1;
	n = b->n;
	m = b->m;
	processed = false;

	if(bdir == 0) {
		inst = new Inst(*b->inst);
	}

	if(bdir == 0) {
		inst->removeNode(v);
	} else if(bdir == 1) {
		inst->f1[v] = true;
		inst->T[v] = true;
	}
}
void BBNode::updateNodeSize()
{
	n = 0;
	nfree = 0;
	for(int i = 0; i < inst->n; i++) {
		if(inst->f0[i]) continue;
		n++;
		if(inst->f1[i]) continue;
		nfree++;
	}

	m = 0;
	for(int ij = 0; ij < inst->m; ij++) {
		if(inst->fe0[ij]) continue;
		m++;
	}
}
