/**
 * \file   bbnode.h
 * \brief  branch-and-bound node
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#ifndef BBNODE_H_
#define BBNODE_H_

#include "inst.h"
#include "ds.h"

class BBNode {
public:
	weight_t lb = 0.0;
	int v = -1, bdir = -1; // branching node and direction (0 or 1)
	int depth = 0;
	int n = 0, m = 0, nfree = 0;
	bool feas = false;
	bool processed = false;
	int state = -1;
	int v2 = -1;
	Inst* inst;

	// queue positions in B&B
	PQMin<weight_t,BBNode*>::handle_type pqposMin;
	PQMax<weight_t,BBNode*>::handle_type pqposMax;

	BBNode(Inst* inst);

	// creating node rooted at given root
	BBNode(Inst* inst, int root, vector<int>& fe0);

	// create a node by branching on variable var (bdir is 0 or 1)
	BBNode(const BBNode* b, int var, int bdir);

	// updates instance graph size
	void updateNodeSize();
};

#endif // BBNODE_H_
