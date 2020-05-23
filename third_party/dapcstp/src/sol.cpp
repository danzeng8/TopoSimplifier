/**
 * \file   sol.cpp
 * \brief  problem solution class
 *
 * \author Martin Luipersbeck 
 * \date   2015-01-03
 */

#include "sol.h"
#include <stack>
#include <iostream>

Sol::~Sol()
{

}

Sol::Sol(const Inst& inst) : 
	nodes(inst.n, false), arcs(inst.m, false), inst(inst)
{

}

Sol::Sol(const Sol& S) :
	nodes(S.nodes), arcs(S.arcs), r(S.r), obj(S.obj), partial(S.partial), inst(S.inst)
{

}

Sol& Sol::operator=(const Sol& S)
{
	if(&S == this) return *this;

	nodes=S.nodes;
	arcs=S.arcs;
	r=S.r;
	obj=S.obj;
	partial=S.partial;

	return *this;
}

void Sol::update(Sol& S) {
	if(S.obj < obj) {
		arcs = S.arcs;
		nodes = S.nodes;
		r = S.r;
		obj = S.obj;
	}
}

int Sol::rootSolution(int r)
{
	int rootOutgoing = 0;
	for(int ij : inst.dout[r]) {
		if(arcs[ij]) {
			rootOutgoing++;
		}
	}
	
	if(rootOutgoing == 1)
		return 0;

	vector<int> visited(inst.n, false);
	visited[r] = true;

	stack<int> stack;
	stack.push(r);
	
	int cnt = 0;
	while(!stack.empty()) {
		const int i = stack.top();
		stack.pop();
		for(int ij : inst.dout[i]) {
			const int j = inst.head[ij];

			if(inst.opposite[ij] == -1) {
				continue;
			}
			const int ji = inst.opposite[ij];

			if(visited[j]) continue;

			if(!arcs[ij] && arcs[ji]) {
				arcs[ij] = 1;
				arcs[ji] = 0;
				cnt++;
			}

			if(arcs[ij]) {
				visited[j] = true;
				stack.push(j);
			}
		}
	}

	return cnt;
}


bool Sol::validate()
{
	vector<int> visited(inst.n, false); stack<int> Q;

	if(r == -1)
		return false;

	if(inst.bigM > 0) {
		int cnt = 0;
		for(int ri : inst.dout[inst.r]) {
			if(arcs[ri]) cnt++;
		}
		if(cnt != 1) {
			return false;
		}
	}

	Q.push(r);
	visited[r] = true;
	weight_t obj1 = inst.offset;
	
	while(!Q.empty()) {
		const int i = Q.top();
		Q.pop();
		for(int ij : inst.dout[i]) {
			const int j = inst.head[ij];
			if(visited[j] || !arcs[ij]) continue;
			obj1 += inst.c[ij];

			visited[j] = true;
			Q.push(j);
		}
	}
	// are all nodes reachable/non-reachable based on specified arcs/nodes?
	for(int i = 0; i < inst.n; i++) {
		if(nodes[i] && !visited[i]) {
			cout << "not reachable " << endl;
			return false;
		}
		if(!nodes[i] && visited[i]) {
			cout << "not reachable 2 " << i << endl;
			return false;
		}
		if(!visited[i]) {
			obj1 += inst.p[i];
		}
	}
	
	// does the specified objective match the actual one?
	if(obj != obj1) {
		return false;
	}
	if(obj < 0) {
		return false;
	}

	return true;
}

weight_t Sol::recomputeObjective()
{
	obj = inst.offset;
	for(int i = 0; i < inst.m; i++) {
		if(arcs[i]) obj += inst.c[i];
	}
	for(int i = 0; i < inst.n; i++) {
		if(!nodes[i]) obj += inst.p[i];
	}
	return obj;
}

