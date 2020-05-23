/**
 * \file   heur.cpp
 * \brief  primal heuristics implementation
 *
 * \author Martin Luipersbeck
 * \date   2013-05-21
 */

using namespace std;

#include "heur.h"
#include "options.h"
#include "ds.h"
#include "util.h"
#include "bounds.h"

#include <list>
#include <stack>

Sol primI(int r, Inst& inst, vector<weight_t>& tw)
{
	Sol sol(inst);
	PQMin<weight_t,int> PQ;
	vector<weight_t> dist(inst.n, WMAX);
	vector<int> pred(inst.n, -1);

	if(r == -1) {
		vector<int> roots;
		for(int i = 0; i < inst.n; i++) {
			if(inst.f0[i]) continue;
			roots.push_back(i);
		}
		r = roots[rand()%roots.size()];
	}

	sol.obj = 0;
	if(inst.bigM > 0) {
		int f1 = 0;
		for(int i = 0; i < inst.n; i++) {
			if(inst.f1[i])
				f1++;
		}
		vector<int> rootArcs;
		for(int ri : inst.dout[inst.r]) {
			if(tw[ri] < WMAX-1) {
				int cnt = inst.countReachableFixed(inst.head[ri]);
				if(cnt == f1)
					rootArcs.push_back(ri);
			}
		}
		if(rootArcs.size() == 0) {
			sol.obj = WMAX;
			return sol;
		}
		
		int ri = rootArcs[rand()%rootArcs.size()];
		sol.arcs[ri] = true;
		
		sol.obj += inst.c[ri];
		sol.nodes[inst.r] = true;
		r = inst.head[ri];
		sol.r = inst.r;
	} else {
		sol.r = r;
	}

	dist[r] = 0;
	sol.nodes[r] = true;
	PQ.push(make_pair(0, r));
	int terms = -1;

	for(int i = 0; i < inst.n; i++) {
		if(!inst.f1[i] && i != r) {
			sol.obj += inst.p[i];
		}
		if(!inst.f0[i] && (inst.T[i] || inst.f1[i]))
			terms++;
	}
	
	while ( !PQ.empty() ) {
		int i = PQ.top().second;
		
		PQ.pop();

		if((!inst.T[i] && !inst.f1[i]) || sol.nodes[i]) {
			for(int ij : inst.dout[i]) {
				const int j = inst.head[ij];
				weight_t d;
				if(dist[i] < WMAX-1 - tw[ij]) // use this to find a valid path even on graphs with WMAX weights
					d = dist[i] + tw[ij];
				else
					d = WMAX-1;

				if(d < dist[j]) {
					dist[j] = d;
					pred[j] = ij;
					PQ.push(make_pair(dist[j], j));
				}
			}
		} else {
			while(!sol.nodes[i]) {
				sol.nodes[i] = true;
				if(!inst.f1[i]) {
					sol.obj -= inst.p[i];
				}

				dist[i] = 0;
				PQ.push(make_pair(dist[i], i));

				int ij = pred[i];
				int j = inst.tail[ij];
				if(i == j) {
					ij = inst.opposite[ij];
					assert(ij != -1);
					j = inst.tail[ij];
				}

				sol.obj += inst.c[ij];
				sol.arcs[ij] = true;
				i = j;
			}
		}
	}
	sol.obj += inst.offset;

	sol.recomputeObjective();
	if(!sol.validate()) {
		sol.obj = WMAX;
		return sol;
	}
	
	sol.recomputeObjective();
	strongprune(inst, sol);
	assert(sol.obj >= 0);

	return sol;
}

void strongprune(Inst& inst, Sol& sol)
{
	vector<int> pred(inst.n, -1);
	vector<weight_t> l(inst.n);
	vector<flag_t> fixed(inst.n, 0);
	vector<flag_t> visited(inst.n, false);
	list<int> bfst;
	stack<int> Q;
	int v, w, u;
	int r = sol.r;

	// fixed means a node in the subtree is fixed to one, so it cannot be pruned
	for(int i = 0; i < inst.n; i++) {
		l[i] = sol.nodes[i] ? inst.p[i] : 0;
		if(inst.f1[i]) fixed[i] = true;
	}

	if(inst.bigM > 0) {
		// prevent artificial root arc from being pruned
		for(int ri : inst.dout[inst.r]) {
			if(sol.arcs[ri]) {
				fixed[inst.head[ri]] = true;
				break;
			}
		}
	}

	Q.push(r);
	visited[r] = true;
	while (!Q.empty()) {
		int i = Q.top();
		Q.pop();

		for(int ij : inst.dout[i]) {
			int j = inst.head[ij];
			if (!sol.arcs[ij] || visited[j]) continue;

			pred[j] = ij;
			
			visited[j] = true;
			Q.push(j);
			bfst.push_front(j);
		}
	}

	for(int i : bfst) {
		int ij = pred[i];
		u = inst.tail[ij];
		
		const weight_t childvalue = l[i] - inst.c[ij];

		if (childvalue > 0 || fixed[i]) {
			// overflow check
			if(l[u] <= WMAX - childvalue)
				l[u] += childvalue;
			else {
				l[u] = WMAX;
			}

			if(fixed[i])
				fixed[u] = true;

		} else {
			sol.nodes[i] = false;
			sol.arcs[ij] = false;
			pred[i] = -1;

			sol.obj += childvalue;

			Q.push(i);
			while (!Q.empty()) {
				i = Q.top();
				Q.pop();

				for(int f : inst.dout[i]) {
					w = inst.head[f];
					if(pred[w] != f) continue;
					
					sol.nodes[w] = false;
					sol.arcs[f] = false;
					pred[w] = -1;
					Q.push(w);
				}
			}
		}
	}
}

