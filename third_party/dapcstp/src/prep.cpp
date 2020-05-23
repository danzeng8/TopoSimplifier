/**
 * \file   prep.cpp
 * \brief  reduction tests implementation
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#include "prep.h"
#include "util.h"
#include "stats.h"
#include "options.h"
//#include "tbb/parallel_for.h"

#include <iostream>
#include <unordered_set>
#include <stack>

void costShift(Inst& inst)
{
	for (int i = 0; i < inst.n; i++) {
		if(inst.f0[i] || inst.r == i || !inst.T[i]) continue;

		weight_t delta = WMAX;
		for(int ij : inst.din[i]) {
			if(inst.c[ij] < delta) { delta = inst.c[ij]; }
		}
		if(inst.din[i].size() == 0 || delta == 0) continue;

		// do not apply shift if the instance is asymmetric unrooted, or if the shift would apply to a potential root
		if(inst.r == -1 && (inst.isAsym || delta < inst.p[i])) {
			continue;
		}

		delta = min(delta, inst.p[i]);

		for(int ij : inst.din[i])
			inst.c[ij] -= delta;
		inst.p[i] -= delta;
		inst.offset += delta;

		if(inst.f1[i]) {
			inst.T[i] = true;
			inst.p[i] = WMAX;
		}
		if(inst.p[i] == 0) {
			if(!inst.f1[i]) {
				inst.T[i] = false;
			}
		}
	}
}

// merges APs to their successor if they just have one successor
int MAcutarc(Inst& inst, vector<flag_t>& ap, vector<int>& lastap)
{
	auto& c = inst.c; auto& p = inst.p; auto& T = inst.T;

	int cnt = 0;
	//tbb::parallel_for(tbb::blocked_range<int>(0, inst.m),
		//[&](tbb::blocked_range<int> r)
	//{
		//for (int ij = r.begin(); ij < r.end(); ++ij)
		//{
			for(int ij = 0; ij < inst.m; ij++) {
			if (inst.fe0[ij]) continue;

			const int i = inst.tail[ij];
			const int j = inst.head[ij];

			if (inst.r == i && inst.bigM > 0) continue; // don't merge in case of artificial root arcs
			if (inst.c[ij] > inst.p[j]) continue;
			if (!ap[i] || !ap[j]) continue;

			unordered_set<int> succ;
			for (int ij2 : inst.dout[i]) {
				const int j = inst.head[ij2];
				if (lastap[j] != i) continue;
				succ.insert(j);
			}

			if (succ.size() == 1 && *succ.begin() == j) {

				inst.updatebmMerge(ij, false);
				inst.merge(ij, i, j);
				cnt++;
			}
			}
		//}
	//});

	stats.ss += cnt;
	return cnt;
}

// if i is an articulation point, and its successor has no cheaper incoming arc, merge
int MAcutnode(Inst& inst, vector<flag_t> ap, vector<int>& lastap)
{
	auto& c = inst.c; auto& p = inst.p; auto& T = inst.T;
	
	int cnt = 0;
	for(int ij = 0; ij < inst.m; ij++) {
		if(inst.fe0[ij]) continue;

		const int i = inst.tail[ij];
		const int j = inst.head[ij];
		if(inst.r == i && inst.bigM > 0) continue; // don't merge in case of artificial root arcs

		if (ap[i] && lastap[j] == i && c[ij] <= p[j] &&
			!inst.hasCheaperIncomingArc(j, ij)) {

			inst.updatebmMerge(ij, false);
			inst.merge(ij, i, j);
			cnt++;
		}
	}

	stats.ms += cnt;
	return cnt;
}

int APfixing(Inst& inst, vector<flag_t> ap, vector<int>& lastap)
{
	int cnt = 0;
	for(int i = 0; i < inst.n; i++) {
		if(inst.f0[i] || !inst.f1[i]) continue;
		if(i == inst.r) continue;

		int j = i;
		while(j != inst.r) {
			j = lastap[j];
			if(j == -1) break;
			
			assert(ap[j]);
			if(!inst.f1[j]) {
				inst.f1[j] = true;
				inst.T[j] = true;
				inst.p[j] = WMAX;
				cnt++;
			}
		}
	}

	return cnt;
}

int MA(Inst& inst)
{
	auto& c = inst.c; auto& p = inst.p; auto& T = inst.T;
	Inst& inst1 = *inst.inst1;

	vector<weight_t> cheapestInc(inst.n, WMAX);
	vector<flag_t> checked(inst.m, false);
	for(int ij = 0; ij < inst.m; ij++) {
		if(inst.fe0[ij]) continue;
		const int j = inst.head[ij];
		if(inst.c[ij] < cheapestInc[j])
			cheapestInc[j] = inst.c[ij];
	}
	
	int cnt = 0;
	for(int ij = 0; ij < inst.m; ij++) {
		if(inst.fe0[ij]) continue;
		const int ji = inst.opposite[ij];

		if(ji == -1 || inst.fe0[ji]) continue;
		if(checked[ij]) continue;

		const int i = inst.tail[ij];
		const int j = inst.head[ij];

		if(inst.r == i || inst.r == j) continue;

		// if unrooted, asymmetric test cannot be applied
		if(inst.r == -1 && abs(c[ij] - c[ji]) > 0) {
			continue;
		}

		checked[ij] = true;
		checked[ji] = true;

		const weight_t pmin = min(p[i], p[j]);
		const weight_t cmax = max(c[ij], c[ji]);

		bool test = c[ij] <= p[j] && c[ji] <= p[i];

		if ( test ) {
			if(cheapestInc[i] < c[ji]) continue;
			if(cheapestInc[j] < c[ij]) continue;

			if(abs(c[ij] - c[ji]) > 0) {
				weight_t diff = c[ij] - c[ji];
				if(diff < 0) inst.increaseRevenue(j, -diff);
				else         inst.increaseRevenue(i, diff);
			}

			assert(c[ij] == c[ji]);

			inst.updatebmMerge(ij, true);
			inst.merge(ij, i, j);
			cnt++;

			for(int ij2 : inst.din[i]) {
				if(inst.c[ij2] < cheapestInc[i])
					cheapestInc[i] = inst.c[ij2];
			}
		}
	}

	stats.ma += cnt;
	return cnt;
}

int ntd1(Inst& inst)
{
	auto& c = inst.c; auto& p = inst.p;

	int cnt = 0;
	stack<int> deg1;
	for(int i = 0; i < inst.n; i++) {
		if(inst.f0[i]) continue;

		if(inst.singleAdjacency(i))
			deg1.push(i);
	}

	while (!deg1.empty()) {
		const int i = deg1.top();
		deg1.pop();

		if(i == inst.r) continue;
		if(!inst.singleAdjacency(i))
			continue;
		
		const int ji = inst.din[i].front();
		const int j = inst.tail[ji];

		if ( inst.c[ji] < p[i] || inst.f1[i] ) {
			inst.contractArc(ji);
		} else {
			inst.removeNode(i);
		}

		cnt++;

		if(inst.singleAdjacency(j))
			deg1.push(j);

	}
	
	stats.d1 += cnt;
	return cnt;
}

int ntd2(Inst& inst)
{
	unordered_set<int> deg2;
	int ij, ji, ki, ik, j, k;

	for(int i = 0; i < inst.n; i++) {
		if(i == inst.r || inst.f0[i] || inst.f1[i] || inst.T[i]) continue;
		deg2.insert(i);
	}

	int cnt = 0;
	while (!deg2.empty()) {
		int i = *deg2.begin();
		deg2.erase(i);
		
		if(i == inst.r || inst.f0[i] || inst.f1[i] || inst.T[i]) continue;

		if(!inst.doubleAdjacency(i, j, k, ki, ij, ji, ik))
			continue;

		int ij2 = -1;
		bool yes = true;
		if(ji != -1 && ik != -1) {
			for(int ij3 : inst.dout[j]) if(inst.head[ij3] == k) ij2 = ij3;
			const weight_t newWeight = inst.c[ji] + inst.c[ik] - inst.p[i];
			if(ij2 != -1) {
				yes = false;
			}
		}
		int ij4 = -1;
		if(ki != -1 && ij != -1) {
			for(int ij3 : inst.dout[k]) if(inst.head[ij3] == j) ij4 = ij3;
			const weight_t newWeight = inst.c[ki] + inst.c[ij] - inst.p[i];
			if(ij4 != -1) {
				yes = false;
			}
		}

		ij2 = -1;
		if(ji != -1 && ik != -1) {
			for(int ij3 : inst.dout[j]) if(inst.head[ij3] == k) ij2 = ij3;
			
			const weight_t newWeight = inst.c[ji] + inst.c[ik] - inst.p[i];
			
			if(ij2 != -1) {
				// triangle
				if(inst.c[ij2] > newWeight) {
					
					assert(inst.c[ij2] > inst.c[ji] + inst.c[ik]);
					inst.c[ij2] = newWeight;

					inst.updatebmNTD2triangle(ij2, ji, ik);
				}

				inst.delArc(ji);
				inst.delArc(ik);
				inst.fe0[ji] = true;
				inst.fe0[ik] = true;

			} else {
				
				inst.c[ji] = newWeight;

				inst.updatebmNTD2(ji, ik);

				inst.moveHead(ji, k);

				inst.delArc(ik);
				inst.fe0[ik] = true;
				ij2 = ji;
			}

			cnt++;
		}

		ij4 = -1;
		if(ki != -1 && ij != -1) {
			for(int ij3 : inst.dout[k]) if(inst.head[ij3] == j) ij4 = ij3;

			const weight_t newWeight = inst.c[ki] + inst.c[ij] - inst.p[i];

			if(ij4 != -1) {
				// triangle
				if(inst.c[ij4] > newWeight) {
					assert(inst.c[ij4] > inst.c[ki] + inst.c[ij]);
					inst.c[ij4] = newWeight;

					inst.updatebmNTD2triangle(ij4, ki, ij);
				}

				// this implies that i is fixed to zero
				inst.delArc(ki);
				inst.delArc(ij);
				inst.fe0[ki] = true;
				inst.fe0[ij] = true;

			} else {
				inst.c[ki] = newWeight;
				
				inst.updatebmNTD2(ki, ij);
				
				inst.moveHead(ki, j);

				inst.delArc(ij);
				inst.fe0[ij] = true;
				ij4 = ki;
			}

			cnt++;
		}

		// update opposite arcs if possible
		if(ij4 != -1 && ij2 != -1) {
			inst.opposite[ij2] = ij4;
			inst.opposite[ij4] = ij2;
		} else {
			if(ij2 != -1) {
				int opp1 = inst.opposite[ij2];
				inst.opposite[ij2] = -1;
				if(opp1 != -1) inst.opposite[opp1] = -1;
			}
			if(ij4 != -1) {
				int opp1 = inst.opposite[ij4];
				inst.opposite[ij4] = -1;
				if(opp1 != -1) inst.opposite[opp1] = -1;
			}
		}

		inst.removeNode(i);
		//inst.bmnn[i].clear();
	}
	
	stats.d2 += cnt;
	return cnt;
}

pair<int,int> bbred(Inst& inst, weight_t lb, weight_t ub, vector<weight_t>& cr, vector<weight_t>& pi)
{
	if(!params.boundbased) return make_pair(0,0);

	vector<weight_t> distR(inst.n, WMAX), distT(inst.n, WMAX);
	//PQMin<weight_t,int> PQ;
	priority_queue<pair<weight_t, int>, std::vector<pair<weight_t, int>>, greater<pair<weight_t, int>>> PQ;
	
	distR[inst.r] = 0;

	PQ.push(make_pair(0.0, inst.r));

	while (!PQ.empty()) {
		int i = PQ.top().second;
		PQ.pop();

		for(int ij : inst.dout[i]) {
			int j = inst.head[ij];
			const weight_t d = distR[i] + cr[ij];
			if (distR[j] > d) {
				distR[j] = d;
				PQ.push(make_pair(d, j));
			}
		}
	}
	//PQ.clear();

	for(int i = 0; i < inst.n; i++) {
		if(inst.f0[i]) continue;
		if(inst.T[i]) {
			distT[i] = 0;
			PQ.push(make_pair(0.0, i));
		}
	}
	while (!PQ.empty()) {
		int i = PQ.top().second;
		PQ.pop();

		for(int ij : inst.din[i]) {
			int j = inst.tail[ij];
			const weight_t d = distT[i] + cr[ij];
			if (distT[j] > d) {
				distT[j] = d;
				PQ.push(make_pair(d, j));
			}
		}
	}
	
	int narcs = 0, nnodes = 0;
	vector<flag_t> arcs2del(inst.m, 0);

	// fix Steiner nodes to zero (remove sets of incident arcs in next loop, since they may overlap)
	for(int i = 0; i < inst.n; i++) {
		if(inst.f0[i] || inst.f1[i] || inst.T[i]) continue;

		const weight_t lhs = lb + distR[i] + distT[i];

		if ( lhs >= ub ) {
			inst.removeNode(i);
			stats.boundbased++;
		}
	}

	// fix potential terminals to one
	for(int i = 0; i < inst.n; i++) {
		if(inst.f0[i] || inst.f1[i] || !inst.T[i]) continue;

		const weight_t lhs = lb + pi[i];

		if( lhs >= ub ) {
			inst.f1[i] = true;
			inst.T[i] = true;
			inst.p[i] = WMAX;
			nnodes++;
			stats.boundbased++;
		}
	}
	
	// fix arcs to zero
	
	for(int ij = 0; ij < inst.m; ij++) {
		if(inst.fe0[ij]) continue;
		
		const int i = inst.tail[ij];
		const int j = inst.head[ij];

		const weight_t lhs = lb + distR[i] + cr[ij] + distT[j];

		if ( lhs >= ub ) {
			inst.delArc(ij);
			inst.fe0[ij] = true;
			narcs++;
			stats.boundbased++;
		}
	}
	
	return make_pair(nnodes, narcs);
}

int lc(Inst& inst)
{
	int cnt = 0;
	vector<weight_t> dist(inst.n, WMAX);
	vector<flag_t> hops(inst.n, 0), processed(inst.n, 0);
	vector<int> marked, todel;

	marked.reserve(inst.n);
	todel.reserve(inst.n);

	//PQMin<weight_t,int> PQ;
	priority_queue<pair<weight_t, int>, std::vector<pair<weight_t, int>>, greater<pair<weight_t, int>>> PQ;
	for(int k = 0; k < inst.n; k++) {
		if(inst.f0[k]) continue;

		weight_t bound = 0;
		for(int ij : inst.dout[k]) {
			if(inst.c[ij] > bound) bound = inst.c[ij];
		}

		dist[k] = 0.0;
		processed[k] = true;
		PQ.push(make_pair(0.0, k));
		marked.push_back(k);
		while(!PQ.empty()) {
			const int i = PQ.top().second;
			processed[i] = true;
			PQ.pop();
			if(dist[i] >= bound || hops[i] > 2) {
				break;
			}

			for(int ij : inst.dout[i]) {
				const int j = inst.head[ij];
				if(processed[j]) continue;

				const weight_t d = dist[i] + inst.c[ij];
				if(d >= bound) continue;
				
				if(d < dist[j]) {
					dist[j] = d;
					hops[j] = hops[i] + 1;
					PQ.push(make_pair(d, j));
					marked.push_back(j);
				} else if(d == dist[j]) {
					dist[j] = d;
					hops[j] = hops[i] + 1;
				}
			}
		}
		//PQ.clear();

		for(int ij : inst.dout[k]) {
			const int j = inst.head[ij];
			
			if ( dist[j] < inst.c[ij] || (dist[j] == inst.c[ij] && hops[j] > 1) ) {
				todel.push_back(ij);
			}
		}
		for(int ij : todel) {
			inst.delArc(ij);
			inst.fe0[ij] = true;
		}
		
		for(int i : marked) {
			dist[i] = WMAX;
			hops[i] = 0;
			processed[i] = false;
		}
		cnt += todel.size();
		marked.clear();
		todel.clear();
	}

	stats.lc += cnt;

	return cnt;
}

