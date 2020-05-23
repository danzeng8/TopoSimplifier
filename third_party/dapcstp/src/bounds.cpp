/**
 * \file   bounds.cpp
 * \brief  procedures for computing lower bounds
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

using namespace std;

#include "bounds.h"
#include "util.h"
#include "options.h"
#include <stdlib.h>
#include <iostream>
//#include "tbb/parallel_for.h"
//#include <tbb/concurrent_priority_queue.h>
//#include <tbb/concurrent_vector.h>

class comparePair {
public:
	bool operator()(const pair<int, int> & u, const pair<int, int> & v) const {
		return u.first > v.first;
	}
};

/**
int getPrio(tbb::concurrent_priority_queue < pair<int, int>, comparePair > & pq) {
	pair<int, int> prioPair;
	pq.try_pop(prioPair);
	pq.push(prioPair);
	return prioPair.first;
}**/



template<typename U> weight_t daR(int r, Inst& inst, vector<U>& c, vector<U>& cr, vector<U>& pi, weight_t ub, double eager, Sol* inc, bool heur)
{
	weight_t lb = inst.offset;
	const int n = inst.n, m = inst.m;

	vector<weight_t>& p = inst.p;
	//tbb::concurrent_vector<weight_t> p(inst.p.begin(), inst.p.end());

	vector<flag_t>& T = inst.T;
	
	int * Q = new int[n];
	int qL1 = 0, qL2 = 0;

	PQMin<int,int> PQ;

	//tbb::concurrent_priority_queue<pair<int, int>, comparePair> PQ;

	bool * cut = new bool[n]; bool * active = new bool[n];

	//tbb::parallel_for(tbb::blocked_range<int>(0, inst.m),
		//[&](tbb::blocked_range<int> r)
	//{
		//for (int i = r.begin(); i < r.end(); ++i)
		//{
		for(int i = 0; i < m; i++) {
			if (inst.fe0[i]) continue;
			cr[i] = c[i];
		}
		//}
	//});
	
	for(int i = 0; i < n; i++) {
		active[i] = false;
		cut[i] = false;
		if(!T[i] || i == r) continue;
		pi[i] = inst.f1[i] ? WMAX : p[i];
		///boost::mutex::scoped_lock lock(mutexLock);
		PQ.push(make_pair(1, i));
		if(inst.f1[i]) active[i] = true;
	}
	
	active[r] = true;
	if(PQ.size() == 0) {
		return lb;
	}
	
	int *L = new int[m];
	int *Ltails = new int[m];
	
	int iter = 0;
	int v = -1;
	bool augmentroot = false;
choose_element:
	//boost::mutex::scoped_lock lock(mutexLock);
	while ( !PQ.empty() ) {
		//pair<int, int> entry;
		//PQ.try_pop(entry);
		//boost::mutex::scoped_lock lock(mutexLock);
		pair<int,int> entry = PQ.top();
		PQ.pop();
		
		int nextpr = (!PQ.empty() ? PQ.top().first : IMAX); //PQ.top().first getPrio(PQ)
		int prio = entry.first;
		v = entry.second;
		//cout << "process entry " << prio << " " << v << endl;

		if(params.lastcomp && PQ.empty() && !augmentroot) {
			augmentroot = true;
			break;
		}

		do {
			int cL = 0, vc = 1;
			int deg = inst.din[v].size();
			for(int j = 0; j < qL1; j++) {
				cut[Q[j]] = false;
			}
			qL1 = 0;
			qL2 = 0;
			Q[qL1++] = v;
			cut[v] = true;

			// identify component
			while ( qL2 != qL1 ) {
				const int w = Q[qL2++];
				for(int ij : inst.din[w]) {
					const int u = inst.tail[ij];
					if ( cr[ij] <= params.dasat ) {
						if(!cut[u]) {
							if ( active[u] ) {
								active[v] = false;
								goto choose_element;
							}

							Q[qL1++] = u;
							cut[u] = true;
							vc++;
							deg += inst.din[u].size();
						}
					} else {
						L[cL] = ij;
						Ltails[cL++] = u;
					}
				}
			}

			int i = 0, j = 0;
			while(i < cL) {
				while(i < cL && cut[Ltails[i]]) {
					i++;
				}
				if(i == cL) break;
				// here !cut[Ltails[i]] holds
				if(j < i) {
					Ltails[j] = Ltails[i];
					L[j] = L[i];
				}
				j++;
				i++;
			}
			cL = j;
			
			// compute component priority
			prio = deg - (vc-1);
			if(params.daguide && inc != nullptr) {
				int nSolArcsInCut = 0;
				for(int i = 0; i < cL; i++) {
					if(inc->arcs[L[i]]) nSolArcsInCut++;
				}
				if(nSolArcsInCut > 1) prio += inst.m*nSolArcsInCut;
			}
			
			iter++;
			
			if(prio > nextpr * eager) {
				break;
			}

			// perform augmentation
			U delta = std::numeric_limits<U>::max();
			for(int i = 0; i < cL; i++) {
				if ( cr[L[i]] < delta ) delta = cr[L[i]];
			}
			delta = min(pi[v], delta);

			for(int i = 0; i < cL; i++)
				cr[L[i]] -= delta;
			pi[v] -= delta;
			lb += delta;

			int cnt = 0;
			for(int i = 0; i < cL; i++) {
				if(cr[L[i]] <= params.dasat) {
					deg += inst.din[Ltails[i]].size();
					vc++;
					if(inc != nullptr && params.daguide && inc->arcs[L[i]] == 1) cnt++;
				}
			}

			if ( pi[v] <= 0 || lb >= ub) {
				break;
			}

			prio = (deg - (vc-1));
			if(cnt > 1) prio += cnt * inst.m;

		} while ( prio <= nextpr );
		
		if (pi[v] != 0) {
			//boost::mutex::scoped_lock lock(mutexLock);
			PQ.push(make_pair(prio, v));
		}
		if (lb >= ub) {
			break;
		}
	}
	 
	// compute last component cuts via shortest path
	if(params.lastcomp && augmentroot && lb < ub) {
		
		// dijkstra
		PQMin<U,int> PQ2;
		vector<U> dist(inst.n, std::numeric_limits<U>::max());
		dist[v] = 0;
		PQ2.push(make_pair(0, v));
		while ( !PQ2.empty() ) {
			const int i = PQ2.top().second;
			if(i == r || dist[i] > pi[v]) break;
			PQ2.pop();

			for(int ij : inst.din[i]) {
				const int j = inst.tail[ij];
				const U d = dist[i] + cr[ij];

				if (d < dist[j]) {
					dist[j] = d;
					PQ2.push(make_pair(d, j));
				}
			}
		}

		if(!inst.f1[v] && dist[r] > pi[v]) {
			dist[r] = pi[v];
		}

		// in case of prune-by-bound we do not need the reduced costs
		if(lb + dist[r] < ub) {

			// compute reduced costs
			for(int i = 0; i < n; i++) {
				if(inst.f0[i]) continue;
				if(dist[i] > dist[r]) dist[i] = dist[r];
			}
			
			for(int i = 0; i < n; i++) {
				if(inst.f0[i]) continue;

				for(int ij : inst.din[i]) {
					const int j = inst.tail[ij];
					if(inst.f0[j]) continue;
					
					const U delta = dist[j] - dist[i];
					if(delta <= 0)
						continue;

					cr[ij] -= delta;
				}
			}

			// only update pi for potential terminals
			if(!inst.f1[v]) {
				pi[v] -= dist[r];
			}
		}

		lb += dist[r];
	}

	// for big M instance, re-adjust bigM 
	if(inst.bigM >= 0 && !heur) {
		// if potential ran out for nodes, none of the bigM arcs may be saturated
		weight_t minSat = WMAX;
		for(int rj : inst.dout[inst.r]) {
			if(cr[rj] < minSat) {
				minSat = cr[rj];
			}
		}
		
		// reduce bigM by minsat -> at least one artificial root arc is saturated
		if(minSat > params.dasat) {
			weight_t delta = minSat;
			inst.bigM -= delta;
			for(int rj : inst.dout[inst.r]) {
				inst.c[rj] -= delta;
				cr[rj] -= delta;
			}
			inst.offset += delta;
			lb += delta;
		}
	}

	delete[] L;
	delete[] Ltails;

	return lb;
}

template weight_t daR<weight_t>(int r, Inst& inst, vector<weight_t>& c, vector<weight_t>& cr, vector<weight_t>& pi, weight_t ub, double rel, Sol* inc, bool heur);
template weight_t daR<double>(int r, Inst& inst, vector<double>& c, vector<double>& cr, vector<double>& pi, weight_t ub, double rel, Sol* inc, bool heur);

