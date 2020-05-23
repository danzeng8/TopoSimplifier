/**
 * \file   inst.cpp
 * \brief  problem instance class
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#include "inst.h"
#include "util.h"

#include <stack>

Inst::Inst()
{

}

Inst::~Inst()
{

}

Inst::Inst(const Inst& src)
{
	din = src.din;
	dout = src.dout;

	tail = src.tail;
	head = src.head;

	opposite = src.opposite;

	pin = src.pin;
	pout = src.pout;

	T = src.T;
	f0 = src.f0;
	f1 = src.f1;
	fe0 = src.fe0;

	p = src.p;
	c = src.c;

	offset = src.offset;
	n = src.n;
	m = src.m;
	t = src.t;
	r = src.r;

	bmna = src.bmna;
	bmaa = src.bmaa;

	isInt = src.isInt;
	isAsym = src.isAsym;
	isMWCS = src.isMWCS;
	inst1 = src.inst1;
	bigM = src.bigM;

	transformation = src.transformation;
}

void Inst::newArc(int i, int j, int ij, int ji, weight_t w)
{
	pout[ij] = (int)dout[i].size();
	pin[ij] = (int)din[j].size();
	dout[i].push_back(ij);
	din[j].push_back(ij);
	tail[ij] = i;
	head[ij] = j;
	opposite[ij] = ji;
	
	c[ij] = (weight_t)w;
}

void Inst::delArc(int ij)
{
	int i = tail[ij];
	int j = head[ij];
	int deg = din[j].size();
	int p1 = pin[ij];
	
	const int ji = opposite[ij];
	if(ji != -1) {
		opposite[ji] = -1;
	}
	opposite[ij] = -1;

	int k;
	for(k = 0; k < din[j].size(); k++) {
		if(din[j][k] == ij) break;
	}
	assert(k < din[j].size());

	if(p1 == deg-1) {
		din[j].pop_back();
	} else {
		assert(k == p1);
		int ij2 = din[j].back();

		din[j].pop_back();
		din[j][p1] = ij2;
		pin[ij2] = p1;
		pin[ij] = -1;
	}

	for(k = 0; k < dout[i].size(); k++) {
		if(dout[i][k] == ij) break;
	}
	assert(k < dout[i].size());
	
	deg = dout[i].size();
	p1 = pout[ij];
	if(p1 == deg-1) {
		dout[i].pop_back();
	} else {
		assert(k == p1);
		int ij2 = dout[i].back();

		dout[i].pop_back();
		dout[i][p1] = ij2;
		pout[ij2] = p1;
		pout[ij] = -1;
	}

	if(bmaa.size() > 0)
		bmaa[ij].clear();
}

void Inst::moveHead(int ij, int k)
{
	int i = tail[ij];
	int j = head[ij];
	int p1 = pin[ij];
	int deg = din[j].size();

	// remove ij from original adj list
	if(p1 == deg-1) {
		din[j].pop_back();
	} else {
		// take the last and replace ij with it
		int ij2 = din[j].back();
		din[j].pop_back();

		din[j][p1] = ij2;
		pin[ij2] = p1;
	}
	
	assert(din[j].size() < deg);

	// place ij in new adj list
	head[ij] = k;
	pin[ij] = din[k].size();
	din[k].push_back(ij);
}

void Inst::moveTail(int ij, int k)
{
	int i = tail[ij];
	int j = head[ij];
	int p1 = pout[ij];
	int deg = dout[i].size();

	// remove ij from original adj list
	if(p1 == deg-1) {
		assert(dout[i][p1] == ij);
		dout[i].pop_back();
	} else {
		// take the last and replace ij with it
		assert(dout[i].size() != 0);
		int ij2 = dout[i].back();
		
		dout[i].pop_back();

		dout[i][p1] = ij2;
		pout[ij2] = p1;
	}

	// place ij in new adj list
	tail[ij] = k;
	pout[ij] = dout[k].size();
	dout[k].push_back(ij);
}

void Inst::merge(int ij, int i, int j)
{
	vector<int> todel, tomove_head, tomove_tail;
	todel.reserve(m);
	tomove_head.reserve(m);
	tomove_tail.reserve(m);

	for(int ij2 : din[j]) {
		if(tail[ij2] == i)
			todel.push_back(ij2);
		else
			tomove_head.push_back(ij2);
	}
	for(int ij2 : dout[j]) {
		if(head[ij2] == i)
			todel.push_back(ij2);
		else
			tomove_tail.push_back(ij2);
	}
	for(int ij2 : todel) {
		delArc(ij2);
		fe0[ij2] = true;
	}
	for(int ij2 : tomove_head) {
		moveHead(ij2, i);
	}
	for(int ij2 : tomove_tail) {
		moveTail(ij2, i);
	}

	assert(din[j].size() == 0);
	assert(dout[j].size() == 0);
	din[j].shrink_to_fit();
	dout[j].shrink_to_fit();

	// seek cheapest incoming/outgoing arcs to neighbors and keep them
	vector<int> ndistin(n, -1);
	for(int ij2 : din[i]) {
		int j2 = tail[ij2];
		if(ndistin[j2] == -1 || c[ij2] < c[ndistin[j2]]) {
			ndistin[j2] = ij2;
		}
	}

	vector<int> ndistout(n, -1);
	for(int ij2 : dout[i]) {
		int j2 = head[ij2];
		if(ndistout[j2] == -1 || c[ij2] < c[ndistout[j2]]) {
			ndistout[j2] = ij2;
		}
	}

	todel.clear();
	for(int ij2 : din[i]) {
		if(ndistin[tail[ij2]] != ij2)
			todel.push_back(ij2);
	}

	for(int ij2 : dout[i]) {
		if(ndistout[head[ij2]] != ij2)
			todel.push_back(ij2);
	}

	for(int ij2 : todel) {
		delArc(ij2);
		fe0[ij2] = true;
	}

	offset += c[ij];

	// overflow
	if(p[i] == WMAX || p[j] == WMAX) {
		p[i] = WMAX;
	} else {
		p[i] = p[i] + p[j] - c[ij];
	}
	
	if(p[i] > 0)
		T[i] = true;
	if(f1[j]) {
		f1[i] = true;
		T[i] = true;
		p[i] = WMAX;
	}

	// eliminate j
	T[j] = false;
	p[j] = 0;
	f0[j] = true;
	f1[j] = false;
}

void Inst::contractArc(int ji)
{
	const int j = tail[ji];
	const int i = head[ji];

	T[j] = true;
	T[i] = false;

	if(f1[i]) {
		f1[i] = false;
		f1[j] = true;
	}

	offset += c[ji];

	// overflow check
	if(p[i] == WMAX || p[j] == WMAX) {
		p[j] = WMAX;
	} else {
		p[j] = p[i] + p[j] - c[ji];
	}
	p[i] = 0;

	updatebmNTD1(ji);

	removeNode(i);
}

void Inst::removeNode(int i)
{
	assert(!f0[i]);
	assert(!f1[i]);

	vector<int> todel;
	todel.reserve(din[i].size()+dout[i].size());
	for(int ij : din[i])  todel.push_back(ij);
	for(int ij : dout[i]) todel.push_back(ij);

	for(int ij : todel) {
		delArc(ij);
		fe0[ij] = true;
	}

	T[i] = false;
	f0[i] = true;
	offset += p[i];
	
	p[i] = 0;
	
	if(bmna.size() > 0)
		bmna[i].clear();
}

void Inst::increaseRevenue(int i, weight_t val)
{
	for(int ji : din[i]) {
		c[ji] += val;
	}
	if(p[i] > WMAX - val) {
		p[i] = WMAX;
	} else {
		p[i] += val;
	}
	offset -= val;
	T[i] = true;
}

weight_t Inst::decreaseRevenue(int i)
{
	weight_t delta = WMAX;
	for(int ji : din[i]) {
		if(c[ji] < delta) {
			delta = c[ji];
		}
	}
	delta = min(delta, p[i]);

	for(int ji : din[i]) {
		c[ji] -= delta;
	}
	p[i] -= delta;
	offset += delta;

	if(f1[i]) {
		T[i] = true;
		p[i] = WMAX;
	}

	if(p[i] == 0) {
		if(!f1[i])
			T[i] = false;
	} else {
		T[i] = true;
	}
	return delta;
}

InstSizeData Inst::countInstSize()
{
	InstSizeData sdata;
	sdata.n = 0;
	sdata.m = 0;
	sdata.t = 0;
	sdata.tr = 0;
	sdata.f1 = 0;

	for(int ij = 0; ij < m; ij++) {
		if(!fe0[ij]) {
			sdata.m++;
		}
	}
	for(int i = 0; i < n; i++) {
		if(!f0[i]) {
			sdata.n++;
			if(f1[i]) {
				sdata.f1++;
			} else {
				if(p[i] > 0) {
					sdata.t++;
				}

				bool realT = false;
				for(int ij : din[i]) {
					if(c[ij] < p[i]) {
						realT = true;
						break;
					}
				}
				if(realT) {
					sdata.tr++;
				}
			}
		}
	}
	return sdata;
}

weight_t Inst::nonreachableRevenue(int start)
{
	stack<int> Q;
	vector<flag_t> visited(n, false);
	Q.push(start);
	visited[start] = true;
	weight_t P = 0;
	for(int i = 0; i < n; i++) {
		if(f0[i] || f1[i] || i == start) continue;
		P += p[i];
	}
	while(!Q.empty()) {
		int i = Q.top();
		Q.pop();
		for(int ij : dout[i]) {
			int j = head[ij];
			if(!visited[j]) {
				if(!f1[j])
					P -= p[j];
				visited[j] = true;
				Q.push(j);
			}
		}
	}
	
	return P;
}

void Inst::printNodeSet(set<int>& s)
{
	printf("nodes = { ");
	for(int i : s)
		printf("p(%d)=%ld ", i, p[i]);
	printf("}\n");
}

void Inst::printArcSet(vector<int>& s)
{
	printf("arcs = { ");
	for(int a : s)
		printf("c(%d %d)=%ld ", tail[a], head[a], c[a]);
	printf("}  ");
}

bool Inst::hasCheaperIncomingArc(int i, int ij)
{
	for(int ij2 : din[i]) {
		if(c[ij2] < c[ij] && ij2 != ij)
			return true;
	}
	return false;
}

bool Inst::singleAdjacency(int i)
{
	return ((din[i].size() == 1 && dout[i].size() <= 1) && 
			(dout[i].size() == 0 || 
			 tail[din[i].front()] == head[dout[i].front()]));
}

bool Inst::doubleAdjacency(int i, int& j, int& k, int& ki, int& ij, int& ji, int& ik)
{
	const int degin = din[i].size();
	const int degout = dout[i].size();

	set<int> neighbors;
	ij = -1; ki = -1; ji = -1; ik = -1, j = -1, k = -1;

	if(degin > 2 || degout > 2) return false;

	for(int a : din[i])  neighbors.insert(tail[a]);
	for(int a : dout[i]) neighbors.insert(head[a]);

	if(neighbors.size() != 2) return false;

	list<int> nlist(neighbors.begin(), neighbors.end());

	j = nlist.front();
	k = nlist.back();

	for(int a : din[i]) {
		if(tail[a] == j) ji = a;
		if(tail[a] == k) ki = a;
	}
	for(int a : dout[i]) {
		if(head[a] == j) ij = a;
		if(head[a] == k) ik = a;
	}

	return true;
}

void Inst::APsearch(int i, vector<flag_t>& ap, vector<vector<int>>& certs)
{
	stack<int> Q; vector<int> lastk(n, -1);
	vector<flag_t> visited(n, false);
	vector<int> parent(n, -1);
	vector<int> disc(n, 0), low(n, 0);
	int depth = 0;

	Q.push(i);

	while(!Q.empty()) {
		const int i = Q.top();
		const int nIncidentArcs = dout[i].size() + din[i].size();

		// post-order action
		if(lastk[i] != -1) {
			int j, k = lastk[i];
			if(k < dout[i].size())
				j = head[dout[i][k]];
			else
				j = tail[din[i][k-dout[i].size()]];

			low[i] = min(low[i], low[j]);
			if(parent[i] != -1 && low[j] >= disc[i]) {
				certs[i].push_back(j);
				ap[i] = true;
			}
		} else {
			visited[i] = true;
			disc[i] = low[i] = ++depth;
		}

		int k = lastk[i]+1;
		while(k < nIncidentArcs) {
			int j;
			if(k < dout[i].size())
				j = head[dout[i][k]];
			else
				j = tail[din[i][k-dout[i].size()]];

			if(!visited[j]) {
				// pre-order action
				parent[j] = i;
				Q.push(j);
				break;
			} else if(j != parent[i]) {
				low[i] = min(low[i], disc[j]);
			}
			k++;
		}

		lastk[i] = k;
		if(k == nIncidentArcs) {
			Q.pop();
		}
	}
}

vector<flag_t> Inst::AP(vector<flag_t>& ap, vector<int>& lastap)
{
	vector<int> disc(n, 0);
	vector<int> low(n, 0);
	vector<int> parent(n, -1);
	vector<flag_t> visited(n, false);
	vector<vector<int>> certs(n);

	ap = vector<flag_t>(n, false);
	lastap = vector<int>(n, -1);
	if(dout[r].size() == 0) return ap;

	APsearch(r, ap, certs);

	ap[r] = true;
	lastap[r] = r;
	for(int ij : dout[r]) {
		int j = head[ij];
		certs[r].push_back(j);
	}

	findAllSubtrees(ap, lastap, certs);

	return ap;
}

void Inst::findAllSubtrees(vector<flag_t>& ap, vector<int>& lastap, vector<vector<int>>& certs)
{
	vector<flag_t> visited(n, false);
	vector<flag_t> processed(n, false);
	stack<int> Q;

	Q.push(r);

	while(!Q.empty()) {
		const int i = Q.top();

		if(!ap[i] || processed[i])
			Q.pop();

		if(!visited[i] || processed[i]) {
			visited[i] = true;

			if(ap[i] && !processed[i])  {
				// pre-order for ap
				for(int j : certs[i]) {
					lastap[j] = i;
					Q.push(j);
				}
				processed[i] = true;
				continue;
			}
			// post-order
			for(int ij : dout[i]) {
				int j = head[ij];
				if(!visited[j]) {
					lastap[j] = lastap[i];
					Q.push(j);
				}
			}
			for(int ij : din[i]) {
				int j = tail[ij];
				if(!visited[j]) {
					lastap[j] = lastap[i];
					Q.push(j);
				}
			}
		}
	}
}

void Inst::updatebmNTD1(int ji)
{
	const int j = tail[ji];
	const int i = head[ji];

	for(int a : bmaa[ji]) { bmna[j].push_back(a); }
	for(int a : bmna[i])  { bmna[j].push_back(a); }
}


void Inst::updatebmNTD2triangle(int ik, int ij, int jk)
{
	bmaa[ik].clear();

	const int i = head[ij];

	// add all implied arcs from ij, j, jk to ik
	for(int a : bmaa[ij]) { bmaa[ik].push_back(a); }
	for(int a : bmaa[jk]) { bmaa[ik].push_back(a); }
	for(int a : bmna[i])  { bmaa[ik].push_back(a); }
}

void Inst::updatebmNTD2(int ij, int jk)
{
	const int i = head[ij];

	// add all implied arcs from j, jk to ij
	for(int a : bmaa[jk]) { bmaa[ij].push_back(a); }
	for(int a : bmna[i])  { bmaa[ij].push_back(a); }
}

void Inst::updatebmMerge(int ij, bool bidirect)
{
	const int i = tail[ij];
	const int j = head[ij];

	// add all implied arcs from j to i
	for(int a : bmna[j])  { bmna[i].push_back(a); }
	for(int a : bmaa[ij]) { bmna[i].push_back(a); }

	if(bidirect) {
		const int ji = opposite[ij];
		if(ji != -1) {

			// add all arcs from ij to i
			for(int a : bmaa[ji]) { bmna[i].push_back(a); }
			bmaa[ji].clear();
		}
	}

	bmaa[ij].clear();
	bmna[j].clear();
}

void Inst::convertMWCS2PCSTP()
{
	// transform to APCSTP
	double minP = WMAX, sumP = 0.0;
	for (int i = 0; i < n; i++) {
		if(p[i] < minP) minP = p[i];
	}
	for (int i = 0; i < m; i++) {
		c[i] -= minP;
	}
	for (int i = 0; i < n; i++) {
		p[i] -= minP;
		sumP += p[i];
	}

	transformation = new Transformation();

	transformation->P = sumP;
	transformation->minP = minP;
}

weight_t Inst::convertPCSTPBound2MWCS(weight_t bound)
{
	return transformation->P - bound + transformation->minP;
}

Inst Inst::createRootedBigMCopy()
{
	Inst inst;
	inst.offset = offset;
	inst.t = t;

	// add artifical root and arcs
	inst.resizeNodes(n+1);
	inst.resizeEdges(m+n);

	inst.isInt = isInt;
	inst.isAsym = isAsym;
	inst.isMWCS = isMWCS;
	inst.inst1 = inst1;

	inst.transformation = transformation;

	for(int ij = 0; ij < m; ij++) {
		const int i = tail[ij];
		const int j = head[ij];
		inst.newArc(i, j, ij, opposite[ij], c[ij]);
	}

	// last node is now artificial root, start artificial arc indexing at end of original arcs
	inst.r = n;
	int ij = m;
	weight_t M = 1;

	for(int i= 0; i < n; i++) {
		inst.p[i] = p[i];
		inst.f0[i] = f0[i];
		inst.f1[i] = f1[i];
		inst.T[i] = T[i];

		M += inst.p[i];
	}

	for(int i = 0; i < n; i++) {
		inst.newArc(inst.r, i, ij, -1, M);
		ij++;
	}

	inst.T[inst.r] = true;
	inst.f1[inst.r] = true;
	inst.offset -= M;
	inst.bigM = M;
	return inst;
}

void Inst::setBigM(weight_t M)
{
	offset += bigM;
	bigM = M;

	for(int rj : dout[r]) {
		c[rj] = bigM;
	}

	offset -= bigM;
}

void Inst::resizeNodes(int _n)
{
	n = _n;
	T.resize(_n, false);
	din.resize(_n);
	dout.resize(_n);
	f0.resize(_n, false);
	f1.resize(_n, false);
	p.resize(_n, 0.0);
}

std::vector<int> Inst::getDinIndex(int index) {
	return din[index];
}

void Inst::resizeEdges(int _m)
{
	m = _m;
	pin.resize(_m);
	pout.resize(_m);
	opposite.resize(_m);
	tail.resize(_m);
	head.resize(_m);
	fe0.resize(_m, false);
	c.resize(_m);
}

int Inst::countReachableFixed(int r)
{
	vector<int> visited(n, false); stack<int> Q;
	Q.push(r);
	visited[r] = true;
	int fixed = 0;
	if(f1[r])
		fixed++;
	
	while(!Q.empty()) {
		const int i = Q.top();
		Q.pop();
		for(int ij : dout[i]) {
			const int j = head[ij];
			if(!visited[j]) {
				if(f1[j])
					fixed++;
				visited[j] = true;
				Q.push(j);
			}
		}
	}
	if(!visited[this->r]) {
		fixed++;
	}
	return fixed;
}
