/**
 * \file   util.cpp
 * \brief  loaders
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

using namespace std;

#include "util.h"
#include "bbtree.h"

#include "timer.h"
#include "options.h"
#include "sol.h"
#include "ds.h"

#include <stdio.h>
#include <boost/filesystem.hpp>
//#include <sys/resource.h>
#include <boost/pending/disjoint_sets.hpp>

#include "stats.h"

Inst load(const char* fn)
{
	Inst inst;
	// all instances are loaded in their APCSTP representation
	if(params.type.compare("nwstp") == 0 || params.type.compare("stp") == 0) {
		inst = loadNWSTP(fn);
	} else if(params.type.compare("mwcs") == 0) {
		inst = loadMWCS(fn);
	} else if(params.type.compare("pcstp") == 0){
		inst = loadPCSTP(fn);
	} else {
		EXIT("error: specified problem type unknown: %s\n", params.type.c_str());
	}

	inst.t = 0;
	for(int i = 0; i < inst.n; i++) {
		if(inst.p[i] > 0 || inst.f1[i]) inst.t++;
		inst.T[i] = (inst.p[i] != 0);
	}

	stats.initial = inst.countInstSize();

	// associates the anti-parallel arc to each arc if it exists (-1 otherwise)
	int antiparallelArcs = 0;
	if(inst.isAsym) {
		for(int ij = 0; ij < inst.m; ij++)
			inst.opposite[ij] = -1;
		for(int ij = 0; ij < inst.m; ij++) {
			if(inst.opposite[ij] != -1) continue;
			const int i = inst.tail[ij];
			const int j = inst.head[ij];
			for(int jk : inst.dout[j]) {
				const int k = inst.head[jk];
				if(k == i) {
					inst.opposite[ij] = jk;
					inst.opposite[jk] = ij;
					break;
				}
			}
		}
		for(int ij = 0; ij < inst.m; ij++) {
			if(inst.opposite[ij] == -1) continue;
			antiparallelArcs++;
		}
	} else {
		antiparallelArcs = inst.m;
	}

	// compute the ratio of bidirected edges / arcs that have an antiparallel counterpart
	stats.bidirect = (double)antiparallelArcs/(2*inst.m-antiparallelArcs);

	return inst;
}

Inst loadMWCS(const char* fn)
{
	FILE *fp;
	char buf[256];
	int n, m, t, v1, v2, r;
	double w, prize;
	int ij = 0;

	if((fp=fopen(fn, "r")) == NULL) {
		EXIT("error: file not found: %s\n", fn);
	}

	Inst inst;
	inst.offset = 0.0;
	inst.r = -1;
	inst.isInt = true;
	inst.isMWCS = true;
	while(fgets(buf, 256, fp) != NULL) {

		if(sscanf(buf, "Nodes %d", &n) == 1) {
			inst.resizeNodes(n);
		}
		if(sscanf(buf, "Edges %d", &m) == 1) {
			m*=2;
			inst.resizeEdges(m);
		}
		if(sscanf(buf, "Arcs %d", &m) == 1) {
			inst.isAsym = true;
			inst.resizeEdges(m);
		}

		if(sscanf(buf, "A %d %d", &v1, &v2) == 2) {
			int i = v1-1, j = v2-1;

			inst.newArc(i, j, ij, -1, 0.0);
			ij++;
		}

		if(sscanf(buf, "E %d %d", &v1, &v2) == 2) {
			int i = v1-1, j = v2-1;

			inst.newArc(i, j, ij, ij+1, 0.0); ij++;
			inst.newArc(j, i, ij, ij-1, 0.0); ij++;
		}

		if(sscanf(buf, "T %d %lf", &v1, &prize) == 2) {
			int i = v1-1;

			inst.isInt &= (floor(prize) == prize);
			
			if(!inst.isInt) {
				prize = floor(prize * params.precision);
			}

			inst.T[i] = true;
			inst.p[i] = (weight_t)prize;
		}
	}

	assert(ij == m);

	fclose(fp);

	inst.convertMWCS2PCSTP();

	return inst;
}

Inst loadNWSTP(const char* fn)
{
	FILE *fp;
	char buf[256];
	int n, m, t, v1, v2, r;
	double w, prize;
	int ij = 0;

	if((fp=fopen(fn, "r")) == NULL) {
		EXIT("error: file not found: %s\n", fn);
	}

	Inst inst;
	inst.offset = 0;
	inst.r = -1;
	inst.isInt = true;

	bool readingT = false;
	bool finishedT = false;
	int nwcounter = 0;
	vector<double> nw;
	while(fgets(buf, 256, fp) != NULL) {

		if(sscanf(buf, "Fixed %lf", &w) == 1) {
			inst.offset = w;
		}

		if(sscanf(buf, "Nodes %d", &n) == 1) {
			inst.n = n;
			inst.resizeNodes(n);
			nw.resize(n, 0.0);
		}

		if(sscanf(buf, "Edges %d", &m) == 1) {
			m*=2;
			inst.resizeEdges(m);
		}

		if(sscanf(buf, "Terminals %d", &t) == 1) {
			inst.t = t;
		}

		if(sscanf(buf, "E %d %d %lf", &v1, &v2, &w) == 3) {

			int i = v1-1, j = v2-1;

			inst.isInt &= (floor(w) == w);

			if(!inst.isInt) {
				w = floor(w * params.precision);
			}

			inst.newArc(i, j, ij, ij+1, (weight_t)w); ij++;
			inst.newArc(j, i, ij, ij-1, (weight_t)w); ij++;
		}

		if(!finishedT && sscanf(buf, "T %d", &v1) == 1) {
			readingT = true;
			int i = v1-1;

			//inst.isInt &= (floor(prize) == prize);
			
			inst.T[i] = true;
			inst.f1[i] = 1;
			inst.p[i] = WMAX;
		}

		if(sscanf(buf, "NW %lf", &w) == 1) {
			int i = nwcounter++;

			if(!inst.isInt) {
				w = floor(w * params.precision);
			}

			nw[i] = w;
			for(int ij : inst.din[i]) {
				inst.c[ij] = inst.c[ij] + w;
			}
		}

		if(readingT && strncmp(buf, "END", 3) == 0) {
			finishedT = true;
		}
	}
	
	for (int i = 0; i < inst.n; i++) {
		if(inst.T[i]) {
			inst.r = i;
			break;
		}
	}
	
	inst.offset += nw[inst.r];

	for (int i = 0; i < inst.n; i++) {
		inst.din[i].shrink_to_fit();
		inst.dout[i].shrink_to_fit();
	}

	assert(ij == m);

	fclose(fp);

	return inst;
}

Inst loadPCSTP(const char* fn)
{
	FILE *fp;
	char buf[256];
	int n, m, t, v1, v2, r;
	double w, prize;
	int ij = 0;

	if((fp=fopen(fn, "r")) == NULL) {
		EXIT("error: file not found: %s\n", fn);
	}

	Inst inst;
	inst.offset = 0.0;
	inst.r = -1;
	inst.isInt = true;

	vector<double> tmpW, tmpP;

	while(fgets(buf, 256, fp) != NULL) {

		if(sscanf(buf, "Nodes %d", &n) == 1) {
			inst.resizeNodes(n);
			tmpP.resize(n);
		}

		if(sscanf(buf, "Edges %d", &m) == 1) {
			m*=2;
			inst.resizeEdges(m);
			tmpW.resize(m);
		}

		if(sscanf(buf, "Arcs %d", &m) == 1) {
			inst.isAsym = true;
			inst.resizeEdges(m);
			tmpW.resize(m);
		}

		if(sscanf(buf, "Terminals %d", &t) == 1) {
			inst.t = t;
		}

		if(sscanf(buf, "E %d %d %lf", &v1, &v2, &w) == 3) {

			int i = v1-1, j = v2-1;

			inst.isInt &= (floor(w) == w);

			inst.newArc(i, j, ij, ij+1, (weight_t)w);
			tmpW[ij] = w;
			ij++;
			inst.newArc(j, i, ij, ij-1, (weight_t)w);
			tmpW[ij] = w;
			ij++;
		}

		if(sscanf(buf, "A %d %d %lf", &v1, &v2, &w) == 3) {

			int i = v1-1, j = v2-1;

			inst.isInt &= (floor(w) == w);

			inst.newArc(i, j, ij, -1, w);
			tmpW[ij] = w;
			ij++;
		}

		if(sscanf(buf, "TP %d %lf", &v1, &prize) == 2) {
			int i = v1-1;

			inst.isInt &= (floor(prize) == prize);

			inst.T[i] = true;
			tmpP[i] = prize;
			
			inst.p[i] = (weight_t)prize;
		}

		if(sscanf(buf, "RootP %d", &v1) == 1) {
			inst.r = v1-1;
			inst.T[inst.r] = true;
		}
	}

	if(!inst.isInt) {
		for (int i = 0; i < inst.n; i++) {
			inst.p[i] = round(tmpP[i] * params.precision);
		}
		for (int i = 0; i < inst.m; i++) {
			inst.c[i] = round(tmpW[i] * params.precision);
		}
	}

	if(inst.r != -1) {
		inst.f1[inst.r] = true;
		inst.p[inst.r] = WMAX;
	}
	
	assert(ij == m);

	fclose(fp);

	return inst;
}

double getBestKnownBound(const char* fn, const char* boundfile)
{
	FILE *fp;
	char buf[1024];
	string name = boost::filesystem::path(params.file).stem().string();
	if((fp=fopen(boundfile, "r")) == NULL) {
		EXIT("error: file not found: %s\n", fn);
	}

	double dBest = -1;
	while(fgets(buf, 1024, fp) != NULL) {
		if(strncmp(buf, name.c_str(), name.length()) == 0) {
			int len = strlen(buf);
			char* sBest = strchr(buf, ' ')+1;
			dBest = strtod(sBest, NULL);
		}
	}

	fclose(fp);

	return dBest;
}

void writeSolution(const char* file, Inst& inst, Sol& sol)
{
	weight_t obj = inst.offset;
	for(int i = 0; i < inst.n; i++) {
		if(!sol.nodes[i]) obj += inst.p[i];
	}
	for(int i = 0; i < inst.m; i++) {
		if(sol.arcs[i]) obj += inst.c[i];
	}

	FILE* fp;
	if((fp=fopen(file, "w")) == NULL) {
		EXIT("error: writing solution file: %s\n", file);
	}

	int nVertices = 0, nEdges = 0;
	for(int i = 0; i < inst.n; i++)    if(sol.nodes[i]) nVertices++;
	for(int ij = 0; ij < inst.m; ij++) if(sol.arcs[ij]) nEdges++;

	fprintf(fp, "SECTION Comment\n");
	fprintf(fp, "Name %s\n", stats.name.c_str());
	fprintf(fp, "Program %s\n", PROGRAM_NAME);
	fprintf(fp, "Version %s\n", PROGRAM_VERSION);
	fprintf(fp, "END\n\n");

	fprintf(fp, "SECTION Solutions\n");
	fprintf(fp, "Solution %.6lf %0.3lf\n", format(sol.obj, inst), Timer::total.elapsed().getSeconds());
	fprintf(fp, "END\n\n");

	fprintf(fp, "SECTION BestSolution\n");
	fprintf(fp, "Vertices %d\n", nVertices);
	for(int i = 0; i < inst.n; i++) {
		if(sol.nodes[i]) fprintf(fp, "V %d\n", i+1);
	}
	fprintf(fp, "Edges %d\n", nEdges);
	for(int ij = 0; ij < inst.m; ij++) {
		if(sol.arcs[ij]) fprintf(fp, "E %d %d\n", inst.tail[ij]+1, inst.head[ij]+1);
	}
	fprintf(fp, "END\n\n");

	fclose(fp);
}

double format(weight_t bound, Inst& inst)
{
	if(inst.isMWCS)
		bound = inst.convertPCSTPBound2MWCS(bound);
	if(!inst.isInt) {
		return (double)bound / params.precision;
	} else {
		return bound;
	}
}

void printBoundPadded(Inst& inst, weight_t bound) {
	if(inst.isInt)
		printf("%15.0lf", (double)format(bound, inst));
	else
		printf("%15.6lf", (double)format(bound, inst));
}

void printBound(Inst& inst, weight_t bound) {
	if(inst.isInt)
		printf("%.0lf", (double)format(bound, inst));
	else
		printf("%.6lf", (double)format(bound, inst));
}

void recoverPartialSol(Sol& solP, Inst& inst1)
{
	vector<int> rnmap, ramap;
	Inst instP = genInst(solP, inst1, ramap, rnmap);

	weight_t storedObj = solP.obj;
	
	BBTree bb(instP);
	bb.setRecover(true);
	bb.setOutput(false);
	bb.solve();

	Sol recoveredSolution = bb.getSol();

	for(int i = 0; i < inst1.n; i++) solP.nodes[i] = false;
	for(int i = 0; i < inst1.m; i++) solP.arcs[i] = false;

	for(int i = 0; i < instP.n; i++) {
		if(!recoveredSolution.nodes[i]) continue;
		solP.nodes[rnmap[i]] = true;
	}

	for(int ij = 0; ij < instP.m; ij++) {
		if(!recoveredSolution.arcs[ij]) continue;
		solP.arcs[ramap[ij]] = true;
	}

	solP.obj = recoveredSolution.obj;
	weight_t obj = solP.recomputeObjective();
	
	assert(obj==solP.obj);
	assert(solP.r == rnmap[recoveredSolution.r]);
	
	if(storedObj < obj) {
		EXIT("%sbackmap error: computed %8.2lf but stored %8.2lf%s\n", RED, format(obj, inst1), format(storedObj, inst1), NORMAL);
	}
}

// converts the backmapping from a solution on a reduced instance into a partial solution on
// the original instance
Sol genPartialSol(Sol& sol, Inst& inst)
{
	Inst& inst1 = *inst.inst1;
	Sol sol1(inst1);
	sol1.obj = sol.obj;
	sol1.r = -1;
	sol1.partial = true;

	int m = 0;
	set<int> roots;
	for(int i = 0; i < inst.n; i++) {
		if(!sol.nodes[i]) continue;
		for(int a : inst.bmna[i]) {
			sol1.nodes[inst1.tail[a]] = 1;
			sol1.nodes[inst1.head[a]] = 1;
			sol1.arcs[a] = 1;
			m++;
			if(i == sol.r) {
				roots.insert(inst1.tail[a]);
				roots.insert(inst1.head[a]);
			}
		}
	}
	for(int ij = 0; ij < inst.m; ij++) {
		if(!sol.arcs[ij]) continue;
		for(int a : inst.bmaa[ij]) {
			sol1.nodes[inst1.tail[a]] = 1;
			sol1.nodes[inst1.head[a]] = 1;
			sol1.arcs[a] = 1;
			m++;
		}
	}
	if(m == 0) {
		for(int i = 0; i < inst.n; i++) {
			if(!sol.nodes[i]) continue;
			sol1.nodes[i] = sol.nodes[i];
		}
	}
	if(roots.empty()) {
		roots.insert(sol.r);
	}

	int n = 0;
	for(int i = 0; i < inst1.n; i++) {
		if(!sol1.nodes[i]) continue;
		n++;
	}
	
	// if root not uniquely specified, infer it by finding a root form which all nodes are reachable
	if(inst1.r != -1) {
		for(int i : roots) {
			if(i == inst1.r) {
				sol1.r = i;
				break;
			}
		}
	} else {
		for(int i : roots) {
			int cnt = cntReachable(i, sol1, inst1);
			
			if(cnt == n) {
				sol1.r = i;
				break;
			}
		}
	}

	assert(sol1.r != -1);

	return sol1;
}

Inst genInst(Sol& sol, Inst& inst1, vector<int>& ramap, vector<int>& rnmap)
{
	Inst inst0;
	int m = 0; int n = 0;

	for(int i = 0; i < inst1.n; i++)    if(sol.nodes[i]) n++;
	for(int ij = 0; ij < inst1.m; ij++) if(sol.arcs[ij]) m++;

	inst0.din.resize(n);
	inst0.dout.resize(n);

	inst0.tail.resize(m);
	inst0.head.resize(m);

	inst0.pin.resize(m);
	inst0.pout.resize(m);
	inst0.opposite.resize(m);

	inst0.p.resize(n, 0);
	inst0.T.resize(n);
	inst0.c.resize(m);

	inst0.f1.resize(n, false);
	inst0.f0.resize(n, false);
	inst0.fe0.resize(m, false);

	inst0.bmaa.resize(m);
	inst0.bmna.resize(n);

	inst0.n = n;
	inst0.m = m;
	inst0.t = 0;
	inst0.isInt = inst1.isInt;
	inst0.isAsym = inst1.isAsym;
	inst0.bigM = inst1.bigM;
	inst0.inst1 = &inst1;

	vector<int> nmap(inst1.n, -1), amap(inst1.m, -1);

	rnmap.resize(n, -1);
	ramap.resize(m, -1);

	int idx = 0;
	for(int i = 0; i < inst1.n; i++) {
		if(!sol.nodes[i]) continue;
		nmap[i] = idx++;
		rnmap[nmap[i]] = i;

		inst0.p[nmap[i]] = inst1.p[i];
		inst0.T[nmap[i]] = inst1.T[i];
		inst0.f1[nmap[i]] = inst1.f1[i];

		if(inst0.T[nmap[i]])  inst0.t++;
	}
	idx = 0;
	for(int ij = 0; ij < inst1.m; ij++) {
		if(!sol.arcs[ij]) continue;
		amap[ij] = idx++;
		ramap[amap[ij]] = ij;

		const int i = inst1.tail[ij];
		const int j = inst1.head[ij];

		inst0.newArc(nmap[i], nmap[j], amap[ij], -1, inst1.c[ij]);
	}
	weight_t penalty = 0;
	for(int i = 0; i < inst1.n; i++) {
		if(!sol.nodes[i]) penalty += inst1.p[i];
	}
	inst0.offset = penalty;
	inst0.r = nmap[sol.r];

	return inst0;
}

int cntReachable(int r, Sol& sol, Inst& inst)
{
	vector<int> visited(inst.n, false);
	stack<int> stack;

	stack.push(r);
	visited[r] = true;
	
	int cnt = 1;
	while(!stack.empty()) {
		const int i = stack.top();
		stack.pop();

		for(int ij : inst.dout[i]) {
			if(!sol.arcs[ij]) continue;

			const int j = inst.head[ij];

			if(visited[j]) continue;

			visited[j] = true;
			cnt++;
			stack.push(j);
		}
	}
	
	return cnt;
}

/**
void enlargeStack()
{
	const rlim_t kStackSize = 128L * 1024L * 1024L;
	struct rlimit rl;
	int result;

	result = getrlimit(RLIMIT_STACK, &rl);
	if(result == 0) {
		if(rl.rlim_cur < kStackSize) {
			rl.rlim_cur = kStackSize;
			result = setrlimit(RLIMIT_STACK, &rl);
			if(result != 0) {
				fprintf(stderr, "setrlimit returned result = %d\n", result);
			}
		}
	}
}**/

Sol loadSol(const char* fn, Inst& inst)
{
	FILE *fp;
	char buf[256];
	int n, m, t, w, v1, v2;
	Sol sol(inst);

	sol.nodes.resize(inst.n, false);
	sol.arcs.resize(inst.m, false);
	
	if((fp=fopen(fn, "r")) == NULL) {
		EXIT("error: file not found: %s\n", fn);
	}

	while(fgets(buf, 256, fp) != NULL) {

		if(sscanf(buf, "V %d", &v1) == 1) {
			sol.nodes[v1-1] = true;
		}
		if(sscanf(buf, "E %d %d", &v1, &v2) == 2) {
			for(int ij : inst.dout[v1-1]) {
				if(inst.head[ij] == v2-1) {
					sol.arcs[ij] = true;
					break;
				}
			}
		}
	}

	// find root
	int roots = 0;
	for(int i = 0; i < inst.n; i++) {
		if(!sol.nodes[i]) continue;
		int indeg = 0;
		for(int ji : inst.din[i]) {
			if(sol.arcs[ji])
				indeg++;
		}
		if(indeg == 0) {
			roots++;
			sol.r = i;
		}
	}
	if(roots != 1) {
		printf("WARNING: invalid solution! #roots=%d", roots);
	}
	printf("loaded solution %10s obj=%lf r=%d\n", params.solfile.c_str(), format(sol.recomputeObjective(), inst), sol.r);

	return sol;
}

// calls external solver for debugging purposes
/**
weight_t oracle(Inst& inst)
{
	FILE* fp;
	char output[1024];
	char cmd[256];
	weight_t response = WMAX;

	writeInstance("oracle.stp", inst);

	sprintf(cmd, "./foo oracle.stp 2>/dev/null | grep ^UB | awk '{print $2}'");
	fp = popen(cmd, "r");

	if(fp == NULL) {
		printf("Oracle: Could not open response file.\n");
		exit(1);
	} else {
		while(fgets(output, sizeof(output)-1, fp) != NULL) {
			response = strtol(output, NULL, 10);
		}
		
		pclose(fp);
		return response;
	}
	
	return response;
}**/

void writeInstance(const char* file, Inst& inst)
{
	FILE* fp;
	if((fp=fopen(file, "w")) == NULL) {
		EXIT("error: writing instance file: %s\n", file);
	}

	vector<int> mapn(inst.n);

	int nVertices = 0;
	int nTerminals = 0;
	for(int i = 0; i < inst.n; i++) {
		if(inst.f0[i]) continue;

		mapn[i] = nVertices++;

		if(inst.p[i] > 0) nTerminals++;
	}

	int nEdges = 0;
	for(int ij = 0; ij < inst.m; ij++) {
		if(inst.fe0[ij]) continue;
		nEdges++;
	}
	
	fprintf(fp, "SECTION Graph\nNodes %d\nEdges %d\n", nVertices, nEdges);
	for(int ij = 0; ij < inst.m; ij++) {
		if(inst.fe0[ij]) continue;
		weight_t w = inst.c[ij];
		if(inst.bigM > 0 && inst.r == inst.tail[ij])
			w = 999999999999;
		fprintf(fp, "E %d %d %ld\n", mapn[inst.tail[ij]]+1, mapn[inst.head[ij]]+1, inst.c[ij]);
	}
	fprintf(fp, "END\n\n");

	fprintf(fp, "SECTION Terminals\n");
	for(int i = 0; i < inst.n; i++) {
		if(inst.f0[i]) continue;
		if(inst.p[i] == WMAX)
			fprintf(fp, "TP %d %ld\n", mapn[i]+1, 999999999999);
		else
			if(inst.p[i] > 0) fprintf(fp, "TP %d %ld\n", mapn[i]+1, inst.p[i]);
	}

	fprintf(fp, "Root %d\n", mapn[inst.r]+1);
	fprintf(fp, "Fixed %ld\n", inst.offset);
	fprintf(fp, "END\n\n");

	fclose(fp);
}

void printHeap(PQMinFib<weight_t,int>& PQ, Inst& inst, int id)
{
	printf("Heap %d = { ", id);
	for(auto it = PQ.begin(); it != PQ.end(); ++it) {
		auto h = PQMinFib<weight_t,int>::s_handle_from_iterator(it);
		auto p = *h;
		printf("%0.2lf ", format(p.first, inst));
	}
	printf("}");
}

Sol dmst(Inst& inst, vector<weight_t>& cr)
{
	vector<int> enter(inst.n, -1);
	stack<int> Q;
	vector<PQMinFib<weight_t,int>> PQ(inst.n);
	vector<int> rankW(inst.n);
	vector<int> parentW(inst.n);
	vector<int> rankS(inst.n);
	vector<int> parentS(inst.n);
	vector<flag_t> inH(inst.m, false);
	vector<flag_t> forbidden(inst.m, false);

	// forest
	vector<int> parent(inst.n*2, -1);
	vector<int> leaf(inst.n, -1);
	vector<int> leaf2(inst.n*2, -1);
	vector<vector<int>> lastCycle(inst.n);
	vector<int> g2f(inst.m);
	vector<int> f2g(inst.n*2);
	vector<vector<int>> children(inst.n*2);
	int fIdx = -1;

	vector<PQMinFib<weight_t,int>::handle_type> PQhandle(inst.m);

	boost::disjoint_sets<int*,int*> S(&rankS[0], &parentS[0]);
	boost::disjoint_sets<int*,int*> W(&rankW[0], &parentW[0]);

	for(int ij = 0; ij < inst.m; ij++) {
		if(inst.fe0[ij]) continue;
		int j = inst.head[ij];
		PQhandle[ij] = PQ[j].push(make_pair(inst.c[ij], ij));
		cr[ij] = inst.c[ij];
	}

	vector<int> nodes;
	for(int i = 0; i < inst.n; i++) {
		if(inst.f0[i]) continue;
		nodes.push_back(i);
		S.make_set(i);
		W.make_set(i);
		Q.push(i);
	}

	int iter = 0;
	while(!Q.empty()) {
		int k = Q.top();
		Q.pop();

		int S_k = S.find_set(k);
		if(PQ[S_k].empty()) {
			continue;
		}

		int ij = PQ[S_k].top().second;
		weight_t val1 = PQ[S_k].top().first;
		PQ[S_k].pop();

		int i = inst.tail[ij];
		int j = inst.head[ij];
		int W_i = W.find_set(i);
		int W_j = W.find_set(j);
		int S_i = S.find_set(i);
		int S_j = S.find_set(j);

		if(S_i == S_k) {
			Q.push(k);
		} else {
			inH[ij] = true;

			fIdx++;
			f2g[fIdx] = ij;
			g2f[ij] = fIdx;
			for(int ij2 : lastCycle[S_k]) {
				int fIdx2 = g2f[ij2];
				parent[fIdx2] = fIdx;
				children[fIdx].push_back(fIdx2);
			}
			if(lastCycle[S_k].empty()) {
				leaf[j] = fIdx;
				leaf2[fIdx] = j;
			}

			if(W_i != W_j) {
				W.link(W_j, W_i);
				enter[S_k] = ij;

			} else {
				
				// detect heaviest arc on cycle
				weight_t val = -1; int S_v;
				int costliest = -1;
				int ij2 = ij;
				lastCycle[S_k].clear();
				
				while(ij2 != -1) {
					
					int i = inst.tail[ij2], j = inst.head[ij2];
					int S_i = S.find_set(i), S_j = S.find_set(j);

					if(inst.c[ij2] > val) {
						val = inst.c[ij2];
						costliest = ij2;
						S_v = S_j;
					}
					lastCycle[S_k].push_back(ij2);
					ij2 = enter[S_i];
				}
				
				vector<PQMinFib<weight_t,int>::handle_type> test;
				for(auto it = PQ[S_k].begin(); it != PQ[S_k].end(); ++it) {
					auto h = PQMinFib<weight_t,int>::s_handle_from_iterator(it);
					test.push_back(h);
				}
				for(auto t : test) {
					auto p = *t;
					PQ[S_k].decrease(t, make_pair(p.first-val1, p.second));
					cr[p.second] -= val1;
				}
				
				// merge
				ij2 = enter[S_i];
				while(ij2 != -1) {
					int i = inst.tail[ij2], j = inst.head[ij2];
					int S_i = S.find_set(i), S_j = S.find_set(j);
					
					test.clear();
					for(auto it = PQ[S_j].begin(); it != PQ[S_j].end(); ++it) {
						auto h = PQMinFib<weight_t,int>::s_handle_from_iterator(it);
						test.push_back(h);
					}
					for(auto t : test) {
						auto p = *t;
						PQ[S_j].decrease(t, make_pair(p.first-cr[ij2], p.second));
						cr[p.second] -= cr[ij2];

					}
					S.link(S_j, S_k);
					int new_k = S.find_set(S_k), to_merge;
					if(new_k == S_k) {
						to_merge = S_j;
					} else {
						to_merge = S_k;
						lastCycle[new_k] = lastCycle[S_k];
						S_k = new_k;
						enter[new_k] = -1;
					}
					PQ[new_k].merge(PQ[to_merge]);
					
					ij2 = enter[S_i];
				}
				enter[S_k] = -1;
				Q.push(S_k);

			}
		}
	}

	// recover solution from support graph H
	stack<int> R;
	for(int ij = 0; ij < inst.m; ij++) {
		if(inst.fe0[ij]) continue;
		if(!inH[ij]) continue;

		int i = inst.tail[ij];
		int j = inst.head[ij];
		int f = g2f[ij];
		if(f != -1) {
			if(parent[f] == -1) {
				R.push(f);
			}
		}
	}
	
	vector<int> B;
	while(!R.empty()) {
		int f = R.top();
		R.pop();

		int ij = f2g[f];
		int i = inst.tail[ij];
		int j = inst.head[ij];
		B.push_back(ij);
		f = leaf[j];
		while(parent[f] != -1) {
			int last = f;
			f = parent[f];
			for(int i : children[f]) {
				parent[i] = -1;
				if(last == i) continue;
				R.push(i);
			}
		}
	}

	Sol sol(inst);
	sol.r = inst.r;
	for(int ij : B) {
		sol.arcs[ij] = true;
		sol.nodes[inst.tail[ij]] = true;
		sol.nodes[inst.head[ij]] = true;
	}
	sol.recomputeObjective();
	return sol;
}

