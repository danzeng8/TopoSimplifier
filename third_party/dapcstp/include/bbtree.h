/**
 * \file   bbtree.h
 * \brief  branch-and-bound-tree
 *
 * \author Martin Luipersbeck 
 * \date   2015-10-03
 */

#ifndef BBTREE_H_
#define BBTREE_H_

#include <random>

#include "inst.h"
#include "sol.h"
#include "bbnode.h"
#include "options.h"
//#include "tbb/concurrent_priority_queue.h";
#include "boost/thread/mutex.hpp"
#include <boost/thread/condition.hpp>

class compareB {
public:
	bool operator()(const pair<weight_t, BBNode*> & u, const pair<weight_t, BBNode*> & v) const {
		return u.first > v.first;
	}
};

class compareBR {
public:
	bool operator()(const pair<weight_t, BBNode*> & u, const pair<weight_t, BBNode*> & v) const {
		return u.first > v.first;
	}
};

class BBTree
{
public:

	enum State     { BB_NONE, BB_OPTIMAL, BB_TIMELIMIT, BB_NODELIMIT, BB_SOLLIMIT, BB_MEMLIMIT };
	enum NodeState { BB_INFEAS, BB_CUTOFF, BB_LEAF, BB_BRANCH };

	BBTree(Inst& inst);
	~BBTree();

	void     processRoots();
	void     initHeur();
	int      initPrep();
	void     solve();
	

	// setters
	void     setSolLim(int i)         { solLim = i; }
	void     setNodeLim(int i)        { nodeLim = i; }
	void     setTimeLim(double d)     { if(d < 0) timeLim = DMAX; else timeLim = d; }

	void     setIncumbent(Sol& sol)   { inc = sol; ub = inc.obj; }
	void     setOutput(bool b)        { bOutput = b; }
	void     setCutUp(weight_t d)     { cutup = d; }
	void     setRecover(bool b)       { bRecover = b; }
	void     setBestKnown(double d)   { bestKnown = d; }

	// getters
	int      getNnodes()          { return nIter; }
	double   getTime()            { return bbTime; }
	double   getTimeBest()        { return timeBestSol; }
	double   getHeurTime()        { return heurTime; }
	double   getRootTime()        { return rootTime; }
	double   getHeurBBTime()      { return heurBBTime; }
	double   getRootGap()         { return gapP(rootlb, rootub); }
	double   getGap()             { return gapP(bestlb, ub); }
	State    getBBState()         { return tState; }
	weight_t getUB()              { return ub; }
	weight_t getLB()              { return bestlb; }
	weight_t getLBM()             { return lbM; }
	weight_t getRootUB()          { return rootub; }
	weight_t getRootLB()          { return rootlb; }
	Sol      getSol()             { return inc; }
	State    getState()           { return tState; }
	int      getNroots()          { return nRoots; }
	int      getNrootsProcessed() { return nRootsProcessed; }
	int      getNrootsOpen()      { return nRootsOpen; }

	Inst&    getInst1()           { return inst1; }
	Sol&     getInc1()            { return inc1; }

	mutable boost::mutex mutexLock;
	//boost::mutex::scoped_lock lock(mutexLock);
private:
	
	State tState;
	Inst& inst;
	Inst  inst1;

	//tbb::concurrent_priority_queue < pair<weight_t, BBNode*>, compareBR > PQmax;
	PQMax<weight_t,BBNode*> PQmax;
	//tbb::concurrent_priority_queue < pair<weight_t, BBNode*>, compareB > PQmin;
	PQMin<weight_t,BBNode*> PQmin;
	PQMin<weight_t, BBNode*>::handle_type pushToPQMin(pair<weight_t, BBNode*> pair);
	PQMax<weight_t, BBNode*>::handle_type pushToPQMax(pair<weight_t, BBNode*> pair);
	bool popPQMin(pair<weight_t, BBNode*> & pair);
	bool popPQMax(pair<weight_t, BBNode*> & pair);
	bool	isPQMinEmpty();
	bool	isPQMaxEmpty();
	

	// incumbent on preprocessed and unpreprocessed graph
	Sol inc, inc1;
	std::vector<Sol*> pool;

	mt19937 rndGen;

	weight_t ub, bestlb;
	weight_t rootlb, rootub;
	double bestKnown = -1;

	weight_t bestSingleNodeSolObj;
	int      bestSingleNodeSolNode;
	weight_t cutup = -1.0;

	double bbTime, timeBestSol, heurTime, heurBBTime, rootTime;
	int    nRoots, nRootsProcessed, nRootsOpen;
	int    nImprovements, nIter;

	weight_t lbM = -1;
	Inst instM;
	vector<weight_t> piM, crM;

	vector<int> prio;
	vector<weight_t> cr, pi;
	vector<double> crf, pif;

	// states
	bool processedRoots = false;
	bool bHeur = false;
	// stores solutions for later expansion (disabled when solving the instance to regain original solution)
	// disables tests that would merge/contract nodes/arcs
	bool bRecover = false;
	// disables all output to terminal
	bool bOutput = true;

	// limits
	int nodeLim, solLim; double timeLim;

	// node operations
	BBNode*           select();
	BBTree::NodeState process(BBNode* b);
	void              add(BBNode* b);
	int               selectBranchVariable(BBNode* b);
	void              branch(BBNode* b);
	void              evalLeaf(BBNode* b);
	BBTree::NodeState strengthenBounds(BBNode* b);
	BBNode*           makeRoot(int r, weight_t lb, vector<int>& fe0);

	// general operations
	void              freeOpenNodes();
	bool              updatePrimal(Inst& inst, Sol& sol);
	int               preprocess(Inst& inst);
	bool              isFeas(Inst& inst, bool bDoNRtest = true);
	vector<weight_t>  setSupportGraph(Inst& inst);
	vector<weight_t>  setSupportGraphf(Inst& inst, vector<double>& cr);
	weight_t          perturbedPrimalHeur(Inst& inst);
	vector<int>       sortedListPotentialRoots();
	void              fixTerm(Inst& inst, int t, vector<int>& fe0);

	// print to console
	void printHeader();
	void printBBLine(BBNode* b, NodeState state, bool bExist);
	void printRootHeader();
	void printRootLine(BBNode* b);
	void printHeurHeader();
	void printHeurLine(int it, weight_t obj, bool bImproved, double time);
	void printHeur1Summary();
	void printHeur2Summary();
	void printRootSummary();
	void printBBSummary();
};

#endif // BBTREE_H_
