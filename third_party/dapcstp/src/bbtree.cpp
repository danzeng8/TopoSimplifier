/**
 * \file   bbtree.cpp
 * \brief  branch-and-bound-tree
 *
 * \author Martin Luipersbeck
 * \date   2015-10-03
 */

#include "bbtree.h"
#include "ds.h"
#include "bounds.h"
#include "prep.h"
#include "heur.h"
#include "util.h"
#include "stats.h"
#include "print.h"
#include "options.h"
#include "timer.h"
#include "procstatus.h"
//#include "tbb/parallel_for.h"

#include <boost/thread.hpp>
//#include <boost/lockfree/spsc_queue.hpp>

#include <stack>
#include <map>

BBTree::BBTree(Inst& inst) : inst(inst), inc(inst), inst1(inst), inc1(inst1)
{
	lbM = 0;
	bestlb = 0;
	rootlb = 0;
	rootub = WMAX;
	ub = WMAX;
	tState = BB_NONE;

	solLim = IMAX;
	nodeLim = IMAX;
	timeLim = DMAX;

	prio.resize(inst.n, 0);

	rndGen = mt19937(params.seed);

	inst.bmaa.resize(inst.m);
	inst.bmna.resize(inst.n);

	// backmapping
	for (int ij = 0; ij < inst.m; ij++) inst.bmaa[ij].push_back(ij);

	inst.inst1 = &inst1;

	cr.resize(inst.m, 0);
	pi.resize(inst.n, 0);
	crf.resize(inst.m, 0);
	pif.resize(inst.n, 0);
	piM.resize(inst.n + 1, 0);
	crM.resize(inst.m + inst.n, 0);

	// executing perturbed primal for eps=0.0 is redundant
	if (params.heureps == 0.0) {
		params.perturbedheur = false;
	}
	else {
		if (params.heureps < 0.0) {
			params.heureps = inst.isInt ? 0.05 : 0.005;
		}
	}

	// when dealing with non-integral weights, use arc saturation epsilon during dual ascent
	if (params.dasat < 0.0) {
		params.dasat = inst.isInt ? 0.0 : 1e-4*params.precision;
	}

	// treat best single node solution separately as it might be eliminated by reduction tests
	bestSingleNodeSolObj = WMAX;
	bestSingleNodeSolNode = -1;
	int f1 = 0;
	for (int i = 0; i < inst.n; i++) {
		if (inst.f1[i]) f1++;
	}
	if (inst.r == -1 && f1 == 0) {
		weight_t Pmax = 0.0, P = 0.0;
		for (int i = 0; i < inst.n; i++) {
			P += inst.p[i];
			if (inst.p[i] > Pmax) {
				Pmax = inst.p[i];
				bestSingleNodeSolNode = i;
			}
		}
		bestSingleNodeSolObj = P - Pmax + inst.offset;
	}
}

BBTree::~BBTree()
{
	for (auto S : pool)
		delete S;
}

bool BBTree::updatePrimal(Inst& inst, Sol& sol)
{
	if (sol.obj < ub) {
		// validate solution before acceptance
		const bool bValid = sol.validate();
		if (!bValid) {
			printf("WARNING: obtained solution infeasible, discarding.\n");
			return false;
		}
		nImprovements++;
		inc.update(sol);

		// if in recovery mode, store solution as partial solution on the unpreprocessed inst1
		if (!bRecover) {
			inc1 = genPartialSol(inc, inst);
			timeBestSol = Timer::total.elapsed().getSeconds();
		}

		ub = sol.obj;

		return true;
	}

	return false;
}

weight_t BBTree::perturbedPrimalHeur(Inst& inst)
{
	// compute perturbed cost
	vector<double> c1(inst.m);
	//tbb::parallel_for(tbb::blocked_range<int>(0, inst.m),
		//[&](tbb::blocked_range<int> r)
	//{
		//for (int i = r.begin(); i < r.end(); ++i)
		//{
			//cout << params.heureps << endl;
			for (int i = 0; i < inst.m; i++) {
				const double eps = inc.arcs[i] ? -params.heureps : params.heureps;
				c1[i] = max(0.0, inst.c[i] * (1 + eps));
			}
				
		//}
	//});

	daR(inst.r, inst, c1, crf, pif, ub, params.daeager, &inc, true);

	// apply primI to support graph
	auto c2 = setSupportGraphf(inst, crf);
	Sol sol = primI(inst.r, inst, c2);
	updatePrimal(inst, sol);

	return sol.obj;
}

void BBTree::initHeur()
{
	Timer tHeur(true);
	bHeur = true;
	BBNode* b = new BBNode(&inst);

	vector<int> roots = sortedListPotentialRoots();
	//printHeurHeader();

	map<int, int> rootArcs;
	if (inst.bigM > 0) {
		roots.clear();
		for (int ri : inst.dout[inst.r]) {
			int i = inst.head[ri];
			roots.push_back(i);
			rootArcs[i] = ri;
		}
		sort(roots.begin(), roots.end(), [&](int i, int j) -> bool { return inst.p[i] > inst.p[j]; });
	}

	int origRoot = inst.r;
	int iter = 0;
	for (int r : roots) {
		int currRoot = r;
		Timer tIteration(true);
		weight_t oldub = ub;

		// choose arc costs for applying dual ascent
		vector<double> c1(b->inst->m);
		if (params.heursupportG && params.heureps > 0.0 && iter > 0) {

			 for(int i = 0; i < inst.m; i++) {
				//tbb::parallel_for(tbb::blocked_range<int>(0, inst.m),
					//[&](tbb::blocked_range<int> r)
				//{
					//for (int i = r.begin(); i < r.end(); ++i)
					//{
						const double eps = inc.arcs[i] ? -params.heureps : params.heureps;
						c1[i] = max(0.0, inst.c[i] * (1 + eps));
					//}
				//});
			}
		}
		else {
			for (int i = 0; i < inst.m; i++) {
				c1[i] = inst.c[i];
			}
		}

		// apply primI to G_S or G by adjusting arc costs
		Sol sol(*b->inst);
		weight_t lb = 0.0;
		if (params.heursupportG) {
			if (inst.bigM > 0) {
				int curr_rootArc = rootArcs[currRoot];
				for (int ri : inst.dout[origRoot]) {
					c1[ri] = (ri == curr_rootArc ? 0 : WMAX - 1);
				}
				lb = daR(origRoot, inst, c1, crf, pif, ub, params.daeager, &inc, true);
				//cout << "lower bound 0.5a " << lb << endl;
			}
			else {
				lb = daR(r, inst, c1, crf, pif, ub, params.daeager, &inc);
				//cout << "lower bound 0.5 " << lb << endl;
			}
			auto c2 = setSupportGraphf(*b->inst, crf);
			sol = primI(r, *b->inst, c2);
		}
		else {
			sol = primI(r, *b->inst, b->inst->c);
		}
		//cout << "lower bound 1 " << lb << endl;

		pool.push_back(new Sol(sol));

		// generate subproblem to apply B&B on G_S
		if (params.heurbb && params.heursupportG) {
			vector<int> fe0;
			fe0.reserve(inst.m);
			for (int ij = 0; ij < inst.m; ij++) {
				if (crf[ij] > params.dasat)
					fe0.push_back(ij);
			}
			/**
			tb::parallel_for(tbb::blocked_range<int>(0, inst.m),
				[&](tbb::blocked_range<int> r)
			{
				for (int ij = r.begin(); ij < r.end(); ++ij)
				{
					if (crf[ij] > params.dasat)
						fe0.push_back(ij);
				}
			});**/

			makeRoot(r, lb, fe0);
		}

		bool bImproved = updatePrimal(*b->inst, sol);
		tIteration.stop();

		// an improve solution may also be produced during makeRoot, which is evaluated at once
		bImproved |= (ub < oldub);

		//printHeurLine(iter, sol.obj, bImproved, tIteration.elapsed().getSeconds());

		// if bound is available on the rooted equivalent instance using big-M arcs,
		// perform bound-based reductions and copy them to the unrooted instance.
		if (origRoot == -1 && params.semiBigM && params.initprep) {
			auto p = bbred(instM, lbM, ub, crM, piM);
			for (int ij = 0; ij < inst.m; ij++) {
				if (inst.fe0[ij]) continue;
				if (instM.fe0[ij]) {
					inst.delArc(ij);
					inst.fe0[ij] = true;
				}
			}
			for (int i = 0; i < inst.n; i++) {
				if (inst.f0[i] || inst.f1[i]) continue;
				if (instM.f0[i]) {
					inst.removeNode(i);
				}
				else if (instM.f1[i]) {
					inst.f1[i] = true;
					inst.T[i] = true;
					inst.p[i] = WMAX;
				}
			}
		}

		if (++iter == params.heurroots) break;
		if (ub - lbM <= params.absgap) {
			break;
		}
	}

	heurTime = tHeur.elapsed().getSeconds();
	rootub = ub;
	//printHeur1Summary();

	b->updateNodeSize();

	if (params.heurbb && pool.size() > 0) {

		// compute the union of all feasible solutions, but bidirect in any case to
		// ensure that there exists a feasible solution if instance is unrooted
		vector<flag_t> bInUnionGraph(inst.m, 0);
		for (int a = 0; a < pool.size(); a++) {
			//tbb::parallel_for(tbb::blocked_range<int>(0, inst.m),
			//	[&](tbb::blocked_range<int> r)
			//{
				//for (int ij = r.begin(); ij < r.end(); ++ij)
				//{
				for (int ij = 0; ij < inst.m; ij++) {
					if (inst.fe0[ij]) continue;
					const int ji = inst.opposite[ij];
					if (pool[a]->arcs[ij]) {
						bInUnionGraph[ij] = true;
					}
					if (ji != -1 && pool[a]->arcs[ji]) {
						bInUnionGraph[ij] = true;
					}
				}
				//}
			//});
		}

		vector<int> fe0;
		fe0.reserve(inst.m);
		//tbb::parallel_for(tbb::blocked_range<int>(0, inst.m),
			//[&](tbb::blocked_range<int> r)
		//{
			//for (int ij = r.begin(); ij < r.end(); ++ij)
			//{
				for (int ij = 0; ij < inst.m; ij++) {
					if (bInUnionGraph[ij]) continue;
					fe0.push_back(ij);
				}
			//}
		//});
		// choose root at random from the set of initial heuristic iterations
		int rndroot = rand() % min(params.heurroots, (int)roots.size());
		int r = roots[rndroot];

		makeRoot(r, 0.1, fe0);

		// process all added B&B nodes to improve the upper bound
		Timer tHeurBB(true);
		processedRoots = true; // skip root processing step
		double origTimelimit = timeLim;
		setTimeLim(params.heurbbtime);
		
		//boost::thread_group threads;
		//for (int i = 0; i < 2; ++i) {
			//boost::thread* threadTree = new boost::thread(&BBTree::solve, this);
			//threads.add_thread(threadTree);
		//}
		//threads.interrupt_all();
		//threads.join_all();
		//boost::thread * threads[10];
		//for (int t = 0; t < 10; t++) {
			//threads[t] = new boost::thread(solve());
		//}
																   //boost::thread t{ this->solve };

			solve();
		//}

		//for (int t = 0; t < 10; t++) {
			//threads[t]->join();
		//	delete threads[t];
		//}

		// if bigM is used, a valid lower bound is available, otherwise it should be reset,
		// as the bestlb of the restricted B&B is not globally valid
		bestlb = params.semiBigM ? lbM : 0;

		timeLim = origTimelimit;
		nIter = 0;
		processedRoots = false;
		heurBBTime = tHeurBB.elapsed().getSeconds();
		rootub = ub;
		rootub = ub;

		//printHeur2Summary();
	}

	inst.r = origRoot;
	bHeur = false;
}

PQMin<weight_t, BBNode*>::handle_type BBTree::pushToPQMin(pair<weight_t, BBNode*> pair)
{
	
	//lock.unlock();

	PQMin<weight_t, BBNode*>::handle_type index = PQmin.push(pair);
	//the_condition_variable.notify_one();
	return index;
}

PQMax<weight_t, BBNode*>::handle_type BBTree::pushToPQMax(pair<weight_t, BBNode*> pair)
{
	
	PQMax<weight_t, BBNode*>::handle_type index = PQmax.push(pair);
	//the_condition_variable.notify_one();
	return index;

}

bool BBTree::popPQMin(pair<weight_t, BBNode*> & pair) {
	//if (PQmin.empty())
	//{
	//	return false;
	//}

	//pair = PQmin.top();

	PQmin.pop();
	return true; 
	//PQmin.pop();
}

bool BBTree::popPQMax(pair<weight_t, BBNode*> & pair) {
	//if (PQmax.empty())
	//{
		//return false;
	//}

	//pair = PQmax.top();
	PQmax.pop();
	return true;
}

bool BBTree::isPQMinEmpty() {

	return PQmin.empty();
}

bool BBTree::isPQMaxEmpty() {

	return PQmax.empty();
}

void BBTree::processRoots()
{
	Timer tRoot(true);

	// during root node processing, nodes are iteratively fixed to zero by
	// settting the incoming arc costs to infinity
	vector<weight_t> c = inst.c, p = inst.p;
	vector<flag_t> f0 = inst.f0, f1 = inst.f1, T = inst.T;
	weight_t offset = inst.offset;

	if (!bRecover) {
		ProgramStats::initRootNodeStats();
	}


//	printRootHeader();
	vector<int> roots = sortedListPotentialRoots();
	vector<int> fe0;
	int processed = 0;

	tState = BB_NONE;
	for (int k : roots) {
		if (params.semiBigM && inst.r == -1 && lbM > 0 && lbM + crM[inst.m + k] >= ub) {
			continue;
		}

		inst.f1[k] = true;
		inst.T[k] = true;
		inc.rootSolution(k);
		weight_t lb = daR(k, inst, inst.c, cr, pi, ub, params.daeager, &inc);

		if (ub - lb > params.absgap) {
			BBNode* b = makeRoot(k, lb, fe0);
			if (b != nullptr) {
				InstSizeData sdata = b->inst->countInstSize();
				if (!bRecover) {
					ProgramStats::addRootNodeStats(sdata);
				}

				//printRootLine(b);
			}
		}

		if (!inst.isAsym) {
			fixTerm(inst, k, fe0);
		}
		else {
			inst.f1[k] = f1[k];
			inst.T[k] = T[k];
		}
		processed++;

		if (ProcStatus::mem() > params.memlimit || tState == BB_MEMLIMIT) {
			tState = BB_MEMLIMIT;
			if (bOutput)
				printf(" --- out of memory during root node processing\n");
			break;
		}
		if (tRoot.elapsed().getSeconds() > timeLim) {
			tState = BB_TIMELIMIT;
			if (bOutput)
				printf(" --- reached time limit during root node processing\n");
			break;
		}
		if (inst.offset >= ub) break;
	}
	inst.c = c;
	inst.T = T;
	inst.p = p;
	inst.f0 = f0;
	inst.f1 = f1;
	inst.offset = offset;

	timeLim = timeLim; //max(0.0, timeLim - tRoot.elapsed().getSeconds());

	if (!bRecover) {
		if (tState != BB_MEMLIMIT && tState != BB_TIMELIMIT) {
			if (PQmin.size() == 0) {
				// solved in root to optimality
				//if (bOutput)
					//printf(" --- no root remaining open\n");
				rootlb = ub;
				bestlb = rootlb;
			}
			else {
				rootlb = PQmin.top().first;
				bestlb = rootlb;
			}
		}
		else {
			if (params.bigM) {
				bestlb = lbM;
			}
		}

		rootub = ub;
		rootTime = tRoot.elapsed().getSeconds();
		nRoots = (int)roots.size();
		nRootsProcessed = processed;
		nRootsOpen = (int)PQmin.size();

		//if (!PQmin.empty()) {
			//ProgramStats::averageRootNodeStats(nRootsOpen);
		//}
	}

	//printRootSummary();
}

void BBTree::solve()
{
	if (cutup >= 0.0 && cutup < ub) {
		ub = cutup;
	}

	bool bSolvedInRoot = false, bOutOfMemInRoot = false;
	if (!processedRoots) {
		processRoots();
		processedRoots = true;
	}
	if (tState == BB_MEMLIMIT || tState == BB_TIMELIMIT) {
		bOutOfMemInRoot = true;
		freeOpenNodes();
		return;
	}
	else {

		//if(isPQMinEmpty()) {
		//if (!PQmin.empty()) {
			//printHeader();
		//}
		//else {
			if (!PQmin.empty()) {
				bSolvedInRoot = true;
				tState = BB_OPTIMAL;
			}
		//}
	}
	nImprovements = 0;
	nIter = 0;
	Timer tBB(true);
	while (!PQmin.empty()) {
	//while(!isPQMinEmpty()) {
		BBNode* b = select();
		if (b == NULL) { break; }
		NodeState state;
		if (!b->processed) {
			state = process(b);
			//cout << "before p lock " << endl;
			//mutexLock.lock();
			//cout << "after p lock " << endl;
			switch (state) {
			case BB_INFEAS:

			case BB_CUTOFF:
				if (b->v != -1) {
					prio[b->v]++;
				}
				//mutexLock.unlock();
				break;
			case BB_LEAF:
				//cout << "leaf " << endl;
				evalLeaf(b);
				//mutexLock.unlock();
				break;
			case BB_BRANCH:
				//cout << "branch " << endl;
				selectBranchVariable(b);
				//mutexLock.unlock();
				branch(b);
				
				break;
			//default:
				//mutexLock.unlock();
			}
			//cout << "before unlock " << endl;
			
			//cout << "after unlock " << endl;
		}
		else {
			//cout << "processed already " << endl;
			state = (NodeState)b->state;
			branch(b);
			//cout << "after branch 2 " << endl;
		}
		//cout << "before pq min empty " << endl;
		//if (!PQmin.empty()) {
		if(!isPQMinEmpty()) {
			bestlb = PQmin.top().first;
		}
		else {
			bestlb = ub;
		}
		//cout << "termination checks " << endl;
		bool bExit = false;
		if (++nIter >= nodeLim) { tState = BB_NODELIMIT; bExit = true; }
		if (nImprovements >= solLim) { tState = BB_SOLLIMIT;  bExit = true; }
		if (tBB.elapsed().getSeconds() > timeLim) { tState = BB_TIMELIMIT; bExit = true; }
		else {
			if (ProcStatus::mem() > params.memlimit) { tState = BB_MEMLIMIT;  bExit = true; }
			else {
				if (PQmin.size() == 0 || PQmin.top().first >= ub) { 
					tState = BB_OPTIMAL;   bExit = true;
				}
			}
		}
		//printBBLine(b, state, bExit);
		
		// case two occurs when both nodes get pruned
		if (state != BB_BRANCH || (state == BB_BRANCH && !b->feas)) {
			delete b->inst;
			b->inst = nullptr;
			delete b;
		}

		if (bExit) break;
	}
	//cout << "after while " << endl;
	const double finishedTime = tBB.elapsed().getSeconds();

	// heuristic should not affect optimality gap
	if (!bHeur)
		if (tState == BB_OPTIMAL) {
			bestlb = ub;
		}

	// adjust bound if a single node solution is optimal, which may get removed by preprocessing
	if (!bHeur && bestSingleNodeSolObj < bestlb) {
		ub = bestSingleNodeSolObj;
		inc1 = Sol(inst1);
		inc1.nodes[bestSingleNodeSolNode] = true;
		inc1.r = bestSingleNodeSolNode;
	}
	//cout << "before free open nodes " << endl;
	freeOpenNodes();
	//cout << "after free nodes " << endl;

	bOutput = true;
	/**if (bOutput) {
		printf(" --- ");
		if (bSolvedInRoot)
			printf("solved to optimality in root");
		else
			switch (tState) {
			case BB_OPTIMAL:
				printf("solved to optimality");
				break;
			case BB_TIMELIMIT:
				printf("time limit reached");
				break;
			case BB_NODELIMIT:
				printf("node limit reached");
				break;
			case BB_SOLLIMIT:
				printf("solution limit reached");
				break;
			case BB_MEMLIMIT:
				printf("memory limit reached");
				break;
			}
		printf(" ( %0.1lf sec. )\n", finishedTime);
	}**/

	// recover a partial solution (which may contain antiparallel arcs)
	if (!bRecover && inc1.partial) {
		recoverPartialSol(inc1, inst1);
	}
	bbTime = finishedTime;
	//if (!bHeur)
		//printBBSummary();
}

void BBTree::add(BBNode* b)
{
	//boost::mutex::scoped_lock lock(mutexLock);
	switch (params.nodeselect) {
	case 0: // worst-bound
		//std::pair<weight_t, BBNode*> maxPair = make_pair(b->lb, b);pushToPQMax(make_pair(b->lb, b));
		b->pqposMax = PQmax.push(make_pair(b->lb, b));
		break;
	case 1:  // dfs
		//std::pair<weight_t, BBNode*> maxPair = make_pair(b->depth, b);
		b->pqposMax =  PQmax.push(make_pair(b->depth, b)); //pushToPQMax(make_pair(b->depth, b));
		break;
	}

	// always add to min-queue for displaying best bound
	b->pqposMin = PQmin.push(make_pair(b->lb, b));
	//b->pqposMin = pushToPQMin(make_pair(b->lb, b));
}

BBNode* BBTree::makeRoot(int r, weight_t lb, vector<int>& fe0)
{
	BBNode* b = new BBNode(&inst, r, fe0);
	b->updateNodeSize();
	NodeState state = process(b);

	b->state = (int)state;

	if (state == BB_INFEAS || state == BB_CUTOFF) {
		delete b->inst;
		delete b;
		return nullptr;
	}
	else {
		switch (state) {
		case BB_LEAF:
			evalLeaf(b);
			break;
		case BB_BRANCH:
			selectBranchVariable(b);
			add(b);
			break;
		}
	}

	return b;
}

vector<weight_t> BBTree::setSupportGraph(Inst& inst)
{
	vector<weight_t> c1 = inst.c;
	int sat = 0;
	for (int ij = 0; ij < inst.m; ij++) {
		if (inst.fe0[ij]) continue;
		if (cr[ij] > params.dasat) {
			c1[ij] = WMAX - 1;
		}
		else {
			sat++;
		}
	}

	return c1;
}

vector<weight_t> BBTree::setSupportGraphf(Inst& inst, vector<double>& crf)
{
	auto c1 = inst.c;
	for (int ij = 0; ij < inst.m; ij++) {
		if (inst.fe0[ij]) continue;
		if (crf[ij] > params.dasat) {
			c1[ij] = WMAX - 1;
		}
	}

	return c1;
}

BBTree::NodeState BBTree::strengthenBounds(BBNode* b)
{
	if (params.daiterations <= 1) return BB_BRANCH;
	if (ub - b->lb <= params.absgap) return BB_CUTOFF;

	// direct arborescence if they come from a different root
	for (int p = 0; p < pool.size(); p++) {
		pool[p]->rootSolution(b->inst->r);
	}

	// execute dual ascent with different guiding solutions
	shuffle(pool.begin(), pool.end(), rndGen);
	weight_t lb;
	int maxsize = min(params.daiterations - 1, (int)pool.size());
	for (int i = 0; i < maxsize; i++) {

		lb = daR(b->inst->r, *b->inst, b->inst->c, cr, pi, ub, params.daeager, pool[i]);
		b->lb = max(b->lb, lb);
		if (ub - b->lb <= params.absgap) {
			return BB_CUTOFF;
		}

		if (b->depth == 0 || !params.redrootonly)
			bbred(*b->inst, lb, ub, cr, pi);
		// test for infeasibility as instance might become infeasible when applying
		// reduction tests using ub <= than the obj. of optimal solution available in this node
		b->feas = isFeas(*b->inst);
		if (!b->feas) {
			return BB_INFEAS;
		}
	}

	return BB_BRANCH;
}

BBTree::NodeState BBTree::process(BBNode* b)
{
	b->processed = true;
	Inst& inst = *b->inst;

	if (b->depth == 0 || !params.redrootonly)
		preprocess(inst);

	b->feas = isFeas(inst);
	if (!b->feas) return BB_INFEAS;

	// improve dual bound
	weight_t lb = daR(b->inst->r, inst, inst.c, cr, pi, ub, params.daeager, &inc);
	b->lb = max(b->lb, lb);
	if (ub - b->lb <= params.absgap) return BB_CUTOFF;

	// after bound-based reductions the node might become infeasible if no improving solution is available
	if (b->depth == 0 || !params.redrootonly)
		bbred(inst, lb, ub, cr, pi);
	b->feas = isFeas(inst);
	if (!b->feas) return BB_INFEAS;

	if (ub - b->lb <= params.absgap) return BB_CUTOFF;

	// try to improve bounds further
	auto state = strengthenBounds(b);
	if (state == BB_CUTOFF || state == BB_INFEAS) return state;

	// try perturbed heuristic
	if (params.perturbedheur) {
		perturbedPrimalHeur(inst);
	}
	else {
		auto c1 = setSupportGraph(inst);
		auto sol = primI(inst.r, inst, c1);
		updatePrimal(inst, sol);
	}
	if (ub - b->lb <= params.absgap) return BB_CUTOFF;

	if (state != BB_CUTOFF && state != BB_INFEAS) {
		b->updateNodeSize();
		if (b->nfree == 0)
			return BB_LEAF;
		else
			return BB_BRANCH;
	}

	return state;
}

void BBTree::branch(BBNode* b)
{
	//cout << "in branch " << endl;
	int v = b->v2;
	if (v != -1) {

		// test feasibility of both branches to avoid unnecessary copying
		bool b0feas, b1feas;
		b->inst->f0[v] = true;
		b0feas = isFeas(*b->inst, false);
		b->inst->f0[v] = false;

		b->inst->f1[v] = true;
		b1feas = isFeas(*b->inst, false);
		b->inst->f1[v] = false;

		if (b0feas && b1feas) {
			// both feasible, need to copy
			//cout << 11 << endl;
			BBNode* b0 = new BBNode(b, v, 0);
			add(b0);

			b->inst->f1[v] = true;
			b->inst->T[v] = true;
			b->inst->p[v] = WMAX;
			b->v = v;
			b->bdir = 1;
			b->depth++;
			b->processed = false;
			//cout << 21 << endl;
			add(b);
			//cout << 31 << endl;

		}
		else if (b0feas) {
			//cout << 12 << endl;
			b->inst->removeNode(v);
			b->v = v;
			b->bdir = 0;
			b->depth++;
			b->processed = false;
			//cout << 22 << endl;
			add(b);
			//cout << 32 << endl;
		}
		else if (b1feas) {
			//cout << 13 << endl;
			b->inst->f1[v] = true;
			b->inst->T[v] = true;
			b->inst->p[v] = WMAX;
			b->v = v;
			b->bdir = 1;
			b->depth++;
			b->processed = false;
			//cout << 23 << endl;
			add(b);
			//cout << 33 << endl;
		}
		else {
			// both pruned, delete later
			b->feas = false;
		}

		// one branch
	}
	else {
		cout << "b->n " << b->n << endl;
		cout << "b->m " << b->m << endl;
		cout << "b->nfree " << b->nfree << endl;
		EXIT("error: no node to branch found\n");
	}
	//cout << "after branch " << endl;
}

int BBTree::selectBranchVariable(BBNode* b)
{
	Inst& inst = (*b->inst);
	int dmax = -1, dmaxS = -1, priomax = -1, v = -1;
	weight_t pmax = 0.0;

//	tbb::parallel_for(tbb::blocked_range<int>(0, inst.n),
		//[&](tbb::blocked_range<int> r)
	//{
		//for (int i = r.begin(); i < r.end(); ++i)
		//{
			for (int i = 0; i < inst.n; i++) {

			if (inst.f1[i] || inst.f0[i]) continue;

			int deg = 0;
			for (int ij : inst.din[i])  if (cr[ij] <= params.dasat) deg++;
			for (int ij : inst.dout[i]) if (cr[ij] <= params.dasat) deg++;

			int degS = 0;
			for (int ij : inst.din[i])  if (inc.arcs[ij]) degS++;
			for (int ij : inst.dout[i]) if (inc.arcs[ij]) degS++;

			const int prio1 = prio[i];

			if (params.branchtype == 0) {
				if ((priomax < prio1) ||
					(priomax == prio1 && deg > dmax) ||
					(priomax == prio1 && deg == dmax && degS > dmaxS)) {
					priomax = prio1;
					dmax = deg;
					dmaxS = degS;
					v = i;
				}
			}
			else if (params.branchtype == 1) {
				if ((deg > dmax) ||
					(deg == dmax && degS > dmaxS)) {
					priomax = prio1;
					dmax = deg;
					dmaxS = degS;
					v = i;
				}
			}
			else if (params.branchtype == 2) {
				if ((deg > dmax)) {
					priomax = prio1;
					dmax = deg;
					dmaxS = degS;
					v = i;
				}
			}
			else if (params.branchtype == 3) {
				if ((degS > dmaxS)) {
					priomax = prio1;
					dmax = deg;
					dmaxS = degS;
					v = i;
				}
			}
			}
		//}
	//});

	b->v2 = v;

	return v;
}

int BBTree::preprocess(Inst& inst)
{
	vector<flag_t> ap;
	vector<int> lastap;
	int riter = 0, removed_total = 0;

	do {
		riter = 0;

		costShift(inst);

		if (inst.r != -1 && !bRecover) {
			if (params.d1) riter += ntd1(inst);
			if (params.d2) riter += ntd2(inst);
		}

		if (!bRecover) {
			if (params.ma) riter += MA(inst);
		}

		if (inst.r != -1) {
			inst.AP(ap, lastap);
			riter += APfixing(inst, ap, lastap);

			if (!bRecover) {
				if (params.ms) riter += MAcutnode(inst, ap, lastap);

				inst.AP(ap, lastap);
				if (params.ss) riter += MAcutarc(inst, ap, lastap);
			}
		}

		if (params.lc) riter += lc(inst);

		removed_total += riter;

	} while (riter > 0);

	return removed_total;
}

void BBTree::freeOpenNodes()
{
	//boost::mutex::scoped_lock lock(mutexLock);
	while (!PQmin.empty()) { // !isPQMinEmpty()

		auto b = PQmin.top().second;
		//std::pair<weight_t, BBNode*> pair;
		//popPQMin(pair);
		PQmin.pop();
		delete b->inst;
		b->inst = 0;
		delete b;
	}
	PQmax.clear();
	PQmin.clear();
}

vector<int> BBTree::sortedListPotentialRoots()
{
	// gather potential roots
	vector<int> roots;
	if (inst.r != -1) {
		roots.push_back(inst.r);
	}
	else {
		for (int i = 0; i < inst.n; i++) {
			if (!inst.isAsym && !inst.T[i]) continue;
			roots.push_back(i);
		}
		sort(roots.begin(), roots.end(), [&](int i, int j) -> bool { return inst.p[i] > inst.p[j]; });
	}
	return roots;
}

void BBTree::fixTerm(Inst& inst, int t, vector<int>& fe0)
{
	for (int ij : inst.din[t]) {
		inst.c[ij] = WMAX;
		fe0.push_back(ij);
	}
	for (int ij : inst.dout[t]) {
		inst.c[ij] = WMAX;
		fe0.push_back(ij);
	}

	inst.f1[t] = false;
	inst.f0[t] = true;
	inst.T[t] = false;

	inst.offset += inst.p[t];
	inst.p[t] = 0;
}

BBNode* BBTree::select()
{
	//boost::mutex::scoped_lock lock(mutexLock);
	BBNode* b = NULL;
	switch (params.nodeselect) {
	case 0: // worst-bound
	{
		b = PQmax.top().second;
		PQmax.erase(b->pqposMax);
		PQmin.erase(b->pqposMin);
		
		break;
	}
	case 1:	// dfs
	{
		b = PQmax.top().second;
		PQmax.erase(b->pqposMax);
		PQmin.erase(b->pqposMin);
		break;
	}
	case 2: // best-bound
	{
		b = PQmin.top().second;
		PQmin.erase(b->pqposMin);
		break;
	}
	}

	return b;
}

void BBTree::evalLeaf(BBNode* b)
{
	Inst& inst = *b->inst;
	if (!isFeas(inst))
		return;

	vector<int> todel;
	for (int ir : inst.din[inst.r]) {
		todel.push_back(ir);
	}
	for (int ir : todel) {
		inst.delArc(ir);
		inst.fe0[ir] = true;
	}

	if (!isFeas(inst))
		return;

	Sol sol = dmst(inst, cr);

	updatePrimal(inst, sol);
}

int BBTree::initPrep()
{
	preprocess(inst);
	if (inst.r == -1 && params.semiBigM) {
		instM = inst.createRootedBigMCopy();
		lbM = daR(instM.r, instM, instM.c, crM, piM, WMAX, params.daeager, nullptr);

		stats.lb = lbM;
		bestlb = lbM;
		rootlb = lbM;
	}
	return 0;
}

bool BBTree::isFeas(Inst& inst, bool bDoNRtest)
{
	int r;
	int f1 = 0;
	for (int i = 0; i < inst.n; i++) {
		if (!inst.f1[i]) continue;
		f1++;
	}
	if (inst.bigM > 0) {
		bool reachAll = false;
		for (int ri : inst.dout[inst.r]) {
			int cnt = inst.countReachableFixed(inst.head[ri]);
			if (cnt == f1) {
				reachAll = true;
				break;
			}
		}

		if (!reachAll) {
			return false;
		}
	}
	if (inst.r != -1) {
		r = inst.r;
	}
	else {
		// if unrooted, start with an arbitrary fixed node
		r = -1;
		for (int i = 0; i < inst.n; i++) {
			if (!inst.f1[i]) continue;
			if (r == -1) {
				r = i;
				break;
			}
		}
		// cannot be infeasible if just one node is fixed to one
		if (f1 <= 1) return true;
	}

	stack<int> Q;
	vector<flag_t> visited(inst.n, false);
	int t = 0;
	for (int i = 0; i < inst.n; i++) {
		if (inst.f1[i]) t++;
	}

	t--;
	visited[r] = true;
	Q.push(r);

	while (!Q.empty()) {
		const int i = Q.top();
		Q.pop();

		for (int ij : inst.dout[i]) {
			const int j = inst.head[ij];
			if (inst.f0[j]) continue; // may have set some nodes to zero for testing, reachability after branching
			if (!visited[j]) {
				visited[j] = true;
				if (inst.f1[j]) t--;
				Q.push(j);
			}
		}
	}
	// apply non-reachability test if instance is feasible
	if (t == 0 && params.nr && inst.r != -1 && bDoNRtest) {
		for (int i = 0; i < inst.n; i++) {
			if (inst.f0[i] || visited[i]) continue;
			inst.removeNode(i);
			stats.nr++;
		}
	}

	return (t == 0);
}

