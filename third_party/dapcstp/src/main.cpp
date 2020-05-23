/**
 * \file   main.cpp
 * \brief  solver for apcstp instances
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#include <stdio.h>
#include <vector>
#include <boost/filesystem.hpp>

#include "stats.h"
#include "options.h"
#include "procstatus.h"
#include "timer.h"
#include "ds.h"
#include "bbtree.h"
#include "prep.h"

Inst load();
void solve(Inst& inst);

double bestKnown = -1;

int main(int argc, char *argv[])
{
	enlargeStack();
	ProgramOptions po(argc, argv);
	ProcStatus::setMemLimit(params.memlimit);
	srand(params.seed);

	Inst inst = load();
	solve(inst);

	return 0;
}

Inst load()
{
	// load instance file
	if(params.file.empty()) {
		EXIT("Input file missing.\n");
	}

	Timer tLoad(true);
	Inst inst = load(params.file.c_str());
	if(params.bigM) {
		if(inst.r == -1) {
			inst = inst.createRootedBigMCopy();
		} else {
			params.bigM = false;
		}
	}

	stats.name = boost::filesystem::path(params.file).stem().string();

	if(!params.boundsfile.empty()) {
		bestKnown = getBestKnownBound(params.file.c_str(), params.boundsfile.c_str());
	}
	
	printf("[ %sload%s   ] [ %s%5.1lf s%s ] ", GREEN, NORMAL, GRAY, tLoad.elapsed().getSeconds(), NORMAL);
	printf("n %5d m %5d t %5d ", inst.n, inst.m, inst.t);
	printf("integer %d asym %d bidir %5.2lf ", inst.isInt, inst.isAsym, stats.bidirect);
	printf("( %s%s%s )", GREEN, stats.name.c_str(), NORMAL);

	printf(" ( best: %s", YELLOWBI);
	if(bestKnown < 0) {
		printf("N/A");
	} else if(inst.isInt)
		printf("%.0lf", bestKnown);
	else
		printf("%.6lf", bestKnown);
	printf("%s )\n\n", NORMAL);

	return inst;
}

void solve(Inst& inst)
{
	BBTree bbtree(inst);

	if(!params.solfile.empty()) {
		Sol start = loadSol(params.solfile.c_str(), inst);
		bbtree.setIncumbent(start);
	}

	if(params.initprep) {
		Timer tPrep(true);
		bbtree.initPrep();
		stats.prep = inst.countInstSize();
		stats.preptime = tPrep.elapsed().getSeconds();
		printf("[ %sprep%s   ] [ %s%5.1lf s%s ] ( %4.1lf %% )", GREEN, NORMAL, GRAY, stats.preptime, NORMAL, stats.prep.m*100.0/inst.m);
		// big-M is only relevant for unrooted instances
		if(params.semiBigM && inst.r == -1) {
			printf(" lbM ");
			if(inst.isInt)
				printf("%13.0lf", format(bbtree.getLBM(), inst));
			else
				printf("%13.6lf", format(bbtree.getLBM(), inst));
		}
		printf("\n\n");
	}

	if(params.cutoff > 0.0)
		bbtree.setCutUp(params.cutoff);
	if(params.cutoffopt) 
		bbtree.setCutUp(bestKnown);
	
	bbtree.setBestKnown(bestKnown);

	// solve
	if(params.heuronly) {
		bbtree.initHeur();
	} else if(params.rootonly) {
		bbtree.initHeur();
		
		if(params.timelimit >= 0)
			bbtree.setTimeLim(max(0.0,params.timelimit-Timer::total.elapsed().getSeconds()));
		if(params.nodelimit >= 0)
			bbtree.setNodeLim(params.nodelimit);
		bbtree.processRoots();
	} else {
		bbtree.initHeur();
		
		if(params.timelimit >= 0)
			bbtree.setTimeLim(max(0.0,params.timelimit-Timer::total.elapsed().getSeconds()));
		if(params.nodelimit >= 0)
			bbtree.setNodeLim(params.nodelimit);
		bbtree.solve();
	}

	stats.bbnodes = bbtree.getNnodes();
	stats.isInt   = inst.isInt;
	stats.isAsym  = inst.isAsym;

	// timing
	stats.bbtime     = bbtree.getTime();
	stats.timeBest   = bbtree.getTimeBest();
	stats.roottime   = bbtree.getRootTime();
	stats.heurtime   = bbtree.getHeurTime();
	stats.heurbbtime = bbtree.getHeurBBTime();

	// bounds
	stats.ub       = format(bbtree.getUB(), inst);
	stats.lb       = format(bbtree.getLB(), inst);
	stats.gap      = gapP(stats.lb, stats.ub);

	// root
	stats.rootub   = format(bbtree.getRootUB(), inst);
	stats.rootlb   = format(bbtree.getRootLB(), inst);
	stats.rootgap  = gapP(stats.rootlb, stats.rootub);
	stats.roots    = bbtree.getNroots();
	stats.oroots   = bbtree.getNrootsOpen();
	stats.proots   = bbtree.getNrootsProcessed();

	if(bbtree.getState() == BBTree::State::BB_MEMLIMIT)
		stats.memout = 1;

	// get running time and backmapped solution
	stats.time = Timer::total.elapsed().getSeconds();
	auto S = bbtree.getInc1();
	weight_t ub = S.obj;
	weight_t lb = bbtree.getLB();

	// check if value from bound file matches computed optimum value
	const bool match = (abs(format(ub, inst)-bestKnown) < 1e-5);
	// validates solution
	stats.valid = S.validate();
	
	//printf("v %6d e %6d t %6d tr %6d bb %6d ", stats.initial.n, stats.initial.m, stats.initial.t, stats.initial.tr, bbtree.getIter());
	printf("bbnodes  %15d\n", bbtree.getNnodes());
	printf("ub       ");
	printBoundPadded(inst, ub);
	printf("\n");
	printf("lb       ");
	printBoundPadded(inst, lb);
	printf("\n");
	printf("rootub   ");
	printBoundPadded(inst, bbtree.getRootUB());
	printf("\n");
	printf("rootlb   ");
	printBoundPadded(inst, bbtree.getRootLB());
	printf("\n");
	printf("gap      %15.3lf\n", bbtree.getGap());
	printf("gapR     %15.3lf\n", bbtree.getRootGap());
	printf("timeBest %15.1lf\n", stats.timeBest);
	printf("time     %15.1lf\n", stats.time);
	printf("matches  %15d\n", match);
	printf("valid    %15d\n", (int)stats.valid);

	printf("lc %5d d1 %5d d2 %5d ma %5d ms %5d ss %5d nr %5d bb %5d\n", stats.lc, stats.d1, stats.d2, stats.ma, stats.ms, stats.ss, stats.nr, stats.boundbased);

	// write output files (solution + stats)
	if(!params.statsfile.empty()) {
		ProgramStats::writeStats(params.statsfile.c_str());
	}

	if(!params.soloutfile.empty()) {
		writeSolution(params.soloutfile.c_str(), bbtree.getInst1(), bbtree.getInc1());
	}

	if(inst.isMWCS) {
		delete inst.transformation;
	}

	if(params.printstatsline)
		printf("STAT;%s;%d;%d;%d;%.6lf;%.6lf;%.6lf;%.3lf;%d;%d;%d\n", stats.name.c_str(), inst.n, inst.m, stats.bbnodes, stats.gap, stats.lb, stats.ub, stats.time, stats.valid, match, stats.memout);
}
