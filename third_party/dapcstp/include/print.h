/**
 * \file   print.h
 * \brief  terminal output for branch-and-bound tree
 *
 * \author Martin Luipersbeck 
 * \date   2015-10-03
 */

#include "procstatus.h"

void BBTree::printHeader()
{
	if(!bOutput)
		return;

	printf("%7s %6s %6s %s", "nodes", "depth", "open", "d");
	printf("%13s %13s %13s ", "nodelb", "lb", "ub");
	printf("%8s %8s ", "n", "m");
	printf("%7s  %8s", "action", "gap");
	printf("\n");
}

void BBTree::printBBLine(BBNode* b, NodeState state, bool bExit)
{
	// print every bbinfofreq-th, first and last line
	if(!bOutput || (nIter % params.bbinfofreq != 0 && !bExit && nIter != 1))
		return;

	weight_t prize = b->v == -1 ? -1 : b->inst->p[b->v];
	if(prize >= WMAX) prize = -1;

	double outUB = ub;
	double outLBNode = b->lb;
	double outLB = bestlb;
	if(b->inst->isMWCS) {
		outUB = inst.convertPCSTPBound2MWCS(outUB);
		outLBNode = inst.convertPCSTPBound2MWCS(outLBNode);
		outLB = inst.convertPCSTPBound2MWCS(outLB);
	}
	if(!b->inst->isInt) {
		outUB /= params.precision;
		outLB /= params.precision;
		outLBNode /= params.precision;
	}

	printf("%7d %6d %6d %d", nIter, b->depth, (int)PQmin.size(), b->bdir);

	if(inst.isInt) printf("%13.0lf %13.0lf %13.0lf ", (double)outLBNode, (double)outLB, (double)outUB);
	else           printf("%13.6lf %13.6lf %13.6lf ", (double)outLBNode, (double)outLB, (double)outUB);
	
	printf("%8d %8d ", b->n, b->m);

	switch(state) {
		case BB_INFEAS: printf("(infeas) "); break;
		case BB_CUTOFF: printf("(cutoff) "); break;
		case BB_LEAF:   printf("(leaf  ) "); break;
		case BB_BRANCH: printf("(branch) "); break;
	}
	printf("%8.5lf %%\n", gapP(bestlb, ub));
}

void BBTree::printRootHeader()
{
	if(!bOutput) return;

	printf("%7s %15s %15s %8s %8s %7s %8s\n", "root", "lb", "ub", "n", "m", "gap", "memory");
}

void BBTree::printRootLine(BBNode* b)
{
	if(!bOutput) return;

	printf("%7d ", b->inst->r);
	printBoundPadded(inst, b->lb);
	printf(" ");
	printBoundPadded(inst, ub);
	printf(" %8d %8d", b->n, b->m);
	printf(" %7.3lf", gapP(format(b->lb, inst), format(ub, inst)));
	printf(" %8d", ProcStatus::mem());
	printf("\n");
}

void BBTree::printHeurHeader()
{
	printf("     ");
	printf(" %15s", "sol");
	if(bestKnown >= 0)
		printf(" %8s    ", "Pgap");
	printf(" %15s", "best");
	if(bestKnown >= 0)
		printf(" %8s    ", "Pgap");
	printf("   %5s    ", "time");
	printf("\n");
}

void BBTree::printHeurLine(int it, weight_t obj, bool bImproved, double time)
{
	printf(" %c", bImproved ? '*' : ' ');
	printf(" %2d:", (it+1));
	printBoundPadded(inst, obj);

	if(bestKnown >= 0)
		printf("   %6.2lf %%  ", gapP(bestKnown, format(obj, inst)));

	printf(" ");
	printBoundPadded(inst, ub);
	if(bestKnown >= 0)
		printf("   %6.2lf %%  ", gapP(bestKnown, format(ub, inst)));

	printf(" [ %s%5.1lf s%s ]", GRAY, time, NORMAL);
	printf("\n");
}

void BBTree::printHeur1Summary()
{
	if(!bOutput)
		return;

	printf("[ %sheur 1%s ] [ %s%5.1lf s%s ] ", GREEN, NORMAL, GRAY, heurTime, NORMAL);

	if(inst.isInt)
		printf("ub %.0lf", format(ub, inst));
	else
		printf("ub %.6lf", format(ub, inst));
	printf("\n\n");
}

void BBTree::printHeur2Summary()
{
	if(!bOutput)
		return;

	printf("[ %sheur 2%s ] [ %s%5.1lf s%s ] ", GREEN, NORMAL, GRAY, heurBBTime, NORMAL);

	if(inst.isInt)
		printf("ub %.0lf", format(ub, inst));
	else
		printf("ub %.6lf", format(ub, inst));
	printf("\n\n");
}

void BBTree::printRootSummary()
{
	if(!bOutput)
		return;

	printf("[ %sroot%s   ] [ %s%5.1lf s%s ] ", GREEN, NORMAL, GRAY, rootTime, NORMAL);
	printf("gap %7.4lf %% %d/%d/%d ", gapP(rootlb, rootub), nRoots, nRootsProcessed, nRootsOpen);
	if(inst.isInt)
		printf("lb %.0lf ub %.0lf ", format(rootlb, inst), format(rootub, inst));
	else
		printf("lb %.6lf ub %.6lf ", format(rootlb, inst), format(rootub, inst));
	printf("avgn %5.1lf avgm %5.1lf\n\n", stats.rootavgn, stats.rootavgm);
}

void BBTree::printBBSummary()
{
	if(!bOutput)
		return;

	printf("[ %sbb%s     ] [ %s%5.1lf s%s ] ", GREEN, NORMAL, GRAY, bbTime, NORMAL);
	printf("gap %7.4lf %% ", gapP(bestlb, ub));
	if(inst.isInt)
		printf("lb %.0lf ub %.0lf ", format(bestlb, inst), format(ub, inst));
	else
		printf("lb %.6lf ub %.6lf ", format(bestlb, inst), format(ub, inst));
	printf("bbnodes %d", nIter);
	printf("\n\n");
}
