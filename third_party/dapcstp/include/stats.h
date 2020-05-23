/**
 * \file   stats.h
 * \brief  stores and writes statistics
 *
 * \author Martin Luipersbeck 
 * \date   2015-10-03
 */

#ifndef STATS_H_
#define STATS_H_

#include <def.h>
#include <inst.h>

class ProgramStats {
public:

	struct Stats {
		// instance data
		string name;
		InstSizeData initial;
		double bidirect = 0.0;

		// initial prep data
		InstSizeData prep;
		double preptime = 0.0;

		// heur
		double heurtime = 0.0, heurbbtime = 0.0;

		// root data
		double rootlb = 0.0, rootub = -1, rootgap = 100.0;
		int roots = -1, proots = -1, oroots = -1;
		double rootavgn = -1, rootavgm = -1, rootavgt = -1, rootavgtr = -1, rootavgf1 = -1;
		double roottime = 0.0;

		// bb data
		int bbnodes = 0;
		double lb = 0.0, ub = -1, bbtime = 0.0, timeBest = 0.0, gap = 100.0;

		// total
		double time = 0.0;

		int d1 = 0, d2 = 0, ma = 0, ms = 0, ss = 0;
		int lc = 0, nr = 0, boundbased = 0;

		int memout = 0;
		bool isInt = false;
		bool isAsym = false;
		bool valid = false;
	};

	static void writeStats(const char* file);
	static void initRootNodeStats();
	static void addRootNodeStats(InstSizeData& sdata);
	static void averageRootNodeStats(int nRootsOpen);

};
extern ProgramStats::Stats stats;

#endif // STATS_H_
