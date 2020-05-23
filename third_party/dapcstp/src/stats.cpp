/**
 * \file   stats.cpp
 * \brief  stores and writes statistics
 *
 * \author Martin Luipersbeck 
 * \date   2015-10-03
 */

#include "stats.h"
#include "options.h"
#include <stdio.h>
#include <boost/filesystem.hpp>

ProgramStats::Stats stats;

void ProgramStats::writeStats(const char* file)
{
	FILE* fp;
	if((fp=fopen(file, "w")) == NULL)
		EXIT("error writing stats: %s\n", file);

	string dirname = boost::filesystem::path(params.file).parent_path().filename().string();

	// instance properties
	fprintf(fp, "%s;", stats.name.c_str());
	fprintf(fp, "%s;", dirname.c_str());
	fprintf(fp, "%d;", (int)stats.isInt);
	fprintf(fp, "%d;", (int)stats.isAsym);
	fprintf(fp, "%d;", stats.initial.n);
	fprintf(fp, "%d;", stats.initial.m);
	fprintf(fp, "%d;", stats.initial.t);
	fprintf(fp, "%d;", stats.initial.tr);
	fprintf(fp, "%d;", stats.initial.f1);

	// initial preprocessing
	fprintf(fp, "%d;", stats.prep.n);
	fprintf(fp, "%d;", stats.prep.m);
	fprintf(fp, "%d;", stats.prep.t);
	fprintf(fp, "%d;", stats.prep.tr);
	fprintf(fp, "%d;", stats.prep.f1);
	fprintf(fp, "%.3lf;", stats.preptime);

	// root node
	fprintf(fp, "%.6lf;", stats.rootlb);
	fprintf(fp, "%.6lf;", stats.rootub);
	fprintf(fp, "%.6lf;", stats.rootgap);
	fprintf(fp, "%d;",    stats.roots);
	fprintf(fp, "%d;",    stats.proots);
	fprintf(fp, "%d;",    stats.oroots);
	fprintf(fp, "%.3lf;", stats.rootavgn);
	fprintf(fp, "%.3lf;", stats.rootavgm);
	fprintf(fp, "%.3lf;", stats.rootavgt);
	fprintf(fp, "%.3lf;", stats.rootavgtr);
	fprintf(fp, "%.3lf;", stats.rootavgf1);
	fprintf(fp, "%.3lf;", stats.roottime);

	// B&B
	fprintf(fp, "%d;",    stats.bbnodes);
	fprintf(fp, "%.6lf;", stats.lb);
	fprintf(fp, "%.6lf;", stats.ub);
	fprintf(fp, "%.6lf;", stats.gap);

	// prep
	fprintf(fp, "%d;", (int)stats.d1);
	fprintf(fp, "%d;", (int)stats.d2);
	fprintf(fp, "%d;", (int)stats.ma);
	fprintf(fp, "%d;", (int)stats.ms);
	fprintf(fp, "%d;", (int)stats.ss);
	fprintf(fp, "%d;", (int)stats.lc);
	fprintf(fp, "%d;", (int)stats.nr);
	fprintf(fp, "%d;", (int)stats.boundbased);

	fprintf(fp, "%.3lf;", stats.heurtime);
	fprintf(fp, "%.3lf;", stats.heurbbtime);
	fprintf(fp, "%.3lf;", stats.timeBest);
	fprintf(fp, "%.3lf;", stats.time);
	
	fprintf(fp, "%d;", stats.memout);

	// paramters
	fprintf(fp, "%d;", (int)params.heursupportG);
	fprintf(fp, "%d;", (int)params.heurbb);
	fprintf(fp, "%.6lf;", params.heureps);
	fprintf(fp, "%d;", (int)params.heurroots);
	fprintf(fp, "%d;", (int)params.heurbbtime);
	fprintf(fp, "%d;", (int)params.perturbedheur);
	fprintf(fp, "%d;", (int)params.daiterations);
	fprintf(fp, "%d;", (int)params.nodeselect);
	fprintf(fp, "%d;", (int)params.branchtype);
	fprintf(fp, "%d;", (int)params.daguide);

	fprintf(fp, "%d;", (int)params.d1);
	fprintf(fp, "%d;", (int)params.d2);
	fprintf(fp, "%d;", (int)params.ma);
	fprintf(fp, "%d;", (int)params.ms);
	fprintf(fp, "%d;", (int)params.ss);
	fprintf(fp, "%d;", (int)params.lc);
	fprintf(fp, "%d;", (int)params.nr);
	fprintf(fp, "%d",  (int)params.boundbased);
	fprintf(fp, "\n");

	fclose(fp);
}

void ProgramStats::initRootNodeStats()
{
	stats.rootavgn = 0; 
	stats.rootavgm = 0;
	stats.rootavgt = 0;
	stats.rootavgtr = 0;
	stats.rootavgf1 = 0;
}

void ProgramStats::addRootNodeStats(InstSizeData& sdata)
{
	stats.rootavgn += sdata.n;
	stats.rootavgm += sdata.m;
	stats.rootavgt += sdata.t;
	stats.rootavgtr += sdata.tr;
	stats.rootavgf1 += sdata.f1;
}

void ProgramStats::averageRootNodeStats(int nRootsOpen)
{
	stats.rootavgn /= nRootsOpen;
	stats.rootavgm /= nRootsOpen;
	stats.rootavgt /= nRootsOpen;
	stats.rootavgtr /= nRootsOpen;
	stats.rootavgf1 /= nRootsOpen;
}
