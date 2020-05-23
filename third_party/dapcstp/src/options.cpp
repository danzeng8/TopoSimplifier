/**
 * \file   options.cpp
 * \brief  handles program options
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#include "options.h"
#include <boost/program_options.hpp>
#include <iostream>

using namespace std;
namespace po = boost::program_options;

ProgramOptions::Parameters params;

ProgramOptions::ProgramOptions(int &argc, char ** &argv)
{
	po::options_description general_options("General options");
	general_options.add_options()
			("help,h", "produce help message")
			("file,f", po::value<string>(&params.file)->default_value(""), "instance file to process")
			("solout,o", po::value<string>(&params.soloutfile)->default_value(""), "solution file for output")
			("stats", po::value<string>(&params.statsfile)->default_value(""), "statistics file for output")
			("sol", po::value<string>(&params.solfile)->default_value(""), "solution file for starting solution")
			("bounds", po::value<string>(&params.boundsfile)->default_value(""), "bounds file for input")
			("precision", po::value<long>(&params.precision)->default_value(-1), "decimal precision read from file (-1: choose automatically 12 for mwcs and 6 for the rest)")
			("printstatsline", po::value<bool>(&params.printstatsline)->default_value(true)->implicit_value(true), "print line containing stats values for quick parsing")
			("type", po::value<string>(&params.type)->default_value("pcstp"), "instance problem type (pcstp|stp|mwcs|nwstp)")
			("seed", po::value<int>(&params.seed)->default_value(0), "random seed")
			("timelimit,t", po::value<double>(&params.timelimit)->default_value(10), "timelimit")
			("memlimit,m", po::value<int>(&params.memlimit)->default_value(10000000000), "memory limit")
			;

			// B&B parameters
	po::options_description bb_options("B&B options");
	bb_options.add_options()
			("bb.cutup", po::value<double>(&params.cutoff)->default_value(-1), "upper bound cutup for branch-and-bound")
			("bb.cutupopt", po::value<bool>(&params.cutoffopt)->default_value(false)->implicit_value(true), "choose cutup from bounds ficcNeighborFill6Connle (specify with --bounds pathtofile)")
			("bb.absgap", po::value<double>(&params.absgap)->default_value(0), "absolute optimality gap")
			("bb.infofreq", po::value<int>(&params.bbinfofreq)->default_value(100), "number of nodes after which the B&B status information is updated")
			("bb.branchtype", po::value<int>(&params.branchtype)->default_value(0), "branching type")
			("bb.nodeselect", po::value<int>(&params.nodeselect)->default_value(0), "node selection strategy")
			("bb.daiterations", po::value<int>(&params.daiterations)->default_value(10), "number of dual ascent iterations per B&B node (minimum: 1)")
			("bb.perturbedheur", po::value<bool>(&params.perturbedheur)->default_value(true)->implicit_value(true), "calls the primal heuristic on the support graph with perturbed cost (deactivated automatically if --heur.eps=0)")
			("bb.nodelimit,m", po::value<int>(&params.nodelimit)->default_value(-1), "node limit")
			;

			// dual ascent parameters
	po::options_description da_options("Dual ascent options");
	da_options.add_options()
			("da.eager", po::value<double>(&params.daeager)->default_value(1.25), "threshold below which an element of the DA priority queue is processed, even if its score is higher than the next element")
			("da.sat", po::value<double>(&params.dasat)->default_value(-1), "threshold below which an arc is viewed as saturated within the dual ascent algorithm")
			("da.guide", po::value<bool>(&params.daguide)->default_value(true)->implicit_value(true), "use guiding solutions")
			("da.lastcomp", po::value<bool>(&params.lastcomp), "lastcomp")
			;

			// enable/disable components
	po::options_description comp_options("Component options");
	comp_options.add_options()
			("rootonly", po::value<bool>(&params.rootonly)->default_value(false)->implicit_value(true), "only process root node of bb")
			("heuronly", po::value<bool>(&params.heuronly)->default_value(false)->implicit_value(true), "only process heuristic without lower bound")
			("initheur", po::value<bool>(&params.initheur)->default_value(true)->implicit_value(true), "perform initialization heuristic to compute starting solution")
			("initprep", po::value<bool>(&params.initprep)->default_value(true)->implicit_value(true), "perform initial preprocessing")
			("redrootonly", po::value<bool>(&params.redrootonly)->default_value(false)->implicit_value(true), "perform reduction tests only in root node")
			("bigM", po::value<bool>(&params.bigM)->implicit_value(true)->default_value(false), "transform instance to rooted one during whole procedure")
			("semiBigM", po::value<bool>(&params.semiBigM)->implicit_value(true)->default_value(false), "transform instance to rooted one during preprocessing by adding artificial root arcs, enables earlier reductions")
			("eagerroot", po::value<bool>(&params.eagerroot)->default_value(false)->implicit_value(true), "process roots eager when low on memory")
			;

			// heuristic parameters
	po::options_description heur_options("Heuristic options");
	heur_options.add_options()
			("heur.eps", po::value<double>(&params.heureps)->default_value(-1), "epsilon parameter used in perturbed construction heuristic (-1: choose automatically)")
			("heur.roots", po::value<int>(&params.heurroots)->default_value(1000), "number of roots for initial heuristics")
			("heur.bb", po::value<bool>(&params.heurbb)->default_value(true)->implicit_value(true), "heuristic that applies B&B on the support graphs created during the initialization heuristic and union of starting solutions")
			("heur.bbtime", po::value<double>(&params.heurbbtime)->default_value(10000), "time limit for b&b heuristic")
			("heur.supportG", po::value<bool>(&params.heursupportG)->default_value(true)->implicit_value(true), "apply shortest path heuristic on support graph computed by dual ascent")
			;

			// reduction tests
	po::options_description red_options("Reduction test options");
	red_options.add_options()
			("red.d1", po::value<bool>(&params.d1)->default_value(true)->implicit_value(true), "degree 1")
			("red.d2", po::value<bool>(&params.d2)->default_value(true)->implicit_value(true), "degree 2")
			("red.ma", po::value<bool>(&params.ma)->default_value(true)->implicit_value(true), "(asymmetric) minimum adjacency")
			("red.ms", po::value<bool>(&params.ms)->default_value(true)->implicit_value(true), "minimum successor")
			("red.ss", po::value<bool>(&params.ss)->default_value(true)->implicit_value(true), "single successor")
			("red.lc", po::value<bool>(&params.lc)->default_value(true)->implicit_value(true), "least cost")
			("red.nr", po::value<bool>(&params.nr)->default_value(true)->implicit_value(true), "non-reachability")
			("red.boundbased", po::value<bool>(&params.boundbased)->default_value(true)->implicit_value(true), "bound-based tests")
			;

	po::options_description all("Allowed parameters");
	all.add(general_options);
	all.add(bb_options);
	all.add(da_options);
	all.add(heur_options);
	all.add(red_options);
	all.add(comp_options);

	po::positional_options_description descPos;
	descPos.add("file", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all).positional(descPos).run(), vm);
	po::notify(vm);

	int p = params.precision;

	// autoselect precision for convenience
	if(p < 0) {
		if(params.type.compare("mwcs") == 0) {
			p = 12;
		} else {
			p = 6;
		}
	}

	params.precision = 1;
	for(int i = 0; i < p; i++) {
		params.precision *= 10;
	}

	if (vm.count("help")) {
	    cout << general_options << endl;
	    cout << bb_options << endl;
	    cout << da_options << endl;
	    cout << comp_options << endl;
	    cout << heur_options << endl;
	    cout << red_options << endl;
	    exit(0);
	}
	else if (!vm.count("file"))
	{
		cout << "No input file given!" << endl;
		exit(0);
	}
}

ProgramOptions::~ProgramOptions()
{

}

