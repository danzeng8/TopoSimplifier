/**
 * \file   options.h
 * \brief  handles program options
 *
 * \author Martin Luipersbeck 
 * \date   2015-05-03
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <string>

class ProgramOptions
{
public:
	struct Parameters
	{
		/// input
		std::string file;
		std::string solfile;
		std::string boundsfile;
		int         seed;

		// output
		std::string soloutfile;
		std::string statsfile;
		bool printstatsline;

		// problem type
		std::string type;
		long        precision;

		double timelimit;
		int    nodelimit;
		int    memlimit;

		// branch-and-bound 
		double cutoff;
		bool   cutoffopt;
		double absgap;
		int bbinfofreq;

		// dual ascent
		double daeager;
		double dasat;
		bool   daguide;

		// heuristic
		int    heurroots;
		bool   heursupportG;
		bool   heurbb;
		double heurbbtime;
		double heureps;

		// enable/disable components
		bool initprep;
		bool initheur;
		bool heuronly;
		bool rootonly;
		bool eagerroot;
		bool bigM;
		bool semiBigM;
		bool redrootonly;

		// b&b
		int  daiterations;
		bool perturbedheur;
		int  nodeselect;
		int  branchtype;
		bool lastcomp;

		// reductions
		bool d1;
		bool d2;
		bool ma;
		bool ms;
		bool ss;
		bool lc;
		bool nr;
		bool boundbased;
	};

	ProgramOptions(int &argc, char ** &argv);
	virtual ~ProgramOptions();

};

extern ProgramOptions::Parameters params;

#endif /* OPTIONS_H_ */
