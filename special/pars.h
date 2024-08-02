/*
code developed by
Chen Peng (whumai@gmail.com, RUC, China)
Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2017-08 - 2017-09-05
*/

#ifndef _PARA_H_
#define _PARA_H_
#include "spcfcs.h"

class Pars {
private:
	Int np;			// number of mpi processes
public:

	Real sigma;				// rbm parameters are randomly set to sigma * Cmplx(random - 0.5, random - 0.5);
	Int seed;				// seed for random number generator used in rbm parameters initialization

	// model related
	Int lx;					// number of unit cell along x-direction
	Int ly;					// number of unit cell along y-direction
	Int nsu;				// number of sites in each unit cell
	Int nf;					// number of all filters
	Int ifnec;				// if ne conserved
	Real hx;				// transverse field hx
	Real J2;
	Bool isPBC;				// TRUE: PBC. FALSE: OBC.

	// derived
	Int nc;					// total number of lattice unitcells
	Int ns;					// total number of lattice sites
	Int ne;					// number of electrons

							//20180516
	Int nf_spec;			// number of filters which need to be optimized

	//gfmc
	Real beta;
	Int	 nP;
	Int  nStep;
	Int  ntotalSample;

	Int Nw;
	Int Nw_per_prog;
	Real gamma;
	Real max_relnpsi;
	Int nbin;
	Bool ifreconfig;

	//mode: 0--sampling according gfmc, 1--reading samples from file.
	Int sampling_mode;
	Int if_savesamplings;

private:
	void set_inert_values();
	void set_values();
	void derive();
	void check_consistency();
public:
	Pars(const MyMpi &mm);
	void print(std::ostream &os = std::cout);
	void print(std::ostream &os, const Str &var, const Str &val, const Str &comment) {
		using namespace std;
		Str true_var = var == "\"\"" ? "" : var;
		os << rght_justify(true_var, 16) << (true_var == "" ? "   " : " = ") << left_justify(val, w_Real)
			<< "    # " + comment << endl;
	}
};

#endif /* _PARA_H_ */
