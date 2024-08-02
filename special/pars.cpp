/*
code developed by
Chen Peng (whumai@gmail.com, RUC, China)
Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2017-08 - 2017-09-05
*/

#include "pars.h"

Pars::Pars(const MyMpi &mm) : np(mm.np())
{
	set_inert_values();
	set_values();
	derive();
	check_consistency();
	if (mm) print();
	// if (mm) { OFS ofs((iox + "info.txt").c_str(), std::ios::app);  print(ofs); }
}

void Pars::set_inert_values()
{
}

void Pars::set_values()
{
	// rbm related
	sigma = 0.05;			// default value: 0.05
	seed = 1;				// default value: 1

	//lattice
	lx = 4;					// default value: 40
	ly = 4;					// default value: 1
	nsu = 1;				// default value: 1
	nf = 2;					// default value: 2
	ifnec = 1;				// default value: 1

	//transverse Ising model
	hx = 1.;				// defualt value: 0.

	//2D Heisenberg
	J2 = 0.50;				// defualt value: 0.
	isPBC = true;

	//mode: 0--sampling according gfmc, 1--reading samples from file.
	sampling_mode = 0;
	if_savesamplings = true;

	//gfmc
	beta = 0.02;				//the evolution of time for the Continuous-Time limit process
	nP = 16;				//the number of reconfig processes. beta * nP should equal to 20~40.
	nStep = 50;			//the number of sample step

	ifreconfig = false;	//if false, then Nw must equal to the number of processes

	Nw = 1;				//the number of walkers
	ntotalSample = Nw * nStep;	//the number of total samples
	nbin = 25;
	
	gamma = 0.0;
	max_relnpsi = 0.;
}

void Pars::derive()
{
	// model related
	nc = lx * ly;				// total number of lattice unitcells
	ns = nsu * nc;				// total number of lattice sites
	ne = ns / 2;				// ne and sgm_z_t must be consistent

	Nw_per_prog = (Nw - 1) / np + 1;
	Nw = Nw_per_prog * np; //make sure that Nw is an integer multiple of np
	
	ntotalSample = Nw * nStep;
}

void Pars::check_consistency()
{
	// model related assertion
	if (lx < 1 || ly < 1 || nc != lx * ly || ne < 0 || ne > ns) WRN(NAV4(lx, ly, ns, ne));
}

void Pars::print(std::ostream &os)
{
#define para_print(var, comment) print(os, NAME(var), STR(var), comment)

	using namespace std;
	os << present() << endl;
	os << endl << "// constant" << endl << endl;
	para_print(np, "number of mpi processes");

	os << endl << "// model" << endl << endl;
	para_print(lx, "number of unit cell along x-direction");
	para_print(ly, "number of unit cell along y-direction");
	para_print(nsu, "number of sites in each unit cell");
	para_print(nf, "number of filters");
	para_print(ifnec, "if ne conserved"); 
	para_print(J2, "the nnn iteraction of 2D-Heisenberg model");
	para_print(isPBC, "Whether the model satisfies periodic boundary conditions");

	os << endl << "// derive" << endl << endl;
	para_print(nc, "total number of lattice unitcells");
	para_print(ns, "number of lattice sites");
	para_print(ne, "number of electrons");

	os << endl << "// gfmc" << endl << endl;
	para_print(sampling_mode, "the mode of sampling");
	para_print(beta, "the evolution of time for the Continuous-Time limit process");
	para_print(Nw, "the number of total walkers");
	para_print(Nw_per_prog, "the number of walkers for each progress");
	para_print(nP, "the number of reconfig processes");
	para_print(ifreconfig, "whether the program executes the reconfiguration process");
	para_print(nStep, "the number of sample step");
	para_print(nbin, "the number of Measured average value");
	para_print(ntotalSample, "the number of total samples");
	para_print(gamma, "the parameter in the fixed-node Hamiltonian");
	para_print(max_relnpsi, "the displacement of ln psi");

	os << endl;

#undef para_print
}
