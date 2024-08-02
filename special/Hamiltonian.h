#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_
#include "spcfcs.h"
#include "pars.h"
#include "Lattice.h"

class Model {
public:
	Int sgmz(const Int &occ) const { return 2 * occ - 1; }			// map from occupancy to pauli matrix sigma_z

	Real Sz(const Int &occ) const { return occ - 0.5; }			    // map from occupancy to spin matrix s_z
	Int occ(const Int &sgm) const { return (1 + sgm) >> 1; }	// inverse map from occupancy to pauli matrix sigma_z

	virtual void find_conn(const OccCfg &cfg, VEC<Real> &elem, VEC<VEC<Int> > &conn) const = 0;
	virtual ~Model() {}
};


class Heisenberg2d : public Model {
public:
	const Int ns;		// number of spins
	const Lattice ltt;	// lattice
	Real J;			// longitudinal coupling constant
	VEC<VEC<Int> > bonds;

	//20180528
	Real Jp;
	VEC<VEC<Int> > next_bonds;

	//20230518
	Real isPBC;


public:
	Heisenberg2d(const Lattice &ltt_i, const Pars &p_i)
		: ltt(ltt_i)
		, ns(ltt_i.ns)
		, J(1.0 * 0.25)
		, Jp(p_i.J2 * 0.25)
		, isPBC(p_i.isPBC)
	{
		get_all_bonds();
	}
	void get_all_bonds();
	//void find_conn(const OccCfg &cfg, VEC<Real> &elem, VEC<VEC<Int> > &conn) const;
	void find_conn(const OccCfg &cfg, VEC<Real> &elem, VEC<VEC<Int> > &conn) const;

};
#endif /* _HAMILTONIAN_H_ */
