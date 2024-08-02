#ifndef _LATTICE_H_
#define _LATTICE_H_
#include <string>
#include <map>
#include "spcfcs.h"

using namespace std;

// lattice vector in real space (or reciprocal space)
class LttVec {
	friend class SiteID;
	friend class Lattice;
private:
	Int a0;		// unit cell coordinate in a0 (or b0) axis
	Int a1;		// unit cell coordinate in a1 (or b1) axis
public:
	LttVec(const Int &a0_i = 0, const Int &a1_i = 0) : a0(a0_i), a1(a1_i) {}
	Int getX() { return a0; }
	Int getY() { return a1; }
};

// site identity
class SiteID {
	friend class LttVec;
	friend class Lattice;
private:
	LttVec lv;		// lattice vector, i.e., unit cell coordinate
	Int i;			// site number in one unit cell
public:
	SiteID(const Int &a0 = 0, const Int &a1 = 0, const Int &i_i = 0) : lv(a0, a1), i(i_i) {}
	explicit SiteID(const LttVec &lv_i, const Int &i_i = 0) : lv(lv_i), i(i_i) {}
	LttVec getlv() { return lv; }
	const LttVec &getlv() const { return lv; }
	Int	getid() const { return i; }
};

class Lattice {
public:
	Int n0;			// number of unitcells in a0 direction
	Int n1;			// number of unitcells in a1 direction
	Int nsc;			// number of sites in one unit cell
	Int nc;			// total number of unitcells
	Int ns;			// total number of sites
	VecReal va0, va1;		// real-space lattice basis vectors
	VecReal vb0, vb1;		// reciprocal lattice basis vectors
private:
	// initialize basis vectors va0, va1, vb0, vb1
	void init_basis_vectors() {
		// add if necessary
	}
public:
	Lattice(Int n0_i = 1, Int n1_i = 1, Int nsu_i = 1) : n0(n0_i), n1(n1_i), nsc(nsu_i), nc(n0 * n1), ns(nsc * nc) {
		if (n0 < 1 || n1 < 1 || nsc < 1 || nc != n0 * n1 || ns != nsc * nc) ERR(NAV5(n0, n1, nsc, nc, ns));
		init_basis_vectors();
	}

	//unit cell coordinate/unit cell index
	LttVec lttvec(const Int &idx) const { return LttVec(idx % n0, idx / n0); }
	Int index(const LttVec &lv) const { return lv.a0 + lv.a1 * n0; }

	LttVec minus(const LttVec &lv) const { return LttVec(-lv.a0, -lv.a1); }
	LttVec add(const LttVec &lv0, const LttVec &lv1) const {
		return LttVec((lv0.a0 + lv1.a0) % n0, (lv0.a1 + lv1.a1) % n1);
	}
	LttVec sub(const LttVec &lv0, const LttVec &lv1) const {
		return LttVec((lv0.a0 - lv1.a0 + n0) % n0, (lv0.a1 - lv1.a1 + n1) % n1);
	}

	//unit cell index/unit cell index
	Int minus(const Int &idx) const { return index(minus(lttvec(idx))); }
	Int add(const Int &idx0, const Int &idx1) const { return index(add(lttvec(idx0), lttvec(idx1))); }
	Int sub(const Int &idx0, const Int &idx1) const { return index(sub(lttvec(idx0), lttvec(idx1))); }

	//site index/site id
	SiteID siteid(const Int &idx) const { return SiteID(lttvec(idx % nc), idx / nc); }
	Int index(const SiteID &sid) const { return index(sid.lv) + sid.i * nc; }

	//site index/site id
	//20180313 linheyu change
	//SiteID siteid(const Int &idx) const { return SiteID(lttvec(idx / nsc), idx % nsc); }
	//Int index(const SiteID &sid) const { return index(sid.lv) * nsc + sid.i; }

	SiteID add(const SiteID &sid, const LttVec &lv) const { return SiteID(add(sid.lv, lv), sid.i); }
	SiteID sub(const SiteID &sid, const LttVec &lv) const { return SiteID(sub(sid.lv, lv), sid.i); }

	Int left(const Int &idx) const { return index(sub(lttvec(idx), lttvec(1))); }
	Int right(const Int &idx) const { return index(add(lttvec(idx), lttvec(1))); };
	Int top(const Int &idx) const { return index(sub(lttvec(idx), lttvec(n0))); };
	Int buttom(const Int &idx) const { return index(add(lttvec(idx), lttvec(n0))); };

	LttVec left(const LttVec &lv) const { return sub(lv, lttvec(1)); }
	LttVec right(const LttVec &lv) const { return add(lv, lttvec(1)); };
	LttVec top(const LttVec &lv) const { return sub(lv, lttvec(n0)); };
	LttVec buttom(const LttVec &lv) const { return add(lv, lttvec(n0)); };

	string vecchar2string(OccCfg cfg) const
	{
		string res;
		res.insert(res.begin(), cfg.begin(), cfg.end());

		return res;
	}
};

#endif /* _LATTICE_H_ */
