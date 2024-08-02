/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _SPCFCS_H_
#define _SPCFCS_H_

#include <vector>
#include <string>

#ifdef _MSC_VER
#include "..\\general\\toolbox.h"
#include "..\\general\\ludcmp.h"
#include "..\\general\\frprmn.h"
#else
#include "../general/toolbox.h"
#include "../general/ludcmp.h"
#include "../general/frprmn.h"
#endif

// directory prefix
extern Str iox;		//        io directory prefix
extern Str bix;		// binary io directory prefix

// directory
#ifdef _MSC_VER
const Str iopfx = "io\\";	//        io directory prefix
const Str bipfx = "bi\\";	// binary io directory prefix
const Str dirslash = "\\";					// directory slash
#else
const Str iopfx = "io/";		//        io directory prefix
const Str bipfx = "bi/";		// binary io directory prefix
const Str dirslash = "/";		// directory slash
#endif
inline void iox_exist()
{
	if (!direxist(iox)) ERR(NAV(iox) + " does not exist!");
}
inline void io_init()
{
	iox_exist();
}

typedef Real Hnf;				// number field for Hamiltonian matrix
typedef Cmplx Wnf;				// number field for wave function, density matrix, etc
typedef Vec<Hnf> VecHnf;
typedef Mat<Hnf> MatHnf;
typedef Vec<Wnf> VecWnf;
typedef Mat<Wnf> MatWnf;
typedef Vec<Char> OccCfg;

inline OccCfg cfg_flip(const OccCfg &cfg_i, const VEC<Int> &flips) {
	OccCfg cfg_new = cfg_i;
	for_Int(i, 0, flips.size()) cfg_new[flips[i]] ^= 1;
	return cfg_new;
}

//20190619
inline Cmplx random_weight(UniformUnitRandom &rng, Real sigma)
{
	return Cmplx(rng() - 0.5, rng() - 0.5) * sigma;
}

inline Real MUL(const VecReal& vec)
{
	Real res = 1.;
	for_Int(i, 0, vec.size())
	{
		res *= vec[i];
	}
	return res;
}

inline std::string OccCfg2String(const OccCfg& cfg)
{
	std::string res;
	for_Int(i, 0, cfg.size())
	{
		res += std::to_string(int(cfg[i]));
	}
	return res;
}

#endif /* _SPCFCS_H_ */
