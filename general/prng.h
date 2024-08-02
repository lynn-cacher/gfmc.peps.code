/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2017-09-05
*/

#ifndef _PRNG_H_
#define _PRNG_H_

#include "stdpfx.h"
#include "mympi.h"
#include "random.h"

// parallel random number engine provider
class PRng : public UniformUnitRandom {
private:
	const MyMpi &mm;
	UniformUnitRandomSFMT rng;
	double uniform_unit_random() { return rng(); }
public:
	PRng(const MyMpi &mm_i) : mm(mm_i) {
		int clck = 0;		// you might want to use this: Int clck = clock();
		int seed = para_seed() + clck;
		rng.init_random(seed);
	}
	void init_random(int seed) { ERR("PRng can not be initialized freely; " + NAV(seed)); }
private:
	inline int para_seed();
};

int PRng::para_seed()
{
	CRandomSFMT sfmt(1);
	int shift = mm.id();
	for_Int(i, 0, shift) {
		std::rand();
		std::rand();
		std::rand();
		sfmt.IRandom(0, 2147483647);
	}
	double rand_a = std::rand() * std::rand() * std::rand() * (pi_Real / 3);
	double rand_b = sfmt.IRandom(0x80000000, 0x7fffffff) * golden_double;
	return int(rand_a) + int(rand_b);
}

#endif /* _PRNG_H_ */
