/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2016 - 2017
*/

#ifndef _RANDOM_H_
#define _RANDOM_H_

#include "stdpfx.h"
#include "vec.h"
#ifdef _MSC_VER
#include <random>
#endif

#ifdef _MSC_VER
#include "..\\randomc\\randomc.h"
#include "..\\randomc\\sfmt.h"
#else
#include "../randomc/randomc.h"
#include "../randomc/sfmt.h"
#endif

// uniform unit random number [0, 1)
class UniformUnitRandom {
private:
	virtual double uniform_unit_random() = 0;
public:
	virtual void init_random(int seed) = 0;
	double operator()() {
		return uniform_unit_random();
	}
	virtual ~UniformUnitRandom() {}
};

// uniform unit random number [0, 1), mine
class UniformUnitRandomMine: public UniformUnitRandom {
private:
	static const double rand_unit_dec;
	static const double inv_rand_unit_max;
	double uniform_unit_random() {
		const Real &d = rand_unit_dec;
		return ((rand() * d + rand()) * d + rand()) * inv_rand_unit_max;
	}
public:
	UniformUnitRandomMine(int seed = 1) {
		init_random(seed);
		if (RAND_MAX < 32767) ERR("RAND_MAX < 32767, RAND_MAX = " + STR(RAND_MAX));
	}
	void init_random(int seed) { srand(seed); }
};

// uniform unit random number [0, 1), Mersenne
class UniformUnitRandomMersenne: public UniformUnitRandom {
private:
	CRandomMersenne rng;
	double uniform_unit_random() { return rng.Random(); }
public:
	UniformUnitRandomMersenne(int seed = 1): rng(seed) {}
	void init_random(int seed) { rng.RandomInit(seed); }
};

// uniform unit random number [0, 1), SFMT
class UniformUnitRandomSFMT: public UniformUnitRandom {
private:
	CRandomSFMT rng;
	double uniform_unit_random() { return rng.Random(); }
public:
	UniformUnitRandomSFMT(int seed = 1): rng(seed) {}
	void init_random(int seed) { rng.RandomInit(seed); }
};

#ifdef _MSC_VER
// uniform unit random number [0, 1), STL
class UniformUnitRandomSTL: public UniformUnitRandom {
private:
	std::default_random_engine generator;
	std::uniform_real_distribution<Real> distribution;
public:
	UniformUnitRandomSTL(int seed = 1): generator(seed), distribution(0., 1.) {}
	void init_random(int seed) { generator.seed(seed); }
	double uniform_unit_random() { return distribution(generator); }
};
#endif

// find a random Real number between [0, 1]
inline Real rand_unit(UniformUnitRandom &rng, const Real &r) { return rng(); }
// find a random Cmplx number in the unit square [0, I] * [0, 1]
inline Cmplx rand_unit(UniformUnitRandom &rng, const Cmplx &c) { return rng() + I * rng(); }

inline Real half_unit(const Real &r)  { return 0.5; }
inline Cmplx half_unit(const Cmplx &c) { return 0.5 + I * 0.5; }

template<typename T>
inline void rand_unit(UniformUnitRandom &rng, Vec<T> &v)
{
	for_Idx (i, 0, v.size()) v[i] = rand_unit(rng, v[i]);
}

void random_number_test_period_rand();
void random_number_test(const Str &name, UniformUnitRandom &rng, Int nrolls, Int nintervals);

#endif /* _RANDOM_H_ */
