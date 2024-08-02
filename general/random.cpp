/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2016 - 2017
*/

#include "random.h"
#include "vec.h"
#include "iofmt.h"

const double UniformUnitRandomMine::rand_unit_dec = RAND_MAX + 1.;
const double UniformUnitRandomMine::inv_rand_unit_max = 1. / (rand_unit_dec * rand_unit_dec * rand_unit_dec);

// find period of rand() in C
// we find that the period of rand() in C is 2147483648
void random_number_test_period_rand()
{
	using namespace std;
	UInt r0 = rand();
	UInt r1 = rand();
	UInt r2 = rand();
	cout << "random number test" << endl;
	cout << "the three consective random numbers are" << endl;
	cout << r0 << endl;
	cout << r1 << endl;
	cout << r2 << endl;
	UInt d0 = r1;
	UInt d1 = r2;
	UInt d2 = rand();
	LLInt period = 1;
	LLInt period_max = 4400000000LL;
	while (period < period_max && (r0 != d0 || r1 != d1 || r2 != d2)) {
		d0 = d1;
		d1 = d2;
		d2 = rand();
		++period;
	}
	if (period < period_max) {
		cout << "the period of rand() in C is " << period << endl;
	} else {
		cout << "the period of rand() in C is larger than " << period_max << endl;
	}
}


// uniform real distribution in [0, 1]
// nrolls is number of experiments
// nintervals is number of intervals
// test rand_unit(Real)
void random_number_test(const Str &name, UniformUnitRandom &rng, Int nrolls, Int nintervals)
{
	using namespace std;
	cout << "uniform random [0, 1) by " << name << ": (nrolls = " << setw(10) << nrolls << ", nintervals = " << setw(4) << nintervals << ")";
	VecInt p(nintervals, 0);
	Real min = +1.E+128;
	Real max = -1.E+128;
	clock_t t = clock();
	Real r0;
	Real r1 = rand_unit(rng, r1);
	Real r2 = rand_unit(rng, r2);
	Int ntriple = 0;
	for_Int (i, 0, nrolls) {
		Real r = rand_unit(rng, r);
		++p[Int(nintervals * r)];
		if (min > r) min = r;
		if (max < r) max = r;
		r0 = r1;
		r1 = r2;
		r2 = r;
		if (r0 < 0.1 && 0.5 < r1 && r1 < 0.6 && 0.9 < r2) ++ntriple;
	}
	cout << fixed;
	cout << "    time = " << setprecision(3) << setw(8) << seconds_since(t) << "s";
	Real ideal_histgram = 1. * nrolls / nintervals;
	VecReal deviation(nintervals);
	for_Int (i, 0, nintervals) deviation[i] = (p[i] - ideal_histgram) / ideal_histgram;
	cout << "    DEV(dev) = " << setprecision(3) << setw(8) << DEV(deviation) * sqrt(Real(nrolls));
	cout << "    DEV(ntriple) = " << setprecision(3) << setw(8) << (ntriple - nrolls * 0.001) / (nrolls * 0.001) * sqrt(Real(nrolls));
	cout << "    log10(min, 1-max) = (" << setprecision(3) << setw(7) << log10(min) << ", " << setw(7) << log10(1 - max) << ")" << endl;
	cout << iofmt("def");
}
