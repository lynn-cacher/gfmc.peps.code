/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _MYOMP_H_
#define _MYOMP_H_

#include "stdpfx.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mymkl.h"

inline void use_omp(Int nthreads = 1)
{
	using namespace std;
#ifdef _OPENMP
	omp_set_num_threads(nthreads);
	int omp_num_threads;
#pragma omp parallel
	omp_num_threads = omp_get_num_threads();
	cout << "OpenMP: number of procs = " << omp_get_num_procs() << endl;
	cout << "OpenMP: omp_num_threads = " << omp_num_threads << endl;
#endif
}

inline void use_mkl(Int nthreads = 1)
{
	using namespace std;
	mkl_set_num_threads(nthreads);
	cout << "MKL: mkl_num_threads = " << mkl_get_max_threads() << endl;
}

#endif /* _MYOMP_H_ */
