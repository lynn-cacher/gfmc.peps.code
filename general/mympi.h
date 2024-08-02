/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _MYMPI_H_
#define _MYMPI_H_

#include "stdpfx.h"
#include "mpi.h"
#include "vec.h"
#include "myomp.h"
#include "safebool.h"

#define MYMPI_REAL MPI_DOUBLE_PRECISION
#define MYMPI_INT MPI_INT

//class VecPartition;

class VecPartition {
private:
	Int np;
	Idx n;
	VecIdx _bgn;
	VecIdx _end;
	VecIdx _len;
	VecIdx _dsp;
public:
	VecPartition(Int np_i, Idx n_i) : np(np_i), n(n_i), _bgn(np), _end(np), _len(np) {
		if (np < 1) ERR(NAV2(np, n));
		Idx k = n / np;
		Idx j = n % np;
		Idx m = 0;
		for_Int(i, 0, np) {
			_bgn[i] = m;
			m += k + (i < j);
			_len[i] = m - _bgn[i];
			_end[i] = m;
		}
		if (m != n) ERR(NAV3(np, n, m));
		_dsp = _bgn;
	}
	Idx bgn(const Int& id) const { return _bgn[id]; }
	Idx end(const Int& id) const { return _end[id]; }
	Idx len(const Int& id) const { return _len[id]; }
	Idx dsp(const Int& id) const { return _dsp[id]; }
	Idx operator[](const Int& id) const { return _len[id]; }
	const Idx* bgn() const { return _bgn.p(); }
	const Idx* end() const { return _end.p(); }
	const Idx* len() const { return _len.p(); }
	const Idx* dsp() const { return _dsp.p(); }
	Idx vec_size() const { return n; }
};

class MyMpi: public safe_bool<MyMpi> {
private:
	typedef std::complex<double> _Cmplx;
	MPI_Comm comm;
	int nprc;
	int myid;
	int mstr;
private:
	// Bcase
	int Bcast(char *buffer, int count) const {
		return MPI_Bcast(buffer, count, MPI_CHAR, mstr, comm);
	}
	int Bcast(int *buffer, int count) const {
		return MPI_Bcast(buffer, count, MPI_INT, mstr, comm);
	}
	int Bcast(double *buffer, int count) const {
		return MPI_Bcast(buffer, count, MPI_DOUBLE, mstr, comm);
	}
	int Bcast(_Cmplx *buffer, int count) const {
		return MPI_Bcast(buffer, count, MPI_DOUBLE_COMPLEX, mstr, comm);
	}
	// Reduce
	int Reduce(const char *sendbuf, char *recvbuf, int sendcount) const {
		return MPI_Reduce(sendbuf, recvbuf, sendcount, MPI_CHAR, MPI_SUM, mstr, comm);
	}
	int Reduce(const int *sendbuf, int *recvbuf, int sendcount) const {
		return MPI_Reduce(sendbuf, recvbuf, sendcount, MPI_INT, MPI_SUM, mstr, comm);
	}
	int Reduce(const double *sendbuf, double *recvbuf, int sendcount) const {
		return MPI_Reduce(sendbuf, recvbuf, sendcount, MPI_DOUBLE, MPI_SUM, mstr, comm);
	}
	int Reduce(const _Cmplx *sendbuf, _Cmplx *recvbuf, int sendcount) const {
		return MPI_Reduce(sendbuf, recvbuf, sendcount, MPI_DOUBLE_COMPLEX, MPI_SUM, mstr, comm);
	}
	// Allreduce
	int Allreduce(const char *sendbuf, char *recvbuf, int sendcount) const {
		return MPI_Allreduce(sendbuf, recvbuf, sendcount, MPI_CHAR, MPI_SUM, comm);
	}
	int Allreduce(const int *sendbuf, int *recvbuf, int sendcount) const {
		return MPI_Allreduce(sendbuf, recvbuf, sendcount, MPI_INT, MPI_SUM, comm);
	}
	int Allreduce(const double *sendbuf, double *recvbuf, int sendcount) const {
		return MPI_Allreduce(sendbuf, recvbuf, sendcount, MPI_DOUBLE, MPI_SUM, comm);
	}
	int Allreduce(const _Cmplx *sendbuf, _Cmplx *recvbuf, int sendcount) const {
		return MPI_Allreduce(sendbuf, recvbuf, sendcount, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
	}
	// Gatherv
	int Gatherv(const double *sendbuf, int sendcount, double *recvbuf, const int *recvcount, const int *displs) const {
		return MPI_Gatherv(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, displs, MPI_DOUBLE, mstr, comm);
	}
	int Gatherv(const _Cmplx *sendbuf, int sendcount, _Cmplx *recvbuf, const int *recvcount, const int *displs) const {
		return MPI_Gatherv(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, displs, MPI_DOUBLE_COMPLEX, mstr, comm);
	}

	// Gather 20191015 linheyu
	int Gather(const double *sendbuf, int sendcount, double *recvbuf, int recvcount) const {
		return MPI_Gather(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, mstr, comm);
	}
	int Gather(const _Cmplx *sendbuf, int sendcount, _Cmplx *recvbuf, int recvcount) const {
		return MPI_Gather(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, MPI_DOUBLE_COMPLEX, mstr, comm);
	}

	//ALL Gather 20190923 lhy
	int Allgather(const double *sendbuf, int sendcount, double *recvbuf, int recvcount) const {
		return MPI_Allgather(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, comm);
	}
	int Allgather(const char *sendbuf, int sendcount, char *recvbuf, int recvcount) const {
		return MPI_Allgather(sendbuf, sendcount, MPI_CHAR, recvbuf, recvcount, MPI_CHAR, comm);
	}
public:
	MyMpi(MPI_Comm comm_i = MPI_COMM_WORLD, int mstr_i = 0) : comm(comm_i), mstr(mstr_i) {
		MPI_Comm_size(comm, &nprc);
		MPI_Comm_rank(comm, &myid);
		if (myid < 0 || myid >= nprc || mstr < 0 || mstr >= nprc) ERR(NAV4(comm, nprc, myid, mstr));
		display();
	}
	MPI_Comm cm() const { return comm; }
	int np() const { return nprc; }
	int id() const { return myid; }
	int ms() const { return mstr; }
	bool boolean_test() const { return myid == mstr; }
	int barrier() const { return MPI_Barrier(comm); }
	void print_process_names(std::ostream &os = std::cout) const;
	void display() const {
		using namespace std;
		if (*this) cout << "MPI: number of processes = " << np() << endl;
		if (*this) cout << "MPI: master process id = " << ms() << endl;
		print_process_names();
	}
public:
	template <typename T>
	int Bcast(T &v) const {
		return Bcast(&v, 1);
	}
	template <typename T>
	int Bcast(Vec<T> &v) const {
		return Bcast(v.p(), (int)v.size());
	}
	template <typename T>
	T Reduce(const T &sendbuf) const {
		T recvbuf;
		Reduce(&sendbuf, &recvbuf, 1);
		return recvbuf;
	}
	template <typename T>
	Vec<T> Reduce(const Vec<T> &sendbuf) const {
		Vec<T> recvbuf(*this ? sendbuf.size() : 0);
		Reduce(sendbuf.p(), recvbuf.p(), (int)sendbuf.size());
		return *this ? recvbuf : sendbuf;
	}
	template <typename T>
	Mat<T> Reduce(const Mat<T> &sendbuf) const {
		Mat<T> recvbuf(*this ? sendbuf.nrows() : 0, *this ? sendbuf.ncols() : 0);
		Reduce(sendbuf.p(), recvbuf.p(), (int)sendbuf.size());
		return *this ? recvbuf : sendbuf;
	}
	template <typename T>
	T Allreduce(const T &sendbuf) const {
		T recvbuf;
		Allreduce(&sendbuf, &recvbuf, 1);
		return recvbuf;
	}
	template <typename T>
	Vec<T> Allreduce(const Vec<T> &sendbuf) const {
		Vec<T> recvbuf(sendbuf.size());
		Allreduce(sendbuf.p(), recvbuf.p(), (int)sendbuf.size());
		return recvbuf;
	}
	template <typename T>
	Vec<T> Gatherv(const T &sendbuf) const {
		Vec<T> recvbuf(*this ? np() : 0);
		VecPartition vp(np());
		Gatherv(&sendbuf, 1, recvbuf.p(), (const int *)vp.len(), (const int *)vp.dsp());
		return recvbuf;
	}
	template <typename T>
	Vec<T> Gatherv(const Vec<T> &sendbuf, const VecPartition &vp) const {
		Vec<T> recvbuf(*this ? vp.vec_size() : 0);
		Gatherv(sendbuf.p(), (int)sendbuf.size(), recvbuf.p(), (const int *)vp.len(), (const int *)vp.dsp());
		return recvbuf;
	}
	
	//Gather 20191015 linheyu
	template <typename T>
	Vec<T> Gather(const T &sendbuf) const {
		Vec<T> recvbuf(np());
		Gather(&sendbuf, 1, recvbuf.p(), 1);
		return recvbuf;
	}
	template <typename T>
	Vec<T> Gather(const Vec<T> &sendbuf) const {
		Vec<T> recvbuf(np()*sendbuf.size());
		Gather(sendbuf.p(), (int)sendbuf.size(), recvbuf.p(), (int)sendbuf.size());
		return recvbuf;
	}

	//ALL Gather 20190923 lhy
	template <typename T>
	Vec<T> Allgather(const T &sendbuf) const {
		Vec<T> recvbuf(np());
		Allgather(&sendbuf, 1, recvbuf.p(), 1);
		return recvbuf;
	}

	template <typename T>
	Vec<T> Allgather(const Vec<T> &sendbuf) const {
		Vec<T> recvbuf(np()*sendbuf.size());
		Allgather(sendbuf.p(), (int)sendbuf.size(), recvbuf.p(), (int)sendbuf.size());
		return recvbuf;
	}
};



inline void use_omp(const MyMpi &mm, Int nthreads = 1)
{
	using namespace std;
#ifdef _OPENMP
	omp_set_num_threads(nthreads);
	int omp_num_threads;
#pragma omp parallel
	omp_num_threads = omp_get_num_threads();
	if (mm) cout << "OpenMP: number of procs = " << omp_get_num_procs() << endl;
	if (mm) cout << "OpenMP: omp_num_threads = " << omp_num_threads << endl;
#endif
}

inline void use_mkl(const MyMpi &mm, Int nthreads = 1)
{
	using namespace std;
	//mkl_set_num_threads(nthreads);
	//mkl_set_dynamic;
	if (mm) cout << "MKL: mkl_num_threads = " << mkl_get_max_threads() << endl;
}

#endif /* _MYMPI_H_ */
