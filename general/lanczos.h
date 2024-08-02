/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2012-12-17
*/

#ifndef _LANCZOS_H_
#define _LANCZOS_H_

#include "toolbox.h"
#include "eigsym.h"

template<typename T, typename F>
void lanczos_one_step(const Int &n, Vec<T> &hp, Vec<T> &v0, Vec<T> &v1, 
					  const Int &s, VecReal &ta, VecReal &tb, Real &beta, F &ho);
template<typename T>
void orthonormalize_one_to_multiple(Vec<T> &v, const Mat<T> &m, Int e);
inline Int line_identify(const VecReal &gse, const Int &kssr);
inline Int line_best(const VecReal &gsr, const Int &left, const Int &kssr);

template<typename T, typename F>
Int lanczos(Idx n, Int nev, VecReal &eval, Mat<T> &evec, Int &nev_eff, F &ho, 
			Int giv = 0, Int nlvl = 0, Real tol = 0., Int dsp = 0, Int ncv = 999)
//	Lanczos algorithm to find several lowest eigenpairs of a Hermitian matrix
//	only matrix-vector multiplication (functor ho) is required, as F ho(hp, p) implied
//	F ho is a callable object with call signature void ho(Vec<T> &hp, const Vec<T> &p), which returns H * |p> in hp
//	typename T can only be Real or Cmplx
//	n      ( input) is the dimension of the matrix
//	nev    ( input) is the expected maximum eigenstates to be found
//	ncv    ( input) is the maximum size of the Krylov space
//	eval   (output) is the eigenvalues
//	evec   (output) is the eigenvectors, one row contains one vector
//	nlvl   ( input) is the number of different eigen levels
//	giv    ( input) how many guessed initial lanczos vectors contained in evec
//	dsp    ( input) is if detailed information will be displayed
//	actual number of eigenpairs is returned
//	nev_eff is the number of eigenvectors corresponding to complete degenerate levels reached
{
	using namespace std;

	Int ncnv;		// number of eigenpairs converged
	Int nlvl_f;		// nlvl found
	Int val_cnv;	// if e-th eval converged
	Int vec_cnv;	// if e-th evec converged
	VecReal evr(ncv);	// eigenvector residuals

	// for tridiagonal matrix
	VecReal ta(ncv);	//     diagonal elements / eigenvalues
	VecReal tb(ncv);	// sub-diagonal elements
	VecReal tr(ncv);	// residual

	// record the history of expanding of Krylov space
	VecReal gse(ncv);	// ground state energies
	VecReal gsr(ncv);	// ground state residuals
	Int left;		// left endpoint of the history
	Int kssm;		// Krylov space size, maximum
	Int kssr;		// Krylov space size, reached
	Int kssu;		// Krylov space size, used
	Real gsel;		// gse obtained from diagonalization in Krylov space of size kssr

	Real havg;		// <psi| H |psi>
	Real gsre;		// gsr exact, = |H * |psi> - havg * |psi>|
	Real beta;		// one minor diagonal value of the tridiagonal matrix
	Int ierr;		// error code for diagonalization subroutines
	T ttmp;

	Vec<T> hp(n);					// H * |psi>
	Vec<T> v0(n);
	Vec<T> v1(n);

	nev_eff = 0;
	nlvl_f = 0;
	if (nev > Int(n)) nev = Int(n);
	// determining each eigenvector on by on, later one is force to be orthogonal to all previous
	for (Int e = 0; e != nev; ++e) {

		// prepare an initial Lanczos vector
		Vec <T> ev(evec[e]);
		if (giv <= e) for_Idx (i, 0, n) ev[i] = rand_unit(ev[i]) - half_unit(ev[i]);
		orthonormalize_one_to_multiple(ev, evec, e);

//		for_Idx (i, 0, n) ev[i] = mc[0][i];

		// run Lanczos steps for the first time
		kssm = ncv;
		if (kssm > Int(n - e)) kssm = Int(n - e);
		v0 = ev;
		v1 = 0.;
		beta = 0.;
		val_cnv = 0;
		vec_cnv = 0;
		for (Int s = 0; s != kssm; ++s) {
			lanczos_one_step(n, hp, v0, v1, s, ta, tb, beta, ho);
			if (beta > 1.E+10 * sqrt(n * min_Real)) {
				orthonormalize_one_to_multiple(v1, evec, e);
			} else {
				val_cnv = 1;
				vec_cnv = 1;
			}

			// check if converged
			kssr = s + 1;
			for_Int (i, 0, s) tr[i] = 0.;
			tr[s] = 1.;
			VecReal va(kssr); va = VecReal(kssr, ta.p());
			VecReal vb(kssr); vb = VecReal(kssr, tb.p());
			ierr = tqlis(VecReal(kssr, tr.p()), va, vb);
			if (ierr >= 0) ERR("tqlis ierr >= 0, " + NAV4(ierr, s, e, n));
			for_Int (i, 0, s) tr[i] = ABS(tr[i] * beta);
			gsr[s] = tr[0];
			gse[s] = va[0];
			Real tol_live = MAX(tol, 1.E+2 * eps_Real * sqrt(1. * kssr));
			if (gsr[s] < tol_live) val_cnv = vec_cnv = 1;
			//if (kssr >= 5) {
			//	if (!val_cnv && gse[kssr - 5] <= MIN(VecReal(4, &gse[kssr - 4]))) val_cnv = 1;
			//	if (val_cnv && !vec_cnv && gsr[kssr - 5] <= tol_live && gsr[kssr - 5] <= MIN(VecReal(4, &gsr[kssr - 4]))) vec_cnv = 1;
			//}
			if (dsp == 1) {
				cout << "[" << e << "] " 
					<< setw(4) << s << " " 
					<< iofmt("sci_def") 
					<< setw(w_Real_def) << tr[0] << " " 
					<< setw(w_Real_def) << tr[s > 0 ? 1 : 0] << " " 
					<< iofmt("sci") 
					<< setw(w_Real) << va[0] << " " 
					<< setw(w_Real) << va[s > 0 ? 1 : 0] << endl 
					<< iofmt("def");
			}
			if (vec_cnv || kssr == n - e) break;

			// swap v0 and v1
			swapall(v0, v1);
		}

		left = line_identify(gse, kssr);
		kssu = line_best(gsr, left, kssr);
		MatReal tu(kssu, kssu, 0.);			// unitary transformation matrix, the cols constitue eigenvectors
		tu.assign_diag(1.);
		VecReal va(kssu); va = VecReal(kssu, ta.p());
		VecReal vb(kssu); vb = VecReal(kssu, tb.p());
		ierr = tqli(tu, va, vb);
		if (ierr >= 0) ERR("tqli ierr >= 0, " + NAV3(ierr, e, n));
		gsel = va[0];		// the lowest eigenvalue

		// run Lanczos steps again to obtain the eigenvector ev for the matrix to diagonalize
		v0 = ev;
		v1 = 0.;
		beta = 0.;
		ev = 0.;
		for (Int s = 0; s != kssu; ++s) {
			for_Idx (i, 0, n) ev[i] = ev[i] + tu[s][0] * v0[i];
			lanczos_one_step(n, hp, v0, v1, s, ta, tb, beta, ho);
			orthonormalize_one_to_multiple(v1, evec, e);
			swapall(v0, v1);
		}
		
		tu.reset();
		orthonormalize_one_to_multiple(ev, evec, e);

		// check if ev is correct and accurate
		ho(hp, ev);
		ttmp = DOT(ev, hp);
		for_Idx (i, 0, n) v0[i] = hp[i] - ev[i] * ttmp;
		ttmp = ttmp + DOT(ev, v0);
		if (ABS(imag(ttmp)) / (1. + ABS(ttmp)) > sqrt(n * 1.E+6) * eps_Real)
			WRN("eigenvalue is not real, matrix may not be Hermitian because " + NAV3(ttmp, e, n));
		havg = real(ttmp);
		for_Idx (i, 0, n) v0[i] = hp[i] - ev[i] * havg;
		gsre = sqrt(real(DOT(v0, v0)));
		cout << iofmt("sci") << "[" << e << "] " 
			<< setw(4) << kssm << " " 
			<< setw(4) << kssr << " " 
			<< setw(4) << kssu << " " 
			<< setw(4) << left << " " 
			<< setw(w_Real) << gsre << " " 
			<< setw(w_Real) << havg << endl;

		Real tol_live = MAX(tol, 1.E+3 * eps_Real * sqrt(1. * n * kssu));
		if (ABS((havg - gsel)) / (ABS(havg) + ABS(gsel) + 1.) > tol_live) {
			std::ostringstream oss;
			oss << "eigenvalues mismatch, havg = " << iofmt("sci") << setw(w_Real) << havg 
				<< ", gsel = " << setw(w_Real) << gsel 
				<< ", (left of >) = " << setw(w_Real) << ABS((havg - gsel)) / (ABS(havg) + ABS(gsel) + 1.) 
				<< ", (right of >) = " << setw(w_Real) << tol_live;
			WRN(oss.str());
		}

		ncnv = e + 1;
		eval[e] = havg;
		evr[e] = gsre;

		if (e && ABS(eval[e] - eval[e - 1]) / (ABS(eval[e]) + ABS(eval[e - 1]) + 1.) > 1.E+3 * eps_Real) {
			nlvl_f = nlvl_f + 1;
			nev_eff = e;
		}
		if (nlvl_f == nlvl) break;
	}

	if (dsp == 1) {
		if (nlvl_f >= nlvl || nlvl >= nev) std::cout << "summary of eigenvalues and residuals" << std::endl;
		if (nlvl_f <  nlvl && nlvl <  nev) std::cout << "summary of eigenvalues and residuals - WRN: nev is too small" << std::endl;
		cout << iofmt("sci");
		for_Int (e, 0, ncnv) cout << e << "\t" << setw(w_Real) << eval[e] << "\t" << setw(w_Real) << evr[e] << endl;
		cout << iofmt("def");
	}

	if (nlvl_f < nlvl && nlvl < nev) WRN("nev is too small, " + NAV3(nev, nlvl, nlvl_f));

	if (nlvl >= nev) nev_eff = ncnv;

	cout << iofmt("def");
	return ncnv;
}
template<typename T, typename F>
void lanczos_one_step(const Int &n, Vec<T> &hp, Vec<T> &v0, Vec<T> &v1, 
					  const Int &s, VecReal &ta, VecReal &tb, Real &beta, F &ho)
//	enlarge the Krylov space by one
//	now it is the s-th step, s goes from 0 to ncv - 1
//	v0 is the s-th Lanczos vector
//	v1 (on  input) is the (s - 1)-st Lanczos vector
//	v1 (on output) is the (s + 1)-st Lanczos vector
//	hp is a temporary vector
//	ta stores the       diagonal elements of the Lanczos tridiagonal matrix
//	tb stores the minor diagonal elements of the Lanczos tridiagonal matrix
//	ta[0, ncv - 1] and tb[0, ncv - 1]
//	tb[0] is irrelavent
//	BGN ---------------- ATTENTION ---------------- BGN
//	on output v1 is not normalized, but on input it is must be normalized
//	this work is left to orthonormalize_one_to_multiple due to efficency reason
//	END ---------------- ATTENTION ---------------- END
{
	using namespace std;
	ho(hp, v0);
	tb[s] = beta;
	if (s) for_Idx (i, 0, n) hp[i] = hp[i] - v1[i] * beta;
	T ttmp = DOT(v0, hp);
	if (ABS(imag(ttmp)) / (1. + ABS(ttmp)) > sqrt(n * 1.E+6) * eps_Real)
		WRN("matrix may not be Hermitian because " + NAV3(ttmp, s, n));
	ta[s] = real(ttmp);
	if (!s) {
		// this "if" find back some significant digits in ta[s] if v0 is close to an eigenvector of matrix
		// v1 is seen as a tmp vector here
		for_Idx (i, 0, n) v1[i] = hp[i] - v0[i] * ta[s];
		ttmp = DOT(v0, v1);
		ta[s] = ta[s] + real(ttmp);
	}

	for_Idx (i, 0, n) v1[i] = hp[i] - v0[i] * ta[s];
	beta = sqrt(real(DOT(v1, v1)));
	if (!s) {
		// this "if" orthogonalize v1 to v0, v1 may not be orthogonal to v0 if |v1| is small, 
		// which is possible if v0 is close to an eigenvector of matrix
		ttmp = DOT(v0, v1);
		for_Idx (i, 0, n) v1[i] = v1[i] - v0[i] * ttmp;
	}
}
template<typename T>
void orthonormalize_one_to_multiple(Vec<T> &v, const Mat<T> &m, Int e)
//	orthonormalize v to m[0] through m[e - 1], which are assumed to be orthonormal
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ1_PS(v.size(), m.ncols(), NAVC3(e, m.nrows(), m.ncols()));
#endif
	Idx n = v.size();
	Real mag0, mag1;
	mag1 = sqrt(real(DOT(v, v)));		// mag1 != 0 is supposed
	mag0 = mag1 * 1.E+16;
	while (mag1 / mag0 < 0.9) {
		for (Idx j = 0; j != e; ++j) {
			T mag2 = DOT(m[j], v);
			for (Idx i = 0; i != n; ++i) v[i] = v[i] - m[j][i] * mag2;
		}
		mag0 = mag1;
		mag1 = sqrt(real(DOT(v, v)));		// mag1 != 0 is supposed
	}
	mag1 = 1. / mag1;
	for (Idx i = 0; i != v.size(); ++i) v[i] = v[i] * mag1;
}
inline Int line_identify(const VecReal &gse, const Int &kssr)
// return the left most index satisfying gse[index] <= avg(gse[left : right_most_index]) + eps;
{
	Int left = kssr - 1;
	Real eps = ABS(gse[left]) * 4 * eps_Real;
	Real last_gse = gse[left];
	Real dlt, tmp = 0.;
	for (Int i = left; i >= 0; --i) {
		dlt = gse[i] - last_gse;
		tmp = tmp + dlt;
		if ((gse[i] - last_gse) * (kssr - i) <= tmp + eps * (kssr - i)) left = i;
	}
	return left;
}
inline Int line_best(const VecReal &gsr, const Int &left, const Int &kssr)
// return index + 1 that index in [left, kssr) that gsr[index] is minimal
{
	Int best = kssr - 1;
	Real tmp = gsr[best];
	for (Int i = best - 1; i >= left; --i) {
		if (tmp >= gsr[i]) {
			tmp = gsr[i];
			best = i;
		}
	}
	return best + 1;
}

#endif /* _LANCZOS_H_ */
