/*
 * Solver_matrix.h
 *
 *  Created on: May 3, 2015
 *      Author: zhou
 */

#ifndef _SOLVER_MATRIX_H_
#define _SOLVER_MATRIX_H_

#include <iostream>
#include <list>
#include "algebra_define.hpp"
#include "array_list.hpp"
#include "matrix_SparCompRow.hpp"
#include "BLAS_Level1.hpp"
#include "preconditioner.hpp"

#ifdef VIENNACL_WITH_OPENCL
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/tools/timer.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
#endif

namespace carpio {
#ifdef VIENNACL_WITH_OPENCL
template<class VALUE>
class Solver_ {
public:
	typedef VALUE vt;
	typedef Solver_<vt> Self;
	typedef Solver_<vt>* pSelf;
	typedef const Solver_<vt>* const_pSelf;

	typedef MatrixSCR_<vt> MatSCR;
	typedef ArrayListV<vt> Arr;
	// type define for viennaCL
	/*
	 * compressed_matrix<T, A> represents a sparse matrix
	 * using a compressed sparse row scheme.
	 * T is the floating point type.
	 * A is the alignment and defaults to 1 at present.
	 */
	typedef viennacl::compressed_matrix<vt> MatSCR_vcl;
	typedef viennacl::vector<vt> Arr_vcl;
	typedef viennacl::vcl_size_t st_vcl;

	static void Copy(MatSCR & cpu_matrix,
			MatSCR_vcl & gpu_matrix,
			st_vcl nonzeros) {
		ASSERT(
				(gpu_matrix.size1() == 0 || cpu_matrix.size1() == gpu_matrix.size1())
				&& bool("Size mismatch"));
		ASSERT(
				(gpu_matrix.size2() == 0 || cpu_matrix.size2() == gpu_matrix.size2())
				&& bool("Size mismatch"));

		viennacl::backend::typesafe_host_array<unsigned int> row_buffer(
				gpu_matrix.handle1(), cpu_matrix.size1() + 1);
		viennacl::backend::typesafe_host_array<unsigned int> col_buffer(
				gpu_matrix.handle2(), nonzeros);
		std::vector<vt> elements(nonzeros);

		viennacl::vcl_size_t row_index = 0;
		viennacl::vcl_size_t data_index = 0;

		//	for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1();
		//		row_it != cpu_matrix.end1(); ++row_it) {
		for (; row_index < cpu_matrix.size1();++row_index) {
			row_buffer.set(row_index, data_index);

			for (viennacl::vcl_size_t ir = cpu_matrix.row_ptr(row_index);
					ir < cpu_matrix.row_ptr(row_index + 1); ++ir) {
				col_buffer.set(data_index, cpu_matrix.col_ind(data_index));
				elements[data_index] = cpu_matrix.val(data_index);
				data_index++;
			}
			data_index = viennacl::tools::align_to_multiple<viennacl::vcl_size_t>(
					data_index, 1); //take care of alignment
		}
		row_buffer.set(row_index, data_index);

		gpu_matrix.set(row_buffer.get(), col_buffer.get(), &elements[0],
				cpu_matrix.size1(), cpu_matrix.size2(), nonzeros);
	}

	static void Solve(MatSCR const& A, Arr& x, Arr const& b, vt& tol, int& max_iter) {
		Arr_vcl bv(A.iLen());
		MatSCR_vcl Av(A.size1(), A.size2(),
				A.NumNonzeros());
		copy(b.begin(), b.end(), bv.begin());
		//std::cout<<bv;
		viennacl::linalg::bicgstab_tag tag(tol, max_iter);
		//compute Jacobi preconditioner:
		viennacl::linalg::jacobi_precond< MatSCR_vcl > vcl_jacobi(Av, viennacl::linalg::jacobi_tag());
		//solve (e.g. using conjugate gradient solver)
		Arr_vcl xv = viennacl::linalg::solve(Av, bv, tag,vcl_jacobi);
		//Arr_vcl xv;
		std::cout<<" av   "<< Av;
		//std::cout<< xv;
		std::vector<vt> stdx(b.size());
		viennacl::copy(xv.begin(), xv.end(), x.begin());
	}

};
#endif

/**
 * solve
 * a1 x + b1 y = c1
 * a2 x + b2 y = c2
 *
 */
template<class VALUE>
int solve(VALUE a1, VALUE b1, VALUE c1, VALUE a2, VALUE b2, VALUE c2, VALUE& x,
		VALUE& y) {
	ASSERT(!(a2 == 0 && b2 == 0));
	ASSERT(!(a1 == 0 && b1 == 0));
	if (a2 == 0) {
		y = double(c2) / double(b2);
	} else {
		y = double(c1 - a1 * c2 / a2) / double(b1 - a1 * b2 / a2 + SMALL);
	}
	x = double(c1 - b1 * y) / double(a1 + SMALL);
}

template<class VALUE>
int solver_gaussian_elimination(  //solver
		MatrixV<VALUE>& A,      //the Matrix     [in]  solver will change matrix
		ArrayListV<VALUE>& b    //the b          [in][out]  solver will change b
		) {
	//Assert
	ASSERT(A.size_i() == A.size_j());   //the matrix must be n x n;
	ASSERT(b.size() == A.size_i());     //the size b=n
	int n = b.size();
	int i, j, max_i;
	VALUE max, dum;
	// for each variable find pivot row and perform forward substitution
	for (i = 0; i < (n - 1); ++i) {
		//  find the pivot row
		max_i = i;
		max = Abs(A[i][i]);
		for (int ii = i + 1; ii < n; ++ii)
			if ((dum = Abs(A[ii][i])) > max) {
				max = dum;
				max_i = ii;
			}
		if (max == 0.0) {
			return -1;                // the matrix A is singular
		}
		// and if it differs from the current row, interchange the two rows.
		if (max_i != i) {
			for (j = 0; j < n; j++) {
				dum = A[i][j];
				A[i][j] = A[max_i][j];
				A[max_i][j] = dum;
			}
			dum = b[i];
			b[i] = b[max_i];
			b[max_i] = dum;
		}
		// Perform forward substitution
		for (int ii = i + 1; ii < n; ++ii) {
			dum = -A[ii][i] / A[i][i];
			for (j = i + 1; j < n; ++j) {
				A[ii][j] += dum * A[i][j];
			}
			b[ii] += dum * b[i];
		}
	}
	// Perform backward substitution
	for (i = n - 1; i >= 0; i--) {
		if (A[i][i] == 0.0)
			return -1;           // matrix is singular
		dum = b[i];
		for (j = i + 1; j < n; j++) {
			dum -= A[i][j] * b[j];
		}
		b[i] = dum / A[i][i];
	}
	return 0;
}

//spares ==================================================
//=========================================================

// Jacobi solver
template<class VALUE>
int Jacobi(const MatrixSCR_<VALUE> &A,    // A  The matrix
		ArrayListV<VALUE> &x,            // x
		const ArrayListV<VALUE>& b,      // b
		int &max_iter,                   // max iter
		Float &tol,                      // Tolerance
		std::list<Float>& lresid)        // list residual
		{
	Float resid;
	typename ArrayListV<VALUE>::size_type n = b.size();
	MatrixSCR_<VALUE> T(A);
	ArrayListV<VALUE> C(n, 1.0 / SMALL);
	ArrayListV<VALUE> newx(n);

	Float normb = nrm2(b);
	ArrayListV<VALUE> r = b - A * x;
	//
	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	// construct T
	St M = T.iLen();
	//st N = T.jLen();

	for (St i = 0; i < M; ++i) {
		// find D
		VALUE dv;
		int flag = -1;
		for (St j = T.row_ptr(i); j < T.row_ptr(i + 1); ++j) {
			if (T.col_ind(j) == i) {
				flag = j;
				dv = T.val(j);
				T.val(j) = 0;
				break;
			}
		}
		for (St j = T.row_ptr(i); j < T.row_ptr(i + 1); ++j) {
			if (j == St(flag)) {
				C(i) = b(T.col_ind(j)) / dv;
			} else {
				T.val(j) = -T.val(j) / dv;
			}
		}
	}

	for (int i = 1; i <= max_iter; ++i) {
		//----
		newx = T * x + C;
		r = b - A * newx;
		//----
		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid <= tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
		x = newx;
	}

	tol = resid;
	return 1;
}

// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
template<class VALUE>
int CG(const MatrixSCR_<VALUE> &A,    // A  The matrix
		ArrayListV<VALUE> &x,        // x
		const ArrayListV<VALUE>& b,  // b
		int &max_iter, Float &tol,   // Tolerance
		std::list<Float>& lresid)        // list residual
		{
	Float resid;
	ArrayListV<VALUE> p(b.size()), z(b.size()), q(b.size());
	//Array alpha(1), beta(1), rho(1), rho_1(1);
	VALUE alpha, beta, rho, rho_1;

	Float normb = nrm2(b);
	ArrayListV<VALUE> r = b - A * x;
	//
	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; ++i) {
		z = r;           //z = M.solve(r);   ===!!!!!
		rho = dot(r, z); //rho(0) = dot(r, z);

		if (i == 1)
			p = z;
		else {
			beta = rho / rho_1; //beta(0) = rho(0) / rho_1(0);
			p = z + beta * p;   //p = z + beta(0) * p;
		}

		q = A * p;
		alpha = rho / dot(p, q);  //alpha(0) = rho(0) / dot(p, q);
		x += alpha * p;           //x += alpha(0) * p;
		r -= alpha * q;           //r -= alpha(0) * q;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid <= tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
		rho_1 = rho;            //rho_1(0) = rho(0);
	}

	tol = resid;
	return 1;
}

template<class VALUE>
int IC_CG(const MatrixSCR_<VALUE> &A,  //
		ArrayListV<VALUE> &x,         //
		const ArrayListV<VALUE>& b,   //
		int &max_iter, Float &tol,    //
		std::list<Float>& lresid)         //
		{
	Float resid;
	ArrayListV<VALUE> p, z, q;
	Float alpha, beta, rho, rho_1;
	ICPre<VALUE> preA(A);
	Float normb = nrm2(b);
	ArrayListV<VALUE> r = b - A * x;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		//if(i<=max_iter){
		z = preA.solve(r);
		//}else{
		//	z = r;                    //z = M.solve(r);   ===!!!!!
		//}
		rho = dot(r, z);             //rho(0) = dot(r, z);

		if (i == 1)
			p = z;
		else {
			beta = rho / rho_1; //beta(0) = rho(0) / rho_1(0);
			p = z + beta * p;   //p = z + beta(0) * p;
		}

		q = A * p;
		alpha = rho / dot(p, q);  //alpha(0) = rho(0) / dot(p, q);
		x += alpha * p;           //x += alpha(0) * p;
		r -= alpha * q;           //r -= alpha(0) * q;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid <= tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}

		rho_1 = rho;            //rho_1(0) = rho(0);
	}

	tol = resid;
	return 1;
}

// CGS solves the unsymmetric linear system Ax = b
// using the Conjugate Gradient Squared method
template<class VALUE>
int IC_CGS(const MatrixSCR_<VALUE> &A,    // A  The matrix
		ArrayListV<VALUE> &x,        // x
		const ArrayListV<VALUE>& b,  // b
		int &max_iter,   //max_iter
		Float &tol,   // Tolerance
		std::list<Float>& lresid) {
	ICPre<VALUE> preA(A);
	Float resid;
	VALUE rho_1, rho_2, alpha, beta;
	ArrayListV<VALUE> p, phat, q, qhat, vhat, u, uhat;

	Float normb = nrm2(b);
	ArrayListV<VALUE> r = b - A * x;
	ArrayListV<VALUE> rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		rho_1 = dot(rtilde, r);
		if (rho_1 == 0) {
			tol = nrm2(r) / normb;
			return 2;
		}
		if (i == 1) {
			u = r;
			p = u;
		} else {
			beta = rho_1 / rho_2;
			u = r + beta * q;
			p = u + beta * (q + beta * p);
		}
		phat = preA.solve(p);
		vhat = A * phat;
		alpha = rho_1 / dot(rtilde, vhat);
		q = u - alpha * vhat;
		uhat = preA.solve(u + q);
		x += alpha * uhat;
		qhat = A * uhat;
		r -= alpha * qhat;
		rho_2 = rho_1;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid < tol) {
			tol = resid;
			max_iter = i;
			return 0; //converged to the desired tolerance tol within maxit iterations.
		}
	}

	tol = resid;
	return 1;         // iterated maxit times but did not converge.
}

//============================================================
// BiCG solves the unsymmetric linear system Ax = b
// using the Preconditioned BiConjugate Gradient method
template<class VALUE>
int IC_BiCG(const MatrixSCR_<VALUE> &A,    // A  The matrix
		ArrayListV<VALUE> &x,        // x
		const ArrayListV<VALUE>& b,  // b
		int &max_iter,   //max_iter
		Float &tol,   // Tolerance
		std::list<Float>& lresid) {
	ICPre<VALUE> preA(A);
	Float resid;
	VALUE rho_1, rho_2, alpha, beta;
	ArrayListV<VALUE> z, ztilde, p, ptilde, q, qtilde;

	Float normb = nrm2(b);
	ArrayListV<VALUE> r = b - A * x;
	ArrayListV<VALUE> rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		z = preA.solve(r);
		ztilde = preA.trans_solve(rtilde);
		rho_1 = dot(z, rtilde);
		if (rho_1 == 0) {
			tol = nrm2(r) / normb;
			max_iter = i;
			return 2;
		}
		if (i == 1) {
			p = z;
			ptilde = ztilde;
		} else {
			beta = rho_1 / rho_2;
			p = z + beta * p;
			ptilde = ztilde + beta * ptilde;
		}
		q = A * p;
		qtilde = A.transMult(ptilde);
		alpha = rho_1 / dot(ptilde, q);
		x += alpha * p;
		r -= alpha * q;
		rtilde -= alpha * qtilde;

		rho_2 = rho_1;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
	}

	tol = resid;
	return 1;
}

//===============================================================
// BiCGSTAB solves the unsymmetric linear system Ax = b
// using the Preconditioned BiConjugate Gradient Stabilized method
template<class VALUE>
int IC_BiCGSTAB( //
		const MatrixSCR_<VALUE> &A,    // A  The matrix
		ArrayListV<VALUE> &x,        // x
		const ArrayListV<VALUE>& b,  // b
		int &max_iter,   //max_iter
		Float &tol,      // Tolerance
		std::list<Float>& lresid) {
	ICPre<VALUE> preA(A);
	Float resid;
	VALUE rho_1, rho_2, alpha, beta, omega;
	ArrayListV<VALUE> p, phat, s, shat, t, v;

	Float normb = nrm2(b);
	ArrayListV<VALUE> r = b - A * x;
	ArrayListV<VALUE> rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		rho_1 = dot(rtilde, r);
		if (rho_1 == 0) {
			tol = nrm2(r) / normb;
			return 2;
		}
		if (i == 1)
			p = r;
		else {
			beta = (rho_1 / rho_2) * (alpha / omega);
			p = r + beta * (p - omega * v);
		}
		phat = preA.solve(p);
		v = A * phat;
		alpha = rho_1 / dot(rtilde, v);
		s = r - alpha * v;
		if ((resid = nrm2(s) / normb) < tol) {
			x += alpha * phat;
			tol = resid;
			return 0;
		}
		shat = preA.solve(s);
		t = A * shat;
		omega = dot(t, s) / dot(t, t);
		x += alpha * phat + omega * shat;
		r = s - omega * t;

		rho_2 = rho_1;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}

		if (omega == 0) {
			tol = nrm2(r) / normb;
			return 3;
		}
	}

	tol = resid;
	return 1;
}

template<class VALUE>
double Residual( //
		const MatrixSCR_<VALUE> &A,    // A  The matrix
		const ArrayListV<VALUE> &x,    // x
		const ArrayListV<VALUE>& b)    // b
		{
	double resid = 0.0;
	double normb = nrm2(b);
	ArrayListV<VALUE> r = b - A * x;
	if (normb == 0.0) {
		normb = 1;
	}
	resid = nrm2(r) / normb;
	return resid;
}

template<class VALUE>
int Dia_BiCGSTAB( //
		const MatrixSCR_<VALUE> &A,  // A  The matrix [in]
		ArrayListV<VALUE> &x,        // x             [in, out]
		const ArrayListV<VALUE>& b,  // b             [in]
		int &max_iter,               //max_iter       [in, out]
		Float &tol,                  // Tolerance     [in, out]
		std::list<Float>& lresid) {
	DiaPre<VALUE> preA(A, 1);
	Float resid;
	VALUE rho_1, rho_2, alpha, beta, omega;
	ArrayListV<VALUE> p, phat, s, shat, t, v;

	Float normb = nrm2(b);
	ArrayListV<VALUE> r = b - A * x;
	ArrayListV<VALUE> rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		lresid.push_back(resid);
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		rho_1 = dot(rtilde, r);                         //dot(v,v)
		if (rho_1 == 0) {
			tol = nrm2(r) / normb;                      //norm(v) / s
			return 2;
		}
		if (i == 1)
			p = r;                                      // v=v
		else {
			beta = (rho_1 / rho_2) * (alpha / omega);   // s
			p = r + beta * (p - omega * v);             // v + s*(v-s*v)
		}
		phat = preA.solve(p);
		v = A * phat;                                   // v = A * phat
		alpha = rho_1 / dot(rtilde, v);
		s = r - alpha * v;
		if ((resid = nrm2(s) / normb) < tol) {
			x += alpha * phat;
			tol = resid;
			lresid.push_back(resid);
			return 0;
		}
		shat = preA.solve(s);
		t = A * shat;
		omega = dot(t, s) / dot(t, t);
		x += alpha * phat + omega * shat;
		r = s - omega * t;

		rho_2 = rho_1;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}

		if (omega == 0) {
			tol = nrm2(r) / normb;
			return 3;
		}
	}

	tol = resid;
	return 1;
}

template<class VALUE>
int BiCGSTAB( //
		const MatrixSCR_<VALUE> &A,    // A  The matrix
		ArrayListV<VALUE> &x,        // x
		const ArrayListV<VALUE>& b,  // b
		int &max_iter,   //max_iter
		Float &tol,      // Tolerance
		std::list<Float>& lresid) {
	DiaPre<VALUE> preA(A, 2);
	Float resid;
	VALUE rho_1, rho_2, alpha, beta, omega;
	ArrayListV<VALUE> p, phat, s, shat, t, v;

	Float normb = nrm2(b);
	ArrayListV<VALUE> r = b - A * x;
	ArrayListV<VALUE> rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		rho_1 = dot(rtilde, r);          //dot(v,v)
		if (rho_1 == 0) {
			tol = nrm2(r) / normb;       //norm(v) / s
			return 2;
		}
		if (i == 1)
			p = r;                       // v=v
		else {
			beta = (rho_1 / rho_2) * (alpha / omega);   // s
			p = r + beta * (p - omega * v);             // v + s*(v-s*v)
		}
		phat = preA.solve(p);
		v = A * phat;
		alpha = rho_1 / dot(rtilde, v);
		s = r - alpha * v;
		if ((resid = nrm2(s) / normb) < tol) {
			x += alpha * phat;
			tol = resid;
			return 0;
		}
		shat = preA.solve(s);
		t = A * shat;
		omega = dot(t, s) / dot(t, t);
		x += alpha * phat + omega * shat;
		r = s - omega * t;

		rho_2 = rho_1;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}

		if (omega == 0) {
			tol = nrm2(r) / normb;
			return 3;
		}
	}

	tol = resid;
	return 1;
}

// abstract class solver

template<class VALUE>
class Solver_ {
public:
	typedef VALUE Vt;
	typedef typename ArrayListV<Vt>::size_type St;
	typedef MatrixSCR_<Vt> MatSCR;
	typedef ArrayListV<Vt> Arr;
	typedef std::list<double> Listr;
protected:
	int _max_iter;   //max_iter
	double _tol;     // Tolerance

	int _num_iter;
	double _resid;
	Listr _lresid; //
public:

	Solver_(int max_iter = 100, Vt tol = 1e-3) {
		_max_iter = max_iter;
		_tol = tol;

		_num_iter = 0;
		_resid = 1e10;
	}

	virtual int solve(const MatSCR &A, // A  The matrix
			Arr &x,          // x
			const Arr& b    // b
			) {
		SHOULD_NOT_REACH;
		return 0;
	}
	;

	int max_iter() const {
		return this->_max_iter;
	}

	int tolerance() const {
		return this->_tol;
	}

	int num_iter() const {
		return this->_num_iter;
	}

	double residual() const {
		return this->_resid;
	}

	Arr get_residual_array() const{
		Arr res(this->_lresid.size());
		St i = 0;
		for(auto val : this->_lresid)
		{
			res[i] = val;
			i++;
		}
		return res;
	}

	virtual ~Solver_() {
	}
	;
protected:
	void _init() {
		this->_lresid;
	}

};

template<class VALUE>
class Solver_Jacobi_: public Solver_<VALUE> {
public:
	typedef VALUE Vt;
	typedef typename ArrayListV<Vt>::size_type St;
	typedef MatrixSCR_<VALUE> MatSCR;
	typedef ArrayListV<VALUE> Arr;
	typedef std::list<double> Listr;

	Solver_Jacobi_(int max_iter = 100, Vt tol = 1e-3) :
			Solver_<VALUE>(max_iter, tol) {
	}

	int solve(const MatSCR &A, // A  The matrix
			Arr &x,          // x
			const Arr& b    // b
			) {
		this->_init();

		return this->_solve(A, x, b);
	}

	int _solve( //
			const MatSCR &A, // A  The matrix
			Arr &x,          // x
			const Arr& b    // b
			) {
		Vt resid;

		St n = b.size();
		MatSCR T(A);
		Arr C(n, 1.0 / SMALL);
		Arr newx(n);

		Vt normb = nrm2(b);
		Arr r = b - A * x;
		//
		if (normb == 0.0)
			normb = 1;

		if ((resid = nrm2(r) / normb) <= this->_tol) {
			this->_resid = resid;
			this->_num_iter = 0;
			return 0;
		}

		// construct T
		St M = T.iLen();
		//st N = T.jLen();

		for (St i = 0; i < M; ++i) {
			// find D
			Vt dv;
			int flag = -1;
			for (St j = T.row_ptr(i); j < T.row_ptr(i + 1); ++j) {
				if (T.col_ind(j) == i) {
					flag = j;
					dv = T.val(j);
					T.val(j) = 0;
					break;
				}
			}
			for (St j = T.row_ptr(i); j < T.row_ptr(i + 1); ++j) {
				if (j == St(flag)) {
					C(i) = b(T.col_ind(j)) / dv;
				} else {
					T.val(j) = -T.val(j) / dv;
				}
			}
		}

		for (int i = 1; i <= this->_max_iter; ++i) {
			//----
			newx = T * x + C;
			r = b - A * newx;
			//----
			resid = nrm2(r) / normb;
			this->_lresid.push_back(resid);
			if (resid <= this->_tol) {
				this->_resid = resid;
				this->_num_iter = i;
				return 0;
			}
			x = newx;
		}

		this->_resid = resid;
		this->_num_iter = this->_max_iter;
		return 1;
	}
};

template<class VALUE>
class Solver_SOR_: public Solver_<VALUE> {
public:
	typedef VALUE Vt;
	typedef typename ArrayListV<Vt>::size_type St;
	typedef MatrixSCR_<VALUE> MatSCR;
	typedef ArrayListV<VALUE> Arr;
	typedef std::list<double> Listr;

protected:
	Vt _omega;
public:
	Solver_SOR_(int max_iter = 100, Vt tol = 1e-3, Vt omega = 1.0) :
			Solver_<VALUE>(max_iter, tol) {
		_omega = omega;
	}

	int solve(const MatSCR &A, // A  The matrix
			Arr &x,          // x
			const Arr& b    // b
			) {
		this->_init();

		return this->_solve(A, x, b);
	}

	int _solve( //
			const MatSCR &A, // A  The matrix
			Arr &x,          // x
			const Arr& b    // b
			) {
		Vt resid;

		St n = b.size();
		MatSCR T(A);
		Arr C(n, 1.0 / SMALL);
		//Arr newx(n);

		Vt normb = nrm2(b);
		Arr r = b - A * x;
		//
		if (normb == 0.0)
			normb = 1;

		if ((resid = nrm2(r) / normb) <= this->_tol) {
			this->_resid = resid;
			this->_num_iter = 0;
			return 0;
		}

		// construct T
		St M = T.iLen();
		//st N = T.jLen();

		for (int i = 1; i <= this->_max_iter; ++i) {
			//----
			for (St i = 0; i < M; ++i) {  //row
				Vt sumkp = 0;
				Vt sumk = 0;
				Vt aa = 1.0;
				for (St idx = T.row_ptr(i); idx < T.row_ptr(i + 1); ++idx) {
					St j = T.col_ind(idx);
					if (j < i) {
						sumkp += T.val(idx) * x[j];
					} else if (j > i) {
						sumk += T.val(idx) * x[j];
					} else {
						aa = T.val(idx);
					}
				}
				x[i] = (1 - this->_omega) * x[i]
						+ (this->_omega / aa) * (b[i] - sumkp - sumk);
			}
			r = b - A * x;
			//----
			resid = nrm2(r) / normb;
			this->_lresid.push_back(resid);
			if (resid <= this->_tol) {
				this->_resid = resid;
				this->_num_iter = i;
				return 0;
			}
		}

		this->_resid = resid;
		this->_num_iter = this->_max_iter;
		return 1;
	}
};

} //end namespace

#endif /* ALGEBRA_SOLVER_MATRIX_H_ */
