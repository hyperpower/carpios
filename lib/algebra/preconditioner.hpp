/************************
 //  \file   ICpre.h
 //  \brief
 // 
 //  \author zhou
 //  \date   12 mai 2014 
 ***********************/
#ifndef PRECONDITIONER_H_
#define PRECONDITIONER_H_

#include "matrix_SparCompRow.hpp"
#include "matrix.hpp"
#include "sort.hpp"
#include "array_list.hpp"
#include "algebra_define.hpp"

namespace carpio
{

template<class VALUE>
void ICFactor(const ArrayListV<St> &pntr,
		const ArrayListV<St> &indx, ArrayListV<VALUE> &val)
{

	int d, g, h, i, j, k, n = pntr.size() - 1;
	Float z;

	for (k = 0; k < n - 1; k++) {
		d = pntr[k];
		z = val[d] = sqrt(val[d]);

		for (i = d + 1; i < pntr[k + 1]; i++)
			val[i] /= z;

		for (i = d + 1; i < pntr[k + 1]; i++) {
			z = val[i];
			h = indx[i];
			g = i;

			for (j = pntr[h]; j < pntr[h + 1]; j++)
				for (; g < pntr[k + 1] && indx[g + 1] <= indx[j]; g++)
					if (indx[g] == indx[j])
						val[j] -= z * val[g];
		}
	}
	d = pntr[n - 1];
	val[d] = sqrt(val[d]);
}

template<class VALUE>
class ICPre
{
public:
	typedef VALUE vt;
private:
	ArrayListV<VALUE> val_;
	ArrayListV<St> pntr_;
	ArrayListV<St> indx_;
	int nz_;
	int dim_[2];

public:
	//ICPre(const CompCol_Mat_double &A);
	ICPre(const Matrix &A) :
			val_(0), pntr_(A.size_i() + 1), indx_(0), nz_(0)
	{
		dim_[0] = A.size_i();
		dim_[1] = A.size_j();
		//Count non-zero number in upper triangular part (including diagonal)
		for (int i = 0; i < dim_[0]; i++) {
			for (int j = i; j < dim_[1]; j++) {
				if (A[i][j] != 0.0)
					nz_++;
			}
		}
		val_.reconstruct(nz_);
		indx_.reconstruct(nz_);
		pntr_[0] = 0;
		for (int i = 0; i < dim_[0]; i++) {
			pntr_[i + 1] = pntr_[i];
			for (int j = i; j < dim_[1]; j++) {
				if (A[i][j] != 0.0) {
					int k = pntr_[i + 1]++;
					val_[k] = A[i][j];
					indx_[k] = j;
				}
			}
		}

		ICFactor(pntr_, indx_, val_);
	}

	ICPre(const MatrixSCR_<VALUE> &A) :
			val_(0), pntr_(A.iLen() + 1), indx_(0), nz_(0)
	{

		dim_[0] = A.iLen();
		dim_[1] = A.jLen();

		int i, j, k;
		//Count non-zero number in upper triangular part (including diagonal)
		for (k = 0; k < dim_[0]; k++) {
			for (j = A.row_ptr(k); j < A.row_ptr(k + 1); j++) {
				if (A.col_ind(j) >= k)
					nz_++;
			}
		}

		val_.reconstruct(nz_);
		indx_.reconstruct(nz_);

		// Copy just triangular part (including diagonal)
		pntr_[0] = 0;
		for (k = 0; k < dim_[0]; k++) {
			pntr_[k + 1] = pntr_[k];
			for (j = A.row_ptr(k); j < A.row_ptr(k + 1); j++) {
				if (A.col_ind(j) >= k) {
					i = pntr_[k + 1]++;
					val_[i] = A.val(j);
					indx_[i] = A.col_ind(j);
				}
			}
		}

		for (i = 0; i < dim_[0]; i++)
			QSort_ascend(indx_, val_, pntr_[i], pntr_[i + 1] - pntr_[i]);
		for (i = 0; i < dim_[0]; i++)
			if (indx_[pntr_(i)] != i) {
				std::cerr << "IC Preconditioner: diagonal not found!" << i
						<< "\n";
				exit(1);
			}
		ICFactor(pntr_, indx_, val_);

	}
	~ICPre()
	{
	}
	;

	ArrayListV<VALUE> solve(const ArrayListV<VALUE> &x) const
	{
		ArrayListV<VALUE> y(x);
		for (int i = 0; i < dim_[1]; i++) {
			Float z = y[i] / val_[pntr_[i]];
			y[i] = z;
			for (int j = pntr_[i] + 1; j < pntr_[i + 1]; j++) {
				y[indx_[j]] -= z * val_[j];
			}
		}
		for (int i = dim_[0] - 1; i >= 0; i--) {
			Float z = 0;
			for (int j = pntr_[i] + 1; j < pntr_[i + 1]; j++) {
				z += y[indx_[j]] * val_[j];
			}
			y[i] = (y[i] - z) / val_[pntr_[i]];
		}
		return y;
	}
	ArrayListV<VALUE> trans_solve(const ArrayListV<VALUE> &x) const
	{
		return this->solve(x);
	}
	MatrixSCR_<VALUE> getPre()
	{
		MatrixSCR_<VALUE> res(dim_[0], dim_[1], nz_, val_, pntr_, indx_);
		return res;
	}
	MatrixSCR_<VALUE> getPreTrans()
	{
		MatrixSCR_<VALUE> res(dim_[0], dim_[1], nz_, val_, pntr_, indx_);
		res.trans();
		return res;
	}
};

// ===========================================================
// ===========================================================

template<class VALUE>
class DiaPre
{
public:

private:
	ArrayListV<VALUE> diag_;
	int inv_;   //if inv_==1  diag_= 1. / M.val(i,j);
				//else        diag_= M.val(i,j);

public:
	//DiaPre (const CompCol_Mat_double &);
	//DiaPre(const Matrix &M, int inv);
	DiaPre(const MatrixSCR_<VALUE> &M, int inv) :
			diag_(M.iLen())
	{
		inv_ = inv;
		int flagi = -1;
		/* Find the diagonal elements */
		for (St i = 0; i < M.iLen(); i++) {
			for (St j = M.row_ptr(i); j < M.row_ptr(i + 1); j++) {
				if (M.col_ind(j) == i) {
					if (M.val(j) == 0.0) {
						flagi = i;
						break;
					} else {
						if (inv_ == 1) {
							diag_[i] = 1. / (M.val(j));
						} else {
							diag_[i] = M.val(j);
						}
						break;
					}
				}
			}
			if (diag_[i] == 0) {
				flagi = i;
				break;
			}
		}

		if (flagi != -1) {
			std::cerr << " >! Diagonal preconditioner failure.";
			std::cerr << " >! Zero detected  \n";
			exit(0);
		}
	}
	~DiaPre()
	{
	}
	;

	ArrayListV<VALUE> solve(const ArrayListV<VALUE> &x) const
	{
		ArrayListV<VALUE> y(x.size());

		for (St i = 0; i < x.size(); i++)
			y[i] = x[i] * diag_[i];

		return y;
	}
	ArrayListV<VALUE> trans_solve(const ArrayListV<VALUE> &x) const
	{
		ArrayListV<VALUE> y(x.Len());

		for (int i = 0; i < x.Len(); i++)
			y[i] = x[i] * diag_[i];

		return y;
	}

	const VALUE& diag(St i) const
	{
		return diag_[i];
	}
	VALUE& diag(St i)
	{
		return diag_[i];
	}

};

// ==============================================
// ==============================================
// ==============================================

template<class VALUE>
class ILUPre
{
public:

private:
	ArrayListV<VALUE> l_val_;
	ArrayListV<St> l_rowptr_;
	ArrayListV<St> l_colind_;
	int l_nz_;

	ArrayListV<VALUE> u_val_;
	ArrayListV<St> u_rowptr_;
	ArrayListV<St> u_colind_;
	int u_nz_;

	int dim_[2];

public:
	ILUPre(const MatrixSCR_<VALUE> &A) :
			l_val_(0), l_rowptr_(A.jLen() + 1), l_colind_(0), l_nz_(0), u_val_(
					0), u_rowptr_(A.jLen() + 1), u_colind_(0), u_nz_(0)
	{
		int i, j, k, pn, qn, rn;
		Float multiplier;

		// Copy
		dim_[0] = A.iLen();
		dim_[1] = A.jLen();

		// Get size of l and u
		for (i = 0; i < dim_[1]; i++) {
			for (j = A.row_ptr(i); j < A.row_ptr(i + 1); j++) {
				if (A.col_ind(j) < i)
					l_nz_++;
				else
					u_nz_++;
			}
		}

		l_val_.reconstruct(l_nz_);
		u_val_.reconstruct(u_nz_);
		l_colind_.reconstruct(l_nz_);
		u_colind_.reconstruct(u_nz_);

		l_rowptr_[0] = u_rowptr_[0] = 0;

		// Split up A into l and u
		for (i = 0; i < dim_[1]; i++) {
			l_rowptr_(i + 1) = l_rowptr_(i);
			u_rowptr_(i + 1) = u_rowptr_(i);

			for (j = A.row_ptr(i); j < A.row_ptr(i + 1); j++)
				if (A.col_ind(j) < i) {
					k = l_rowptr_[i + 1]++;
					l_val_[k] = A.val(j);
					l_colind_[k] = A.col_ind(j);
				} else if (A.col_ind(j) >= i) {
					k = u_rowptr_[i + 1]++;
					u_val_[k] = A.val(j);
					u_colind_[k] = A.col_ind(j);
				}
		}

		for (i = 0; i < dim_[1]; i++) {
			QSort(l_colind_, l_val_, l_rowptr_[i],
					l_rowptr_[i + 1] - l_rowptr_[i]);
			QSort(u_colind_, u_val_, u_rowptr_[i],
					u_rowptr_[i + 1] - u_rowptr_[i]);
		}

		// Factor matrix
		for (i = 1; i < dim_[0]; i++) {
			for (j = l_rowptr_[i]; j < l_rowptr_[i + 1]; j++) {
				pn = u_rowptr_[l_colind_[j]];
				multiplier = (l_val_[j] /= u_val_[pn]);

				qn = j + 1;
				rn = u_rowptr_[i];

				for (pn++;
						pn < u_rowptr_[l_colind_[j] + 1] && u_colind_[pn] < i;
						pn++) {
					while (qn < l_rowptr_[i + 1]
							&& l_colind_[qn] < u_colind_[pn])
						qn++;
					if (qn < l_rowptr_[i + 1] && u_colind_[pn] == l_colind_[qn])
						l_val_[qn] -= multiplier * u_val_[pn];
				}
				for (; pn < u_rowptr_[l_colind_[j] + 1]; pn++) {
					while (rn < u_rowptr_[i + 1]
							&& u_colind_[rn] < u_colind_[pn])
						rn++;
					if (rn < u_rowptr_[i + 1] && u_colind_[pn] == u_colind_[rn])
						u_val_[rn] -= multiplier * u_val_[pn];
				}
			}
		}
	}

	~ILUPre(void)
	{
	}
	;

	//Array solve(const Array &x) const;
	//Array trans_solve(const Array &x) const;
	MatrixSCR_<VALUE> getPreL()
	{
		MatrixSCR_<VALUE> res(dim_[0], dim_[1], l_nz_, l_val_, l_rowptr_, l_colind_);
		return res;
	}
	MatrixSCR_<VALUE> getPreU()
	{
		MatrixSCR_<VALUE> res(dim_[0], dim_[1], u_nz_, u_val_, u_rowptr_, u_colind_);
		return res;
	}
};

}

#endif /* ICPRE_H_ */
