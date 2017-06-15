/************************
 //  \file   MatrixSCR.h
 //  \brief
 // 
 //  \author zhou
 //  \date   25 avr. 2014 
 ***********************/
#ifndef MATRIXSPARCOMPROW_H_
#define MATRIXSPARCOMPROW_H_

#include "../type_define.hpp"
#include "array_list.hpp"
#include "matrix.hpp"
#include "arithmetic.hpp"

namespace carpio {

template<class VALUE> class MatrixSCO_;
template<class VALUE> class MatrixSCC_;
/*
 * Example:
 * row_prt() 0        3        6        9   10      12
 * val()     1  2  3, 4  5  6, 7  8  9, 10, 11  12,
 * col_ind() 0  1  4  0  1  2  1  2  4  3   0   4
 */
template<class VALUE>
class MatrixSCR_ {
public:
	typedef VALUE Vt;
private:
	ArrayListV<Vt> val_;       // data values (nz_ elements)
	ArrayListV<St> rowptr_;    // row_ptr     (dim_[0]+1 elements)
	ArrayListV<St> colind_;    // col_ind     (nz_ elements)

	St nz_;                   // number of nonzeros
	St dim_[2];               // number of rows, cols

public:
	MatrixSCR_(void) :
			val_(0), rowptr_(0), colind_(0), nz_(0) {
		dim_[0] = 0;
		dim_[1] = 0;
	}
	//MatrixSCC(const Matrix &M);
	MatrixSCR_(const MatrixSCR_<Vt> &S) :
			val_(S.val_), rowptr_(S.rowptr_), colind_(S.colind_), nz_(S.nz_) {
		dim_[0] = S.dim_[0];
		dim_[1] = S.dim_[1];
	}
	MatrixSCR_(const MatrixSCC_<Vt> &C) :
			val_(C.NumNonzeros()),  //
			rowptr_(C.iLen() + 1),  //
			colind_(C.NumNonzeros()),  //
			nz_(C.NumNonzeros())  //
	{
		dim_[0] = C.getiLen();
		dim_[1] = C.getjLen();

		St i, j;

		ArrayListV<St> tally(C.iLen() + 1, 0);
		//      First pass through nonzeros.  Tally entries in each row.
		//      And calculate rowptr array.
		for (i = 0; i < nz_; i++) {
			tally[C.row_ind(i)]++;
		}
		rowptr_[0] = 0;
		for (j = 0; j < dim_[0]; j++)
			rowptr_(j + 1) = rowptr_(j) + tally(j);
		//      Make copy of rowptr for use in second pass.
		tally = rowptr_;
		//      Second pass through nonzeros.   Fill in index and value entries.
		St count = 0;
		for (i = 1; i <= dim_[1]; i++) {
			for (j = count; j < C.col_ptr(i); j++) {
				val_[tally(C.row_ind(j))] = C.val(j);
				colind_[tally(C.row_ind(j))] = i - 1;
				tally[C.row_ind(count)]++;
				count++;
			}
		}

	}

	MatrixSCR_(const MatrixSCO_<Vt> &CO) :
			val_(CO.NumNonzeros()), rowptr_(CO.iLen() + 1), colind_(
					CO.NumNonzeros()), nz_(CO.NumNonzeros()) {
		dim_[0] = CO.iLen();
		dim_[1] = CO.jLen();

		St i;
		ArrayListV<St> tally(CO.iLen() + 1, 0);
		//      First pass through nonzeros.  Tally entries in each row.
		//      And calculate rowptr array.
		for (i = 0; i < nz_; i++) {
			tally[CO.row_ind(i)]++;
		}
		rowptr_(0) = 0;
		for (i = 0; i < dim_[0]; i++) {
			rowptr_(i + 1) = rowptr_(i) + tally(i);
		}
		//      Make copy of rowptr for use in second pass.
		tally = rowptr_;
		//      Second pass through nonzeros.   Fill in index and value entries.
		for (i = 0; i < nz_; i++) {
			val_[tally(CO.row_ind(i))] = CO.val(i);
			colind_[tally(CO.row_ind(i))] = CO.col_ind(i);
			tally[CO.row_ind(i)]++;
		}

	}

	MatrixSCR_(St M, St N, St nz, Vt *val, St *r, St *c) :
			val_(val, nz), rowptr_(*r, M + 1), colind_(*c, nz), nz_(nz) {
		dim_[0] = M;
		dim_[1] = N;
	}

	MatrixSCR_(St M, St N, St nz, const ArrayListV<Vt> &val,
			const ArrayListV<St> &r, const ArrayListV<St> &c) :
			val_(val), rowptr_(r), colind_(c), nz_(nz) {
		dim_[0] = M;
		dim_[1] = N;
	}

	Vt* get_nz_pointer() {
		return this->val_.getPointer();
	}

	Vt& val(St i) {
		return val_(i);
	}

	St& col_ind(St i) {
		return colind_(i);
	}

	St& row_ptr(St i) {
		return rowptr_(i);
	}

	const Vt& val(St i) const {
		return val_(i);
	}

	const St& col_ind(St i) const {
		return colind_(i);
	}

	const St& row_ptr(St i) const {
		return rowptr_(i);
	}

	St size() const {
		return dim_[0] * dim_[1];
	}

	St iLen() const {
		return dim_[0];
	}
	St jLen() const {
		return dim_[1];
	}
	/*
	 * make it looks like ublas
	 */
	St size1() const {
		return dim_[0];
	}
	St size2() const {
		return dim_[1];
	}

	St NumNonzeros() const {
		return nz_;
	}

	Vt max() const {
		return val_.findMax();
	}

	Vt min() const {
		return val_.findMin();
	}

	Vt sum() const {
		return val_.sum();
	}
	/*
	 * returns the sum along dimension dim.
	 * For example, if A is a matrix,
	 * then sum_row() is a column vector
	 * containing the sum of each row.
	 */
	ArrayListV<Vt> sum_row() const {
		ArrayListV<Vt>  res(this->iLen());
		for (St i = 0; i < this->iLen(); ++i) {
			res[i] = 0;
			for (St j = this->row_ptr(i); j < this->row_ptr(i + 1); ++j) {
				res[i] += this->val(j);
			}
		}
		return res;
	}

	MatrixSCR_<Vt>& operator=(const MatrixSCR_<Vt> &R) {
		dim_[0] = R.dim_[0];
		dim_[1] = R.dim_[1];
		nz_ = R.nz_;
		val_ = R.val_;
		rowptr_ = R.rowptr_;
		colind_ = R.colind_;
		return *this;
	}

	MatrixSCR_<Vt>& newsize(St M, St N, St nz) {
		dim_[0] = M;
		dim_[1] = N;
		nz_ = nz;
		val_.reconstruct(nz);
		rowptr_.reconstruct(M + 1);
		colind_.reconstruct(nz);
		return *this;
	}

	Vt operator()(St i, St j) const {
		ASSERT(i >= 0 && i < dim_[0]);
		ASSERT(j >= 0 && j < dim_[1]);
		for (St t = rowptr_(i); t < rowptr_(i + 1); t++) {
			if (colind_(t) == j)
				return val_(t);
		}
		return 0.0;
	}

	ArrayListV<Vt> operator*(const ArrayListV<Vt> &x) const {
		St M = dim_[0];
		St N = dim_[1];
		//  Check for compatible dimensions:
		ASSERT(x.size() == N);

		ArrayListV<Vt> res(M);
		for (St i = 0; i < M; ++i) {
			for (St j = rowptr_[i]; j < rowptr_[i + 1]; ++j) {
				res[i] += x[colind_[j]] * val_[j];
			}
		}
		return res;
	}

	ArrayListV<Vt> transMult(const ArrayListV<Vt> &x) const {
		St Mt = dim_[1];
		St Nt = dim_[0];
		//  Check for compatible dimensions:
		ASSERT(x.size() == Nt);
		ArrayListV<Vt> res(Mt);
		for (St i = 0; i < Mt; ++i) {
			for (St j = rowptr_[i]; j < rowptr_[i + 1]; ++j) {
				res[i] += x[colind_[j]] * val_[j];
			}
		}
		return res;
	}

	void trans() {
		MatrixSCC_<Vt> C(dim_[1], dim_[0], nz_, val_, colind_, rowptr_);
		dim_[0] = C.getiLen();
		dim_[1] = C.getjLen();

		St i, j;

		ArrayListV<St> tally(C.getiLen() + 1, 0);
		//      First pass through nonzeros.  Tally entries in each row.
		//      And calculate rowptr array.
		for (i = 0; i < nz_; i++) {
			tally[C.row_ind(i)]++;
		}
		rowptr_[0] = 0;
		for (j = 0; j < dim_[0]; j++)
			rowptr_(j + 1) = rowptr_(j) + tally(j);
		//      Make copy of rowptr for use in second pass.
		tally = rowptr_;
		//      Second pass through nonzeros.   Fill in index and vt entries.
		St count = 0;
		for (i = 1; i <= dim_[1]; i++) {
			for (j = count; j < C.col_ptr(i); j++) {
				val_[tally(C.row_ind(j))] = C.val(j);
				colind_[tally(C.row_ind(j))] = i - 1;
				tally[C.row_ind(count)]++;
				count++;
			}
		}
	}
	/*
	 *  get a new matrix. It is the Transpose of this matrix
	 *  Expensive function
	 */
	MatrixSCR_<Vt> get_trans() const {
		MatrixSCR_ res(dim_[0], dim_[1], nz_, val_, rowptr_, colind_);
		res.trans();
		return res;
	}

	/*
	 *  check this is a diagonally dominant matrix or not
	 */
	bool is_diagonally_dominant() const {
		St Mt = dim_[1];
		//  Check for compatible dimensions:
		for (St i = 0; i < Mt; ++i) {
			Vt sum_ = 0;
			Vt vdiag = 0;
			for (St j = rowptr_[i]; j < rowptr_[i + 1]; ++j) {
				St col_idx = this->colind_(j);
				if (i == col_idx) {
					vdiag = Abs(this->val_[j]);
				} else {
					sum_ += Abs(this->val_[j]);
				}
			}
			if (sum_ > vdiag) {
				return false;
			}
		}
		return true;
	}

	void show(St a) const {
		if (a == 0) {
			std::cout << "RowPtr " << "ColInd " << "vt " << std::endl;
			for (St i = 0; i < dim_[0]; i++) {
				for (St ii = rowptr_[i]; ii < rowptr_[i + 1]; ii++) {
					std::cout << std::scientific << i << "  ";
					std::cout << std::scientific << colind_[ii] << "  ";
					std::cout << std::scientific << val_[ii] << "\n";
				}
			}
		} else {
			for (St i = 0; i < dim_[0]; i++) {
				for (St j = 0; j < dim_[1]; j++) {
					bool flag = 0;
					for (St ii = rowptr_[i]; ii < rowptr_[i + 1]; ii++) {
						if (colind_[ii] == j) {
							std::cout << std::scientific << val_[ii] << "  ";
							flag = 1;
						}
					}
					if (flag == 0) {
						std::cout << std::scientific << 0.0 << "  ";
					}
				}
				std::cout << std::endl;
			}
		}
	}
};

}
// This is the end of namespace

#endif /* MATRIXSPARCOMPROW_H_ */
