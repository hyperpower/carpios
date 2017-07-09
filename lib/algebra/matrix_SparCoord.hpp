/************************
 //  \file   MatrixSCO.h
 //  \brief
 // 
 //  \author zhou
 //  \date   16 avr. 2014 
 ***********************/
#ifndef MATRIXSPARCOORD_H_
#define MATRIXSPARCOORD_H_

#include "algebra_define.hpp"
#include "array_list.hpp"
#include "matrix.hpp"

namespace carpio {

template<class VALUE> class MatrixSCC_;
template<class VALUE> class MatrixSCR_;

template<class VALUE>
class MatrixSCO_ {
public:
	typedef VALUE vt;
	typedef MatrixV<VALUE> Matrix;
private:
	ArrayListV<vt> val_;    // data values (nz_ elements)
	ArrayListV<St> rowind_;    // row_ind (nz_ elements)
	ArrayListV<St> colind_;    // col_ind (nz_ elements)

	St nz_;                   // number of nonzeros
	St dim_[2];               // number of rows, cols
public:

	MatrixSCO_(void) :
			val_(0), rowind_(0), colind_(0), nz_(0) {
		dim_[0] = 0;
		dim_[1] = 0;
	}

	MatrixSCO_(const Matrix &M) {
		dim_[0] = M.size_i();
		dim_[1] = M.size_j();
		// count non-zero
		St countnz = 0;
		for (St i = 0; i < M.size_i(); i++) {
			for (St j = 0; j < M.size_j(); j++) {
				if (M[i][j] != 0) {
					countnz++;
				}
			}
		}
		nz_ = countnz;
		val_.reconstruct(nz_);
		rowind_.reconstruct(nz_);
		colind_.reconstruct(nz_);
		St idx = 0;
		for (St i = 0; i < M.size_i(); i++) {
			for (St j = 0; j < M.size_j(); j++) {
				if (M[i][j] != 0) {
					val_[idx] = M[i][j];
					rowind_[idx] = i;
					colind_[idx] = j;
					idx++;
				}
			}
		}
	}

	MatrixSCO_(const MatrixSCO_ &S) :
			val_(S.val_), rowind_(S.rowind_), colind_(S.colind_), nz_(S.nz_) {
		dim_[0] = S.dim_[0];
		dim_[1] = S.dim_[1];
	}

	MatrixSCO_(St M, St N, St nz, vt *val, St *r, St *c) :
			val_(val, nz), rowind_((*r), nz), colind_(*c, nz), nz_(nz) {
		dim_[0] = M;
		dim_[1] = N;
	}

	MatrixSCO_(const MatrixSCR_<vt> &R) :
			val_(R.NumNonzeros()), rowind_(R.NumNonzeros()), colind_(
					R.NumNonzeros()), nz_(R.NumNonzeros()) {
		dim_[0] = R.iLen();
		dim_[1] = R.jLen();
		St count = 0;
		St i, j;
//  Loop through rows...
		for (i = 1; i <= dim_[0]; i++) {
			for (j = count; j < R.row_ptr(i); j++) {
				val_[count] = R.val(count);
				colind_[count] = R.col_ind(count);
				rowind_[count] = i - 1;
				count++;
			}
		}
	}

	MatrixSCO_(const MatrixSCC_<vt> &C) :
			val_(C.NumNonzeros()), rowind_(C.NumNonzeros()), colind_(
					C.NumNonzeros()), nz_(C.NumNonzeros()) {
		dim_[0] = C.getiLen();
		dim_[1] = C.getjLen();

		int count = 0;
		int i, j;
//  Loop through columns...
		for (j = 1; j <= dim_[1]; j++) {
			for (i = count; i < C.col_ptr(j); i++) {
				val_[count] = C.val(count);
				rowind_[count] = C.row_ind(count);
				colind_[count] = j - 1;
				count++;
			}
		}
	}

	MatrixSCO_<vt>& operator=(const MatrixSCO_<vt> &C) {
		dim_[0] = C.dim_[0];
		dim_[1] = C.dim_[1];
		nz_ = C.nz_;
		val_ = C.val_;
		rowind_ = C.rowind_;
		colind_ = C.colind_;
		return *this;
	}

	MatrixSCO_<vt>& newsize(St M, St N, St nz) {
		dim_[0] = M;
		dim_[1] = N;
		nz_ = nz;
		val_.reconstruct(nz);
		rowind_.reconstruct(nz);
		colind_.reconstruct(nz);
		return *this;
	}

	MatrixSCO_<vt>& resize(St M, St N, St nz) {
		dim_[0] = M;
		dim_[1] = N;
		nz_ = nz;
		val_.resize(nz);
		rowind_.resize(nz);
		colind_.resize(nz);
		return *this;
	}

//slow---

	vt operator()(St i, St j) const {
		ASSERT(i >= 0 && i < dim_[0]);
		ASSERT(j >= 0 && j < dim_[1]);
		for (St t = 0; t < nz_; t++) {
			if (rowind_(t) == i && colind_(t) == j) {
				return val_(t);
			}
		}
		return 0.0;
	}

	ArrayListV<vt> operator*(const ArrayListV<vt> &x) const {
		St M = dim_[0];
		St N = dim_[1];
//  Check for compatible dimensions:
		ASSERT(x.size() == N);
		ArrayListV<vt> res(M);
		for (St i = 0; i < nz_; i++) {
			res[rowind_[i]] += x[colind_[i]] * val_[i];
		}
		return res;
	}

	ArrayListV<vt> transMult(const ArrayListV<vt> &x) const {
		St tM = dim_[1];
		St tN = dim_[0];
		ASSERT(!(x.Len() == tN));
		ArrayListV<vt> res(tM);
		for (St i = 0; i < nz_; i++) {
			res[colind_[i]] += x[rowind_[i]] * val_[i];
		}
		return res;
	}

	vt max() const {
		return val_.findMax();
	}

	vt min() const {
		return val_.findMin();
	}

	void show(St a) const {
		//if a==0 show matrix in Coordinate formate
		if (a == 0) {
			std::cout << "RowIdx " << "ColIdx " << "Value " << std::endl;
			for (St i = 0; i < nz_; i++) {
				std::cout << std::scientific << rowind_[i] << "  ";
				std::cout << std::scientific << colind_[i] << "  ";
				std::cout << std::scientific << val_[i] << "\n";
			}
		} else {
			for (St i = 0; i < this->iLen(); i++) {
				for (St j = 0; j < this->jLen(); j++) {
					bool isnz = 0;
					for (St t = 0; t < nz_; t++) {
						if (rowind_(t) == i && colind_(t) == j) {
							std::cout << std::scientific << val_(t) << "  ";
							isnz = 1;
						}
					}
					if (isnz == 0) {
						std::cout << std::scientific << 0.0 << "  ";
					}
				}
				std::cout << std::endl;
			}
		}

	}

	vt& val(St i) {
		return val_(i);
	}

	St& row_ind(St i) {
		return rowind_(i);
	}

	St& col_ind(St i) {
		return colind_(i);
	}

	const vt& val(St i) const {
		return val_(i);
	}

	const St& row_ind(St i) const {
		return rowind_(i);
	}

	const St& col_ind(St i) const {
		return colind_(i);
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

	St NumNonzeros() const {
		return nz_;
	}

};

}
// This is the end of namespace

#endif /* MATRIXSPARCOORD_H_ */
