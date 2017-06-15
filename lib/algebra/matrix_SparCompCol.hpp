/************************
 //  \file   MatrixSCC.h
 //  \brief
 //
 //  \author zhou
 //  \date   17 avr. 2014
 ***********************/
#ifndef MATRIXSPARCOMPCOL_H_
#define MATRIXSPARCOMPCOL_H_

#include "../type_define.hpp"
#include "array_list.hpp"
#include "matrix.hpp"

namespace carpio {

template<class VALUE> class MatrixSCO_;
template<class VALUE> class MatrixSCR_;

template<class VALUE>
class MatrixSCC_ {
public:
    typedef VALUE vt;
private:
    ArrayListV<VALUE> val_;             // data values (nz_ elements)
    ArrayListV<St> rowind_;    // row_ind     (nz_ elements)
    ArrayListV<St> colptr_;    // col_ptr     (dim_[1]+1 elements)

    St nz_;                   // number of nonzeros
    St dim_[2];               // number of rows, cols

public:

    MatrixSCC_(void) :
        val_(0), rowind_(0), colptr_(0), nz_(0) {
        dim_[0] = 0;
        dim_[1] = 0;
    }

    MatrixSCC_(const MatrixSCC_<vt> &S) :
        val_(S.val_), rowind_(S.rowind_), colptr_(S.colptr_), nz_(S.nz_) {
        dim_[0] = S.dim_[0];
        dim_[1] = S.dim_[1];
    }

    MatrixSCC_( St M,
                St N,  St nz,
                vt *val,  St *r,
                St *c) :
        val_(val, nz), rowind_(*r, nz), colptr_(*c, N + 1), nz_(nz) {
        dim_[0] = M;
        dim_[1] = N;
    }

    MatrixSCC_(
        St M, //
        St N, //
        St nz, //
        const ArrayListV<vt> &val,
        const ArrayListV<St> &rid,
        const ArrayListV<St> &cp) :
        val_(val), rowind_(rid), colptr_(cp), nz_(nz) {
        dim_[0] = M;
        dim_[1] = N;
    }

    MatrixSCC_(const MatrixSCR_<vt> &R) :
        val_(R.NumNonzeros()), rowind_(R.NumNonzeros()), colptr_(
            R.getjLen() + 1), nz_(R.NumNonzeros()) {
        dim_[0] = R.getiLen();
        dim_[1] = R.getjLen();

        int i, j;
        ArrayListV<int> tally(R.getjLen() + 1, 0);
        //      First pass through nonzeros.  Tally entries in each column.
        //      And calculate colptr array.
        for (i = 0; i < nz_; i++) {
            tally[R.col_ind(i)]++;
        }
        colptr_[0] = 0;
        for (j = 0; j < dim_[1]; j++)
            colptr_(j + 1) = colptr_(j) + tally(j);
        //      Make copy of colptr for use in second pass.
        tally = colptr_;
        //      Second pass through nonzeros.  Fill in index and value entries.
        int count = 0;
        for (i = 1; i <= dim_[0]; i++) {
            for (j = count; j < R.row_ptr(i); j++) {
                val_[tally(R.col_ind(j))] = R.val(j);
                rowind_[tally(R.col_ind(j))] = i - 1;
                tally[R.col_ind(count)]++;
                count++;
            }
        }
    }

    MatrixSCC_(const MatrixSCO_<vt> &CO) :
        val_(CO.NumNonzeros()), rowind_(CO.NumNonzeros()), colptr_(
            CO.getjLen() + 1), nz_(CO.NumNonzeros()) {

        dim_[0] = CO.getiLen();
        dim_[1] = CO.getjLen();

        int i, j;
        ArrayListV<int> tally(CO.getjLen() + 1, 0);
//  First pass through nonzeros.  Tally entries in each column.
//  And calculate colptr array.
        for (i = 0; i < nz_; i++) {
            tally[CO.col_ind(i)]++;
        }
        colptr_[0] = 0;
        for (j = 0; j < dim_[1]; j++) {
            colptr_[j + 1] = colptr_[j] + tally[j];
        }
//  Make copy of colptr for use in second pass.
        tally = colptr_;
//  Second pass through nonzeros.   Fill in index and value entries.
        for (i = 0; i < nz_; i++) {
            val_[tally[CO.col_ind(i)]] = CO.val(i);
            rowind_[tally[CO.col_ind(i)]] = CO.row_ind(i);
            tally[CO.col_ind(i)]++;
        }
    }

    vt& val( St i) {
        return val_(i);
    }

    St& row_ind(
        St i) {
        return rowind_(i);
    }

    St& col_ptr(
        St i) {
        return colptr_(i);
    }

    const vt& val(
        St i) const {
        return val_(i);
    }

    const  St& row_ind(
        St i) const {
        return rowind_(i);
    }

    const  St& col_ptr(
        St i) const {
        return colptr_(i);
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

    MatrixSCC_<vt>& operator=(const MatrixSCC_ &C) {
        dim_[0] = C.dim_[0];
        dim_[1] = C.dim_[1];
        nz_ = C.nz_;
        val_ = C.val_;
        rowind_ = C.rowind_;
        colptr_ = C.colptr_;
        return *this;
    }

    MatrixSCC_<vt>& newsize(
        St M,  St N,
        St nz) {
        dim_[0] = M;
        dim_[1] = N;

        nz_ = nz;
        val_.reconstruct(nz);
        rowind_.reconstruct(nz);
        colptr_.reconstruct(N + 1);
        return *this;
    }

    vt operator()( St i,
                   St j) const {
        ASSERT(i >= 0 && i < dim_[0]);
        ASSERT(j >= 0 && j < dim_[1]);
        for (int t = colptr_(j); t < colptr_(j + 1); t++) {
            if (rowind_(t) == i) {
                return val_(t);
            }
        }
        return 0.0;
    }

    ArrayListV<vt> operator*(
        const ArrayListV<vt> &x) const {
        int M = dim_[0];
        int N = dim_[1];
        //  Check for compatible dimensions:
        ASSERT(x.Len() == N);
        ArrayListV < vt > res(M);
        for (int i = 0; i < N; i++) {
            for (int j = colptr_[i]; j < colptr_[i + 1]; j++) {
                res[rowind_[j]] += x[i] * val_[j];
            }
        }
        return res;
    }

    ArrayListV<vt> transMult(
        const ArrayListV<vt> &x) const {
        int Mt = dim_[1];
        int Nt = dim_[0];
        //  Check for compatible dimensions:
        ASSERT(x.Len() == Nt);
        ArrayListT < vt > res(Mt);
        for (int i = 0; i < Nt; i++) {
            for (int j = colptr_[i]; j < colptr_[i + 1]; j++) {
                res[rowind_[j]] += x[i] * val_[j];
            }
        }
        return res;
    }

    void show( St a) const {
        if (a == 0) {
            std::cout << "RowIdx " << "ColPtr " << "Value " << "\n";
            for (int i = 0; i < dim_[1]; i++) {
                for (int ii = colptr_[i]; ii < colptr_[i + 1]; ii++) {
                    std::cout << std::scientific << rowind_[ii] << "  ";
                    std::cout << std::scientific << i << "  ";
                    std::cout << std::scientific << val_[ii] << "\n";
                }
            }
        } else {
            for (int i = 0; i < dim_[0]; i++) {
                for (int j = 0; j < dim_[1]; j++) {
                    bool flag = 0;
                    for (int ii = colptr_[j]; ii < colptr_[j + 1]; ii++) {
                        if (rowind_[ii] == i) {
                            std::cout << std::scientific << val_[ii] << "  ";
                            flag = 1;
                        }
                    }
                    if (flag == 0) {
                        std::cout << std::scientific << 0.0 << "  ";
                    }
                }
                std::cout << "\n";
            }
        }
    }
};

}
// This is the end of namespace

#endif /* MATRIXSPARCOMPCOL_H_ */
