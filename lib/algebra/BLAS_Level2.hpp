/************************
 //  \file   BLAS_Level2.h
 //  \brief
 //
 //  \author czhou
 //  \date   16 avr. 2015
 ***********************/
#ifndef BLAS_LEVEL2_H_
#define BLAS_LEVEL2_H_

#include "../type_define.hpp"
#include "arithmetic.hpp"
#include "array_list.hpp"

namespace Larus {
/**
 * \brief   GBMV  performs one of the matrix-vector operations \n
 *          y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
 *          where alpha and beta are scalars, x and y are vectors and A is an
 *          m by n band matrix, with kl sub-diagonals and ku super-diagonals.
 *
 * \param   array
 *
 * \return  int
 */
template<class VALUE>
int gbmv( //
		char trans,      //[in]
		int m,           //[in]
		int n,           //[in]
		int kl,          //[in]
		int ku,          //[in]
		VALUE alpha,     //[in]
		const VALUE* a, //[in]
		int lda,         //[in]
		const VALUE* x, //[in]
		int incx,        //[in]
		VALUE beta,      //[in]
		VALUE* y,        //[in,out]
		int incy         //[in]
		) {
	VALUE one = 1.0e+0,zero = 0.0e+0;
	//Test the input parameters.
	int info = 0;
	if (trans != 'N' && trans != 'T' && trans != 'C') {
		info = 1;
	} else if (m < 0) {
		info = 2;
	} else if (n < 0) {
		info = 3;
	} else if (kl < 0) {
		info = 4;
	} else if (ku < 0) {
		info = 5;
	} else if (lda < (kl + ku + 1)) {
		info = 8;
	} else if (incx == 0) {
		info = 10;
	} else if (incy == 0) {
		info = 13;
	}
	if (info != 0) {
		//error msg
		return -1;
	}
	if ((m == 0) || (n == 0) || ((alpha == zero) && (beta == one)))
		return;

	// Set LENX and LENY, the lengths of the vectors x and y, and set
	// up the start points in X and Y.
	int lenx,leny,kx,ky;
	if (trans == 'N') {
		lenx = n;
		leny = m;
	} else {
		lenx = m;
		leny = n;
	}
	if (incx > 0) {
		kx = 1;
	} else {
		kx = (lenx - 1) * incx;
	}
	if (incy > 0) {
		ky = 1;
	} else {
		ky = (leny - 1) * incy;
	}
	// First form y := beta*y.
	if (beta != one) {
		if (incy == 1) {
			if (beta == zero) {
				for (int i = 0; i < leny; ++i) {
					y[i] = zero;
				}
			} else {
				for (int i = 0; i < leny; ++i) {
					y[i] = beta * y[i];
				}
			}
		} else {
			int iy = ky;
			if (beta == zero) {
				for (int i = 0; i < leny; ++i) {
					y[iy] = zero;
					iy = iy + incy;
				}
			} else {
				for (int i = 0; i < leny; ++i) {
					y[iy] = beta * y[iy];
					iy = iy + incy;
				}
			}
		}
	}
	if (alpha == zero)
		return;
	int kup1 = ku + 1;
	if (trans == 'N') {
		// Form y := alpha*A*x + y.

	} else {
		// Form y := alpha*A**T*x + y.

	}
}

}

#endif /* BLAS_LEVEL2_H_ */
