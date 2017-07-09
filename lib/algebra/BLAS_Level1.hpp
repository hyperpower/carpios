/*
 * BLAS_Level1.h
 *
 *  Created on: Apr 14, 2015
 *      Author: zhou
 */

#ifndef _BLAS_LEVEL1_H_
#define _BLAS_LEVEL1_H_

#include "algebra_define.hpp"
#include "arithmetic.hpp"
#include "array_list.hpp"
#include <cmath>

namespace carpio {
//Function ROTG =================================
//construct givens plane rotation.
template<class VALUE>
int rotg(          //
		VALUE& sa, //
		VALUE& sb, //
		VALUE& c,  //
		VALUE& s)  //
		{
	VALUE r,roe,scale,z;
	roe = sb;
	if (Abs(sa) > Abs(sb)) {
		roe = sa;
	}
	scale = Abs(sa) + Abs(sb);
	if (scale == 0.0) {
		c = 1.0;
		s = 0.0;
		r = 0.0;
		z = 0.0;
	} else {
		r = scale
				* sqrt(
						(sa / scale) * (sa / scale)
								+ (sb / scale) * (sb / scale));
		r = sign2(roe) * r;
		c = sa / r;
		s = sb / r;
		z = 1.0;
		if (Abs(sa) > Abs(sb)) {
			z = s;
		}
		if (Abs(sb) > Abs(sa) && c != 0.0) {
			z = 1.0 / c;
		}
	}
	sa = r;
	sb = z;
	return 1;
} // << -------------------------------------------

//Function ROTG =================================
//construct givens plane rotation.
//SPARAM is array, dimension 5
//SPARAM(0)=SFLAG
//SPARAM(1)=SH11
//SPARAM(2)=SH21
//SPARAM(3)=SH12
//SPARAM(4)=SH22
template<class VALUE>
int srotmg( //
		VALUE& sd1,    //[in,out]
		VALUE& sd2,    //[in,out]
		VALUE& sx1,    //[in,out]
		VALUE& sy1,    //[in]
		VALUE* sparam) //[in,out] is array size=5
		{
	// .. Local Scalars ..
	VALUE sflag,sh11,sh12,sh21,sh22,sp1,sp2,sq1,sq2,stemp,su;
	// .. Data statements ..
	Float zero = 0.0,one = 1.0,two = 2.0;
	Float gam = 4096.e0,gamsq = 1.67772e7,rgamsq = 5.96046e-8;

	if (sd1 < zero) {
		// GO ZERO-H-D-AND-SX1..
		sflag = -one;
		sh11 = zero;
		sh12 = zero;
		sh21 = zero;
		sh22 = zero;
		//
		sd1 = zero;
		sd2 = zero;
		sx1 = zero;
	} else {
		// CASE - SD1 - NONNEGATIVE
		sp2 = sd2 * sy1;
		if (sp2 == zero) {
			sflag = -two;
			sparam[1] = sflag;
			return 0;
		}
		//REGULAR-CASE..
		sp1 = sd1 * sx1;
		sq2 = sp2 * sy1;
		sq1 = sp1 * sx1;
		//
		if (Abs(sq1) > Abs(sq2)) {
			sh21 = -sy1 / sx1;
			sh12 = sp2 / sp1;
			//
			su = one - sh12 * sh21;
			//
			if (su > zero) {
				sflag = zero;
				sd1 = sd1 / su;
				sd2 = sd2 / su;
				sx1 = sx1 * su;
			}
		} else {
			if (sq2 < zero) {
				// GO ZERO-H-D-AND-SX1..
				sflag = -one;
				sh11 = zero;
				sh12 = zero;
				sh21 = zero;
				sh22 = zero;
				//
				sd1 = zero;
				sd2 = zero;
				sx1 = zero;
			} else {
				sflag = one;
				sh11 = sp1 / sp2;
				sh22 = sx1 / sy1;
				su = one + sh11 * sh22;
				stemp = sd2 / su;
				sd2 = sd1 / su;
				sd1 = stemp;
				sx1 = sy1 * su;
			}
		}
		// PROCESURE..SCALE-CHECK
		if (sd1 != zero) {
			while ((sd1 <= rgamsq) || (sd1 >= gamsq)) {
				if (sflag == zero) {
					sh11 = one;
					sh22 = one;
					sflag = -one;
				} else {
					sh21 = -one;
					sh12 = one;
					sflag = -one;
				}
				if (sd1 <= rgamsq) {
					sd1 = sd1 * gam * gam;
					sx1 = sx1 / gam;
					sh11 = sh11 / gam;
					sh12 = sh12 / gam;
				} else {
					sd1 = sd1 / (gam * gam);
					sx1 = sx1 * gam;
					sh11 = sh11 * gam;
					sh12 = sh12 * gam;
				}
			}
		}

		if (sd2 != zero) {
			while ((Abs(sd2) <= rgamsq) || (Abs(sd2) >= gamsq))
				if (sflag == zero) {
					sh11 = one;
					sh22 = one;
					sflag = -one;
				} else {
					sh21 = -one;
					sh12 = one;
					sflag = -one;
				}
			if (Abs(sd2) <= rgamsq) {
				sd2 = sd2 * gam * gam;
				sh21 = sh21 / gam;
				sh22 = sh22 / gam;
			} else {
				sd2 = sd2 / (gam * gam);
				sh21 = sh21 * gam;
				sh22 = sh22 * gam;
			}
		}
	}
	if (sflag < zero) {
		sparam[1] = sh11;
		sparam[2] = sh21;
		sparam[3] = sh12;
		sparam[4] = sh22;
	} else if (sflag == zero) {
		sparam[2] = sh21;
		sparam[3] = sh12;
	} else {
		sparam[1] = sh11;
		sparam[4] = sh22;
	}
	sparam[0] = sflag;
	return 1;
}	// << ------------------------------------------

//Function ROT =================================
//applies a plane rotation.
template<class VALUE>
int rot(	//
		int n,     //
		VALUE* sx, //
		int incx,  //
		VALUE* sy, //
		int incy,  //
		const VALUE& c,  //
		const VALUE& s)  //
		{
//.. Local Scalars ..
	VALUE stemp;
//
	if (n <= 0)
		return 0;
	if (incx == 1 && incy == 1) {
		//code for both increments equal to 1
		for (int i = 0; i < n; ++i) {
			stemp = c * sx[i] + s * sy[i];
			sy[i] = c * sy[i] - s * sx[i];
			sx[i] = stemp;
		}
	} else {
		// code for unequal increments or equal increments not equal
		// to 1
		int ix,iy;
		ix = 0;
		iy = 0;
		if (incx < 0)
			ix = (1 - n) * incx;
		if (incy < 0)
			iy = (1 - n) * incy;
		for (int i = 0; i < n; ++i) {
			stemp = c * sx[ix] + s * sy[iy];
			sy[iy] = c * sy[iy] - s * sx[ix];
			sx[ix] = stemp;
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return 1;
}		// << ------------------------------------------

/**
 * \brief   interchanges two vectors.
 *
 * \param   array
 *
 * \return  int
 */
template<class VALUE>
int swap( //
		int n,        //size of the array, sx.size==sy.size
		VALUE* sx,    //
		int incx, //
		VALUE* sy,    //
		int incy  //
		) {
	VALUE stemp;
	int ix,iy,m,mp1;
	if (n <= 0)
		return 0;
	if (incx == 0 || incy == 0)
		return 0;

	if (incx == 1 && incy == 1) {
		//increase = 1
		m = (n - 1) % 3;
		for (int i = 0; i <= m; ++i) {
			stemp = sx[i];
			sx[i] = sy[i];
			sy[i] = stemp;
		}
		if (n <= 3) {
			return 1;
		}
		mp1 = m + 1;
		for (int i = mp1; i < n; i += 3) {  //unrolled loop
			stemp = sx[i];
			sx[i] = sy[i];
			sy[i] = stemp;
			stemp = sx[i + 1];
			sx[i + 1] = sy[i + 1];
			sy[i + 1] = stemp;
			stemp = sx[i + 2];
			sx[i + 2] = sy[i + 2];
			sy[i + 2] = stemp;
		}
	} else {
		ix = 0;
		iy = 0;
		if (incx < 0)
			ix = (1 - n) * incx;
		if (incy < 0)
			iy = (1 - n) * incy;
		for (int i = 0; i < n; ++i) {
			stemp = sx[ix];
			sx[ix] = sy[iy];
			sy[iy] = stemp;
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return 1;
}  // << ------------------------------------------
//

/**
 * \brief   sx = sa * sx \n
 *          scales a vector by a constant.
 *
 * \param   array
 *
 * \return  int
 */
template<class VALUE>
int scal( //
		int n,        //size of the array, sx.size==sy.size
		VALUE sa,     //
		VALUE* sx,    //
		int incx     //
		) {
	if (n < 0 || incx <= 0) {
		return 0;
	}
	if (incx == 1) {
		//increase = 1
		int LN = 5;
		int m = (n - 1) % LN;
		for (int i = 0; i <= m; ++i) {
			sx[i] = sa * sx[i];
		}
		if (n <= LN) {
			return 1;
		}
		int mp1 = m + 1;
		for (int i = mp1; i < n; i += LN) {  //unrolled loop
			sx[i] = sa * sx[i];
			sx[i + 1] = sa * sx[i + 1];
			sx[i + 2] = sa * sx[i + 2];
			sx[i + 3] = sa * sx[i + 3];
			sx[i + 4] = sa * sx[i + 4];
		}
	} else {
		int nincx = n * incx;
		for (int i = 0; i < nincx; i += incx) {
			sx[i] = sa * sx[i];
		}
	}
}

/**
 * \brief   sum (abs())
 *
 * \param   array
 *
 * \return  int
 */
template<class VALUE>
VALUE asum( //
		int n,        //size of the array, sx.size==sy.size
		const VALUE* sx,    //
		int incx     //
		) {
	VALUE asum = 0.0e0;
	if (n < 0 || incx <= 0) {
		return asum;
	}
	if (incx == 1) {
		//increase = 1
		int LN = 6;
		int m = (n - 1) % LN;
		for (int i = 0; i <= m; ++i) {
			asum = asum + Abs(sx[i]);
		}
		if (n <= LN) {
			return asum;
		}
		int mp1 = m + 1;
		for (int i = mp1; i < n; i += LN) {  //unrolled loop
			asum = asum + Abs(sx[i]) + Abs(sx[i + 1]) + Abs(sx[i + 2])
					+ Abs(sx[i + 3]) + Abs(sx[i + 4]) + Abs(sx[i + 5]);
		}
	} else {
		int nincx = n * incx;
		for (int i = 0; i < nincx; i += incx) {
			asum = asum + Abs(sx[i]);
		}
	}
	return asum;
}  // << ------------------------------------------

/**
 * \brief   max(abs())
 *
 * \param   array
 *
 * \return  int
 */
template<class VALUE>
VALUE amax( //
		int n,        //size of the array, sx.size
		const VALUE* sx,    //
		int incx     //
		) {
	VALUE max = 0.0e0;
	if (n < 0 || incx <= 0) {
		return max;
	}
	if (incx == 1) {
		//increase = 1
		max = Abs(sx[0]);
		for (int i = 1; i < n; ++i) {
			if (Abs(sx[i]) > max) {
				max = Abs(sx[i]);
			}
		}
	} else {
		int nincx = n * incx;
		max = Abs(sx[0]);
		for (int i = 0; i < nincx; i += incx) {
			if (Abs(sx[i]) > max) {
				max = Abs(sx[i]);
			}
		}
	}
	return max;
}  // << ------------------------------------------

/**
 * \brief   copies a vector, x, to a vector, y.
 *
 * \param   array
 *
 * \return  int
 */
template<class VALUE>
int copy( //
		int n,        //size of the array, sx.size==sy.size
		VALUE* sx,     //
		int incx,      //
		VALUE* sy,    //
		int incy      //
		) {
	int ix,iy,m,mp1;
	if (n <= 0)
		return 0;
	if (incx == 0 || incy == 0)
		return 0;

	if (incx == 1 && incy == 1) {
		//increase = 1
		int LN = 7;
		m = (n - 1) % LN;
		for (int i = 0; i <= m; ++i) {
			sy[i] = sx[i];
		}
		if (n <= LN) {
			return 1;
		}
		mp1 = m + 1;
		for (int i = mp1; i < n; i += LN) {  //unrolled loop
			sy[i] = sx[i];
			sy[i + 1] = sx[i + 1];
			sy[i + 2] = sx[i + 2];
			sy[i + 3] = sx[i + 3];
			sy[i + 4] = sx[i + 4];
			sy[i + 5] = sx[i + 5];
			sy[i + 6] = sx[i + 6];
		}
	} else {
		ix = 0;
		iy = 0;
		if (incx < 0)
			ix = (1 - n) * incx;
		if (incy < 0)
			iy = (1 - n) * incy;
		for (int i = 0; i < n; ++i) {
			sy[i] = sx[i];
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return 1;
}  // << ------------------------------------------
//

/**
 * \brief   copies a vector, x, to a vector, y.
 *
 * \param   array
 *
 * \return  int
 */
template<class VALUE>
VALUE dot( //
		int n,        //size of the array, sx.size==sy.size
		const VALUE* sx,    //
		int incx,     //
		const VALUE* sy,    //
		int incy      //
		) {
	VALUE dot = 0.0e0;
	int ix,iy,m,mp1;
	if (n <= 0)
		return dot;
	if (incx == 0 || incy == 0)
		return dot;

	if (incx == 1 && incy == 1) {
		//increase = 1
		int LN = 5;
		m = (n - 1) % LN;
		for (int i = 0; i <= m; ++i) {
			dot = dot + sx[i] * sy[i];
		}
		if (n <= LN) {
			return dot;
		}
		mp1 = m + 1;
		for (int i = mp1; i < n; i += LN) {  //unrolled loop
			dot = dot + sx[i] * sy[i] + sx[i + 1] * sy[i + 1]
					+ sx[i + 2] * sy[i + 2] + sx[i + 3] * sy[i + 3]
					+ sx[i + 4] * sy[i + 4];
		}
	} else {
		ix = 0;
		iy = 0;
		if (incx < 0)
			ix = (1 - n) * incx;
		if (incy < 0)
			iy = (1 - n) * incy;
		for (int i = 0; i < n; ++i) {
			dot = dot + sx[i] * sy[i];
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return dot;
}  // << ------------------------------------------
//

/**
 * \brief   SNRM2 returns the euclidean norm of a vector via the function
 *          name, so that
 *          SNRM2 := sqrt( x'*x ).
 *
 * \param   array
 *
 * \return  int
 */
template<class VALUE>
VALUE nrm2( //
		int n,        //size of the array, sx.size==sy.size
		const VALUE* x,    //
		int incx    //
		) {
	VALUE zero = 0.0e+0;
	VALUE norm = zero;
	if (n <= 0 || incx <= 0) {
		return norm;
	} else if (n == 1) {
		norm = Abs(x[0]);
	} else {
		for (int i = 0; i < (n - 1) * incx; i += incx) {
			if (x[i] != zero) {
				norm += x[i] * x[i];
			}
		}
		norm = sqrt(norm);
	}
	return norm;
}  // << ------------------------------------------
//

/**
 * \brief   NRM returns the p norm of a vector via the function
 *          name, so that
 *          SNRM2 := ( x ).
 *
 * \param   array
 *
 * \return  int
 */
template<class VALUE, class VALUE2>
VALUE nrmp( //
		int n,        //size of the array, sx.size==sy.size
		const VALUE* x,    //
		VALUE2 p,         //
		int incx       //
		) {
	if (p == 1.0) {
		return asum(n, x, incx);
	}
	if (p == 2.0) {
		return nrm2(n, x, incx);
	}
	VALUE zero = 0.0e+0;
	VALUE norm = zero;
	if (n <= 0 || incx <= 0) {
		return norm;
	} else if (n == 1) {
		norm = Abs(x[0]);
	} else {
		for (int i = 0; i < (n - 1) * incx; i += incx) {
			if (x[i] != zero) {
				norm += pow(Abs(x[i]), p);
			}
		}
		norm = pow(norm, 1.0 / double(p));
	}
	return norm;
}  // << ------------------------------------------
//

//==============================================
//==============================================
template<typename TYPE>
int swap( //
		ArrayListV<TYPE>& asx, //[in,out]
		ArrayListV<TYPE>& asy  //[in,out]
		) {
	if (asx.size() != asy.size()) {
		return 0;
	}
	int n = asx.size();
	TYPE* sx = asx.getPointer();
	TYPE* sy = asy.getPointer();
	return swap(n, sx, 1, sy, 1);
}

template<typename TYPE>
int copy( //
		ArrayListV<TYPE>& asx, //[in,out]
		ArrayListV<TYPE>& asy  //[in,out]
		) {
	if (asx.size() != asy.size()) {
		return 0;
	}
	int n = asx.size();
	TYPE* sx = asx.getPointer();
	TYPE* sy = asy.getPointer();
	return copy(n, sx, 1, sy, 1);
}

template<typename TYPE>
TYPE dot( //
		const ArrayListV<TYPE>& asx, //[in,out]
		const ArrayListV<TYPE>& asy  //[in,out]
		) {
	if (asx.size() != asy.size()) {
		return 0;
	}
	int n = asx.size();
	const TYPE* sx = asx.getPointer();
	const TYPE* sy = asy.getPointer();
	return dot(n, sx, 1, sy, 1);
}

template<typename TYPE>
int scal( //
		const TYPE& sa,       //[in]
		ArrayListV<TYPE>& asx //[in,out]
		) {
	int n = asx.size();
	TYPE* sx = asx.getPointer();
	return scal(n, sa, sx, 1);
}

template<class TYPE>
TYPE nrm2( //
		const ArrayListV<TYPE>& ax  //[in]
		) {
	int n = ax.size();
	const TYPE* x = ax.getPointer();
	return nrm2(n, x, 1);
}

template<class TYPE>
TYPE nrm1( //
		const ArrayListV<TYPE>& ax  //[in]
		) {
	int n = ax.size();
	const TYPE* x = ax.getPointer();
	return asum(n, x, 1);
}

template<class TYPE>
TYPE nrmp( //
		const ArrayListV<TYPE>& ax,  //[in]
		Float p //
		) {
	int n = ax.size();
	const TYPE* x = ax.getPointer();
	return nrmp(n, x, p, 1);
}

template<class TYPE>
TYPE nrminf( //
		const ArrayListV<TYPE>& ax  //[in]
		) {
	int n = ax.size();
	const TYPE* x = ax.getPointer();
	return amax(n, x, 1);
}

}

#endif /* ALGEBRA_BLAS_LEVEL1_H_ */
