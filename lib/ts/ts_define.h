/*
 * ts_define.h
 *
 *  Created on: May 30, 2015
 *      Author: zhou
 */

#ifndef TS_TS_DEFINE_H_
#define TS_TS_DEFINE_H_

#include "stdio.h"
#include <list>
#include <set>
#include <array>
#include <vector>
#include <set>
#include <map>

#include <string>
#include <limits>
#include <memory>
#include "math.h"
#include "../type_define.hpp"


namespace TS {

// data type define  =============================
typedef double Float;
typedef int Int;
typedef unsigned int uInt;
typedef char Char;
typedef std::size_t st;
typedef std::string String;
typedef void* utPointer;

// data structure define =========================
// c++11 support ----alias template--------------
template<typename T>
using List = std::list<T>;
template<typename T, st S>
using Array = std::array<T, S>;
template<typename T>
using Vector = std::vector<T>;
template<typename T>
using Set = std::set<T>;
template<typename KEY, typename VALUE>
using Map = std::map<KEY, VALUE>;
template<typename T1, typename T2>
using Pair = std::pair<T1, T2>;


const static double PI = 3.141592653589793238463;
const double PI_3div4 = 3 * PI / 4;
const double PI_div2 = 1.57079632679489661923;
const double EPSILON = 1e-12;


//
enum Intersect {
	TS_OUT = -1, TS_ON = 0, TS_IN = 1
};

enum ERROR_CODE {
	NO_ERROR = 1, //
	ERR_OTHER = 0, //
	ERR_NULL_POINTER = -1, //
	ERR_DEGERATE = -2, //
};

enum OBJ_TYPE {
	EMPTY = 0,
	POINT = 1,
	VERTEX = 2,
	SEGMENT = 3,
	EDGE = 4,
	TRIANGLE = 5,
	FACE = 6,
	SURFACE = 7
};
enum Aix {
	_X = 0, _Y = 1, _Z = 2
};

enum Plane {
	_XY = 0, _YZ = 1, _ZX = 2
};

inline Aix ToAix(const st& i) {
	ASSERT(i >= 0 && i < 3);
	switch (i) {
	case 0:
		return _X;
	case 1:
		return _Y;
	case 2:
		return _Z;
	default:
		ASSERT_MSG(false, "Error input i");
	}
	SHOULD_NOT_REACH;
	return _X;
}
enum Location {
	_M = 0, _C = 1, _P = 2,
};
#define _return_val_if_fail(expr,val)  {       \
		if (!(expr))                           \
			return (val);                      \
	    }

#define _return_if_fail(expr)  {               \
		if (!(expr))                           \
			return ;                           \
	    }

template<class V>
inline int SIGN(const V& x) {
	return ((x) > 0. ? 1 : -1);
}
template<class V>
inline int ORIENT1D(const V& a, const V& b) {
	return ((a) > (b) ? 1 : (a) < (b) ? -1 : 0);
}
template<class V>
inline V MAX(const V& a, const V& b) {
	return a > b ? a : b;
}

template<class V>
inline V MAX(const V& a, const V& b, const V& c) {
	V tmp = (a > b ? a : b);
	return tmp > c ? tmp : c;
}

template<class V>
inline V MIN(const V& a, const V& b) {
	return a < b ? a : b;
}

template<class V>
inline V MIN(const V& a, const V& b, const V& c) {
	V tmp = (a < b ? a : b);
	return tmp < c ? tmp : c;
}
template<class V>
inline V ABS(const V& a) {
	return a < 0 ? -a : a;
}

template<class V>
static int sortp(V p[], std::size_t n) {
	int sign = 1;
	std::size_t i, j;

	for (i = 0; i < n - 1; ++i) {
		for (j = 0; j < n - 1 - i; ++j) {
			if (long(p[j + 1]) < long(p[j])) {
				V tmp = p[j];

				p[j] = p[j + 1];
				p[j + 1] = tmp;
				sign = -sign;
			}
		}
	}
	return sign;
}

class Statistic {

protected:
	void copy(const Statistic& src, Statistic& dst) const {
		dst.min = src.min;
		dst.max = src.max;
		dst.sum = src.sum;
		dst.sum2 = src.sum2;
		dst.mean = src.mean;
		dst.stddev = src.stddev;
		dst.n = src.n;
	}
public:
	Float min, max, sum, sum2, mean, stddev;
	uInt n;

	Statistic() {
		reset();
	}

	Statistic(const Statistic& src) {
		min = src.min;
		max = src.max;
		sum = src.sum;
		sum2 = src.sum2;
		mean = src.mean;
		stddev = src.stddev;
		n = src.n;
	}

	Statistic& operator=(const Statistic &a) {
		if (this == &a) {
			return *this;
		}
		copy(a, (*this));
		return *this;
	}

	void reset() {
		max = std::numeric_limits<Float>::max();
		min = std::numeric_limits<Float>::min();
		sum = 0.0;
		sum2 = 0.0;
		mean = 0.0;
		stddev = 0.0;
		n = 0;
	}

	void add_value(Float val) {
		if (n == 0) {
			max = val;
			min = val;
			sum = val;
			sum2 = val * val;
			mean = val;
			n = 1;
		} else {
			if (val < min)
				min = val;
			if (val > max)
				max = val;
			sum += val;
			sum2 += val * val;
			n++;
			mean = sum / n;
			stddev = sqrt((sum2 - sum * sum / (Float) n) / (Float) n);
		}
	}

	void show() const {
		std::cout << "min:    " << min << " -> max: " << max << "\n";
		std::cout << "mean:   " << mean << "\n";
		std::cout << "sum:    " << sum << "\n";
		std::cout << "stddev: " << stddev << "\n";
	}
}
;

}

#endif /* TS_TS_DEFINE_H_ */
