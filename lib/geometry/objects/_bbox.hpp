#ifndef _BBOX_HPP_
#define _BBOX_HPP_

#include "../geometry_define.hpp"
#include "_point.hpp"
#include "../operations/_operation.hpp"
#include <array>

namespace carpio {

template<typename TYPE, St DIM>
class Operation_;

template<typename TYPE, St DIM>
class Bbox_ {
public:
	static const St Dim = DIM;
	typedef Point_<TYPE, DIM> Point;
	typedef Bbox_<TYPE, DIM> Self;
	typedef St size_type;
	typedef TYPE Vt;
	typedef TYPE& reference;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef const TYPE& const_reference;

	typedef Operation_<TYPE, DIM> Op;
private:
	Point _min,_max;
	//TYPE _xmin,_ymin,_xmax,_ymax;
public:
	Bbox_() {
	}
	Bbox_(Vt min, Vt max) {
		for (St i = 0; i < Dim; i++) {
			_min[i] = min;
			_max[i] = max;
		}
	}
	Bbox_(Vt x_min, Vt y_min, Vt x_max, Vt y_max) {
		ASSERT(Dim == 2);
		_min[0] = x_min;
		_min[1] = y_min;

		_max[0] = x_max;
		_max[1] = y_max;
	}

	Bbox_(const Point& point) {
		_min = point;
		_max = point;
	}
	Bbox_(const Point& pmin, const Point& pmax) {
		_min = pmin;
		_max = pmax;
	}
	Vt min(int a) const {
		ASSERT(a < Dim);
		return _min[a];
	}
	Vt max(int a) const {
		ASSERT(a < Dim);
		return _max[a];
	}
	Vt xmin() const {
		return _min[_X_];
	}
	Vt xmax() const {
		return _max[_X_];;
	}
	Vt ymin() const {
		return _min[_Y_];;
	}
	Vt ymax() const {
		return _max[_Y_];;
	}

	Self& operator=(const Self& other) const {
		if (this != &other) { // self-assignment check expected
			_min = other._min;
			_max = other._max;
		}
		return *this;
	}

	Self operator+(const Self& b) const {
		Point min = Op::Min(_min, b._min);
		Point max = Op::Max(_max, b._max);
		return Self(min, max);
	}

};

}

#endif

