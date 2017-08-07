#ifndef _BBOX_HPP_
#define _BBOX_HPP_

#include "../geometry_define.hpp"
#include "_point.hpp"
#include "../operations/_operation.hpp"
#include <array>
#include "../../utility/any.hpp"

#include "_triface.hpp"

namespace carpio {

template<typename TYPE, St DIM>
class Operation_;

struct TagBox: public TagGeometry {
	TagBox() {
	}
};

template<typename TYPE, St DIM>
class Box_ {
public:
	static const St Dim = DIM;
	typedef Point_<TYPE, DIM> Point;
	typedef Box_<TYPE, DIM> Self;
	typedef St size_type;
	typedef TYPE Vt;
	typedef TagBox Tag;
	typedef TYPE& reference;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef const TYPE& const_reference;

	typedef Operation_<TYPE, DIM> Op;
	protected:
	Point _min, _max;
	//TYPE _xmin,_ymin,_xmax,_ymax;
public:
	Box_() {
	}
	Box_(Vt min, Vt max) {
		for (St i = 0; i < Dim; i++) {
			_min[i] = min;
			_max[i] = max;
		}
	}
	Box_(Vt x_min, Vt y_min, Vt x_max, Vt y_max) {
		ASSERT(Dim == 2);
		_min[0] = x_min;
		_min[1] = y_min;

		_max[0] = x_max;
		_max[1] = y_max;
	}

	Box_(const Point& point) {
		_min = point;
		_max = point;
	}
	Box_(const Point& pmin, const Point& pmax) {
		_min = pmin;
		_max = pmax;
	}

	Box_(const Self& a, const Self& b) { // union to one box without data
		Self& self = (*this);
		self._min = Op::Min(a._min, b._min);
		self._max = Op::Max(a._max, b._max);
	}

	Vt min(int a) const {
		ASSERT(a < Dim);
		return _min[a];
	}
	Vt max(int a) const {
		ASSERT(a < Dim);
		return _max[a];
	}
	Point& min() {
		return _min;
	}
	const Point& min() const {
		return _min;
	}
	Point& max() {
		return _max;
	}
	const Point& max() const {
		return _max;
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

	Vt center(int a) const {
		ASSERT(a < Dim);
		return (_min[a] + _max[a]) * 0.5;
	}

	Point center() const {
		Point res;
		for (St i = 0; i < Dim; i++) {
			res[i] = center(i);
		}
		return res;
	}

	Self& operator=(const Self& other) {
		if (this != &other) { // self-assignment check expected
			this->_min = other._min;
			this->_max = other._max;
		}
		return *this;
	}

	bool operator<(const Self& rhs) const {
		//compare center point x -> y -> z -> point address
		const Self& self = (*this);
		if (this->center(_X_) < rhs.center(_X_)) {
			return true;
		} else if (this->center(_X_) == rhs.center(_X_)) {
			if (DIM < 2)
				return (this < &rhs);
			if (this->center(_Y_) < rhs.center(_Y_)) {
				return true;
			} else if (this->center(_Y_) == rhs.center(_Y_)) {
				if (DIM < 3)
					return (this < &rhs);
				if (this->center(_Z_) < rhs.center(_Z_)) {
					return true;
				}
			}
		}
		return (this < &rhs);
	}

	Self operator+(const Self& b) const {
		Point min = Op::Min(_min, b._min);
		Point max = Op::Max(_max, b._max);
		return Self(min, max);
	}

	Vt get(Axes aix, Orientation loc) const {
		if (aix == _Z_ && DIM == 2) {
			return 0;
		}
		Vt xv = 0;
		const Self& self = (*this);
		switch (loc) {
		case _M_: {
			xv = _min[aix];
			break;
		}
		case _C_: {
			xv = (_min[aix] + _max[aix]) * 0.5;
			break;
		}
		case _P_: {
			xv = _max[aix];
			break;
		}
		default: {
			break;
		}
		}
		return xv;
	}

	Vt get(Orientation loc, Axes aix) const {
		return get(aix, loc);
	}

	Vt get_d(Axes aix) const {
		const Self& self = (*this);
		return _max[aix] - _min[aix];
	}

	Vt get_dh(Axes aix) const {
		const Self& self = (*this);
		return (_max[aix] - _min[aix]) / 2.0;
	}

	Point get_point(Orientation lx, Orientation ly, Orientation lz) const {
		return Poi(get(_X_, lx), get(_Y_, ly), get(_Z_, lz));
	}

};

struct TagBBox: public TagGeometry {
	TagBBox() {
	}
};

template<class TYPE, St DIM>
class TriSurface_;

template<class TYPE, St DIM, class SURFACE>
class TriFace_;

template<typename TYPE, St DIM>
class BBox_: public Box_<TYPE, DIM> {
public:
	static const St Dim = DIM;
	typedef Point_<TYPE, DIM> Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef BBox_<TYPE, DIM> Self;
	typedef Box_<TYPE, DIM> Box;
	typedef St size_type;
	typedef TagBBox Tag;
	typedef TYPE Vt;
	typedef TYPE& reference;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef const TYPE& const_reference;

	typedef Operation_<TYPE, DIM> Op;

	typedef TriSurface_<TYPE, DIM> TriSurface;
	typedef TriFace_<TYPE, DIM, TriSurface> TriFace;

protected:
	Any _obj;

public:
	BBox_() {
		_obj = nullptr;
	}
	BBox_(const Any& o) {
		_obj = o;
		_set_box();
	}

	BBox_(const Point& min, const Point& max): Box(min, max){
		_obj = nullptr;
	}


	BBox_(const Self& a, const Self& b) :
			Box(a, b) { // union to one box without data
		_obj = nullptr;
	}

	Any& get_obj() {
		return _obj;
	}

	const Any& get_obj() const {
		return _obj;
	}

protected:
	void _set_box() {
		if (_obj.type() == typeid(Segment)) {
			Segment s = any_cast<Segment>(_obj);
			Box bs = s.box();
			this->_min = bs.min();
			this->_max = bs.max();
			return;
		}
		if ((_obj.type() == typeid(Segment*))) {
			Segment* s = any_cast<Segment*>(_obj);
			Box bs = s->box();
			this->_min = bs.min();
			this->_max = bs.max();
			return;
		}
		if ((_obj.type() == typeid(TriFace))) {
			TriFace s = any_cast<TriFace>(_obj);
			Box bs = s.box();
			this->_min = bs.min();
			this->_max = bs.max();
			return;
		}
		if ((_obj.type() == typeid(TriFace*))) {
			TriFace* s = any_cast<TriFace*>(_obj);
			Box bs = s->box();
			this->_min = bs.min();
			this->_max = bs.max();
			return;
		}

	}

}
;

}

#endif

