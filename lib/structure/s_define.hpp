#ifndef _S_DEFINE_HPP
#define _S_DEFINE_HPP

#include "type_define.hpp"

#include "algebra/algebra.hpp"
#include "geometry/geometry.hpp"
#include <map>
#include <memory>
#include <iomanip>
#include <utility>
#include <unordered_map>
#include <cmath>

namespace structure {

typedef carpio::St St;
typedef int Idx;
typedef double Vt;

#define FOR_EACH_DIM for(St d=0;d<DIM;++d)

enum Axes {
	_X_ = 0, //
	_Y_ = 1, //
	_Z_ = 2, //
};

static const St _NORMAL[3][3] = { { 0, _Z_, _Y_ }, { _Z_, 0, _X_ }, { _Y_, _X_,
		0 } };

inline St Normal(St dim1, St dim2) {
	ASSERT(dim1 != dim2);
	return _NORMAL[dim1][dim2];
}

static const St _NORMAL2[3][3] = { { _Y_, _Z_ }, { _Z_, _X_ }, { _X_, _Y_ } };
inline void Normal(St dim, St& dim1, St& dim2) {
	dim1 = _NORMAL2[dim][0];
	dim2 = _NORMAL2[dim][1];
}

enum Orientation {
	_M_ = 0, //
	_P_ = 1, //
	_C_ = 2, //
};

typedef carpio::MatrixSCR_<Vt> Mat;
typedef carpio::ArrayListV<Vt> Arr;
typedef carpio::ArrayListV<St> Arrst;
typedef carpio::ArrayListV<Idx> Arridx;

template<St Dim>
using Poi_ = carpio::Point_<Vt, Dim>;

// Idx is one component of Index
//
// Idx
//  |
// Index
//  |
// Ijk
//

template<St DIM>
class Index_: public std::array<Idx, DIM> {
public:
	static const St Dim = DIM;
	typedef Index_<DIM> Self;
public:
	Index_() {
		this->fill(0);
	}
	Index_(Idx a, Idx b = 0, Idx c = 0) {
		this->at(0) = a;
		if (Dim >= 2) {
			this->at(1) = b;
		}
		if (Dim >= 3) {
			this->at(2) = c;
		}
	}
	Index_(const Index_<DIM>& other) {
		for (St d = 0; d < Dim; ++d) {
			this->at(d) = other[d];
		}
	}
	Idx i() const {
		return this->at(0);
	}
	Idx j() const {
		return (Dim >= 2) ? this->at(1) : 0;
	}
	Idx k() const {
		return (Dim >= 3) ? this->at(2) : 0;
	}

	// plus 1
	Self p(St dim) {
		Self res(*this);
		if (dim == 0) {
			res[0] += 1;
		}
		if (dim == 1) {
			res[1] += 1;
		}
		if (dim == 2) {
			res[2] += 1;
		}
		return res;
	}
	Self p(St dim) const {
		Self res(*this);
		if (dim == 0) {
			res[0] += 1;
		}
		if (dim == 1) {
			res[1] += 1;
		}
		if (dim == 2) {
			res[2] += 1;
		}
		return res;
	}
	Self m(St dim) {
		Self res(*this);
		if (dim == 0) {
			res[0] -= 1;
		}
		if (dim == 1) {
			res[1] -= 1;
		}
		if (dim == 2) {
			res[2] -= 1;
		}
		return res;
	}
	Self m(St dim) const {
		Self res(*this);
		if (dim == 0) {
			res[0] -= 1;
		}
		if (dim == 1) {
			res[1] -= 1;
		}
		if (dim == 2) {
			res[2] -= 1;
		}
		return res;
	}
	// high level shift
	Self shift(St dim, St ori) {
		// dim 0 1 2
		//     x y z
		// ori 0 1 2
		//     m p c
		ASSERT(ori == 0 || ori == 1);
		if (ori == 0) {
			return this->m(dim);
		} else {
			return this->p(dim);
		}
	}
	Self shift(St dim, St ori) const {
		// dim 0 1 2
		//     x y z
		// ori 0 1 2
		//     m p c
		ASSERT(ori == 0 || ori == 1 || ori == 2);
		if (ori == 0) {
			return this->m(dim);
		} else if (ori == 1) {
			return this->p(dim);
		} else {
			return (*this);
		}
	}

	Idx value(St i) const {
		if (i < Dim) {
			return this->at(i);
		} else {
			return 0;
		}
	}

	void set(St i, Idx val) {
		if (i < Dim) {
			this->at(i) = val;
		}
	}

	bool operator==(const Index_<DIM>& other) const {
		bool res = true;
		for (St d = 0; d < Dim; ++d) {
			res = res && this->at(d) == other[d];
			if (!res) {
				return res;
			}
		}
		return res;
	}

	Idx operator()(const St& dim) const {
		return dim < DIM ? this->at(dim) : 0;
	}

	void show() const {
		std::cout << "(";
		for (St d = 0; d < Dim; ++d) {
			std::cout << this->at(d);
			if (d != Dim - 1) {
				std::cout << ", ";
			}
		}
		std::cout << ")";
	}
};
template<St DIM>
struct Index_compare_ {
	typedef Index_<DIM> Index;
	bool operator()(const Index& a, const Index& b) const {
		for (St d = 0; d < DIM; d++) {
			if (a[d] < b[d]) {
				return true;
			} else if (a[d] == b[d]) {
				continue;
			} else {
				return false;
			}
		}
		return false;
	}
};

template<St DIM>
struct Index_hash_ {
	typedef Index_<DIM> Index;
	std::size_t operator()(const Index& a) const {
		std::size_t res[DIM];
		for (St d = 0; d < DIM; d++) {
			res[d] = std::hash<St> { }(a[d]);
		}
		switch (DIM) {
		case 1: {
			return res[0];
			break;
		}
		case 2: {
			return res[0] ^ (res[1] << 1);
			break;
		}
		case 3: {
			return res[0] ^ ((res[1] << 1) >> 1) ^ (res[2] << 1);
			break;
		}
		}
		return false;
	}
};

template<St DIM>
std::ostream& operator<<(std::ostream& stream, const Index_<DIM>& index) {
	stream << "(";
	for (St d = 0; d < DIM; ++d) {
		stream << index[d];
		if (d != DIM - 1) {
			stream << ", ";
		}
	}
	stream << ")";
	return stream;
}

template<St DIM>
class Ijk_ {
public:
	static const St Dim = DIM;
	typedef Index_<DIM> Index;
protected:
	typedef Ijk_<DIM> Self;
	Index _current;
	Index _end;
public:
	Ijk_(Idx i, Idx in, Idx j = 0, Idx jn = 0, Idx k = 0, Idx kn = 0) {
		Idx a[] = { i, j, k };
		Idx an[] = { in, jn, kn };
		for (St d = 0; d < Dim; ++d) {
			_current[d] = a[d];
			_end[d] = an[d];
		}
	}

	Ijk_(const Index& c, const Index& end) :
			_current(c), _end(end) {
	}

	Ijk_(const Self& other) {
		_current = other._current;
		_end = other._end;
	}

	void reset_end(Index end) {
		this->_end = end;
	}

	Index& current() {
		return this->_current;
	}
	const Index& current() const {
		return this->_current;
	}

	Self& operator=(const Self& other) {
		_current = other._current;
		_end = other._end;
		return *this;
	}

	const Idx& operator()(Idx i) const {
		return _current[i];
	}
	const Idx& operator[](Idx i) const {
		return _current[i];
	}

	bool operator==(const Self& other) {
		return _current == other._current && _end == other._end;
	}

	Self& operator++() {
		return _incr();
	}
	Self& operator--() {
		return _decr();
	}

	Idx i() const {
		return _current[0];
	}
	Idx j() const {
		return (Dim >= 2) ? _current[1] : 0;
	}
	Idx k() const {
		return (Dim >= 3) ? _current[2] : 0;
	}

	void show() const {
		_current.show();
		std::cout << "<";
		_end.show();
		std::cout << ">";
	}

	bool is_end() const {
		St ld = Dim - 1;
		bool res = true;
		for (St d = 0; d < ld; ++d) {
			res = res && (_current[d] == 0);
		}
		return res && (_current[ld] == _end[ld]);
	}

	bool is_start() const {
		bool res = true;
		for (St d = 0; d < Dim; ++d) {
			res = res && (_current[d] == 0);
		}
		return res;
	}

	void to_end() {
		St ld = Dim - 1;
		for (St d = 0; d < ld; ++d) {
			_current[d] == 0;
		}
		_current[ld] = _end[ld];
	}

	void to_last() {
		for (St d = 0; d < Dim; ++d) {
			_current[d] = _end[d] - 1;
		}
	}

protected:
	Self& _incr() {
		if (is_end()) {
			return *this;
		}
		for (St d = 0; d < Dim; ++d) {
			Idx add = _current[d] + 1;
			if (add < _end[d]) {
				_current[d] = add;
				return *this;
			} else if (d != Dim - 1) {
				_current[d] = 0;
			} else {  // force add
				_current[d] = add;
			}
		}
		return *this;
	}

	Self& _decr() {
		if (is_start()) {
			return *this;
		}
		for (St d = 0; d < Dim; ++d) {
			Idx sub = _current[d] - 1;
			if (sub > 0) {
				_current[d] = sub;
				return *this;
			} else if (d != 0) {
				_current[d] = 0;
			} else {
				_current[d] = _end[d] - 1;
			}
		}
		return *this;
	}

};

//#define __DEBUG__

#ifdef __DEBUG__
Index_<2> index2(0, 9);
Index_<3> index3(0, 0, 0);

#endif

}

#endif
