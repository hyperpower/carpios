#ifndef CELL_H_
#define CELL_H_

#include "domain_define.hpp"
#include "geometry/geometry.hpp"

namespace carpio {

template<typename VALUE, St DIM>
class Cell_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);

	typedef VALUE vt;

	typedef Cell_<VALUE, DIM> Self;

	typedef void (*pfunction)(Self *, utPointer);

protected:
	vt _center[Dim];
	vt _hd[Dim];
public:
	/*
	 *  constructor
	 */
	Cell_() {
		for (St i = 0; i < Dim; ++i) {
			_center[i] = 0.0;
			_hd[i] = 0.0;
		}
	}

	Cell_(const vt &x, const vt &dhx, //
			const vt &y = 0.0, const vt &dhy = 0.0, //
			const vt &z = 0.0, const vt &dhz = 0.0) {
		for (St i = 0; i < Dim; ++i) {
			if (i == 0) {
				_center[i] = x;
				ASSERT(dhx > 0.0);
				_hd[i] = dhx;
			} else if (i == 1) {
				_center[i] = y;
				ASSERT(dhy > 0.0);
				_hd[i] = dhy;
			} else if (i == 2) {
				_center[i] = z;
				ASSERT(dhz > 0.0);
				_hd[i] = dhz;
			}
		}
	}
	Cell_(const Self& s) {
		for (int i = 0; i < Dim; i++) {
			_center[i] = s._center[i];
			_hd[i] = s._hd[i];
		}
	}

	/*
	 *  get
	 */
	inline vt get(const Orientation &ori, const Axes &axes) const {
		vt res = 0.0;
		if (axes >= Dim) {
			return 0.0;
		}
		switch (ori) {
		case _M_: {
			res = _center[axes] - _hd[axes];
			break;
		}
		case _C_: {
			res = _center[axes];
			break;
		}
		case _P_: {
			res = _center[axes] + _hd[axes];
			break;
		}
		default: {
			break;
		}
		}
		return res;
	}

	inline vt get_d(const Axes &axes) const {
		ASSERT(axes < Dim);
		return 2.0 * _hd[axes];
	}

	inline vt get_hd(const Axes &axes) const {
		return _hd[axes];
	}

	inline vt volume() const {
		vt res = 1.0;
		for (St i = 0; i < Dim; ++i) {
			res *= 2.0 * _hd[i];
		}
		return res;
	}

protected:
	static const int INSIDE = 0; // 000000
	static const int LEFT = 1;   // 000001
	static const int RIGHT = 2;  // 000010
	static const int BOTTOM = 4; // 000100
	static const int TOP = 8;    // 001000
	static const int BACK = 16;  // 010000
	static const int FRONT = 32; // 100000

	int _out_code(vt x, vt y, vt z = 0) {
		vt xmin = this->get(_M_, _X_);
		vt ymin = this->get(_M_, _Y_);
		vt xmax = this->get(_P_, _X_);
		vt ymax = this->get(_P_, _Y_);

		vt zmin = this->get(_M_, _Z_);
		vt zmax = this->get(_P_, _Z_);

		int code;

		code = INSIDE;          // initialised as being inside of clip cell

		if (x < xmin)           // to the left of clip cell
			code |= LEFT;
		else if (x > xmax)      // to the right of clip cell
			code |= RIGHT;

		if (y < ymin)           // below the clip cell
			code |= BOTTOM;
		else if (y > ymax)      // above the clip cell
			code |= TOP;

		if (Dim >= 3) {
			if (z < zmin)           // below the clip cell
				code |= BACK;
			else if (z > zmax)      // above the clip cell
				code |= FRONT;
		}

		return code;
	}

public:
	bool is_in_on(const vt& x, const vt& y = 0, const vt& z = 0) const {
		return (IsInRange(this->get(_M_, _X_), x, this->get(_P_, _X_), _cc_)
				&& ((Dim >= 2) ?
						IsInRange(this->get(_M_, _Y_), y, this->get(_P_, _Y_),
								_cc_) :
						true)
				&& ((Dim == 3) ?
						IsInRange(this->get(_M_, _Z_), z, this->get(_P_, _Z_),
								_cc_) :
						true));
	}
	bool is_in_on(const Point_<vt, Dim>& p) const{
		if (Dim >= 2) {
			return this->is_in_on(p.x(), p.y());
		} else if (Dim >= 3) {
			return this->is_in_on(p.x(), p.y(), p.z());
		} else {
			return false;
		}
	}
	bool is_in_on(const Axes& axes, const vt& cor, Range r = _cc_) const {
		ASSERT(axes < Dim);
		vt m = this->get(_M_, axes);
		vt p = this->get(_P_, axes);
		return IsInRange(m, cor, p, r);
	}
	bool is_in_on(const Segment_<vt, Dim>& seg) {
		if (Dim >= 2) {
			return is_in_on(seg.psx(), seg.pex(), seg.psy(), seg.pey());
		} else if (Dim >= 3) {
			return is_in_on(seg.psx(), seg.pex(), seg.psy(), seg.pey(),
					seg.psz(), seg.pez());
		} else {
			return false;
		}
	}

	bool is_in_on( //
			const vt& xs, const vt& xe, //
			const vt& ys, const vt& ye, //
			const vt& zs = 0, const vt& ze = 0 //
			) const {
		int ocs = _out_code(xs, ys, zs);
		int oce = _out_code(xe, ye, ze);

		if (!(ocs | oce)) {
			// Bitwise OR is 0. If one point inon cell, the segment is inon
			return true;
		}
		if ((ocs & oce) > 0) {
			// Bitwise AND large than 0. the two points are in the same region.
			return false;
		}
		// The poits are out of cell, but the segment may across the cell.
		// #include "geometry.hpp"
		if (Dim == 2) {
			Segment_<vt, Dim> seg(xs, xe, ys, ye);
			int res = 0;
			Point_<vt, Dim> p(this->get(_M_, _X_), this->get(_M_, _Y_));
			res += OnWhichSide3(seg, p);
			p.reconstruct(this->get(_P_, _X_), this->get(_M_, _Y_));
			res += OnWhichSide3(seg, p);
			p.reconstruct(this->get(_P_, _X_), this->get(_P_, _Y_));
			res += OnWhichSide3(seg, p);
			p.reconstruct(this->get(_M_, _X_), this->get(_P_, _Y_));
			res += OnWhichSide3(seg, p);
			if (res == 4 || res == -4) {
				return false;
			} else {
				return true;
			}
		} else if (Dim == 3) {
			ASSERT_MSG(false, "unfinish 3d");
		} else {
			SHOULD_NOT_REACH;
		}
		return false; // make complier happy;
	}

	void show(pfunction fun = nullptr, utPointer utp = nullptr) const {
		fun(this, utp);
	}

	inline void transfer(vt dx, vt dy = 0.0, vt dz = 0.0) {
		_center[0] += dx;
		if (Dim >= 2) {
			_center[1] += dy;
		}
		if (Dim >= 3) {
			_center[2] += dz;
		}
	}
}
;

}

#endif /* CELL_H_ */
