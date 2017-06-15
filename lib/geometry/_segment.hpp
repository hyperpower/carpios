#ifndef _SEGMENT_HPP_
#define _SEGMENT_HPP_

#include "../type_define.hpp"
#include "geometry_define.hpp"
#include "_point.hpp"
#include <array>
#include "math.h"
namespace carpio {

template<typename TYPE, St DIM>
class Segment_: public std::array<Point_<TYPE, DIM>, 2> {
public:
	static const St Dim = DIM;
	typedef TYPE Vt;
	typedef Vt& ref_Vt;
	typedef const Vt& const_ref_Vt;
	typedef Segment_<Vt, DIM> Self;
	typedef Segment_<Vt, DIM>& ref_Self;
	typedef const Segment_<Vt, DIM>& const_ref_Self;
	typedef Point_<Vt, DIM> Point;
	typedef Point_<Vt, DIM>* pPoint;
	typedef Point_<Vt, DIM>& ref_Point;
	typedef const Point_<Vt, DIM>& const_ref_Point;
public:
	Segment_() :
			std::array<Point_<Vt, DIM>, 2>() {
		_set_empty();
	}
	Segment_(const Point& s, const Point& e) {
		ASSERT(s != e);
		this->ps() = s;
		this->pe() = e;
	}
	Segment_(const Vt& ax, const Vt& bx, //
			const Vt& ay, const Vt& by,  //
			const Vt& az = 0, const Vt& bz = 0) {
		Point s(ax, ay, az);
		Point e(bx, by, bz);
		reconstruct(s, e);
	}
	Segment_(const_ref_Self rhs) {
		this->ps() = rhs.ps();
		this->pe() = rhs.pe();
	}
	void reconstruct(const_ref_Self rhs) {
		this->ps() = rhs.ps();
		this->pe() = rhs.pe();
	}
	void reconstruct(const Point& s, const Point& e) {
		ASSERT(s != e);
		this->ps() = s;
		this->pe() = e;
	}
	void reconstruct(const Vt& ax, const Vt& bx, //
			const Vt& ay, const Vt& by,  //
			const Vt& az = 0, const Vt& bz = 0) {
		Point s(ax, ay, az);
		Point e(bx, by, bz);
		reconstruct(s, e);
	}

	bool operator==(const_ref_Self rhs) const {
		return (this->ps() == rhs.ps() && this->pe() == rhs.pe()) ? true : false;
	}

	ref_Point ps() {
		return this->at(0);
	}
	const_ref_Point ps() const {
		return this->at(0);
	}
	ref_Point pe() {
		return this->at(1);
	}
	const_ref_Point pe() const {
		return this->at(1);
	}
	Point pc() const {
		return Point((pex() + psx()) * 0.5, (pey() + psy()) * 0.5,
				(Dim == 3) ? ((pez() + psz()) * 0.5) : 0);
	}

	const_ref_Vt psx() const {
		return this->ps().x();
	}
	const_ref_Vt pex() const {
		return this->pe().x();
	}
	const_ref_Vt psy() const {
		ASSERT(Dim >= 2);
		return this->ps().y();
	}
	const_ref_Vt pey() const {
		ASSERT(Dim >= 2);
		return this->pe().y();
	}
	const_ref_Vt psz() const {
		ASSERT(Dim >= 3);
		return this->ps().z();
	}
	const_ref_Vt pez() const {
		ASSERT(Dim >= 3);
		return this->pe().z();
	}
	ref_Vt psx() {
		return this->ps().x();
	}
	ref_Vt pex() {
		return this->pe().x();
	}
	ref_Vt psy() {
		ASSERT(Dim >= 2);
		return this->ps().y();
	}
	ref_Vt pey() {
		ASSERT(Dim >= 2);
		return this->pe().y();
	}
	ref_Vt psz() {
		ASSERT(Dim >= 3);
		return this->ps().z();
	}
	ref_Vt pez() {
		ASSERT(Dim >= 3);
		return this->pe().z();
	}
	Vt dx() const {
		return pex() - psx();
	}
	Vt dy() const {
		return pey() - psy();
	}
	Vt dz() const {
		return pez() - psz();
	}
	Vt length() const {
		Vt len = 0.0;
		len = sqrt(
				double(
						(psx() - pex()) * (psx() - pex())
								+ (psy() - pey()) * (psy() - pey())));
		return len;
	}
	Vt slope() const {
		ASSERT(Dim == 2);
		return (pey() - psy()) / (pex() - psx() + SMALL);
	}

	void scale(Vt xfactor, Vt yfactor, Vt zfactor = 1) {
		this->ps().x() = psx() * xfactor;
		this->ps().y() = psy() * yfactor;
		if (Dim == 3) {
			psz() = psz() * zfactor;
		}
		pex() = pex() * xfactor;
		pey() = pey() * yfactor;
		if (Dim == 3) {
			pez() = pez() * zfactor;
		}
		if (pe() == ps()) {
			_set_empty();
		}
	}
	void transfer(Vt dx, Vt dy, Vt dz = 0.0) {
		if (!empty()) {
			psx() = psx() + dx;
			psy() = psy() + dy;
			pex() = pex() + dx;
			pey() = pey() + dy;
		}
		if (Dim == 3) {
			psz() = psz() + dz;
			pez() = pez() + dz;
		}
	}

	/** Set the beginning point */
	void set_s(const Point& p) {
		this->ps() = p;
	}
	/** Set the end point */
	void set_e(const Point& p) {
		this->pe() = p;
	}

	/** Change the segment orientation */
	Self& change_orientation() {
		Point tmp = this->ps();
		this->ps() = this->pe();
		this->pe() = tmp;
		return *this;
	}

	bool empty() const {
		if (psx() == 0.0 && psy() == 0.0 && pex() == 0.0 && pey() == 0.0
				&& ((Dim == 3) ? (psz() == 0.0 && pez() == 0.0) : true)) {
			return true;
		} else {
			return false;
		}
	}
	void show() const {
		std::cout.precision(4);
		std::cout << "( " << psx() << ", " << psy();
		if (Dim == 3) {
			std::cout << ", " << psz();
		} else {
			std::cout << "";
		}
		std::cout << " )--->(" << this->pex() << ", " << pey();
		if (Dim == 3) {
			std::cout << ", " << pez();
		} else {
			std::cout << "";
		}
		std::cout << " )\n";
	}

	/*
	 *  compare
	 */

	bool is_gt_x(const Vt& v) const {    //>
		return (this->pex() > v && this->psx() > v);
	}
	bool is_gt_y(const Vt& v) const {    //>
		return (this->pey() > v && this->psy() > v);
	}
	bool is_ge_x(const Vt& v) const {    //>=
		return (this->pex() >= v && this->psx() >= v);
	}
	bool is_ge_y(const Vt& v) const {    //>=
		return (this->pey() >= v && this->psy() >= v);
	}

	bool is_lt_x(const Vt& v) const {    //<
		return (this->pex() < v && this->psx() < v);
	}
	bool is_lt_y(const Vt& v) const {    //<
		return (this->pey() < v && this->psy() < v);
	}
	bool is_lt_z(const Vt& v) const {    //<
		ASSERT(Dim == 3);
		return (this->pez() < v && this->psz() < v);
	}
	bool is_le_x(const Vt& v) const {    //<=
		return (this->pex() <= v && this->psx() <= v);
	}
	bool is_le_y(const Vt& v) const {    //<=
		return (this->pey() <= v && this->psy() <= v);
	}

	bool is_vertical() const {
		ASSERT(!empty());
		return psx() == pex();
	}

	bool is_horizontal() const {
		ASSERT(!empty());
		return psy() == pey();
	}
	bool is_in_box(const Point& pt) const {
		ASSERT(!empty());
		if (is_horizontal()) {
			return (((psx() <= pt.x) && (pt.x <= pex()))
					|| ((pex() <= pt.x) && (pt.x <= psx())));
		}
		if (is_vertical()) {
			return (((psy() <= pt.y) && (pt.y <= pey()))
					|| ((pey() <= pt.y) && (pt.y <= psy())));
		}
		return (((psx() <= pt.x) && (pt.x <= pex()))
				|| ((pex() <= pt.x) && (pt.x <= psx())))
				&& (((psy() <= pt.y) && (pt.y <= pey()))
						|| ((pey() <= pt.y) && (pt.y <= psy())));
	}

protected:
	void _set_empty() {
		ps().x() = 0.0;
		pe().x() = 0.0;
		ps().y() = 0.0;
		pe().y() = 0.0;
		if (Dim == 3) {
			ps().z() = 0.0;
			pe().z() = 0.0;
		}
	}
};

template<typename TYPE, St DIM>
inline std::ostream& operator<<(std::ostream& o, const Segment_<TYPE, DIM>& p) {
	return o << p.ps() << "-" << p.pe();
}

}
#endif
