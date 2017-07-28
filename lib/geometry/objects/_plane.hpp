#ifndef _PLANE_HPP_
#define _PLANE_HPP_

#include "../geometry_define.hpp"
#include "_point.hpp"
#include <array>

namespace carpio {
template<typename TYPE>
class Plane_: public std::array<TYPE, 4> {
	//The Line function defined as ax+by+cz=alpha
public:
	typedef St size_type;
	typedef TYPE vt;
	typedef TYPE& reference;
	typedef const TYPE& const_reference;
public:
	Plane_() :
			std::array<vt, 4>() {
	}
	Plane_(const vt& a, const vt& b, const vt& c, const vt& alpha) {
		this->at(0) = a;
		this->at(1) = b;
		this->at(2) = c;
		this->at(3) = alpha;
	}
	Plane_(const vt& x, const vt& y, const vt& z, const vt& nx, const vt& ny,
			const vt& nz) {
		//assert(!isEqual(ax, bx) || !isEqual(ay,by));
		this->a() = nx;
		this->b() = ny;
		this->c() = nz;
		this->alpha() = (nx * x + ny * y + nz * z);
	}
	void reconstruct(const vt& a, const vt& b, const vt& c, const vt& alpha) {
		if (a == 0.0 && b == 0.0) {
			a = SMALL;
			b = SMALL;
		}
		this->a() = a;
		this->b() = b;
		this->c() = c;
		this->alpha() = alpha;
	}
	inline reference a() {
		return this->at(0);
	}
	inline reference b() {
		return this->at(1);
	}
	inline reference c() {
		return this->at(2);
	}
	inline reference alpha() {
		return this->at(3);
	}
	inline const_reference a() const {
		return this->at(0);
	}
	inline const_reference b() const {
		return this->at(1);
	}
	inline const_reference c() const {
		return this->at(2);
	}
	inline const_reference alpha() const {
		return this->at(2);
	}
	vt cal_x(const vt& y, const vt& z) const {
		if (this->a() == 0.0) {
			return (this->alpha() - this->b() * y - this->c() * z) / SMALL;
		} else {
			return (this->alpha() - this->b() * y - this->c() * z) / this->a();
		}
	}
	vt cal_y(const vt& x, const vt& z) const {
		if (this->b() == 0.0) {
			return (this->alpha() - this->a() * x - this->c() * z) / SMALL;
		} else {
			return (this->alpha() - this->a() * x - this->c() * z) / this->b();
		}
	}
	vt cal_z(const vt& x, const vt& y) const {
		if (this->c() == 0.0) {
			return (this->alpha() - this->a() * x - this->b() * y) / SMALL;
		} else {
			return (this->alpha() - this->a() * x - this->b() * y) / this->c();
		}
	}
	vt intersept_x() const {
		return this->cal_x(0, 0);
	}
	vt intersept_y() const {
		return this->cal_y(0, 0);
	}
	vt intersept_z() const {
		return this->cal_z(0, 0);
	}
	vt norm_x() const {
		return this->a();
	}
	vt norm_y() const {
		return this->b();
	}
	vt norm_z() const {
		return this->c();
	}
	bool empty() const {
		if (this->a() != 0.0 && this->b() != 0.0 && this->c() != 0.0) {
			return true;
		} else {
			return false;
		}
	}
	void show() const {
		std::cout << this->a() << " X + " << this->b() << " Y + " << this->c()
				<< " Z = " << this->alpha() << "\n";
	}
};

}

#endif

