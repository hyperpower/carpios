#ifndef _POINT_HPP_
#define _POINT_HPP_

#include "geometry_define.hpp"
#include <array>
#include <sstream>

namespace carpio {
//Point T ====================================
template<typename TYPE, St DIM>
class Point_: public std::array<TYPE, DIM> {
public:
	static const St Dim = DIM;
	//typedef point__tag self_tag;
	typedef St size_type;
	typedef TYPE Vt;
	typedef TYPE& reference;
	typedef TYPE* pointer;
	typedef const TYPE* const_pointer;
	typedef const TYPE& const_reference;

	//constructor
	Point_() :
			std::array<TYPE, Dim>() {
	}

	Point_(const Vt& a, const Vt& b = 0, const Vt& c = 0) :
			std::array<Vt, Dim>() {
		this->at(0) = a;
		if (Dim >= 2) {
			this->at(1) = b;
		}
		if (Dim == 3) {
			this->at(2) = c;
		}
	}

	const_reference val(Axes axi) const {
		switch (axi) {
		case _X_: {
			return this->at(0);
		}
		case _Y_: {
			ASSERT(Dim >= 2);
			return this->at(1);
		}
		case _Z_: {
			ASSERT(Dim >= 3);
			return this->at(2);
		}
		default: {
			SHOULD_NOT_REACH;
		}
		}
		return this->at(0); //make compile happy;
	}
	reference val(Axes axi) {
		switch (axi) {
		case _X_: {
			return this->at(0);
		}
		case _Y_: {
			ASSERT(Dim >= 2);
			return this->at(1);
		}
		case _Z_: {
			ASSERT(Dim >= 3);
			return this->at(2);
		}
		default: {
			SHOULD_NOT_REACH;
		}
		}
		SHOULD_NOT_REACH;
		return this->at(0); //make compile happy;
	}

	Vt value(St idx) const {
		Vt res;
		if (idx < Dim) {
			return val(ToAxes(idx));
		} else {
			return 0;
		}
	}

	const_reference x() const {
		return this->at(0);
	}

	reference x() {
		return this->at(0);
	}

	const_reference y() const {
		return this->at(1);
	}

	reference y() {
		return this->at(1);
	}

	const_reference z() const {
		ASSERT(Dim == 3);
		return this->at(2);
	}

	reference z() {
		ASSERT(Dim == 3);
		return this->at(2);
	}

	void reconstruct(const Vt& a, const Vt& b, const Vt& c = 0) {
		this->at(0) = a;
		this->at(1) = b;
		if (Dim == 3) {
			this->at(2) = c;
		}
	}

	bool operator==(const Point_<Vt, Dim> &a) const {
		if (Dim == 2) {
			return (this->at(0) == a[0] && this->at(1) == a[1]) ? true : false;
		} else {
			return (this->at(0) == a[0] && this->at(1) == a[1]
					&& this->at(2) == a[2]) ? true : false;
		}
	}
	bool operator!=(const Point_<Vt, Dim> &a) const {
		if (Dim == 2) {
			return !((this->at(0) == a[0] && this->at(1) == a[1]) ? true : false);
		} else {
			return !(
					(this->at(0) == a[0] && this->at(1) == a[1]
							&& this->at(2) == a[2]) ? true : false);
		}
	}
	void show() const {
		std::cout << std::scientific << "( " << this->at(0);
		if (Dim >= 2) {
			std::cout << " , " << this->at(1);
		} else if (Dim == 3) {
			std::cout << " , " << this->at(2);
		}
		std::cout << " )\n";
	}
	std::string to_string() const {
		std::stringstream sstr;
		sstr.precision(4);
		sstr << std::scientific << "( " << this->at(0);
		if (Dim >= 2) {
			sstr << " , " << this->at(1);
		} else if (Dim == 3) {
			sstr << " , " << this->at(2) << " )\n";
		} else {
			sstr << " )\n";
		}
		return sstr.str();
	}
	template<typename T>
	void transfer(const T&dx, const T&dy, const T&dz) {
		this->at(0) = this->at(0) + Vt(dx);
		this->at(1) = this->at(1) + Vt(dy);
		if (Dim == 3) {
			this->at(2) = this->at(2) + Vt(dz);
		}
	}
	template<typename T>
	void scale(const T&dx, const T&dy, const T&dz) {
		this->at(0) = this->at(0) * Vt(dx);
		this->at(1) = this->at(1) * Vt(dy);
		if (Dim == 3) {
			this->at(2) = this->at(2) * Vt(dz);
		}
	}
	inline size_type size() const {
		return size_type(Dim);
	}
	inline pointer ptr() {
		return this->data();
	}
	inline const_pointer ptr() const {
		return this->data();
	}
};

template<typename TYPE, St DIM>
std::ostream& operator<<(std::ostream& stream, const Point_<TYPE, DIM>& point) {
	stream << "(";
	for (St d = 0; d < DIM; ++d) {
		stream << point[d];
		if (d != DIM - 1) {
			stream << ", ";
		}
	}
	stream << ")";
	return stream;
}

} //end namespace

#endif /* POINT_H_ */
