/************************
 //  \file   Point.h
 //  \brief
 // 
 //  \author czhou
 //  \date   21 mai 2015 
 ***********************/
#ifndef TS_POINT_H_
#define TS_POINT_H_

#include <array>
#include <assert.h>
#include "ts_define.h"
//#include "ts_predicates.h"

namespace TS {

template<class TYPE, st DIM> class Segment;
template<class TYPE, st DIM> class Surface;
template<class TYPE, st DIM> class Edge;
template<class TYPE, st DIM> class Vertex;
template<class TYPE, st DIM> class Triangle;

template<class TYPE, st DIM>
class Point: public Array<TYPE, DIM> {
public:
	typedef Array<TYPE, DIM> base_class;
	typedef Point<TYPE, DIM> self_class;

	typedef typename base_class::value_type value_type;
	typedef typename base_class::pointer pointer;
	typedef typename base_class::const_pointer const_pointer;
	typedef typename base_class::reference reference;
	typedef typename base_class::const_reference const_reference;
	typedef typename base_class::iterator iterator;
	typedef typename base_class::const_iterator const_iterator;
	typedef typename base_class::size_type size_type;
	typedef typename base_class::difference_type difference_type;
	typedef typename base_class::reverse_iterator reverse_iterator;
	typedef typename base_class::const_reverse_iterator const_reverse_iterator;

	//constructor
	Point() :
			base_class() {
	}

	Point(const TYPE& a, const TYPE& b, const TYPE& c = 0) :
			base_class() {
		this->at(0) = a;
		this->at(1) = b;
		if (DIM == 3) {
			this->at(2) = c;
		}
	}

	void set(const TYPE& a, const TYPE& b, const TYPE& c = 0) {
		this->at(0) = a;
		this->at(1) = b;
		if (DIM == 3) {
			this->at(2) = c;
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
		assert(DIM == 3);
		return this->at(2);
	}

	reference z() {
		assert(DIM == 3);
		return this->at(2);
	}

	const_reference val(Aix a) const {
		switch (a) {
		case _X: {
			return this->x();
			break;
		}
		case _Y: {
			return this->y();
			break;
		}
		case _Z: {
			assert(DIM == 3);
			return this->z();
			break;
		}
		}
		SHOULD_NOT_REACH;
		return this->x();
	}

	reference val(Aix a) {
		switch (a) {
		case _X: {
			return this->x();
			break;
		}
		case _Y: {
			return this->y();
			break;
		}
		case _Z: {
			assert(DIM == 3);
			return this->z();
			break;
		}
		}
		SHOULD_NOT_REACH;
		return this->x();
	}

	void reconstruct(const TYPE& a, const TYPE& b, const TYPE& c = 0) {
		this->at(0) = a;
		this->at(1) = b;
		if (DIM == 3) {
			this->at(2) = c;
		}
	}

	bool operator==(const Point<TYPE, DIM> &a) const {
		if (DIM == 2) {
			return (this->at(0) == a[0] && this->at(1) == a[1]) ? true : false;
		} else {
			return (this->at(0) == a[0] && this->at(1) == a[1]
					&& this->at(2) == a[2]) ? true : false;
		}
	}
	bool operator!=(const Point<TYPE, DIM> &a) const {
		if (DIM == 2) {
			return !((this->at(0) == a[0] && this->at(1) == a[1]) ? true : false);
		} else {
			return !(
					(this->at(0) == a[0] && this->at(1) == a[1]
							&& this->at(2) == a[2]) ? true : false);
		}
	}
	void show() const {
		std::cout << std::scientific << "( " << this->at(0) << " , "
				<< this->at(1);
		if (DIM == 3) {
			std::cout << " , " << this->at(2) << " )\n";
		} else {
			std::cout << " )\n";
		}
	}

	template<typename T>
	void transfer(const T&dx, const T&dy, const T&dz) {
		this->at(0) = this->at(0) + TYPE(dx);
		this->at(1) = this->at(1) + TYPE(dy);
		if (DIM == 3) {
			this->at(2) = this->at(2) + TYPE(dz);
		}
	}

	template<typename T>
	void scale(const T&dx, const T&dy, const T&dz) {
		this->at(0) = this->at(0) * TYPE(dx);
		this->at(1) = this->at(1) * TYPE(dy);
		if (DIM == 3) {
			this->at(2) = this->at(2) * TYPE(dz);
		}
	}

	inline size_type size() const {
		return size_type(DIM);
	}

};

template<typename TYPE, st DIM>
std::ostream& operator<<(std::ostream& stream, const Point<TYPE, DIM>& point) {
	stream << "(";
	for (st d = 0; d < DIM; ++d) {
		stream << point[d];
		if (d != DIM - 1) {
			stream << ", ";
		}
	}
	stream << ")";
	return stream;
}


template<class POINT>
bool point_is_in_on_rectangle(const POINT& p, const POINT& p1,
		const POINT& p2) {
	bool res = true;
	for (typename POINT::size_type i = 0; i < POINT::Dim; ++i) {
		res = res && p[i] >= p1[i] && p[i] <= p2[i];
	}
	return res;
}
/**
 * gts_point_orientation:
 * @p1: a #GtsPoint.
 * @p2: a #GtsPoint.
 * @p3: a #GtsPoint.
 *
 * Checks for orientation of the projection of three points on the
 * (x,y) plane. The result is also an approximation of twice the
 * signed area of the triangle defined by the three points. This
 * function uses adaptive floating point arithmetic and is
 * consequently geometrically robust.
 *
 * Returns: a positive value if @p1, @p2 and @p3 appear in
 * counterclockwise order, a negative value if they appear in
 * clockwise order and zero if they are colinear.
 */
//template<class POINT>
//Float point_orientation(const POINT& p1, const POINT& p2, const POINT& p3) {
//	return orient2d((double *) p1.data(), (double *) p2.data(),
//			(double *) p3.data());
//}
/**
 * point_orientation_3d:
 * @p1: a Point.
 * @p2: a Point.
 * @p3: a Point.
 * @p4: a Point.
 *
 * Checks if @p4 lies above, below or on the plane passing through the
 * points @p1, @p2 and @p3. Below is defined so that @p1, @p2 and @p3
 * appear in counterclockwise order when viewed from above the
 * plane. The returned value is an approximation of six times the
 * signed volume of the tetrahedron defined by the four points. This
 * function uses adaptive floating point arithmetic and is
 * consequently geometrically robust.
 *
 * Returns:
 *        a positive value if @p4 lies below,
 *        a negative value if @p4 lies above the plane,
 *        zero             if the four points are coplanar.
 */
//template<class POINT>
//Float point_orientation_3d(const POINT& p1, const POINT& p2, const POINT& p3,
//		const POINT& p4) {
//	assert(POINT::Dim == 3);
//	return orient3d((double *) p1.data(), (double *) p2.data(),
//			(double *) p3.data(), (double *) p4.data());
//}
/**
 * point_is_in_triangle:
 * @p: a #GtsPoint.
 * @t: a #GtsTriangle.
 *
 * Tests if the planar projection (x, y) of @p is inside, outside or
 * on the boundary of the planar projection of @t.  This function is
 * geometrically robust.
 *
 * Returns: %GTS_IN if @p is inside @t, %GTS_ON if @p is on the boundary of
 * @t, %GTS_OUT otherwise.
 */
template<class TYPE, st DIM>
Intersect point_is_in_triangle(const Point<TYPE, DIM>& p,
		const Triangle<TYPE, DIM>& t) {
	Vertex<TYPE, DIM> v1, v2, v3;
	Float d1, d2, d3;

	triangle_vertices(t, &v1, &v2, &v3);  //

	d1 = point_orientation((v1), (v2), p);
	if (d1 < 0.0)
		return TS_OUT;
	d2 = point_orientation((v2), (v3), p);
	if (d2 < 0.0)
		return TS_OUT;
	d3 = point_orientation((v3), (v1), p);
	if (d3 < 0.0)
		return TS_OUT;
	if (d1 == 0.0 || d2 == 0.0 || d3 == 0.0)
		return TS_ON;
	return TS_IN;
}

/**
 * gts_point_in_triangle_circle:
 * @p: a #GtsPoint.
 * @t: a #GtsTriangle.
 *
 * Tests if the planar projection (x, y) of @p is inside or outside
 * the circumcircle of the planar projection of @t. This function is
 * geometrically robust.
 *
 * Returns: a positive number if @p lies inside,
 * a negative number if @p lies outside and zero if @p lies on
 * the circumcircle of @t.
 */
//template<class TYPE, st DIM>
//Float point_in_triangle_circle(     //
//		const Point<TYPE, DIM>& p,  //
//		const Triangle<TYPE, DIM>& t) {
//	Vertex<TYPE, DIM> v1, v2, v3;
//	triangle_vertices(t, v1, v2, v3);
//	return incircle((double *) v1.data(), (double *) v2.data(),
//			(double *) v3.data(), (double *) &p.data());
//}
/**
 * point_in_circle:
 * @p: a #GtsPoint.
 * @p1: a #GtsPoint.
 * @p2: a #GtsPoint.
 * @p3: a #GtsPoint.
 *
 * Tests if the planar projection (x, y) of @p is inside or outside the
 * circle defined by the planar projection of @p1, @p2 and @p3.
 *
 * Returns: a positive number if @p lies inside,
 * a negative number if @p lies outside and zero if @p lies on
 * the circle.
 */
//template<class TYPE, st DIM>
//Float gts_point_in_circle(const Point<TYPE, DIM>& p,   //
//		const Point<TYPE, DIM>& p1,  //
//		const Point<TYPE, DIM>& p2,  //
//		const Point<TYPE, DIM>& p3) {
//	return incircle((double *) p1.data(), (double *) p2.data(),
//			(double *) p3.data(), (double *) p.data());
//}
//template<class TYPE, st DIM>
//Float gts_point_in_sphere(  //
//		const Point<TYPE, DIM>& p,  //
//		const Point<TYPE, DIM>& p1, //
//		const Point<TYPE, DIM>& p2,  //
//		const Point<TYPE, DIM>& p3, //
//		const Point<TYPE, DIM>& p4) {
//	return insphere((double *) &p1.data(), (double *) &p2.data(),
//			(double *) &p3.data(), (double *) &p4.data(), (double *) &p.data());
//}




/**
 * gts_point_segment_closest:
 * @p: a #GtsPoint.
 * @s: a #GtsSegment.
 * @closest: a #GtsPoint.
 *
 * Set the coordinates of @closest to the coordinates of the point belonging
 * to @s closest to @p.
 */
template<class TYPE, st DIM>
void gts_point_segment_closest(const Point<TYPE, DIM>& p,
		const Segment<TYPE, DIM>& s, Point<TYPE, DIM>& closest) {
	TYPE t, ns2;
	Point<TYPE, DIM>* p1, *p2;

	p1 = s->v1;
	p2 = s->v2;
	ns2 = point_distance2(p1, p2);

	if (ns2 == 0.0) {
		closest = (*p1);
		return;
	}

	t = ((p2->x() - p1->x()) * (p->x() - p1->x())
			+ (p2->y() - p1->y()) * (p->y() - p1->y())
			+ (p2->z() - p1->z()) * (p->z() - p1->z())) / ns2;

	if (t > 1.0)
		closest = (*p2);
	else if (t < 0.0)
		closest = (*p1);
	else
		closest.set( //
				(1. - t) * p1->x() + t * p2->x(), //
				(1. - t) * p1->y() + t * p2->y(), //
				(1. - t) * p1->z() + t * p2->z());
}

/**
 * gts_point_orientation_sos:
 * @p1: a #GtsPoint.
 * @p2: a #GtsPoint.
 * @p3: a #GtsPoint.
 *
 * Checks for orientation of the projection of three points on the
 * (x,y) plane.
 *
 * Simulation of Simplicity (SoS) is used to break ties when the
 * orientation is degenerate (i.e. @p3 lies on the line defined by
 * @p1 and @p2).
 *
 * Returns: a positive value if @p1, @p2 and @p3 appear in
 * counterclockwise order or a negative value if they appear in
 * clockwise order.
 */

template<class TYPE, st DIM>
int point_orientation_sos( //
		const Point<TYPE, DIM>& p1, //
		const Point<TYPE, DIM>& p2,  //
		const Point<TYPE, DIM>& p3 //
		) {
	Float o;

//	o = orient2d((double *) p1.data(), (double *) p2.data(),
//			(double *) p3.data());
	if (o != 0.)
		return SIGN(o);
	else {
		const Point<TYPE, DIM>* p[3];
		int sign;

		p[0] = &p1;
		p[1] = &p2;
		p[2] = &p3;

		sign = sortp(p, 3);

		/* epsilon^1/4 */
		o = ORIENT1D(p[1]->x(), p[2]->x());
		if (o != 0.)
			return -SIGN(o) * sign;

		/* epsilon^1/2 */
		o = ORIENT1D(p[1]->y(), p[2]->y());
		if (o != 0.)
			return SIGN(o) * sign;

		/* epsilon */
		o = ORIENT1D(p[0]->x(), p[2]->x());
		if (o != 0.)
			return SIGN(o) * sign;

		/* epsilon^3/2 */
		return sign;
	}
}

} //end of namespace

#endif /* POINT_H_ */
