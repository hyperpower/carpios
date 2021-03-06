#ifndef _OPERATION_HPP_
#define _OPERATION_HPP_

#include <array>
#include <cmath>
#include <memory>
#include "../geometry_define.hpp"
#include "../objects/_objects.hpp"
#include "_intersection.hpp"

namespace carpio {

enum BooleanOpType {
	INTERSECTION, UNION, DIFFERENCE, XOR
};

enum BooleanObjectType {
	SUBJECT, CLIPPING
};

template<class TYPE>
class Contour_;

template<class TYPE, St DIM>
class PointChain_;

template<typename TYPE, St DIM>
class Segment_;

template<typename TYPE, St DIM>
class Intersection_;

template<typename TYPE, St DIM>
class Operation_ {
public:
	static const St Dim = DIM;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Point;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Segment& ref_Segment;
	typedef const Segment& const_ref_Segment;
	typedef PointChain_<TYPE, DIM> PointChain;
	typedef Intersection_<TYPE, DIM> Isc;
	typedef Box_<TYPE, DIM> Box;

	typedef Point_<TYPE, 2> Point2;

	template<class CA, class CB>
	static double Distance(const CA& a, const CB& b) {
		_Distance(a, b, typename CA::Tag(), typename CB::Tag());
	}

	/*===============================================
	 * Calculate the distance between 2 Points
	 // for Point2D
	 // @param   p1 Point1
	 // @param   p2 Point1
	 // @return     Distance
	 */
	static double Distance2(const Point& p1, const Point& p2) {
		if (Dim == 1) {
			double dis = p1[0] - p2[0];
			return std::abs(dis);
		}
		if (Dim == 2) {
			return double(
					(p1.x() - p2.x()) * (p1.x() - p2.x())
							+ (p1.y() - p2.y()) * (p1.y() - p2.y()));
		}
		if (Dim == 3) {
			return double(
					(p1.x() - p2.x()) * (p1.x() - p2.x())
							+ (p1.y() - p2.y()) * (p1.y() - p2.y())
							+ (p1.z() - p2.z()) * (p1.z() - p2.z()));
		}
		SHOULD_NOT_REACH;
		return 0.0;
	}

//static double Distance(const Point& p1, const Point& p2) {
//		return std::sqrt(Distance2(p1, p2));
//	}
	static double _Distance(const Point& p1, const Point& p2, TagPoint,
			TagPoint) {
		return std::sqrt(Distance2(p1, p2));
	}

	//===============================================
	// Dot product (v1-v3) . (v2-v3)
	// for Point
	// @param    p1 Point1
	// @param    p2 Point1
	// @return   the result of dot product
	//-----------------------------------------------
	static double Dot(const Point &p1, const Point &p2, const Point &p3) {
		double sum = 0.0;
		for (St d = 0; d < Dim; d++) {
			double ds = p1[d] - p3[d];
			double de = p2[d] - p3[d];
			sum += ds * de;
		}
		return sum;
	}

	//===============================================
	// Cross product (v1 - v4) . ((v2-v4) x (v3-v4))
	// for Point
	// @param    p1 Point1
	// @param    p2 Point1
	// @return      the resualt of cross multiply
	//              The scalar triple product
	//              (also called the mixed product, box product)
	//-----------------------------------------------
	static double TripleScalar(const Point &v1, const Point&v2, const Point &v3,
			const Point &v4 = Point()) {
		if (Dim == 2) {
			/// d1    = v1 - v3
			/// d2    = v2 - v3
			/// cross = d1x * d2y - d1y * d2x
			return ((v1.x() - v3.x()) * (v2.y() - v3.y())
					- (v2.x() - v3.x()) * (v1.y() - v3.y()));
		}

		if (Dim == 3) {
			double a[3][3];
			for (short i = 0; i != 3; ++i) {
				a[0][i] = v1[i] - v4[i];
				a[1][i] = v2[i] - v4[i];
				a[2][i] = v3[i] - v4[i];
			}

			return a[0][0] * a[1][1] * a[2][2]
					+ a[0][1] * a[1][2] * a[2][0]
					+ a[0][2] * a[1][0] * a[2][1]
					- a[0][2] * a[1][1] * a[2][0]
					- a[0][1] * a[1][0] * a[2][2]
					- a[0][0] * a[1][2] * a[2][1];
		}
		SHOULD_NOT_REACH;
		return 0.0;
	}

	/** Signed area of the triangle (p1, p2, p3) */
	static inline double SignedArea( //
			const Point& v1, const Point& v2, const Point& v3) {
		ASSERT(Dim == 2);
		return TripleScalar(v1, v2, v3);
	}

	static bool IsCCW(const Point& p0, const Point& p1, const Point& p2) {
		double tmp = SignedArea(p0, p1, p2);
		if (tmp > 0)
			return true;
		else
			return false;
	}

	static bool IsCW(const Point& p0, const Point& p1, const Point& p2) {
		double tmp = SignedArea(p0, p1, p2);
		if (tmp < 0)
			return true;
		else
			return false;
	}

	static int OnWhichSide3(const Point& p0, const Point& p1, const Point& p2) {
		double tmp = SignedArea(p0, p1, p2);
		if (tmp > 0) {
			return 1;
		} else if (tmp < 0) {
			return -1;
		} else {
			return 0;
		}
	}

	static int OnWhichSide3(
			const Point& p0,
			const Point& p1,
			const Point& p2,
			const Point& p3) {
		double tmp = TripleScalar(p0, p1, p2, p3);
		if (tmp > 0) {
			return 1;
		} else if (tmp < 0) {
			return -1;
		} else {
			return 0;
		}
	}

	/** Signed area of the triangle ( (0,0), p1, p2) */
	static inline double SignedArea(const Point& p1, const Point& p2) {
		ASSERT(Dim == 2);
		return -p2.x() * (p1.y() - p2.y()) - -p2.y() * (p1.x() - p2.x());
	}

	/** Sign of triangle (p1, p2, o) */
	static inline int Sign(const Point& p1, const Point& p2, const Point& o) {
		ASSERT(Dim == 2);
		double det = (p1.x() - o.x()) * (p2.y() - o.y())
				- (p2.x() - o.x()) * (p1.y() - o.y());
		return (det < 0 ? -1 : (det > 0 ? +1 : 0));
	}

//	inline bool PointInTriangle(const Segment& s, Point& o, Point& p) {
//		ASSERT(Dim == 2);
//		int x = Sign(s.ps(), s.pe(), p);
//		return ((x == Sign(s.pe(), o, p)) && (x == Sign(o, s.ps(), p)));
//	}

	/**
	 * Get minimun loaction of
	 *
	 * for example:
	 * a   = ( 1, 2, 3)
	 *            ^
	 * b   = ( 0, 4, 2)
	 *         ^     ^
	 * res = Min(a, b);
	 * res = ( 0, 2, 2)
	 */
	static Point Min(const Point& a, const Point& b) {
		Point res;
		for (St i = 0; i < Dim; i++) {
			res[i] = std::min(a[i], b[i]);
		}
		return res;
	}
	/**
	 * Get max location of
	 */
	static Point Max(const Point& a, const Point& b) {
		Point res;
		for (St i = 0; i < Dim; i++) {
			res[i] = std::max(a[i], b[i]);
		}
		return res;
	}

	static bool CircumCircle(
			Vt xp, Vt yp,  // p
			Vt x1, Vt y1,  // p1
			Vt x2, Vt y2,  // p2
			Vt x3, Vt y3,  // p3
			Vt &xc, Vt &yc, Vt &r) {
		double m1, m2, mx1, mx2, my1, my2;
		double dx, dy, rsqr, drsqr;

		/* Check for coincident points */
		if (std::abs(y1 - y2) < SMALL && std::abs(y2 - y3) < SMALL)
			return (false);
		if (std::abs(y2 - y1) < SMALL) {
			m2 = -(x3 - x2) / (y3 - y2);
			mx2 = (x2 + x3) / 2.0;
			my2 = (y2 + y3) / 2.0;
			xc = (x2 + x1) / 2.0;
			yc = m2 * (xc - mx2) + my2;
		} else if (std::abs(y3 - y2) < SMALL) {
			m1 = -(x2 - x1) / (y2 - y1);
			mx1 = (x1 + x2) / 2.0;
			my1 = (y1 + y2) / 2.0;
			xc = (x3 + x2) / 2.0;
			yc = m1 * (xc - mx1) + my1;
		} else {
			m1 = -(x2 - x1) / (y2 - y1);
			m2 = -(x3 - x2) / (y3 - y2);
			mx1 = (x1 + x2) / 2.0;
			mx2 = (x2 + x3) / 2.0;
			my1 = (y1 + y2) / 2.0;
			my2 = (y2 + y3) / 2.0;
			xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
			yc = m1 * (xc - mx1) + my1;
		}
		dx = x2 - xc;
		dy = y2 - yc;
		rsqr = dx * dx + dy * dy;
		r = std::sqrt(rsqr);
		dx = xp - xc;
		dy = yp - yc;
		drsqr = dx * dx + dy * dy;
		return ((drsqr <= rsqr) ? true : false);
	}

	static bool CircumCircle(
			const Point2& p,
			const Point2& p1,  // p1
			const Point2& p2,  // p2
			const Point2& p3,  // p3
			Point2& pc, Vt &r) {
		return CircumCircle(
				p[0], p[1],
				p1[0], p1[1],
				p2[0], p2[1],
				p3[0], p3[1],
				pc[0], pc[1], r);
	}

	static bool IsInOnCircumCircle(
			const Point2& p,
			const Point2& p1,  // p1
			const Point2& p2,  // p2
			const Point2& p3  // p3
			) {
		Point pc;
		Vt r = 0;
		return CircumCircle(
				p[0], p[1],
				p1[0], p1[1],
				p2[0], p2[1],
				p3[0], p3[1],
				pc[0], pc[1], r);
	}

	/**
	 * p in on Triangle (p1, p2, p3)
	 */
	static bool IsInOn(
			const Point& p1,
			const Point& p2,
			const Point& p3,
			const Point &p) {
		if (IsCCW(p1, p, p2))
			return false;
		if (IsCCW(p2, p, p3))
			return false;
		if (IsCCW(p3, p, p1))
			return false;
		return true;
	}

	/**
	 * centroid  (barycenter)
	 *
	 * returns center point as a point
	 */
	static Point Centroid(const Point& v1, const Point& v2, const Point& v3) {
		Vt x = (v1.x() + v2.x() + v3.x()) / 3.0;
		Vt y = (v1.y() + v2.y() + v3.y()) / 3.0;
		Vt z = (DIM == 3) ? (v1.z() + v2.z() + v3.z()) / 3.0 : 0.0;
		return Point(x, y, z);
	}

	/**
	 * triangle_normal:
	 * @t: a #GtsTriangle.
	 * @x: the x coordinate of the normal.
	 * @y: the y coordinate of the normal.
	 * @z: the z coordinate of the normal.
	 *
	 * Computes the coordinates of the oriented normal of @t as the
	 * cross-product of two edges, using the left-hand rule. The normal is
	 * not normalized.  If this triangle is part of a closed and oriented
	 * surface, the normal points to the outside of the surface.
	 */
	static Point Normal(const Point& v1, const Point& v2, const Point& v3) {
		Vt x1 = v2.x() - v1.x();
		Vt y1 = v2.y() - v1.y();
		Vt z1 = (DIM == 3) ? v2.z() - v1.z() : 0.0;

		Vt x2 = v3.x() - v1.x();
		Vt y2 = v3.y() - v1.y();
		Vt z2 = (DIM == 3) ? v3.z() - v1.z() : 0.0;

		Vt x = y1 * z2 - z1 * y2;
		Vt y = z1 * x2 - x1 * z2;
		Vt z = (DIM == 3) ? x1 * y2 - y1 * x2 : 0.0;

		return Point(x, y, z);
	}

	// Find the circumcenter of three 2-D points by Cramer's Rule to find
	// the intersection of two perpendicular bisectors of the triangle's
	// edges.
	// http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
	//
	// Return true if successful; return false if points are collinear
	static bool Circumcenter(Vt x0, Vt y0,
			Vt x1, Vt y1,
			Vt x2, Vt y2,
			Vt& centerx, Vt& centery)
			{
		Vt D;
		Vt x0m2, y1m2, x1m2, y0m2;
		Vt x0p2, y1p2, x1p2, y0p2;
		x0m2 = x0 - x2;
		y1m2 = y1 - y2;
		x1m2 = x1 - x2;
		y0m2 = y0 - y2;
		x0p2 = x0 + x2;
		y1p2 = y1 + y2;
		x1p2 = x1 + x2;
		y0p2 = y0 + y2;

		D = x0m2 * y1m2 - x1m2 * y0m2;
		if ((D < SMALL) && (D > -SMALL))
			return false;

		centerx = (((x0m2 * x0p2 + y0m2 * y0p2) / 2 * y1m2)
				- (x1m2 * x1p2 + y1m2 * y1p2) / 2 * y0m2) / D;
		centery = (((x1m2 * x1p2 + y1m2 * y1p2) / 2 * x0m2)
				- (x0m2 * x0p2 + y0m2 * y0p2) / 2 * x1m2) / D;

		return true;
	}

	static Box BoundingBox(const Point& v1, const Point& v2, const Point& v3) {
		Point max = Max(v1, v2);
		max = Max(max, v3);
		Point min = Min(v1, v2);
		min = Min(min, v3);
		return Box(min, max);
	}
	static Box BoundingBox(const Box& box, const Point& pnew) {
		Point max = Max(box.max(), pnew);
		Point min = Min(box.min(), pnew);
		return Box(min, max);
	}

	template<class Container>
	static Box BoundingBox(const Container& con) {
		typename Container::value_type dummy;
		return _BoundingBox(con, dummy);
	}

	template<class Container>
	static Box _BoundingBox(const Container& con, Point dummy) {
		Point min, max;

		auto iter = con.begin();
		const Point& p = (*iter);
		max = p;
		min = p;
		for (; iter != con.end(); ++iter) {
			const Point& p = (*iter);
			for (St d = 0; d < Dim; d++) {
				if (p[d] > max[d]) {
					max[d] = p[d];
				}
				if (p[d] < min[d]) {
					min[d] = p[d];
				}
			}
		}
		return Box(min, max);
	}

	template<class Container>
	static Box _BoundingBox(const Container& con, Point* dummy) {
		Point min, max;

		auto iter = con.begin();
		const Point& p = *(*iter);
		max = p;
		min = p;
		for (; iter != con.end(); ++iter) {
			const Point& p = *(*iter);
			for (St d = 0; d < Dim; d++) {
				if (p[d] > max[d]) {
					max[d] = p[d];
				}
				if (p[d] < min[d]) {
					min[d] = p[d];
				}
			}
		}
		return Box(min, max);
	}

	template<class Container>
	static Box _BoundingBox(const Container& con,
			std::shared_ptr<Point> dummy) {
		Point min, max;

		auto iter = con.begin();
		const Point& p = *(*iter);
		max = p;
		min = p;
		for (; iter != con.end(); ++iter) {
			const Point& p = *(*iter);
			for (St d = 0; d < Dim; d++) {
				if (p[d] > max[d]) {
					max[d] = p[d];
				}
				if (p[d] < min[d]) {
					min[d] = p[d];
				}
			}
		}
		return Box(min, max);
	}
	/// normalize itself
	static Point Normalize(const Point &p) {
		Point r;
		r.normalize();
		return r;
	}

	static void Normalize(Point &p, const Point& dis, const Point& min) {
		for (St d = 0; d < Dim; d++) {
			p[d] = (p[d] - min[d]) / dis[d];
		}
	}

	static void UnNormalize(Point &p, const Point& dis, const Point& min) {
		for (St d = 0; d < Dim; d++) {
			p[d] = p[d] * dis[d] + min[d];
		}
	}

	template<class Container>
	static void Normalize(Container& con) {
		typename Container::value_type dummy;
		Box bb = _BoundingBox(con, dummy);

		Point dis = bb.d();
		Point min = bb.min();

		_Normalize(con, dummy, dis, min);
	}

	template<class Container>
	static void _Normalize(Container& con,
			std::shared_ptr<Point> dummy,
			const Point& dis, const Point& min) {
		for (auto& p : con) {
			Normalize(*p, dis, min);
		}
	}

	template<class Container>
	static void _Normalize(Container& con,
			Point* dummy,
			const Point& dis, const Point& min) {
		for (auto& p : con) {
			Normalize(*p, dis, min);
		}
	}

	static bool IsInOn(
			const Vt& xmin, const Vt& xmax,
			const Vt& ymin, const Vt& ymax,
			const Vt& x, const Vt& y) {
		ASSERT(ymin <= ymax);
		ASSERT(xmin <= xmax);
		if (ymin == ymax) {
			return (xmin <= x) && (x <= xmax);
		}
		if (xmin == xmax) {
			return (ymin <= y) && (y <= ymax);
		}
		return ((xmin <= x) && (x <= xmax)) && ((ymin <= y) && (y <= ymax));
	}
	/**
	 * p in on Box (min, max, p)
	 */
	static bool IsInOn(
			const Point& pmin,
			const Point& pmax,
			const Point& p) {
		if (Dim == 2) {
			return IsInOn(pmin[0], pmax[0], pmin[1], pmax[1], p[0], p[1]);
		} else {
			SHOULD_NOT_REACH;
			return false;
		}
	}

	static bool IsInBox(const Segment &s, const Point &pt) {
		ASSERT(!s.empty());
		if (s.is_horizontal()) {
			return (((s.psx() <= pt.x()) && (pt.x() <= s.pex()))
					|| ((s.pex() <= pt.x()) && (pt.x() <= s.psx())));
		}
		if (s.is_vertical()) {
			return (((s.psy() <= pt.y()) && (pt.y() <= s.pey()))
					|| ((s.pey() <= pt.y()) && (pt.y() <= s.psy())));
		}
		return (((s.psx() <= pt.x()) && (pt.x() <= s.pex()))
				|| ((s.pex() <= pt.x()) && (pt.x() <= s.psx())))
				&& (((s.psy() <= pt.y()) && (pt.y() <= s.pey()))
						|| ((s.pey() <= pt.y()) && (pt.y() <= s.psy())));
	}

// Find intersetion on 1d
//                ----------------
//         -------|---------     |
// --------|------|--------|-----|---------
//        u0     v0        u1   v1
//               w[0]      w[1]
// should be move to intersection class
	static int FindIntersection(Vt u0, Vt u1, Vt v0, Vt v1, Vt w[2]) {
		if ((u1 < v0) || (u0 > v1))
			return 0;
		if (u1 > v0) {
			if (u0 < v1) {
				w[0] = (u0 < v0) ? v0 : u0;
				w[1] = (u1 > v1) ? v1 : u1;
				return 2;
			} else {
				// u0 == v1
				w[0] = u0;
				return 1;
			}
		} else {
			// u1 == v0
			w[0] = u1;
			return 1;
		}
	}

// Find intersection on 2d
// should be move to intersection class
	static int FindIntersection(const Segment& seg0, const Segment& seg1,
			Point& pi0, Point& pi1) {
		ASSERT(Dim == 2);
		const Point& p0 = seg0.ps();
		Point d0(seg0.pex() - p0.x(), seg0.pey() - p0.y());
		const Point& p1 = seg1.ps();
		Point d1(seg1.pex() - p1.x(), seg1.pey() - p1.y());
		Vt sqrEpsilon = 1e-8; // it was 0.001 before
		Point E(p1.x() - p0.x(), p1.y() - p0.y());
		Vt kross = d0.x() * d1.y() - d0.y() * d1.x();
		Vt sqrKross = kross * kross;
		Vt sqrLen0 = d0.x() * d0.x() + d0.y() * d0.y();
		Vt sqrLen1 = d1.x() * d1.x() + d1.y() * d1.y();

		if (sqrKross > sqrEpsilon * sqrLen0 * sqrLen1) {
			// lines of the segments are not parallel
			Vt s = (E.x() * d1.y() - E.y() * d1.x()) / kross;
			if ((s < 0) || (s > 1)) {
				return 0;
			}
			Vt t = (E.x() * d0.y() - E.y() * d0.x()) / kross;
			if ((t < 0) || (t > 1)) {
				return 0;
			}
			// intersection of lines is a point an each segment
			pi0.x() = p0.x() + s * d0.x();
			pi0.y() = p0.y() + s * d0.y();
			if (pi0.dist(seg0.ps()) < 1e-9)
				pi0 = seg0.ps();
			if (pi0.dist(seg0.pe()) < 1e-9)
				pi0 = seg0.pe();
			if (pi0.dist(seg1.ps()) < 1e-9)
				pi0 = seg1.ps();
			if (pi0.dist(seg1.pe()) < 1e-9)
				pi0 = seg1.pe();
			return 1;
		}

		// lines of the segments are parallel
		Vt sqrLenE = E.x() * E.x() + E.y() * E.y();
		kross = E.x() * d0.y() - E.y() * d0.x();
		sqrKross = kross * kross;
		if (sqrKross > sqrEpsilon * sqrLen0 * sqrLenE) {
			// lines of the segment are different
			return 0;
		}

		// Lines of the segments are the same. Need to test for overlap of segments.
		// s0 = Dot (D0, E) * sqrLen0
		Vt s0 = (d0.x() * E.x() + d0.y() * E.y()) / sqrLen0;
		// s1 = s0 + Dot (D0, D1) * sqrLen0
		Vt s1 = s0 + (d0.x() * d1.x() + d0.y() * d1.y()) / sqrLen0;
		Vt smin = std::min(s0, s1);
		Vt smax = std::max(s0, s1);
		Vt w[2];
		int imax = FindIntersection(0.0, 1.0, smin, smax, w);

		if (imax > 0) {
			pi0.x() = p0.x() + w[0] * d0.x();
			pi0.y() = p0.y() + w[0] * d0.y();
			if (pi0.dist(seg0.ps()) < 1e-9)
				pi0 = seg0.ps();
			if (pi0.dist(seg0.pe()) < 1e-9)
				pi0 = seg0.pe();
			if (pi0.dist(seg1.ps()) < 1e-9)
				pi0 = seg1.ps();
			if (pi0.dist(seg1.pe()) < 1e-9)
				pi0 = seg1.pe();
			if (imax > 1) {
				pi1.x() = p0.x() + w[1] * d0.x();
				pi1.y() = p0.y() + w[1] * d0.y();
			}
		}

		return imax;
	}

	template<typename _ForwardIterator>
	static bool IsSimple(_ForwardIterator begin, _ForwardIterator end,
			bool isclose = false) {
		bool sametype = std::is_same<typename _ForwardIterator::value_type,
				Point>::value;
		ASSERT(sametype == true);

		typedef _ForwardIterator iterator;

		long count = 0;
		for (iterator iter = begin; iter != end; ++iter) {
			count++;
			if (count > 3) {
				break;
			}
		}
		if (count <= 3) {
			return true;
		}
		//
		iterator iter_end = end;
		std::advance(iter_end, -1);
		iterator iter_end2 = end;
		std::advance(iter_end2, -1);
		for (iterator iter0 = begin;
				iter0 != iter_end; ++iter0) {
			iterator iter1 = std::next(iter0);

			iterator iterm0 = std::prev(iter0);
			for (iterator iter2 = begin;
					iter2 != iter_end; ++iter2) {
				iterator iter3 = std::next(iter2);
				if (iter2 == iterm0 || iter2 == iter0 || iter2 == iter1) {
					continue;
				}
				//
//				std::cout << "P1 = " << (*iter0);
//				std::cout << " P2 = " << (*iter1);
//				std::cout << " P3 = " << (*iter2);
//				std::cout << " P4 = " << (*iter3);

				bool res = Isc::Check_asSegment(
						*iter0, *iter1, *iter2, *iter3,
						INTERSECT_NORMAL
								| INTERSECT_POINT_POINT
								| INTERSECT_POINT_POINT_2
								| INTERSECT_POINT_SEGMENT
								| INTERSECT_POINT_SEGMENT_2);
//				std::cout << " inter = " << res << "\n";
				if (res == true) {
					return false;
				}
			}
		}
		if (isclose) {
			iter_end = end;
			std::advance(iter_end, -2);
			iterator iter0 = std::prev(end, 1);
			iterator iter1 = begin;
			for (iterator iter2 = std::next(iter1);
					iter2 != iter_end; ++iter2) {
				iterator iter3 = std::next(iter2);
				//
//				std::cout << "P1 = " << (*iter0);
//				std::cout << " P2 = " << (*iter1);
//				std::cout << " P3 = " << (*iter2);
//				std::cout << " P4 = " << (*iter3);

				bool res = Isc::Check_asSegment(
						*iter0, *iter1, *iter2, *iter3,
						INTERSECT_NORMAL
								| INTERSECT_POINT_POINT
								| INTERSECT_POINT_POINT_2
								| INTERSECT_POINT_SEGMENT
								| INTERSECT_POINT_SEGMENT_2);
//				std::cout << " inter = " << res << "\n";
				if (res == true) {
					return false;
				}
			}
		}
		return true;
	}

};

//template<typename TYPE, St DIM>
//bool IsBoxCross(const Segment_<TYPE, DIM> &s1, const Segment_<TYPE, DIM> &s2) {
//	return IsInBox(s1, s2.ps()) || IsInBox(s1, s2.pe()) || IsInBox(s2, s1.ps())
//			|| IsInBox(s2, s1.pe());
//}
//template<typename TYPE>
//int OnWhichSide3(const Segment_<TYPE, 2> &s, const Point_<TYPE, 2> &pt) {
//	Float rcro = Cro(s.pe(), pt, s.ps());
//	if (rcro == 0.0) {
//		return 0;
//	} else if (rcro < 0) {
//		return -1;
//	} else {
//		return 1;
//	}
//}
/*
 * Segment -------------------  Segment
 */

//template<typename TYPE>
//int IntersectType(const Segment_<TYPE, 2> &s1, const Segment_<TYPE, 2> &s2) {
//	if (!IsBoxCross(s1, s2)) {
//		return NO_INTERSECT;
//	} else {
//		//step 1
//		int s12s = OnWhichSide3(s1, s2.ps());
//		int s12e = OnWhichSide3(s1, s2.pe());
//		if (s12s == s12e) { //ignore the both equal to 0, overlap is not intersect
//			return NO_INTERSECT;
//		}
//		int s21s = OnWhichSide3(s2, s1.ps());
//		int s21e = OnWhichSide3(s2, s1.pe());
//		if (s21s == s21e) { //ignore the both equal to 0, overlap is not intersect
//			return NO_INTERSECT;
//		}
//		if ((s12s + s12e) == 0 && (s21s + s21e) == 0) {
//			return INTERSECT;
//		}
//		int res = INTERSECT;
//		if (s12s == 0)
//			res = res | START_2;
//		if (s12e == 0)
//			res = res | END_2;
//		if (s21s == 0)
//			res = res | START_1;
//		if (s21e == 0)
//			res = res | END_1;
//		return res;
//	}
//}
////template<typename TYPE, St DIM>
////bool IsIntersect(const Segment_<TYPE, DIM> &s1, const Segment_<TYPE, DIM> &s2) {
////	int type = IntersectType(s1, s2);
////	return (type | INTERSECT) == type ? true : false;
////}
////template<typename TYPE, St DIM>
////bool IsIntersect(const Point_<TYPE, 2>& s1s, const Point_<TYPE, 2>& s1e,
////		const Point_<TYPE, 2>& s2s, const Point_<TYPE, 2>& s2e) {
////	Segment_<TYPE, DIM> s1(s1s, s1e);
////	Segment_<TYPE, DIM> s2(s2s, s2e);
////	return IsIntersect(s1, s2);
////}
////template<typename TYPE>
////Point_<TYPE, 2> CalIntersect(const Segment_<TYPE, 2> &s1,
////		const Segment_<TYPE, 2> &s2) {
////	ASSERT(IsIntersect(s1, s2));
////	Float x1, x2, y1, y2;
////	Float x3, x4, y3, y4;
////	Float resxx;
////	Float resyy;
////	x1 = Float(s1.psx());
////	y1 = Float(s1.psy());
////	x2 = Float(s1.pex());
////	y2 = Float(s1.pey());
////	x3 = Float(s2.psx());
////	y3 = Float(s2.psy());
////	x4 = Float(s2.pex());
////	y4 = Float(s2.pey());
////	Float b1 = (y2 - y1) * x1 + (x1 - x2) * y1;
////
////	Float b2 = (y4 - y3) * x3 + (x3 - x4) * y3;
////	if (s1.is_horizontal() && s2.is_vertical()) {
////		resxx = x3;
////		resyy = y1;
////	} else if (s2.is_horizontal() && s1.is_vertical()) {
////		resxx = x1;
////		resyy = y3;
////	} else if (s1.is_horizontal()) {
////		resxx = (b2 + (x4 - x3) * y1) / (y4 - y3);
////		resyy = y1;
////	} else if (s1.is_vertical()) {
////		resxx = x1;
////		resyy = (b2 - (y4 - y3) * x1) / (x3 - x4);
////	} else if (s2.is_horizontal()) {
////		resxx = (b1 + (x2 - x1) * y3) / (y2 - y1);
////		resyy = y3;
////	} else if (s2.is_vertical()) {
////		resxx = x3;
////		resyy = (b1 - (y2 - y1) * x3) / (x1 - x2);
////	} else {
////		Float d = (x2 - x1) * (y4 - y3) - (x4 - x3) * (y2 - y1);
////		Float d1 = b2 * (x2 - x1) - b1 * (x4 - x3);
////		Float d2 = b2 * (y2 - y1) - b1 * (y4 - y3);
////		ASSERT(d != 0);
////		resxx = d1 / d;
////		resyy = d2 / d;
////	}
////	return Point_<TYPE, 2>(resxx, resyy);
////}
/////*
//// * Segment -------- Line
//// */
////template<typename TYPE, St DIM>
////bool IsIntersect(const Segment_<TYPE, DIM> &s1, //
////		const Point_<TYPE, DIM>& ps, //
////		const Point_<TYPE, DIM>& pe) {
////	// ps --> pe are two points on line
////	Segment_<TYPE, DIM> s(ps, pe);
////	int s12s = OnWhichSide3(s, s1.ps());
////	int s12e = OnWhichSide3(s, s1.pe());
////	int sum = s12s + s12e;
////	// co-linear is not intersect
////	if (sum == 0 || sum == -2 || sum == 2) {
////		return false;
////	} else {
////		return true;
////	}
////}
/////*
//// * Point   -------- Polygon
//// */
////template<typename TYPE>
////double _WindingNumber(const Point_<TYPE, 2>& ref, const Point_<TYPE, 2>& vi,
////		const Point_<TYPE, 2>& vip) {
////	Point_<TYPE, 2> refh(ref.x() + 1.0, ref.y());
////	Float wn = 0;
////	Float a = Cro(refh, vi, ref);
////	Float b = Cro(refh, vip, ref);
////	if (a == 0 && b == 0) {
////		return wn;
////	}
////	if (a * b < 0) {
////		//vi vi+1 crosses the x
////		Float c = Cro(vip, ref, vi);
////		if ((c > 0 && a < 0) || (c < 0 && a > 0)) {
////			//vi vi+1 crosses the positive x
////			if (a < 0) {
////				wn++;
////			} else {
////				wn--;
////			}
////		}
////	} else if (a == 0 && (vi.x() > ref.x())) {
////		if (b > 0) {
////			wn = wn + 0.5;
////		} else {
////			wn = wn - 0.5;
////		}
////	} else if (b == 0 && (vip.x() > ref.x())) {
////		if (a < 0) {
////			wn = wn + 0.5;
////		} else {
////			wn = wn - 0.5;
////		}
//	}
//	return wn;
//}
//template<typename TYPE>
//Float WindingNumber(const Polygon_<TYPE>& poly, const Point_<TYPE, 2>& ref) {
//	Float wn = 0;    // the  winding number counter
//	// loop through all edges of the polygon
//	for (int i = 0; i < poly.size_vertexs() - 1; i++) { // edge from V[i] to  V[i+1]
//		wn += _WindingNumber(ref, poly.v(i), poly.v(i + 1));
//	}
//	wn += _WindingNumber(ref, poly.v(poly.size_vertexs() - 1), poly.v(0));
//	return wn;
//}

//template<typename TYPE>
//bool IsOut(const Polygon_<TYPE>& poly, const Point_<TYPE, 2>& ref) {
//	return (0 == WindingNumber(poly, ref)) ? true : false;
//}
/*
 *  Is the subject inside of the polygon
 */
//template<typename TYPE>
//bool IsIn(const Polygon_<TYPE>& poly, const Polygon_<TYPE>& sub) {
//	for (St i = 0; i < poly.size_vertexs(); ++i) {
//		Float wn = WindingNumber(poly, sub.v(i));
//		if (0 == wn || -1 == wn) {
//			return false;
//		}
//	}
//	return true;
//}
/*
 * Segment -------- Polygon
 */
//template<typename TYPE>
//ArrayListT<Segment_<TYPE, 2> > ToArraySegment(const Polygon_<TYPE>& p) {
//	St n = p.size_vertexs();
//	ArrayListT<Segment_<TYPE, 2> > as(n);
//	for (St i = 0; i < n - 1; i++) {
//		as[i].reconstruct(p.v(i), p.v(i + 1));
//	}
//	as[n - 1].reconstruct(p.v(n - 1), p.v(0));
//	return as;
//}
//
//template<typename TYPE>
//Float Perimeter(const Polygon_<TYPE>& ap) {
//	ASSERT(!ap.empty());
//	Float res = 0;
//	for (int i = 1; i < ap.size_() - 1; i++) {
//		res += Distance(ap.v(i), ap.v(i + 1));
//	}
//	res += Distance(ap.v(ap.size() - 1), ap.v(0));
//	return res;
//}
}

#endif
