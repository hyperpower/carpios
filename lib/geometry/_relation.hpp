#ifndef _RELATION_HPP_
#define _RELATION_HPP_

#include "geometry_define.hpp"
#include <array>
#include "_point.hpp"
#include "_segment.hpp"
#include "_polygon.hpp"

namespace carpio {

template<typename TYPE, St DIM>
class Operation_ {
public:
	static const St Dim = DIM;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Point;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Segment_<TYPE, DIM>& ref_Segment;
	typedef const Segment_<TYPE, DIM>& const_ref_Segment;

	/*===============================================
	 * Calculate the distance between 2 Points
	 // for Point2D
	 // @param   p1 Point1
	 // @param   p2 Point1
	 // @return     Distance
	 */
	static double Distance(const Point& p1, const Point& p2) {
		if (Dim == 1) {
			return double(std::abs(p1.x() - p2.x()));
		}
		if (Dim == 2) {
			return sqrt(
					double(
							(p1.x() - p2.x()) * (p1.x() - p2.x())
									+ (p1.y() - p2.y()) * (p1.y() - p2.y())));
		}
		if (Dim == 3) {
			return sqrt(
					double(
							(p1.x() - p2.x()) * (p1.x() - p2.x())
									+ (p1.y() - p2.y()) * (p1.y() - p2.y())
									+ (p1.z() - p2.z()) * (p1.z() - p2.z())));
		}
		SHOULD_NOT_REACH;
		return 0.0;
	}

	//===============================================
	// Dot multiply (sp-op).(ep-op)
	// for Point2D
	// @param    p1 Point1
	// @param    p2 Point1
	// @return      the resualt of dot multiply
	//-----------------------------------------------
	static double Dot(const Point &sp, const Point &ep, const Point &op) {
		double sum = 0.0;
		for (St d = 0; d < Dim; d++) {
			double ds = sp[d] - op[d];
			double de = ep[d] - op[d];
			sum += ds * de;
		}
		return sum;
	}

	/** Signed area of the triangle (p0, p1, p2) */
	static inline double SignedArea( //
			const Point& p0, const Point& p1, const Point& p2) {
		ASSERT(Dim == 2);
		return (p0.x() - p2.x()) * (p1.y() - p2.y())
				- (p1.x() - p2.x()) * (p0.y() - p2.y());
	}

	/** Signed area of the triangle ( (0,0), p1, p2) */
	static inline float SignedArea(const Point& p1, const Point& p2) {
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

	inline bool PointInTriangle(const Segment& s, Point& o, Point& p) {
		ASSERT(Dim == 2);
		int x = Sign(s.ps(), s.pe(), p);
		return ((x == Sign(s.pe(), o, p)) && (x == Sign(o, s.ps(), p)));
	}

	//===============================================
	// Cross multiply (sp-op).(ep-op)
	// for Point2D
	// @param    p1 Point1
	// @param    p2 Point1
	// @return      the resualt of dot multiply
	//-----------------------------------------------
	static double Cro(const Point &v1, const Point&v2, const Point &v3,
			const Point &v4) {
		if (Dim == 2) {
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

			return a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0]
					+ a[0][2] * a[1][0] * a[2][1] - a[0][2] * a[1][1] * a[2][0]
					- a[0][1] * a[1][0] * a[2][2] - a[0][0] * a[1][2] * a[2][1];
		}
		SHOULD_NOT_REACH;
		return 0.0;
	}

	static bool IsInBox(const Vt& xmin, const Vt& xmax, const Vt& ymin, const Vt& ymax,
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

};

template<typename TYPE, St DIM>
bool IsBoxCross(const Segment_<TYPE, DIM> &s1, const Segment_<TYPE, DIM> &s2) {
	return IsInBox(s1, s2.ps()) || IsInBox(s1, s2.pe()) || IsInBox(s2, s1.ps())
			|| IsInBox(s2, s1.pe());
}
template<typename TYPE>
int OnWhichSide3(const Segment_<TYPE, 2> &s, const Point_<TYPE, 2> &pt) {
	Float rcro = Cro(s.pe(), pt, s.ps());
	if (rcro == 0.0) {
		return 0;
	} else if (rcro < 0) {
		return -1;
	} else {
		return 1;
	}
}
/*
 * Segment -------------------  Segment
 */

enum IntersectType {
	NO_INTERSECT = 0x10000,
	INTERSECT = 0x1,
	START_1 = 0x10,
	END_1 = 0x20,
	START_2 = 0x100,
	END_2 = 0x200,
};

template<typename TYPE>
int IntersectType(const Segment_<TYPE, 2> &s1, const Segment_<TYPE, 2> &s2) {
	if (!IsBoxCross(s1, s2)) {
		return NO_INTERSECT;
	} else {
		//step 1
		int s12s = OnWhichSide3(s1, s2.ps());
		int s12e = OnWhichSide3(s1, s2.pe());
		if (s12s == s12e) { //ignore the both equal to 0, overlap is not intersect
			return NO_INTERSECT;
		}
		int s21s = OnWhichSide3(s2, s1.ps());
		int s21e = OnWhichSide3(s2, s1.pe());
		if (s21s == s21e) { //ignore the both equal to 0, overlap is not intersect
			return NO_INTERSECT;
		}
		if ((s12s + s12e) == 0 && (s21s + s21e) == 0) {
			return INTERSECT;
		}
		int res = INTERSECT;
		if (s12s == 0)
			res = res | START_2;
		if (s12e == 0)
			res = res | END_2;
		if (s21s == 0)
			res = res | START_1;
		if (s21e == 0)
			res = res | END_1;
		return res;
	}
}

template<typename TYPE, St DIM>
bool IsIntersect(const Segment_<TYPE, DIM> &s1, const Segment_<TYPE, DIM> &s2) {
	int type = IntersectType(s1, s2);
	return (type | INTERSECT) == type ? true : false;
}
template<typename TYPE, St DIM>
bool IsIntersect(const Point_<TYPE, 2>& s1s, const Point_<TYPE, 2>& s1e,
		const Point_<TYPE, 2>& s2s, const Point_<TYPE, 2>& s2e) {
	Segment_<TYPE, DIM> s1(s1s, s1e);
	Segment_<TYPE, DIM> s2(s2s, s2e);
	return IsIntersect(s1, s2);
}
template<typename TYPE>
Point_<TYPE, 2> CalIntersect(const Segment_<TYPE, 2> &s1,
		const Segment_<TYPE, 2> &s2) {
	ASSERT(IsIntersect(s1, s2));
	Float x1, x2, y1, y2;
	Float x3, x4, y3, y4;
	Float resxx;
	Float resyy;
	x1 = Float(s1.psx());
	y1 = Float(s1.psy());
	x2 = Float(s1.pex());
	y2 = Float(s1.pey());
	x3 = Float(s2.psx());
	y3 = Float(s2.psy());
	x4 = Float(s2.pex());
	y4 = Float(s2.pey());
	Float b1 = (y2 - y1) * x1 + (x1 - x2) * y1;

	Float b2 = (y4 - y3) * x3 + (x3 - x4) * y3;
	if (s1.is_horizontal() && s2.is_vertical()) {
		resxx = x3;
		resyy = y1;
	} else if (s2.is_horizontal() && s1.is_vertical()) {
		resxx = x1;
		resyy = y3;
	} else if (s1.is_horizontal()) {
		resxx = (b2 + (x4 - x3) * y1) / (y4 - y3);
		resyy = y1;
	} else if (s1.is_vertical()) {
		resxx = x1;
		resyy = (b2 - (y4 - y3) * x1) / (x3 - x4);
	} else if (s2.is_horizontal()) {
		resxx = (b1 + (x2 - x1) * y3) / (y2 - y1);
		resyy = y3;
	} else if (s2.is_vertical()) {
		resxx = x3;
		resyy = (b1 - (y2 - y1) * x3) / (x1 - x2);
	} else {
		Float d = (x2 - x1) * (y4 - y3) - (x4 - x3) * (y2 - y1);
		Float d1 = b2 * (x2 - x1) - b1 * (x4 - x3);
		Float d2 = b2 * (y2 - y1) - b1 * (y4 - y3);
		ASSERT(d != 0);
		resxx = d1 / d;
		resyy = d2 / d;
	}
	return Point_<TYPE, 2>(resxx, resyy);
}
/*
 * Segment -------- Line
 */
template<typename TYPE, St DIM>
bool IsIntersect(const Segment_<TYPE, DIM> &s1, //
		const Point_<TYPE, DIM>& ps, //
		const Point_<TYPE, DIM>& pe) {
	// ps --> pe are two points on line
	Segment_<TYPE, DIM> s(ps, pe);
	int s12s = OnWhichSide3(s, s1.ps());
	int s12e = OnWhichSide3(s, s1.pe());
	int sum = s12s + s12e;
	// co-linear is not intersect
	if (sum == 0 || sum == -2 || sum == 2) {
		return false;
	} else {
		return true;
	}
}
/*
 * Point   -------- Polygon
 */
template<typename TYPE>
double _WindingNumber(const Point_<TYPE, 2>& ref, const Point_<TYPE, 2>& vi,
		const Point_<TYPE, 2>& vip) {
	Point_<TYPE, 2> refh(ref.x() + 1.0, ref.y());
	Float wn = 0;
	Float a = Cro(refh, vi, ref);
	Float b = Cro(refh, vip, ref);
	if (a == 0 && b == 0) {
		return wn;
	}
	if (a * b < 0) {
		//vi vi+1 crosses the x
		Float c = Cro(vip, ref, vi);
		if ((c > 0 && a < 0) || (c < 0 && a > 0)) {
			//vi vi+1 crosses the positive x
			if (a < 0) {
				wn++;
			} else {
				wn--;
			}
		}
	} else if (a == 0 && (vi.x() > ref.x())) {
		if (b > 0) {
			wn = wn + 0.5;
		} else {
			wn = wn - 0.5;
		}
	} else if (b == 0 && (vip.x() > ref.x())) {
		if (a < 0) {
			wn = wn + 0.5;
		} else {
			wn = wn - 0.5;
		}
	}
	return wn;
}
template<typename TYPE>
Float WindingNumber(const Polygon_<TYPE>& poly, const Point_<TYPE, 2>& ref) {
	Float wn = 0;    // the  winding number counter
	// loop through all edges of the polygon
	for (int i = 0; i < poly.size_vertexs() - 1; i++) { // edge from V[i] to  V[i+1]
		wn += _WindingNumber(ref, poly.v(i), poly.v(i + 1));
	}
	wn += _WindingNumber(ref, poly.v(poly.size_vertexs() - 1), poly.v(0));
	return wn;
}

template<typename TYPE>
bool IsOut(const Polygon_<TYPE>& poly, const Point_<TYPE, 2>& ref) {
	return (0 == WindingNumber(poly, ref)) ? true : false;
}
/*
 *  Is the subject inside of the polygon
 */
template<typename TYPE>
bool IsIn(const Polygon_<TYPE>& poly, const Polygon_<TYPE>& sub) {
	for (St i = 0; i < poly.size_vertexs(); ++i) {
		Float wn = WindingNumber(poly, sub.v(i));
		if (0 == wn || -1 == wn) {
			return false;
		}
	}
	return true;
}

/*
 * Segment -------- Polygon
 */
template<typename TYPE>
ArrayListT<Segment_<TYPE, 2> > ToArraySegment(const Polygon_<TYPE>& p) {
	St n = p.size_vertexs();
	ArrayListT<Segment_<TYPE, 2> > as(n);
	for (St i = 0; i < n - 1; i++) {
		as[i].reconstruct(p.v(i), p.v(i + 1));
	}
	as[n - 1].reconstruct(p.v(n - 1), p.v(0));
	return as;
}
template<typename TYPE>
bool IsSimple(const Polygon_<TYPE>& ap) {
	int nap = ap.size_vertexs();
	//ASSERT(nap >= 3);
	if (nap == 3) {
		return true;
	}
	int i = 0;
	for (int j = i + 2; j < nap - 1; j++) {
		if (IsIntersect(ap.v(i), ap.v(i + 1), ap.v(j), ap.v(j + 1))) {
			return false;
		}
	}
	for (i = 1; i < nap - 2; i++) {
		for (int j = i + 2; j < nap; j++) {
			if (IsIntersect(ap.v(i), ap.v(i + 1), ap.v(j),
					ap.v((j + 1) % nap))) {
				return false;
			}
		}
	}
	return true;
}

template<typename TYPE>
Float Perimeter(const Polygon_<TYPE>& ap) {
	ASSERT(!ap.empty());
	Float res = 0;
	for (int i = 1; i < ap.size_() - 1; i++) {
		res += Distance(ap.v(i), ap.v(i + 1));
	}
	res += Distance(ap.v(ap.size() - 1), ap.v(0));
	return res;
}
}

#endif
