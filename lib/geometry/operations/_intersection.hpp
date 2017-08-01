#ifndef _INTERSECTION_HPP_
#define _INTERSECTION_HPP_

#include <cmath>
#include "../geometry_define.hpp"
#include "../objects/_objects.hpp"
#include "_operation.hpp"
#include "../../utility/any.hpp"

namespace carpio {

enum SegmentIntersectType {
	NO_INTERSECT = 0,
	INTERSECT = 1 << 1,
	START_1 = 1 << 2,
	END_1 = 1 << 3,
	START_2 = 1 << 4,
	END_2 = 1 << 5,
	OVERLAP = 1 << 6,
};

enum SegmentIntersectCase {
	INTERSECT_NORMAL = 1 << 1,
	INTERSECT_POINT_SEGMENT = 1 << 2,
	INTERSECT_POINT_POINT = 1 << 3,
	INTERSECT_POINT_SEGMENT_2 = 1 << 4, //OVERLAP
	INTERSECT_POINT_POINT_2 = 1 << 5,   //SAME
};

template<typename TYPE, St DIM>
class Box_;

template<typename TYPE, St DIM>
class BBox_;

template<typename TYPE, St DIM>
class Segment_;

template<typename TYPE, St DIM>
class PointChain_;

template<typename TYPE, St DIM>
class Operation_;

template<typename TYPE, St DIM>
class Intersection_ {
public:
	static const St Dim = DIM;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Point;
	typedef Box_<TYPE, DIM> Box;
	typedef BBox_<TYPE, DIM> BBox;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Segment_<TYPE, DIM>& ref_Segment;
	typedef const Segment_<TYPE, DIM>& const_ref_Segment;
	typedef PointChain_<TYPE, DIM> PointChain;

	typedef Operation_<TYPE, Dim> Op;

public:
	template<class CA, class CB>
	static bool Check(const CA& a, const CB& b) {
		return _Check(a, b, typename CA::Tag(), typename CB::Tag());
	}

	template<class CA, class CB, class ARG>
	static bool Check(const CA& a, const CB& b, const ARG& arg) {
		return _Check(a, b, arg, typename CA::Tag(), typename CB::Tag());
	}

	template<class CA, class CB, typename ARG>
	static bool Find(const CA& a, const CB& b, const ARG& arg = ARG()) {
		_Find(a, b, arg, typename CA::Tag(), typename CB::Tag());
	}

	static bool Check_asBox(
			const Point& pmin1,
			const Point& pmax1,
			const Point& pmin2,
			const Point& pmax2) {
		for (St d = 0; d < Dim; d++) {
			if (pmax1[d] < pmin2[d]) {
				return false;
			}
			if (pmin1[d] > pmax2[d]) {
				return false;
			}
		}
		return true;
	}

	static bool Check_asSegmentBox(
			const Point& p11,
			const Point& p12,
			const Point& p21,
			const Point& p22) {
		Point min1 = Op::Min(p11, p12);
		Point max1 = Op::Max(p11, p12);
		Point min2 = Op::Min(p21, p22);
		Point max2 = Op::Max(p21, p22);
		return Check_asBox(min1, max1, min2, max2);
	}

	/// ================================================
	/// protected part
	template<class CA, class CB>
	static bool _Check(
			const CA& a, const CB& b,
			TagPoint, TagPoint) {
		return a == b;
	}

	static bool _Check(
			const Point& a,
			const Point& b,
			const typename Point::Vt& arg,
			TagPoint,
			TagPoint) {
		return Op::Distance(a, b) < arg;
	}

	static bool _Check(
			const Box& b1, const Box& b2,
			typename Box::Tag, typename Box::Tag) {
		return Check_asBox(b1.min(), b1.max(), b2.min(), b2.max());
	}
	static bool _Check(
			const BBox& b1, const BBox& b2,
			typename BBox::Tag, typename BBox::Tag) {
		return Check_asBox(b1.min(), b1.max(), b2.min(), b2.max());
	}

	static bool _Check(
			const BBox& b1, const BBox& b2, bool inside,
			typename BBox::Tag, typename BBox::Tag) {
		if (inside) {
			const Any& aobj1 = b1.get_obj();
			const Any& aobj2 = b2.get_obj();
			if (aobj1.type() == typeid(Segment)
					&& aobj2.type() == typeid(Segment)) {
				return Check(any_cast<Segment>(aobj1), any_cast<Segment>(aobj2));
			}
			SHOULD_NOT_REACH;
			return false;
		} else {
			return Check(b1.min(), b1.max(), b2.min(), b2.max());
		}
	}

	static int CheckType_asSegment(
			const Point& p11, const Point& p12,
			const Point& p21, const Point& p22) {
		int type = 0;
		if (!Check_asSegmentBox(p11, p12, p21, p22)) {
			return NO_INTERSECT;
		} else {
			//step 1
			int s12s = Op::OnWhichSide3(p11, p12, p21);
			int s12e = Op::OnWhichSide3(p11, p12, p22);
			if (s12s == s12e) { //ignore the both equal to 0, overlap is not intersect
				if (s12s == 0) {
					// overlap colinear
					if (p11 == p21) {
						type = type | START_1 | START_2;
					}
					if (p12 == p21) {
						type = type | END_1 | START_2;
					}
					if (p11 == p22) {
						type = type | START_1 | END_2;
					}
					if (p12 == p22) {
						type = type | END_1 | END_2;
					}
					if (type == NO_INTERSECT) {
						return OVERLAP;
					} else {
						return type;
					}
				}
				return NO_INTERSECT;
			}
			int s21s = Op::OnWhichSide3(p21, p22, p11);
			int s21e = Op::OnWhichSide3(p21, p22, p12);
			if (s21s == s21e) {
				return NO_INTERSECT;
			}
			if ((s12s + s12e) == 0 && (s21s + s21e) == 0) {
				return INTERSECT;
			}
			int res = NO_INTERSECT;
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

	static SegmentIntersectCase _ToSegmentIntersectCase(int type) {
		if (type == INTERSECT) {
			return INTERSECT_NORMAL;
		}
		int flagp1 = 0, flagp2 = 0;
		if (type & START_1) {
			flagp1++;
		}
		if (type & END_1) {
			flagp1++;
		}
		if (type & START_2) {
			flagp2++;
		}
		if (type & END_2) {
			flagp2++;
		}
		if (type & OVERLAP) {
			return INTERSECT_POINT_SEGMENT_2;
		}
		if (flagp1 == 2 && flagp2 == 2) {
			return INTERSECT_POINT_POINT_2;
		} else if (flagp1 == 1 && flagp2 == 1) {
			return INTERSECT_POINT_POINT;
		} else {
			return INTERSECT_POINT_SEGMENT;
		}
	}

	static bool Check_asSegment(
			const Point& p11, const Point& p12,
			const Point& p21, const Point& p22,
			int typeinclude = INTERSECT_NORMAL) {
		int type = CheckType_asSegment(p11, p12, p21, p22);
		if (type == NO_INTERSECT) {
			return false;
		} else {
			int typec = _ToSegmentIntersectCase(type);
			if (typec & typeinclude) {
				return true;
			} else {
				return false;
			}
		}
	}
};

inline static std::string ToString_SegmentIntersectType(
		const int& type) {
	std::stringstream stream;
	if (type == NO_INTERSECT) {
		stream << "NO_INTERSECT ";
		return stream.str();
	}
	if (type & INTERSECT) {
		stream << "INTERSECT ";
	}
	if (type & OVERLAP) {
		stream << "OVERLAP ";
	}
	if (type & START_1) {
		stream << "START_1 ";
	}
	if (type & START_2) {
		stream << "START_2 ";
	}
	if (type & END_1) {
		stream << "END_1 ";
	}
	if (type & END_2) {
		stream << "END_2 ";
	}
	return stream.str();
}
}
#endif
