#ifndef _INTERSECTION_HPP_
#define _INTERSECTION_HPP_

#include <cmath>
#include "../geometry_define.hpp"
#include "../objects/_objects.hpp"
#include "_operation.hpp"
#include "../../utility/any.hpp"
#include "_tri_tri_intersect_02.hpp"

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

template<typename BOX>
class BBTree_;

template<typename TYPE, St DIM>
class Segment_;

template<typename TYPE, St DIM>
class PointChain_;

template<typename TYPE, St DIM>
class Operation_;

template<class TYPE, St DIM>
class TriSurface_;

template<typename TYPE, St DIM>
class Intersection_ {
public:
	static const St Dim = DIM;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Point;
	typedef Box_<TYPE, DIM> Box;
	typedef BBox_<TYPE, DIM> BBox;
	typedef BBTree_<BBox> BBTree;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Segment& ref_Segment;
	typedef const Segment& const_ref_Segment;
	typedef PointChain_<TYPE, DIM> PointChain;

	typedef TriSurface_<TYPE, DIM> TriSurface;
	typedef typename TriSurface::Fac TriFace;

	typedef typename TriFace::Tag TagTriFace;
	typedef typename TriSurface::Edg Edge;
	typedef typename TriSurface::Ver Vertex;

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

	template<class CA, class CB, class RES, typename ARG = int>
	static bool Find(const CA& a, const CB& b, RES& res, const ARG& arg = 0) {
		_Find(a, b, res, arg, typename CA::Tag(), typename CB::Tag());
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
			const TriFace& a, const TriFace& b,
			TagTriFace, TagTriFace) {
		const Point& v0 = *(a.vertex(0));
		const Point& v1 = *(a.vertex(1));
		const Point& v2 = *(a.vertex(2));
		const Point& u0 = *(b.vertex(0));
		const Point& u1 = *(b.vertex(1));
		const Point& u2 = *(b.vertex(2));

		typename Point::Vt V0[] = { v0[0], v0[1], v0[2] };
		typename Point::Vt V1[] = { v1[0], v1[1], v1[2] };
		typename Point::Vt V2[] = { v2[0], v2[1], v2[2] };
		typename Point::Vt U0[] = { u0[0], u0[1], u0[2] };
		typename Point::Vt U1[] = { u1[0], u1[1], u1[2] };
		typename Point::Vt U2[] = { u2[0], u2[1], u2[2] };
		return TriTriIsect02::tri_tri_overlap_test_3d(V0, V1, V2, U0, U1, U2);
	}
	static bool _Check(
			const TriFace* a, const TriFace* b,
			TagTriFace ta, TagTriFace tb) {
		return _Check(*a, *b, ta, tb);
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
	template<class ARG>
	static bool _Check(
			const Segment& b1, const Segment& b2, const ARG& arg,
			typename Segment::Tag, typename Segment::Tag) {
		return Check_asSegment(
				b1.ps(), b1.pe(),
				b2.ps(), b2.ps(),
				arg);
	}

	template<class ARG>
	static bool _Check(
			const BBox& b1, const BBox& b2, const ARG& arg,
			typename BBox::Tag, typename BBox::Tag) {
		bool resb = Check(b1, b2);
		_RETURN_VAL_IF_FAIL(resb, false);
		/// box is intersect, then check the objects
		const Any& aobj1 = b1.get_obj();
		const Any& aobj2 = b2.get_obj();
		if (aobj1.type() == typeid(Segment)
				&& aobj2.type() == typeid(Segment)) {
			return Check(any_cast<Segment>(aobj1), any_cast<Segment>(aobj2),
					arg);
		}
		if (aobj1.type() == typeid(TriFace)
				&& aobj2.type() == typeid(TriFace)) {
			return Check(any_cast<TriFace>(aobj1), any_cast<TriFace>(aobj2));
		}
		if (aobj1.type() == typeid(TriFace*)
				&& aobj2.type() == typeid(TriFace*)) {
			TriFace* a = any_cast<TriFace*>(aobj1);
			TriFace* b = any_cast<TriFace*>(aobj2);
			return Check(*a, *b);
		}

		SHOULD_NOT_REACH;
		return false;
	}

	static bool _Check(
			const BBTree& tree, const BBox& box,
			typename BBTree::Tag, typename BBox::Tag) {
		// if one leaf box a overlap with b, return true
		bool res = false;
		typename BBTree::Func_flag fun = [&box, &res](
				typename BBTree::const_pNode pn, bool& flag) {
			if(res == true) {
				return;
			}
			if (!(Check(pn->box(), box))) {
				flag = false;
			} else {
				flag = true;
			}
			if(pn->is_leaf() && Check(pn->box(), box)) {
				res = true;
			}
		};
		tree.PreOrder_flag(fun);
		return res;
	}
	template<class ARG>
	static bool _Check(
			const BBTree& tree, const BBox& box, ARG arg,
			typename BBTree::Tag, typename BBox::Tag) {
		// if one leaf box a overlap with b, return true
		bool res = false;
		typename BBTree::Func_flag fun = [&box, &res, &arg](
				typename BBTree::const_pNode pn, bool& flag) {
			if(res == true) {
				return;
			}
			if (!(Check(pn->box(), box))) {
				flag = false;
			} else {
				flag = true;
			}
			if(pn->is_leaf() && Check(pn->box(), box, arg)) {
				res = true;
			}
		};
		tree.PreOrder_flag(fun);
		return res;
	}

	static bool _Check(
			const BBTree& ta, const BBTree& tb,
			typename BBTree::Tag, typename BBTree::Tag) {
		// if one leaf box a overlap with b, return true
		if (!Check(ta.root_box(), tb.root_box())) {
			return false;
		}
		bool res = false;
		for (typename BBTree::const_iterator iter = tb.begin();
				iter != tb.end();
				++iter) {
			if (iter->is_leaf()) {
				const BBox& box = (iter->box());
				res = Check(ta, box);
				if (res == true) {
					return res;
				}
			}
		}
		return res;
	}

	template<class ARG>
	static bool _Check(
			const BBTree& ta, const BBTree& tb, const ARG& arg,
			typename BBTree::Tag, typename BBTree::Tag) {
		// if one leaf box a overlap with b, return true
		if (!Check(ta.root_box(), tb.root_box())) {
			return false;
		}
		bool res = false;
		for (typename BBTree::const_iterator iter = tb.begin();
				iter != tb.end();
				++iter) {
			if (iter->is_leaf()) {
				const BBox& box = (iter->box());
				res = Check(ta, box, arg);
				if (res == true) {
					return res;
				}
			}
		}
		return res;
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

	// --------------------Find -----------------------
	template<class CA, class CB, typename ARG>
	static bool _Find(const CA& a, const CB& b,
			std::tuple<bool, bool, Point, Point>& res, const ARG& arg,
			typename TriFace::Tag, typename TriFace::Tag) {
		ASSERT(DIM == 3);
		const Point& v0 = *(a.vertex(0));
		const Point& v1 = *(a.vertex(1));
		const Point& v2 = *(a.vertex(2));
		const Point& u0 = *(b.vertex(0));
		const Point& u1 = *(b.vertex(1));
		const Point& u2 = *(b.vertex(2));

		typename Point::Vt V0[] = { v0[0], v0[1], v0[2] };
		typename Point::Vt V1[] = { v1[0], v1[1], v1[2] };
		typename Point::Vt V2[] = { v2[0], v2[1], v2[2] };
		typename Point::Vt U0[] = { u0[0], u0[1], u0[2] };
		typename Point::Vt U1[] = { u1[0], u1[1], u1[2] };
		typename Point::Vt U2[] = { u2[0], u2[1], u2[2] };
		bool coplane;
		typename Point::Vt start[3], end[3];
		std::get<0>(res) =
				TriTriIsect02::tri_tri_intersection_test_3d(
						V0, V1, V2, U0, U1, U2,
						coplane, start, end);
		std::get<1>(res) = coplane;
		Point& ps = std::get<2>(res);
		Point& pe = std::get<3>(res);
		for (int i = 0; i < 3; i++) {
			ps[i] = start[i];
			pe[i] = end[i];
		}
		return std::get<0>(res);
	}

	template<class CA, class CB, class RES, typename ARG>
	static bool _Find(const CA& b1, const CB& b2,
			RES& res, const ARG& arg,
			typename BBox::Tag, typename BBox::Tag) {
		ASSERT(DIM == 3);
		bool resb = Check(b1, b2);
		_RETURN_VAL_IF_FAIL(resb, false);
		/// box is intersect, then check the objects
		const Any& aobj1 = b1.get_obj();
		const Any& aobj2 = b2.get_obj();
		if (aobj1.type() == typeid(TriFace)
				&& aobj2.type() == typeid(TriFace)) {
			return Find(any_cast<TriFace>(aobj1), any_cast<TriFace>(aobj2),
					res, arg);
		}
		if (aobj1.type() == typeid(TriFace*)
				&& aobj2.type() == typeid(TriFace*)) {
			TriFace* a = any_cast<TriFace*>(aobj1);
			TriFace* b = any_cast<TriFace*>(aobj2);
			return Find(*a, *b, res, arg);
		}
		SHOULD_NOT_REACH;
		return false;
	}

	template<class CA, class CB, class RES, typename ARG>
	static bool _Find(const CA& tree, const CB& box,
			RES& res, const ARG& arg,
			typename BBTree::Tag, typename BBox::Tag) {
		ASSERT(DIM == 3);
		bool check_res = Check(tree, box);
		_RETURN_VAL_IF_FAIL(check_res, false);

		/// RES is container
		typedef typename RES::value_type ResOne;
		//ResOne resone;
		/// box is intersect, then check the objects
		typename BBTree::Func_flag fun = [&box, &arg, &res](
				typename BBTree::const_pNode pn, bool& flag) {
			if (!(Check(pn->box(), box))) {
				flag = false;
			} else {
				flag = true;
			}
			if(pn->is_leaf() && Check(pn->box(), box)) {
				ResOne resone;
				Find(pn->box(), box, resone, arg);
				res.push_back(resone);
			}
		};
		tree.PreOrder_flag(fun);

		return false;
	}

}
;

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
