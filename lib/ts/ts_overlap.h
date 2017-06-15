#ifndef TS_OVERLAP_H_
#define TS_OVERLAP_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_triangle.h"
#include "ts_edge.h"
#include "ts_AABBox.h"
#include "ts_BBTree.h"
#include "ts_tt.h"

namespace TS {

template<class TYPE, st DIM>
class Overlap {
public:
	typedef TYPE vt;
	typedef st size_type;
	typedef Overlap<TYPE, DIM> Self;
	typedef AABBox<TYPE, DIM> Box;
	typedef const AABBox<TYPE, DIM> const_Box;
	typedef Point<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef std::shared_ptr<Poi> spPoi;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Tri* pTri;
	typedef std::shared_ptr<Tri> spTri;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Ver* pVer;
	typedef std::shared_ptr<Ver> spVer;
	typedef Segment<TYPE, DIM> Seg;
	typedef Seg* pSeg;
	typedef std::shared_ptr<Seg> spSeg;

	typedef BBNode<Box> Node;
	typedef const BBNode<Box> const_Node;
	typedef Node* pNode;
	typedef const Node* const_pNode;
	typedef BBTree<Box> Tree;
	typedef const BBTree<Box> const_Tree;

	static const st Dim = DIM;
public:
	static bool DetectBox(const Box& a, const Box& b) {
		return a.are_overlapping(b);
	}

	static bool DetectBox(const Node& a, const Node& b) {
		return Self::DetectBox(a.box, b.box);
	}

	static bool DetectBox(const Tree& a, const Box& b) {
		// if one leaf box a overlap with b, return true
		bool res = false;
		typename Tree::Func_flag fun = [&b, &res](const_pNode pn, bool& flag) {
			if(res == true) {
				return;
			}
			if (!pn->box().are_overlapping(b)) {
				flag = false;
			} else {
				flag = true;
			}
			if(pn->is_leaf() && pn->box().are_overlapping(b)) {
				res = true;
			}
		};
		a.PreOrder_flag(fun);
		return res;
	}

	static bool DetectBox(const Tree& ta, const Tree& tb) {
		// if one leaf box a overlap with b, return true
		bool res = false;
		for (typename Tree::const_iterator iter = tb.begin();
				iter != tb.end() && res == false; ++iter) {
			if (iter->is_leaf()) {
				const Box& b = (iter->box());
				typename Tree::Func_flag fun =
						[&b, &res](const_pNode pn, bool& flag) {
							if(res == true) {
								return;
							}
							if (!pn->box().are_overlapping(b)) {
								flag = false;
							} else {
								flag = true;
							}
							if(pn->is_leaf() && pn->box().are_overlapping(b)) {
								res = true;
							}
						};
				ta.PreOrder_flag(fun);
			}
		}
		return res;
	}

	static bool DetectObj(const Box& a, const Box& b) {
		if (!a.are_overlapping(b)) {
			return false;
		}
		const utPointer uta = a.utp_obj();
		OBJ_TYPE ta = a.type;
		const utPointer utb = b.utp_obj();
		OBJ_TYPE tb = b.type;
		return _Detect(uta, utb, ta, tb);
	}

	static bool DetectObj(const Tree& ta, const Tree& tb) {
		// if one leaf box a overlap with b, return true
		bool res = false;
		for (typename Tree::const_iterator iter = tb.begin();
				iter != tb.end() && res == false; ++iter) {
			if (iter->is_leaf()) {
				const Box& b = (iter->box());
				typename Tree::Func_flag fun =
						[&b, &res](const_pNode pn, bool& flag) {
							if(res == true) {
								return;
							}
							if (!pn->box().are_overlapping(b)) {
								flag = false;
							} else {
								flag = true;
							}
							if(pn->is_leaf()) {
								res = DetectObj(pn->box(), b);
							}
						};
				ta.PreOrder_flag(fun);
			}
		}
		return res;
	}

	static bool DetectIsect(const Box& a, const Box& b, Poi& pa, Poi& pb) {
		if (!a.are_overlapping(b)) {
			return false;
		}
		const utPointer uta = a.utp_obj();
		OBJ_TYPE ta = a.type;
		const utPointer utb = b.utp_obj();
		OBJ_TYPE tb = b.type;
		//Poi pa, pb;
		return _Detect(uta, utb, ta, tb, pa, pb);
	}

	static bool DetectObj(const Tree& ta, const Tree& tb, std::vector<Box>& av, //
			std::vector<Box>& bv //
			) {
		// if one leaf box a overlap with b, return true
		av.clear();
		bv.clear();
		for (typename Tree::const_iterator iter = tb.begin();
				iter != tb.end(); ++iter) {
			if (iter->is_leaf()) {
				const Box& b = (iter->box());
				typename Tree::Func_flag fun =
						[&b, &av, &bv](const_pNode pn, bool& flag) {
							if (!pn->box().are_overlapping(b)) {
								flag = false;
							} else {
								flag = true;
							}
							if(pn->is_leaf()) {
								bool res = DetectObj(pn->box(), b);
								if(res) {
									av.push_back(pn->box());
									bv.push_back(b);
								}
							}
						};
				ta.PreOrder_flag(fun);
			}
		}
		return av.size() > 0;
	}

	static bool Detect(pTri a, pTri b) {
		int res = TriTriIsect(*a, *b);
		return res == 1 ? true : false;
	}

	static bool _Detect(utPointer a, utPointer b, OBJ_TYPE ta, OBJ_TYPE tb) {
		if ((ta == TRIANGLE && tb == TRIANGLE) || (ta == FACE && tb == FACE)) {
			pTri ra = CAST(pTri, a);
			pTri rb = CAST(pTri, b);
			return Detect(ra, rb);
		}
		ASSERT(false);
		return false; //make complier happy
	}

	static bool _Detect(utPointer a, utPointer b, OBJ_TYPE ta, OBJ_TYPE tb,
			Poi& pa, Poi& pb) {
		if ((ta == TRIANGLE && tb == TRIANGLE) || (ta == FACE && tb == FACE)) {
			pTri ra = CAST(pTri, a);
			pTri rb = CAST(pTri, b);
			return Detect(ra, rb, pa, pb);
		}
		ASSERT(false);
		return false; //make complier happy
	}

	static bool Detect(pTri a, pTri b, Poi& ps, Poi& pe) {
		bool co = false;
		int res = TriTriIsect(*a, *b, co, ps, pe);
		bool onep = ((ps) == (pe));
		return (res == 1 && co == false && onep == false) ? true : false;
	}

};

}

#endif
