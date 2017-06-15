#ifndef _INTERPOLATE_H_
#define _INTERPOLATE_H_


#include <functional>
#include <math.h>

#include "../carpio_define.hpp"
#include "domain/domain.hpp"
#include "expression.hpp"


namespace carpio {
/*
 *  static class
 *  template all the function
 */

template<typename COO_VALUE, typename VALUE, int DIM>
class Expression_;

template<typename COO_VALUE, typename VALUE, int DIM>
class Interpolate_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Domain_<COO_VALUE, VALUE, DIM> Domain;
	typedef Domain_<COO_VALUE, VALUE, DIM>& ref_Domain;
	typedef const Domain_<COO_VALUE, VALUE, DIM>& const_ref_Domain;
	typedef Domain_<COO_VALUE, VALUE, DIM>* pDomain;
	typedef const Domain_<COO_VALUE, VALUE, DIM>* const_pDomain;

	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef const Grid_<COO_VALUE, VALUE, DIM> const_Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM>& ref_Grid;
	typedef const Grid_<COO_VALUE, VALUE, DIM>& const_ref_Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pGrid;

	typedef Ghost_<COO_VALUE, VALUE, DIM> Ghost;
	typedef const Ghost_<COO_VALUE, VALUE, DIM> const_Ghost;
	typedef Ghost_<COO_VALUE, VALUE, DIM>& ref_Ghost;
	typedef const Ghost_<COO_VALUE, VALUE, DIM>& ref_const_Ghost;
	typedef Ghost_<COO_VALUE, VALUE, DIM>* pGhost;
	typedef const Ghost_<COO_VALUE, VALUE, DIM>* const_pGhost;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef PData_<cvt, vt, Dim> PData;
	typedef PData_<cvt, vt, Dim>& ref_PData;
	typedef const PData_<cvt, vt, Dim>& const_ref_PData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef const Node_<COO_VALUE, VALUE, DIM> *const_pNode;

	typedef Stencil_<cvt, vt, Dim, 1> Stencil_1;
	typedef Stencil_<cvt, vt, Dim, 1>& ref_Stencil_1;
	typedef const Stencil_<cvt, vt, Dim, 1>& const_ref_Stencil_1;
	typedef Stencil_<cvt, vt, Dim, 2> Stencil_2;
	typedef Stencil_<cvt, vt, Dim, 2>& ref_Stencil_2;
	typedef const Stencil_<cvt, vt, Dim, 2>& const_ref_Stencil_2;

	typedef Expression_<COO_VALUE, VALUE, DIM> Exp;
	typedef Expression_<COO_VALUE, VALUE, DIM>* pExp;
	typedef Expression_<COO_VALUE, VALUE, DIM>& ref_Exp;
	//typedef typename Exp::Term Term;
	typedef typename Exp::iterator ExpIter;
	typedef typename Exp::const_iterator const_ExpIter;

	typedef Point_<COO_VALUE, DIM> Point;

	typedef Vector_<COO_VALUE, DIM> Vector;

	class PExp: public std::pair<Point, Exp> {
	public:
		PExp() {

		}
		PExp(const Point& p, const Exp& e) {
			this->first = p;
			this->second = e;
		}
		PExp(const_pNode pn, const Exp& e) {
			Point p(pn->cp(_X_), pn->cp(_Y_), pn->cp(_Z_));
			this->first = p;
			this->second = e;
		}
		void set(const Point& p, const Exp& e) {
			this->first = p;
			this->second = e;
		}
		const Point& P() const {
			return this->first;
		}
		Point& P() {
			return this->first;
		}
		const Exp& E() const {
			return this->second;
		}
		Exp& E() {
			return this->second;
		}
	};
	/*
	 * static function
	 * ----------------------------------------------------
	 */

	/*
	 *  There is only one node in the stencil
	 *  The point (x,y) is in that node, and
	 *  the value is equal to the center point.
	 *  Stencil:  ---X---C---X---->
	 */
	static int _1Node(vt& res, const St& idx, const Stencil_1& stc, cvt x,
			cvt y, cvt z = 0) {
		ASSERT(stc.non_null_nodes() == 1);
		// the non null nodes must be center node
		const_pNode pc = stc.center_pnode();
		ASSERT(pc->is_in_on(x, y, z));
		// ---
		res = pc->cd(idx);
		return _SUCCESS;
	}

	static int _1Node(ref_PData res, const_ref_Stencil_1 stc) {
		//assert
		ASSERT(stc.non_null_nodes() == 1);
		typename Stencil_2D1::const_pNode pc = stc.center_pnode();
		ASSERT(pc->is_in_on(res.x(), res.y(), res.z()));
		// ---
		_AverangeValueFromCenterLeaf(res, pc);
		return _SUCCESS;
	}

	static Exp _1Node(ref_Stencil_1 stc) {
		//assert
		ASSERT(stc.non_null_nodes() == 1);
		pNode pc = stc.center_pnode();
		//ASSERT(pc->is_in_on(res.x(), res.y(), res.z()));
		// ---
		PExp res = _AverangeExpFromCenterLeaf(pc);
		return res.E();
	}
	static Exp _1Node(ref_Stencil_2 stc) {
		//assert
		ASSERT(stc.non_null_nodes() == 1);
		pNode pc = stc.center_pnode();
		//ASSERT(pc->is_in_on(res.x(), res.y(), res.z()));
		// ---
		PExp res = _AverangeExpFromCenterLeaf(pc);
		return res.E();
	}

	/*
	 *  There are two nodes in the stencil, x is
	 *  the location, somewhere on the stencil axes.
	 *
	 *  Stencil: ---X---C---O---->  or
	 *           ---O---C---X---->
	 *                ^
	 *                x
	 */
	static int _2NodeOnAxes(ref_PData res, ref_Stencil_1 stc) {
		// assert
		ASSERT(stc.non_null_nodes() == 2);
		pNode pc = stc.center_pnode();
		PData_2D pdc(res.arr_idx(), pc->cp(_X_), pc->cp(_Y_));
		_AverangeValueFromCenterLeaf(pdc, pc);

		// another pnode is close to center
		pNode pn = stc.forward_pnode(1, stc.axes());
		if (pn == nullptr) {
			pn = stc.backward_pnode(1, stc.axes());
		}
		PData pdn(res.arr_idx(), pn->cp(_X_), pn->cp(_Y_));
		switch (__2NodeRelation(pc, pn)) {
		case _Equal_: {
			_AverangeValueFromCenterLeaf(pdn, pn);
			break;
		}
		case _FineCoarse_: {
			//unfinish ==================
			SHOULD_NOT_REACH;
			break;
		}
		case _CoarseFine_: {
			_AverangeValueFromCenterLeaf(pdn, pn);
			break;
		}
		default: {
			SHOULD_NOT_REACH;
			break;
		}
		}

		// ---
		for (St i = 0; i < res.size(); ++i) {
			cvt x = res.p(stc.axes(0));
			cvt x1 = pdc.p(stc.axes(0));
			vt y1 = pdc.val(i);
			cvt x2 = pdn.p(stc.axes(0));
			vt y2 = pdn.val(i);
			res.flag(i) = PData::Flag_Center;
			res.val(i) = __LinearInterpolation(x, x1, y1, x2, y2);
		}

		return _SUCCESS;
	}

	/*
	 *  There are two nodes in the stencil, x is
	 *  the location, somewhere on the stencil axes.
	 *
	 *  Stencil: ---X---C---O---->  or
	 *           ---O---C---X---->
	 *                ^
	 *                x
	 */
	static Exp _2NodeOnAxes(cvt x, //
			ref_Stencil_1 stc,  //
			const Direction& dir, //previous direction
			const int& order // previous order
			) {
		// assert
		ASSERT(stc.non_null_nodes() == 2);
		Axes axes = FaceDirectionToAxes(dir);
		pNode pc = stc.center_pnode();
		PExp&& cexp = _AverangeExpFromCenterLeaf(pc);
		// another pnode is close to center
		PExp nexp;
		pNode pn = stc.forward_pnode(1, stc.axes());
		if (pn == nullptr) {
			pn = stc.backward_pnode(1, stc.axes());
			nexp = _2NodeRelation(pc, pn, _M_, axes, order);
		} else {
			nexp = _2NodeRelation(pc, pn, _P_, axes, order);
		}
		// ---
		cvt x1 = pc->cp(stc.axes(0));
		cvt x2 = nexp.P().val(stc.axes(0));
		return __LinearInterpolation(x, x1, cexp.E(), x2, nexp.E());
	}

	/*
	 *  There are two nodes in the stencil, x is
	 *  the location, somewhere on the stencil axes.
	 *
	 *  Stencil: ---C---O---O---->  or
	 *           ---O---C---O---->  or
	 *           ---O---O---C---->
	 *                ^
	 *                x
	 */
	static int _3NodeOnAxes(ref_PData res, ref_Stencil_1 stc) {
		// assert
		ASSERT(stc.non_null_nodes() == 3);
		const_pNode pc = stc.center_pnode();
		PData pdc(res.arr_idx(), pc->cp(_X_), pc->cp(_Y_));
		_AverangeValueFromCenterLeaf(pdc, pc);
		// forward pnode is close to center
		const_pNode pf = stc.forward_pnode(1, stc.axes());
		ASSERT(pf != nullptr);
		PData pdf(res.arr_idx(), pf->cp(_X_), pf->cp(_Y_));
		switch (__2NodeRelation(pc, pf)) {
		case _Equal_: {
			_AverangeValueFromCenterLeaf(pdf, pf);
			break;
		}
		case _FineCoarse_: {
			//unfinish ==================
			SHOULD_NOT_REACH;
			break;
		}
		case _CoarseFine_: {
			_AverangeValueFromCenterLeaf(pdf, pf);
			break;
		}
		default: {
			SHOULD_NOT_REACH;
			break;
		}
		}
		// forward pnode is close to center
		const_pNode pb = stc.backward_pnode(1, stc.axes());
		ASSERT(pb != nullptr);
		PData pdb(res.arr_idx(), pb->cp(_X_), pb->cp(_Y_));
		switch (__2NodeRelation(pc, pb)) {
		case _Equal_: {
			_AverangeValueFromCenterLeaf(pdb, pb);
			break;
		}
		case _FineCoarse_: {
			//unfinish ==================
			SHOULD_NOT_REACH;
			break;
		}
		case _CoarseFine_: {
			_AverangeValueFromCenterLeaf(pdb, pb);
			break;
		}
		}
		// ---
		for (St i = 0; i < res.size(); ++i) {
			cvt x = res.p(stc.axes(0));
			cvt x1 = pdf.p(stc.axes(0));
			vt y1 = pdf.val(i);
			cvt x2 = pdc.p(stc.axes(0));
			vt y2 = pdc.val(i);
			cvt x3 = pdb.p(stc.axes(0));
			vt y3 = pdb.val(i);
			res.flag(i) = PData::Flag_Center;
			res.val(i) = __2OrderInterpolation(x, x1, y1, x2, y2, x3, y3);
		}
		return _SUCCESS;
	}

	/*
	 *  There are two nodes in the stencil, x is
	 *  the location, somewhere on the stencil axes.
	 *
	 *  Stencil: ---C---O---O---->  or
	 *           ---O---C---O---->  or
	 *           ---O---O---C---->
	 *                ^
	 *                x
	 */
	static Exp _3NodeOnAxes(cvt x, ref_Stencil_1 stc, //
			const Direction& dir, //previous direction
			const int& order      // previous order
			) {
		// assert
		ASSERT(stc.non_null_nodes() == 3);
		Axes axes = FaceDirectionToAxes(dir);
		pNode pc = stc.center_pnode();
		PExp expc = _AverangeExpFromCenterLeaf(pc);
		// forward pnode is close to center

		pNode pf = stc.forward_pnode(1, stc.axes());
		ASSERT(pf != nullptr);
		PExp expf = _2NodeRelation(pc, pf, _P_, axes, order);
		// backward pnode is close to center
		pNode pb = stc.backward_pnode(1, stc.axes());
		ASSERT(pb != nullptr);
		PExp expb = _2NodeRelation(pc, pb, _M_, axes, order);
		// ---
		cvt x1 = expf.P().val(stc.axes(0));
		cvt x2 = expc.P().val(stc.axes(0));
		cvt x3 = expb.P().val(stc.axes(0));
		return __2OrderInterpolation(x, x1, expf.E(), x2, expc.E(), x3,
				expb.E());

	}

	static void OnFace_1Order( // 2D QuadTree Node
			pNode pn,                      //node
			const Direction& dir,                //face
			const ArrayListV<St>& arridx,            //data index
			ArrayListV<vt>& arrres                //data res
			) {
		// assert
		ASSERT(pn != nullptr);
		ASSERT(IsFaceDirection(dir));
		// get PData
		cvt x = pn->p(dir, _X_);
		cvt y = pn->p(dir, _Y_);
		cvt z = pn->p(dir, _Z_);
		PData pdata(arridx, x, y, z);
		// stencil
		Orientation ori;
		Axes axe;
		FaceDirectionToOrientationAndAxes(dir, ori, axe);
		Stencil_1 sten(pn, axe, (IsP(ori) ? 1 : 0), (IsM(ori) ? 1 : 0));
		if (sten.non_null_nodes() == 1) {
			_1Node(pdata, sten);
		} else {
			_2NodeOnAxes(pdata, sten);
		}
		arrres = pdata.arr_val();
	}

	static Exp OnFace_1Order(  //
			pNode pn,  //
			const Direction& dir  //
			) {
		// assert
		ASSERT(pn != nullptr);
		ASSERT(IsFaceDirection(dir));
		//
		Orientation ori;
		Axes axe;
		FaceDirectionToOrientationAndAxes(dir, ori, axe);
		Stencil_1 sten(pn, axe, (IsP(ori) ? 1 : 0), (IsM(ori) ? 1 : 0));
		if (sten.non_null_nodes() == 1) {
			return _1Node(sten);
		} else {
			return _2NodeOnAxes(pn->p(ori, axe), sten, dir, 1);
		}
	}

	static Exp OnAxes_1Order(  //
			pNode pn,      //
			const Direction& dir,  //
			const cvt& dis  // relative distance to the center node
			) {
		// assert
		ASSERT(pn != nullptr);
		ASSERT(IsFaceDirection(dir));
		//
		Orientation ori;
		Axes axe;
		FaceDirectionToOrientationAndAxes(dir, ori, axe);
		Stencil_1 sten(pn, axe, (IsP(ori) ? 1 : 0), (IsM(ori) ? 1 : 0));
		if (sten.non_null_nodes() == 1) {
			return _1Node(sten);
		} else {
			return _2NodeOnAxes(pn->cp(axe) + dis, sten, dir, 1);
		}
	}

	static Exp OnCorner_1Order(  //
			pNode pn,      //
			const Direction& dir  //
			) {
		// assert
		ASSERT(pn != nullptr);
		ASSERT(IsCornerDirection(dir));
		//
		Orientation o1, o2;
		Axes a1, a2;
		CornerDirectionToOrientationAndAxes(dir, o1, a1, o2, a2);
		Stencil_2 sten(pn, a1, (IsP(o1) ? 1 : 0), (IsM(o1) ? 1 : 0), a2,
				(IsP(o2) ? 1 : 0), (IsM(o2) ? 1 : 0));
		St non_null = sten.non_null_nodes();
		switch (non_null) {
		case 1: {
			return _1Node(sten);
			break;
		}
		case 2: {
			return _2Node_Corner(sten, o1, a1, o2, a2);
			break;
		}
		case 3: {
			return _3Node_Corner(sten, o1, a1, o2, a2);
			break;

		}
		case 4: {
			return _4Node_Corner(sten, o1, a1, o2, a2);
			break;
		}
		}
		SHOULD_NOT_REACH;
		return Exp();
	}

	static Exp _2Node_Corner(Stencil_2 stc,  //
			const Orientation&o1, const Axes&a1,  //
			const Orientation&o2, const Axes& a2) {
		// assert
		ASSERT(stc.non_null_nodes() == 2);
		pNode pc = stc.center_pnode();
		Direction dir_c = ToCornerDirection(o1, a1, o2, a2);
		pNode pcc = Node::GetChild_CornerLeaf(pc, dir_c);
		if (pcc == nullptr) {
			pcc = pc;
		}
		PExp cexp = _AverangeExpFromCenterLeaf(pcc);
		// another pnode is close to center
		PExp nexp;
		pNode pn = stc.get_pnode(IsP(o1) ? 1 : -1, 0);
		pn = stc.get_pnode(IsP(o1) ? 1 : -1, 0);
		Direction dir_n = ToCornerDirection(Opposite(o1), a1, o2, a2);
		if (pn == nullptr) {
			pn = stc.get_pnode(0, IsP(o2) ? 1 : -1);
			dir_n = ToCornerDirection(o1, a1, Opposite(o2), a2);
		}
		if (pn == nullptr) {
			pn = stc.get_pnode(IsP(o1) ? 1 : -1, IsP(o2) ? 1 : -1);
			dir_n = ToCornerDirection(Opposite(o1), a1, Opposite(o2), a2);
		}
		ASSERT(pn != nullptr);
		pNode pnc = Node::GetChild_CornerLeaf(pn, dir_n);
		if (pnc == nullptr) {
			pnc = pn;
		}
		// ---
		cvt xc = pcc->p(dir_c, a1);
		cvt yc = pcc->p(dir_c, a2);
		cvt x0 = pcc->cp(a1);
		cvt y0 = pcc->cp(a2);
		cvt x1 = pnc->cp(a1);
		cvt y1 = pnc->cp(a2);
		return __LinearInterpolation(xc, yc, x0, y0, cexp.E(), x1, y1, nexp.E());
	}

	static Exp _3Node_Corner(Stencil_2 stc,  //
			const Orientation&o1, const Axes&a1,  //
			const Orientation&o2, const Axes& a2) {
		/*
		 *    pc2   pn    dir_c2  dir_n
		 *    pc    pc1   dir_c   dir_c1
		 */
		// assert
		ASSERT(stc.non_null_nodes() == 3);
		// center node
		pNode pc = stc.center_pnode();
		Direction dir_c = ToCornerDirection(o1, a1, o2, a2);
		pNode pcc = Node::GetChild_CornerLeaf(pc, dir_c);
		if (pcc == nullptr) {
			pcc = pc;
		}
		PExp pexp_c = _AverangeExpFromCenterLeaf(pcc);
		// pn
		pNode pn = stc.get_pnode(IsP(o1) ? 1 : -1, IsP(o2) ? 1 : -1);
		Direction dir_n = ToCornerDirection(Opposite(o1), a1, Opposite(o2), a2);
		pNode pnc = nullptr;
		if (pn != nullptr) {
			pnc = Node::GetChild_CornerLeaf(pn, dir_n);
		}
		if (pnc == nullptr) {
			pnc = pn;
		}
		if (pnc != nullptr) {
			PExp pexp_n = _AverangeExpFromCenterLeaf(pnc);
			cvt xc = pcc->p(dir_c, a1);
			cvt yc = pcc->p(dir_c, a2);
			cvt x0 = pcc->cp(a1);
			cvt y0 = pcc->cp(a2);
			cvt x1 = pnc->cp(a1);
			cvt y1 = pnc->cp(a2);
			return __LinearInterpolation(xc, yc, x0, y0, pexp_c.E(), x1, y1,
					pexp_n.E());
		}
		// ------------------
		// pc1
		pNode pc1 = stc.get_pnode(IsP(o1) ? 1 : -1, 0);
		Direction dir_c1 = ToCornerDirection(Opposite(o1), a1, o2, a2);
		pNode pc1c = Node::GetChild_CornerLeaf(pc1, dir_c1);
		if (pc1c == nullptr) {
			pc1c = pc1;
		}
		PExp pexp_c1 = _AverangeExpFromCenterLeaf(pc1c);
		// pc2
		pNode pc2 = stc.get_pnode(0, IsP(o2) ? 1 : -1);
		Direction dir_c2 = ToCornerDirection(o1, a1, Opposite(o2), a2);
		pNode pc2c = Node::GetChild_CornerLeaf(pc2, dir_c2);
		if (pc2c == nullptr) {
			pc2c = pc2;
		}
		PExp pexp_c2 = _AverangeExpFromCenterLeaf(pc2c);
		cvt xc = pcc->p(dir_c, a1);
		cvt yc = pcc->p(dir_c, a2);
		cvt x0 = pc1c->cp(a1);
		cvt y0 = pc1c->cp(a2);
		cvt x1 = pc2c->cp(a1);
		cvt y1 = pc2c->cp(a2);
		return __LinearInterpolation(xc, yc, x0, y0, pexp_c1.E(), x1, y1,
				pexp_c2.E());
	}

	static Exp _4Node_Corner(Stencil_2 stc,  //
			const Orientation&o1, const Axes&a1,  //
			const Orientation&o2, const Axes& a2) {
		/*
		 *    pc2   pn    dir_c2  dir_n
		 *    pc    pc1   dir_c   dir_c1
		 */
		// assert
		ASSERT(stc.non_null_nodes() == 4);
		// center node
		pNode pc = stc.center_pnode();
		Direction dir_c = ToCornerDirection(o1, a1, o2, a2);
		pNode pcc = Node::GetChild_CornerLeaf(pc, dir_c);
		if (pcc == nullptr) {
			pcc = pc;
		}
		PExp pexp_c = _AverangeExpFromCenterLeaf(pcc);
		// pn
		pNode pn = stc.get_pnode(IsP(o1) ? 1 : -1, IsP(o2) ? 1 : -1);
		Direction dir_n = ToCornerDirection(Opposite(o1), a1, Opposite(o2), a2);
		pNode pnc = Node::GetChild_CornerLeaf(pn, dir_n);
		if (pnc == nullptr) {
			pnc = pn;
		}
		PExp pexp_n = _AverangeExpFromCenterLeaf(pnc);
		cvt xc = pcc->p(dir_c, a1);
		cvt yc = pcc->p(dir_c, a2);
		cvt x0 = pcc->cp(a1);
		cvt y0 = pcc->cp(a2);
		cvt x1 = pnc->cp(a1);
		cvt y1 = pnc->cp(a2);
		Exp pexp1 = __LinearInterpolation(xc, yc, x0, y0, pexp_c.E(), x1, y1,
				pexp_n.E());
		// ------------------
		// pc1
		pNode pc1 = stc.get_pnode((IsP(o1) ? 1 : -1), 0);
		Direction dir_c1 = ToCornerDirection(Opposite(o1), a1, o2, a2);
		pNode pc1c = Node::GetChild_CornerLeaf(pc1, dir_c1);
		if (pc1c == nullptr) {
			pc1c = pc1;
		}
		PExp pexp_c1 = _AverangeExpFromCenterLeaf(pc1c);
		// pc2
		pNode pc2 = stc.get_pnode(0, (IsP(o2) ? 1 : -1));
		Direction dir_c2 = ToCornerDirection(o1, a1, Opposite(o2), a2);
		pNode pc2c = Node::GetChild_CornerLeaf(pc2, dir_c2);
		if (pc2c == nullptr) {
			pc2c = pc2;
		}
		PExp pexp_c2 = _AverangeExpFromCenterLeaf(pc2c);
		x0 = pc1c->cp(a1);
		y0 = pc1c->cp(a2);
		x1 = pc2c->cp(a1);
		y1 = pc2c->cp(a2);
		Exp pexp = __LinearInterpolation(xc, yc, x0, y0, pexp_c1.E(), x1, y1,
				pexp_c2.E());
		// ------
		pexp.plus(pexp1);
		pexp.times(0.5);
		return pexp;
	}

	static void OnFace_1Order( // 2D QuadTree Node
			pNode pn,                      //node
			const Direction& dir,                //face
			const St& idx,            //data index
			vt& res                //data res
			) {
		ArrayListV<St> arridx(1);            //data index
		ArrayListV<Vt> arrres(1);
		arridx[0] = idx;
		OnFace_1Order(pn, dir, arridx, arrres);
		res = arrres[0];
	}

	static void OnFace_2Order( // 2D QuadTree Node
			pNode pn,                      //node
			const Direction& dir,                //face
			const ArrayListV<St>& arridx,            //data index
			ArrayListV<vt>& arrres                //data res
			) {
		// assert
		ASSERT(pn != nullptr);
		ASSERT(IsFaceDirection(dir));
		// get PData
		cvt x = pn->p(dir, _X_);
		cvt y = pn->p(dir, _Y_);
		PData pdata(arridx, x, y);
		// stencil
		Orientation ori;
		Axes axe;
		FaceDirectionToOrientationAndAxes(dir, ori, axe);
		Stencil_1 sten(pn, axe, 1, 1);
		if (sten.non_null_nodes() == 1) {
			_1Node(pdata, sten);
		} else if (sten.non_null_nodes() == 2) {
			_2NodeOnAxes(pdata, sten);
		} else {
			_3NodeOnAxes(pdata, sten);
		}
		arrres = pdata.arr_val();
	}

	static Exp OnFace_2Order(  //
			pNode pn,  //
			const Direction& dir  //
			) {
		// assert
		ASSERT(pn != nullptr);
		ASSERT(IsFaceDirection(dir));
		// stencil
		Orientation ori;
		Axes axe;
		FaceDirectionToOrientationAndAxes(dir, ori, axe);
		Stencil_1 sten(pn, axe, 1, 1);
		cvt x = pn->p(ori, axe);
		if (sten.non_null_nodes() == 1) {
			return _1Node(sten);
		} else if (sten.non_null_nodes() == 2) {
			return _2NodeOnAxes(x, sten, dir, 2);
		} else {
			return _3NodeOnAxes(x, sten, dir, 2);
		}
	}

	static Exp OnAxes_2Order(  //
			pNode pn,  //
			const Direction& dir,  //
			const cvt& dis) {
		// assert
		ASSERT(pn != nullptr);
		ASSERT(IsFaceDirection(dir));
		// stencil
		Orientation ori;
		Axes axe;
		FaceDirectionToOrientationAndAxes(dir, ori, axe);
		Stencil_1 sten(pn, axe, 1, 1);
		cvt x = pn->cp(axe);
		if (sten.non_null_nodes() == 1) {
			return _1Node(sten);
		} else if (sten.non_null_nodes() == 2) {
			return _2NodeOnAxes(x + dis, sten, dir, 2);
		} else {
			return _3NodeOnAxes(x + dis, sten, dir, 2);
		}
	}

	static Exp OnAxes(  //
			pNode pn,  //
			const Direction& dir,  //
			const cvt& dis, //
			const int& order) {
		switch (order) {
		case 1: {
			return OnAxes_1Order(pn, dir, dis);
		}
		case 2: {
			return OnAxes_2Order(pn, dir, dis);
		}
		default: {
			SHOULD_NOT_REACH;
			return Exp();
		}
		}
	}

	static Exp OnFace(  //
			pNode pn,  //
			const Direction& dir,  //
			const int& order) {
		switch (order) {
		case 1: {
			return OnFace_1Order(pn, dir);
		}
		case 2: {
			return OnFace_2Order(pn, dir);
		}
		default: {
			SHOULD_NOT_REACH;
			return Exp();
		}
		}
	}
	static void OnFace_2Order( // 2D QuadTree Node
			pNode_2D pn,                      //node
			const Direction& dir,             //face
			const St& idx,            //data index
			Vt& res                //data res
			) {
		ArrayListV<St> arridx(1);            //data index
		ArrayListV<vt> arrres(1);
		arridx[0] = idx;
		OnFace_2Order(pn, dir, arridx, arrres);
		res = arrres[0];
	}
	/*
	 * Interpolate on the random point in the Node
	 *
	 */
	static Exp OnPlane(   //
			pNode pn,     //
			Axes a1, Axes a2,  // the axes
			cvt c1, cvt c2,    // the coordinate
			int order)         //
			{

		switch (order) {
		case 1: {
			return OnPlane_1Order(pn, a1, a2, c1, c2);
		}
		case 2: {
			SHOULD_NOT_REACH;
			return Exp();
		}
		default: {
			SHOULD_NOT_REACH;
			return Exp();
		}
		}
	}
	static Exp OnPlane_1Order(         //
			pNode pn,     //
			Axes a1, Axes a2,  // the axes
			cvt c1, cvt c2    // the coordinate
			)// the coordinate
			 //
			{
		// assert the point is in the node
		ASSERT(a1 != a2);
		ArrayListV<cvt> point(3);
		point[St(a1)] = c1;
		point[St(a2)] = c2;
		ASSERT(pn->is_in_on(point[0], point[1], point[2]));

		// get two face direction
		Orientation o1 = ((c1 - pn->cp(a1)) >= 0) ? _P_ : _M_;
		Direction d1 = ToFaceDirection(o1, a1);
		Orientation o2 = ((c2 - pn->cp(a2)) >= 0) ? _P_ : _M_;
		Direction d2 = ToFaceDirection(o2, a2);
		// a corner direction
		Direction dc = d1 | d2;
		// get three expression
		Exp&& exp1 = OnFace_1Order(pn, d1);
		Exp&& exp2 = OnFace_1Order(pn, d2);
		Exp&& exp3 = OnCorner_1Order(pn, dc);
		Exp exp0;
		exp0.plus(1.0, pn, 1);
		// get coordinate
		cvt x1 = pn->cp(a1);
		cvt x2 = pn->p(o1,a1);
		cvt y1 = pn->cp(a2);
		cvt y2 = pn->p(o2,a2);

		return __BilinearInterpolation(c1, c2, x1, x2, y1, y2, exp0, exp1, exp3, exp2);
	}
	/*
	 * gradient
	 */
	static Exp OnFace_Gradient(pNode pn, const Direction& dir, int order) {
		switch (order) {
		case 1: {
			return OnFace_1Order_Gradient(pn, dir);
		}
		case 2: {
			return OnFace_2Order_Gradient(pn, dir);
		}
		default: {
			SHOULD_NOT_REACH;
			return Exp();
		}
		}
	}

	static Exp OnFace_1Order_Gradient(pNode pn,  //
			const Direction& dir) {
		// assert
		ASSERT(pn != nullptr);
		ASSERT(IsFaceDirection(dir));
		//
		Orientation ori;
		Axes axe;
		FaceDirectionToOrientationAndAxes(dir, ori, axe);
		Stencil_1 sten(pn, axe, (IsP(ori) ? 1 : 0), (IsM(ori) ? 1 : 0));
		if (sten.non_null_nodes() == 1) {
			return Exp(); // 1 node gradient
		} else {
			return _2NodeGradient(sten, dir, 1);
		}
	}
	static Exp OnFace_2Order_Gradient(pNode pn,  //
			const Direction& dir) {
		// assert
		ASSERT(pn != nullptr);
		ASSERT(IsFaceDirection(dir));
		// stencil
		Orientation ori;
		Axes axe;
		FaceDirectionToOrientationAndAxes(dir, ori, axe);
		Stencil_1 sten(pn, axe, 1, 1);
		cvt x = pn->p(ori, axe);
		if (sten.non_null_nodes() == 1) {
			return Exp();
		} else if (sten.non_null_nodes() == 2) {
			return _2NodeGradient(sten, dir, 2);
		} else {
			return _3NodeGradient(x, sten, dir, 2);
		}
	}
	static Exp _2NodeGradient(ref_Stencil_1 stc, //
			const Direction& dir, //previous direction
			const int& order // previous order
			) {
		// assert
		ASSERT(stc.non_null_nodes() == 2);
		Axes axes = FaceDirectionToAxes(dir);
		pNode pc = stc.center_pnode();
		PExp cexp = _AverangeExpFromCenterLeaf(pc);
		// another pnode is close to center
		PExp nexp;
		pNode pn = stc.forward_pnode(1, stc.axes());
		if (pn == nullptr) {
			pn = stc.backward_pnode(1, stc.axes());
			nexp = _2NodeRelation(pc, pn, _M_, axes, order);
		} else {
			nexp = _2NodeRelation(pc, pn, _P_, axes, order);
		}
		// ---
		cvt x1 = nexp.P().val(axes);
		cvt x2 = pc->cp(axes);
		//  pc - pn / ( dc - dn);
		return __LinearGradient(x1, nexp.E(), x2, cexp.E());
	}
	static Exp _3NodeGradient(cvt x, ref_Stencil_1 stc, //
			const Direction& dir, //previous direction
			const int& order // previous order
			) {
		// assert
		ASSERT(stc.non_null_nodes() == 3);
		Axes axes = FaceDirectionToAxes(dir);
		pNode pc = stc.center_pnode();
		PExp expc = _AverangeExpFromCenterLeaf(pc);
		// forward pnode is close to center
		pNode pf = stc.forward_pnode(1, stc.axes());
		ASSERT(pf != nullptr);
		PExp expf = _2NodeRelation(pc, pf, _P_, axes, order);
		// backward pnode is close to center
		pNode pb = stc.backward_pnode(1, stc.axes());
		ASSERT(pb != nullptr);
		PExp expb = _2NodeRelation(pc, pb, _M_, axes, order);
		// ---
		cvt x1 = expf.P().val(axes);
		cvt x2 = expc.P().val(axes);
		cvt x3 = expb.P().val(axes);
		return __2OrderGradient(x, x1, expf.E(), x2, expc.E(), x3, expb.E());
	}

protected:
	/*
	 *  the idx in res is defined
	 */
	static int _GetPDataFromPNodeCenter(ref_PData res, const_pNode pc) {
		ASSERT(pc != nullptr);
		res.set_point(pc->cp(_X_), pc->cp(_Y_), pc->cp(_Z_)); //set point as the center of pNode
		for (St i = 0; i < res.size(); ++i) {
			res.flag(i) = PData::Flag_Center;     //set flag as Center
			res.val(i) = pc->cd(res.idx(i));     //get val
		}
		return _SUCCESS;
	}

	static void _AverangeValueFromCenterLeaf(        //
			PData& res,        //
			const_pNode pn) {
		// res is just initialized
		//improve for special case
		if (pn->is_leaf()) {
			_GetPDataFromPNodeCenter(res, pn);
			return;
		}
		//========================
		res.set_all_center();
		cvt volume = 0;
		std::function<void(const_pNode, int)> fun =
				[&res, &volume](const_pNode pnode, int dummy) {
					if (pnode->is_leaf()) {
						vt v = pnode->volume();
						for(St i = 0; i<res.size();++i) {
							res.val(i) = res.val(i)
							+ pnode->cd(res.idx(i)) * v;
						}
						volume = volume + v;
					}
				};
		int dummy = 1;
		pn->traversal(fun, dummy);
		for (St i = 0; i < res.size(); ++i) {
			res.val(i) /= volume;
		}
	}
	static void _CF_1(        //
			PData& res,        //
			const_pNode pn, //
			Direction dir) {
		// direction from pc to pn
		// direction is face direction
		// Get PExp on pn
		// One level down
		Direction odir = Opposite(dir);
		Axes aix;
		Orientation ori;
		FaceDirectionToOrientationAndAxes(odir, ori, aix);
		// point
		// Only one axes should be shift
		res.p(aix) = (pn->p(ori, aix) + pn->cp(aix)) * 0.5;
		int count = 0;
		Exp exp;
		for (St i = 0; i < NumFaces; i++) {
			if (pn->has_child(i) && is_on_direction(i, odir)) {
				count++;
				PExp pexp = _AverangeExpFromCenterLeaf(pn->child[i]);
				exp.plus(pexp.E());
			}
		}
		exp.times(1.0 / vt(count));
		//return PExp(poi, exp);
	}
	static PExp _CF_Exp_1(        //Coarse to fine  method 1
			const_pNode pn, const Direction& dir) {
		// direction from pc to pn
		// direction is face direction
		// Get PExp on pn
		// One level down
		Direction odir = Opposite(dir);
		Axes aix;
		Orientation ori;
		FaceDirectionToOrientationAndAxes(odir, ori, aix);
		Point poi(pn->cp(_X_), pn->cp(_Y_), pn->cp(_Z_));
		poi.val(aix) = (pn->p(ori, aix) + pn->cp(aix)) * 0.5;
		int count = 0;
		Exp exp;
		for (St i = 0; i < NumFaces; i++) {
			if (pn->has_child(i) && is_on_direction(i, odir)) {
				count++;
				PExp pexp = _AverangeExpFromCenterLeaf(pn->child[i]);
				exp.plus(pexp.E());
			}
		}
		exp.times(1.0 / vt(count));
		return PExp(poi, exp);
	}

	static PExp _AverangeExpFromCenterLeaf(        //
			const_pNode pn) {
		//improve for special case
		PExp res(pn, Exp());
		if (pn->is_leaf()) {
			res.E().plus(1.0, pn, 1.0);
			return res;
		}
		//========================
		// weight averange by volume
		vt volume = 0;
		std::function<void(const_pNode, int)> fun =
				[&res, &volume](const_pNode pnode, int dummy) {
					if (pnode->is_leaf()) {
						vt v = pnode->volume();
						res.E().plus(v, pnode, 1.0);
						volume = volume + v;
					}
				};
		int dummy = 1;
		pn->traversal(fun, dummy);
		res.E().times(1.0 / volume);
		return res;
	}

	inline static PExp _2NodeRelation(  //
			pNode pc, pNode pn, const Orientation& o, //
			const Axes& a, //
			const int& order) {
		switch (__2NodeRelation(pc, pn)) {
		case _Equal_: {
			return _AverangeExpFromCenterLeaf(pn);
			break;
		}
		case _FineCoarse_: {
			if (Dim == 2) {
				PExp nexp;
				Axes va = VerticalAxes2D(a);
				vt dis = pc->cp(va) - pn->cp(va);
				Orientation ori = dis > 0 ? _P_ : _M_;
				Direction nd = ToFaceDirection(ori, va);
				Exp exp = OnAxes(pn, nd, dis, order);
				Point p;
				p.val(a) = pn->cp(a);
				p.val(va) = pn->cp(va) + dis;
				nexp.set(p, exp);
				return nexp;
			} else { //3d
				SHOULD_NOT_REACH;
			}
			break;
		}
		case _CoarseFine_: {
			Direction dir = ToFaceDirection(o, a);
			return _CF_Exp_1(pn, dir);
			break;
		}
		default: {
			SHOULD_NOT_REACH;
			break;
		}
		}
		SHOULD_NOT_REACH;
		return PExp();
	}

	inline static int __2NodeRelation(const_pNode pc, const_pNode p) {
		St pcl = pc->get_level();
		St pl = p->get_level();
		if (pcl == pl && p->has_child() == false) {
			//case 1  on the same level
			return _Equal_;
		}
		if (pcl > pl) {
			//case 1  neighbor is coarse
			return _FineCoarse_;
		}
		if (pcl == pl) {
			//case 3  neighbor is fine
			return _CoarseFine_;
		}
		SHOULD_NOT_REACH;
		return _Null_;
	}

	inline static vt __LinearInterpolation( //
			const cvt& x, //
			const cvt& x0, const vt& y0, //
			const cvt& x1, const vt& y1 //
			) {
		return (y1 - y0) / (x1 - x0) * (x - x0) + y0;
	}

	inline static Exp __LinearInterpolation(  //
			const cvt& x,  //
			const cvt& x0, const Exp& y0,  //
			const cvt& x1, const Exp& y1) {
		Exp res(y1);
		res.minus(y0);
		res.times((x - x0) / (x1 - x0));
		res.plus(y0);
		return res;
	}
	/*
	 * linear interpolation in 2d space
	 * the (x,y) will project on the line of (x0, y0) to (x1, y1)
	 * the (x0,y0) set as an original point
	 * the project point define as (xp, yp)
	 *
	 */
	inline static Exp __LinearInterpolation(  //
			const cvt& x, const cvt& y, //
			const cvt& x0, const cvt& y0, const Exp& v0,  //
			const cvt& x1, const cvt& y1, const Exp& v1) {
		Vector vo((x1 - x0), (y1 - y0));  //(x0, y0)-->(x1, y1)
		Vector vp((x - x0), (y - y0));    //(x0, y0)-->(x, y)
		cvt r = dot(vo, vp) / dot(vo, vo);
		cvt xx0 = 0.0;
		cvt xx1 = vo.len();
		cvt xx = r * vo.len();
		return __LinearInterpolation(xx, xx0, v0, xx1, v1);
	}

	inline static Exp __BilinearInterpolation( //
			const cvt& x, const cvt& y, //
			const cvt& x1, const cvt& x2, //
			const cvt& y1, const cvt& y2,  //
			const Exp& v11, const Exp& v21, //
			const Exp& v22, const Exp& v12
			) {
		Exp res(v11);
		res.times((x2 - x) * (y2 - y));
		Exp tv21(v21);
		tv21.times((x - x1) * (y2 - y));
		res.plus(tv21);
		Exp tv22(v22);
		tv22.times((x - x1) * (y - y1));
		res.plus(tv22);
		Exp tv12(v12);
		tv12.times((x2 - x) * (y - y1));
		res.plus(tv12);

		res.times(1.0 / (x2 - x1) / (y2 - y1));
		return res;
	}

	inline static vt __2OrderInterpolation( //
			const cvt& x, //
			const cvt& x1, const vt& y1, //
			const cvt& x2, const vt& y2, //
			const cvt& x3, const vt& y3 //
			) {
		return y1 * (x - x2) * (x - x3) / (x1 - x2) / (x1 - x3)
				+ y2 * (x - x1) * (x - x3) / (x2 - x1) / (x2 - x3)
				+ y3 * (x - x1) * (x - x2) / (x3 - x1) / (x3 - x2);
	}

	inline static Exp __2OrderInterpolation( //
			const cvt& x, //
			const cvt& x1, const Exp& y1, //
			const cvt& x2, const Exp& y2, //
			const cvt& x3, const Exp& y3 //
			) {
		Exp res(y1);
		res.times((x - x2) * (x - x3) / (x1 - x2) / (x1 - x3));
		Exp t2(y2);
		t2.times((x - x1) * (x - x3) / (x2 - x1) / (x2 - x3));
		Exp t3(y3);
		t3.times((x - x1) * (x - x2) / (x3 - x1) / (x3 - x2));
		res.plus(t2);
		res.plus(t3);
		return res;
	}

	inline static Exp __LinearGradient( //
			const cvt& x0, const Exp& y0, //
			const cvt& x1, const Exp& y1 //
			) {
		Exp res(y1);
		res.minus(y0);
		res.times(1.0 / (x1 - x0));
		return res;
	}

	inline static Exp __2OrderGradient( //
			const cvt& x,  //
			const cvt& x1, const Exp& y1, //
			const cvt& x2, const Exp& y2, //
			const cvt& x3, const Exp& y3  //
			) {
		Exp res(y1);
		res.times((2.0 * x - x2 - x3) / (x1 - x2) / (x1 - x3));
		Exp res2(y2);
		res2.times((2.0 * x - x1 - x3) / (x2 - x1) / (x2 - x3));
		Exp res3(y3);
		res3.times((2.0 * x - x1 - x2) / (x3 - x1) / (x3 - x2));
		res.plus(res2);
		res.plus(res3);
		return res;
	}
};

}

#endif
