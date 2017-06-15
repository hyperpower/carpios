#include "interpolate.h"

namespace carpio {
/*
 *  There is only one node in the stencil
 *  The point (x,y) is in that node, and
 *  the value is equal to the center point.
 *  Stencil:  ---X---C---X---->
 */
int _1Node(Float& res, const St& idx, const Stencil_2D1& stc, Float x,
		Float y) {
	ASSERT(stc.non_null_nodes() == 1);
	// the non null nodes must be center node
	typename Stencil_2D1::const_pNode pc = stc.center_pnode();
	ASSERT(pc->is_in_on(x, y));
	// ---
	res = pc->cd(idx);
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
int _2NodeOnAxes(Float& res, const St& idx, const Stencil_2D1& stc, Float x) {
	// assert
	ASSERT(stc.non_null_nodes() == 2);
	typename Stencil_2D1::const_pNode pc = stc.center_pnode();
	// another pnode is close to center
	int flag = 0;
	typename Stencil_2D1::const_pNode po = stc.forward_pnode(1, stc.axes());
	if (po == nullptr) {
		flag = 1;
		po = stc.backward_pnode(1, stc.axes());
	}
	// ---
	typename Stencil_2D1::const_pNode ps = nullptr, pb = nullptr;
	//                                 s small       b big
	if (flag == 0) {
		ps = pc;
		pb = po;
	} else {
		ps = po;
		pb = pc;
	}

}

/*
 *  the idx in res is defined
 */
int _GetPDataFromPNodeCenter(PData_2D& res, const_pNode_2D pc) {
	ASSERT(pc != nullptr);
	res.set_point(pc->cp(_X_), pc->cp(_Y_));  //set point as the center of pNode
	for (St i = 0; i < res.size(); ++i) {
		res.flag(i) = PData_2D::Flag_Center;  //set flag as Center
		res.val(i) = pc->cd(res.idx(i));      //get val
	}
	return _SUCCESS;
}

void _AverangeValueFromCenterLeaf(        //
		PData_2D& res,        //
		const_pNode_2D pn) {
//improve for special case
	if (pn->is_leaf()) {
		_GetPDataFromPNodeCenter(res, pn);
		return;
	}

//========================
	res.set_all_center();
	St num = 0;
	std::function<void(const_pNode_2D, int)> fun =
			[&res, &num](const_pNode_2D pnode, int dummy) {
				if (pnode->is_leaf()) {
					for(St i = 0; i<res.size();++i) {
						res.val(i) = res.val(i) + pnode->cd(res.idx(i));
					}
					num++;
				}
			};
	int dummy = 1;
	pn->traversal(fun, dummy);
	for (St i = 0; i < res.size(); ++i) {
		res.val(i) /= num;
	}
}
template<class Node>
inline int __2NodeRelation(const Node* pc, const Node* p) {
	St pcl = pc->get_level();
	St pl = p->get_level();
	if (pcl == pl && p->has_child() == false) {
		//case 1  on the same level
		return _E_;
	}
	if (pcl > pl) {
		//case 1  neighbor is coarse
		return _F_C_;
	}
	if (pcl == pl) {
		//case 3  neighbor is fine
		return _C_F_;
	}
	SHOULD_NOT_REACH;
	return -1;
}
int _1Node(PData_2D& res, const Stencil_2D1& stc) {
	//assert
	ASSERT(stc.non_null_nodes() == 1);
	typename Stencil_2D1::const_pNode pc = stc.center_pnode();
	ASSERT(pc->is_in_on(res.x(), res.y()));
	// ---
	_AverangeValueFromCenterLeaf(res, pc);
	return _SUCCESS;
}
/*
 * Linear interpolation
 */
template<typename TYPE>
inline Float __LinearInterpolation( //
		const TYPE& x, //
		const TYPE& x0, const TYPE& y0, //
		const TYPE& x1, const TYPE& y1  //
		) {
	return (y1 - y0) / (x1 - x0) * (x - x0) + y0;
}
/*
 * Linear interpolation
 */
template<typename CVT, typename VT>
inline VT __2OrderInterpolation( //
		const CVT& x, //
		const CVT& x1, const VT& y1, //
		const CVT& x2, const VT& y2, //
		const CVT& x3, const VT& y3  //
		) {
	return y1 * (x - x2) * (x - x3) / (x1 - x2) / (x1 - x3)
			+ y2 * (x - x1) * (x - x3) / (x2 - x1) / (x2 - x3)
			+ y3 * (x - x1) * (x - x2) / (x3 - x1) / (x3 - x2);
}

int _2NodeOnAxes(PData_2D& res, const Stencil_2D1& stc) {
	// assert
	ASSERT(stc.non_null_nodes() == 2);
	typename Stencil_2D1::const_pNode pc = stc.center_pnode();
	PData_2D pdc(res.arr_idx(), pc->cp(_X_), pc->cp(_Y_));
	_AverangeValueFromCenterLeaf(pdc, pc);

	// another pnode is close to center
	typename Stencil_2D1::const_pNode pn = stc.forward_pnode(1, stc.axes());
	if (pn == nullptr) {
		pn = stc.backward_pnode(1, stc.axes());
	}
	PData_2D pdn(res.arr_idx(), pn->cp(_X_), pn->cp(_Y_));
	switch (__2NodeRelation(pc, pn)) {
	case _E_: {
		_AverangeValueFromCenterLeaf(pdn, pn);
		break;
	}
	case _F_C_: {
		//unfinish ==================
		SHOULD_NOT_REACH;
		break;
	}
	case _C_F_: {
		_AverangeValueFromCenterLeaf(pdn, pn);
		break;
	}
	}

	// ---
	for (St i = 0; i < res.size(); ++i) {
		Float x = res.p(stc.axes(0));
		Float x1 = pdc.p(stc.axes(0));
		Float y1 = pdc.val(i);
		Float x2 = pdn.p(stc.axes(0));
		Float y2 = pdn.val(i);
		res.flag(i) = PData_2D::Flag_Center;
		res.val(i) = __LinearInterpolation(x, x1, y1, x2, y2);
	}

	return _SUCCESS;
}

int _3NodeOnAxes(PData_2D& res, const Stencil_2D1& stc) {
	// assert
	ASSERT(stc.non_null_nodes() == 3);
	typename Stencil_2D1::const_pNode pc = stc.center_pnode();
	PData_2D pdc(res.arr_idx(), pc->cp(_X_), pc->cp(_Y_));
	_AverangeValueFromCenterLeaf(pdc, pc);
	// forward pnode is close to center
	typename Stencil_2D1::const_pNode pf = stc.forward_pnode(1, stc.axes());
	ASSERT(pf != nullptr);
	PData_2D pdf(res.arr_idx(), pf->cp(_X_), pf->cp(_Y_));
	switch (__2NodeRelation(pc, pf)) {
	case _E_: {
		_AverangeValueFromCenterLeaf(pdf, pf);
		break;
	}
	case _F_C_: {
		//unfinish ==================
		SHOULD_NOT_REACH;
		break;
	}
	case _C_F_: {
		_AverangeValueFromCenterLeaf(pdf, pf);
		break;
	}
	}
	// forward pnode is close to center
	typename Stencil_2D1::const_pNode pb = stc.backward_pnode(1, stc.axes());
	ASSERT(pb != nullptr);
	PData_2D pdb(res.arr_idx(), pb->cp(_X_), pb->cp(_Y_));
	switch (__2NodeRelation(pc, pb)) {
	case _E_: {
		_AverangeValueFromCenterLeaf(pdb, pb);
		break;
	}
	case _F_C_: {
		//unfinish ==================
		SHOULD_NOT_REACH;
		break;
	}
	case _C_F_: {
		_AverangeValueFromCenterLeaf(pdb, pb);
		break;
	}
	}
	// ---
	for (St i = 0; i < res.size(); ++i) {
		Cvt x = res.p(stc.axes(0));
		Cvt x1 = pdf.p(stc.axes(0));
		Vt y1 = pdf.val(i);
		Cvt x2 = pdc.p(stc.axes(0));
		Vt y2 = pdc.val(i);
		Cvt x3 = pdb.p(stc.axes(0));
		Vt y3 = pdb.val(i);
		res.flag(i) = PData_2D::Flag_Center;
		res.val(i) = __2OrderInterpolation(x, x1, y1, x2, y2, x3, y3);
	}
	return _SUCCESS;
}

void InterpolateOnFace_1Order( // 2D QuadTree Node
		pNode_2D pn,                      //node
		const Direction& dir,                //face
		const ArrayListV<St>& arridx,            //data index
		ArrayListV<Vt>& arrres                //data res
		) {
	// assert
	ASSERT(pn != nullptr);
	ASSERT(IsFaceDirection(dir));
	// get PData
	Cvt x = pn->p(dir, _X_);
	Cvt y = pn->p(dir, _Y_);
	PData_2D pdata(arridx, x, y);
	// stencil
	Orientation ori;
	Axes axe;
	FaceDirectionToOrientationAndAxes(dir, ori, axe);
	Stencil_2D1 sten(pn, axe, (IsP(ori) ? 1 : 0), (IsM(ori) ? 1 : 0));
	if (sten.non_null_nodes() == 1) {
		_1Node(pdata, sten);
	} else {
		_2NodeOnAxes(pdata, sten);
	}
	arrres = pdata.arr_val();
}

void InterpolateOnFace_1Order( // 2D QuadTree Node
		pNode_2D pn,                      //node
		const Direction& dir,                //face
		const St& idx,            //data index
		Vt& res                //data res
		) {
	ArrayListV<St> arridx(1);            //data index
	ArrayListV<Vt> arrres(1);
	arridx[0] = idx;
	InterpolateOnFace_1Order(pn, dir, arridx, arrres);
	res = arrres[0];
}

void InterpolateOnFace_2Order( // 2D QuadTree Node
		pNode_2D pn,                      //node
		const Direction& dir,                //face
		const ArrayListV<St>& arridx,            //data index
		ArrayListV<Vt>& arrres                //data res
		) {
	// assert
	ASSERT(pn != nullptr);
	ASSERT(IsFaceDirection(dir));
	// get PData
	Cvt x = pn->p(dir, _X_);
	Cvt y = pn->p(dir, _Y_);
	PData_2D pdata(arridx, x, y);
	// stencil
	Orientation ori;
	Axes axe;
	FaceDirectionToOrientationAndAxes(dir, ori, axe);
	Stencil_2D1 sten(pn, axe, 1, 1);
	if (sten.non_null_nodes() == 1) {
		_1Node(pdata, sten);
	} else if (sten.non_null_nodes() == 2) {
		_2NodeOnAxes(pdata, sten);
	} else {
		_3NodeOnAxes(pdata, sten);
	}
	arrres = pdata.arr_val();
}

void InterpolateOnFace_2Order( // 2D QuadTree Node
		pNode_2D pn,                      //node
		const Direction& dir,             //face
		const St& idx,            //data index
		Vt& res                //data res
		) {
	ArrayListV<St> arridx(1);            //data index
	ArrayListV<Vt> arrres(1);
	arridx[0] = idx;
	InterpolateOnFace_2Order(pn, dir, arridx, arrres);
	res = arrres[0];
}

void Interpolate_1Order(        // 2D QuadTree Node
		pNode_2D pn,                  //node
		const Cvt& x, const Cvt& y,              //point
		const ArrayListV<St>& arridx,            //data index
		ArrayListV<Vt>& arrres                //data res
		) {
	// assert
	ASSERT(pn != nullptr);
	ASSERT(pn->is_in_on(x,y));
	// get PData

	// stencil

}

}
