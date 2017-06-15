/************************
 //  \file   ts_instersect.h
 //  \brief
 // 
 //  \author czhou
 //  \date   19 juin 2015 
 ***********************/
#ifndef TS_INTERSECT_H_
#define TS_INTERSECT_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_triangle.h"
#include "ts_edge.h"
#include "ts_AABBox.h"
#include "ts_BBTree.h"

namespace TS {

template<class BOX>
inline bool do_intersect(BOX* pb,
		Triangle<typename BOX::value_type, BOX::Dim>* pt) {
	return pb->are_overlapping(pt);
}

template<class BOX>
inline bool do_intersect(Triangle<typename BOX::value_type, BOX::Dim>* tri1,
		Triangle<typename BOX::vt, BOX::Dim>* tri2) {

}

template<class BOX>
bool do_intersect_box_obj(BOX* pb, utPointer upt, OBJ_TYPE type) {
	if (type == EMPTY) {
		return false;
	}
	if (type == TRIANGLE) {
		Triangle<typename BOX::value_type, BOX::Dim>* pt = (Triangle<
				typename BOX::value_type, BOX::Dim>*) upt;
		return do_intersect(pb, pt);
	}
	return false;
}

template<typename BOX>
void cb_intersect_box(Int& flag, BBNode<BOX>* pn, utPointer utp) {
	Array<utPointer, 2>& arr = (*(Array<utPointer, 2>*) utp);
	BOX& boxin = (*(BOX*) arr[0]);
	if (!pn->box.are_overlapping(&boxin)) {
		flag = -1;
	} else {
		flag = 1;
	}
	if (pn->is_leaf() && (flag != -1)) {
		List<BOX*>& lres = (*(List<BOX*>*) arr[1]);
		lres.push_back(&(pn->box));
	}
}

template<typename BOX>
void cb_intersect_box_obj(Int& flag, BBNode<BOX>* pn, utPointer utp) {
	Array<utPointer, 2>& arr = (*(Array<utPointer, 2>*) utp);
	BOX& boxin = (*(BOX*) arr[0]);
	if (!pn->box.are_overlapping(&boxin)) {
		flag = -1;
	} else {
		flag = 1;
	}
	if (pn->is_leaf() && (flag != -1)) {
		if (boxin.bounded != nullptr) {
			if (do_intersect_box_obj(&(pn->box), boxin.bounded, boxin.type)) {
				List<BOX*>& lres = (*(List<BOX*>*) arr[1]);
				lres.push_back(&(pn->box));
			}
		}
	}
}

template<class BOX>
bool do_intersect_box_box(BBTree<BOX>* pt, BOX* pb) {
	Int flag = 0;
	List<BOX*> lres;
	Array<utPointer, 2> arr;
	arr[0] = pb;
	arr[1] = &lres;
	pt->PreOrder_flag(flag, cb_intersect_box, &arr);
	if (!lres.empty()) {
		return true;
	} else {
		return false;
	}
}

template<class BOX>
bool do_intersect_box_box(BBTree<BOX>* pt, BOX* pb, List<BOX*>& lres) {
	Int flag = 0;
	Array<utPointer, 2> arr;
	arr[0] = pb;
	arr[1] = &lres;
	pt->PreOrder_flag(flag, cb_intersect_box, &arr);
	if (!lres.empty()) {
		return true;
	} else {
		return false;
	}
}
// this function will return the list boxes in pt which intersect with
// the object in box pb
template<class BOX>
bool do_intersect_box_obj(BBTree<BOX>* pt, BOX* pb, List<BOX*>& lres) {
	Int flag = 0;
	Array<utPointer, 2> arr;
	arr[0] = pb;
	arr[1] = &lres;
	pt->PreOrder_flag(flag, cb_intersect_box_obj, &arr);
	if (!lres.empty()) {
		return true;
	} else {
		return false;
	}
}

template<class BOX>
bool do_intersect_obj_obj(BBTree<BOX>* pt, BOX* pb, List<BOX*>& lres) {
	Int flag = 0;
	Array<utPointer, 2> arr;
	arr[0] = pb;
	arr[1] = &lres;
	pt->PreOrder_flag(flag, cb_intersect_box_obj, &arr);
	if (lres.empty()) {
		return false;
	} else {
		// there are some boxes intersect with the obj in pb

	}
}

template<class BOX>
bool DoIntersect(BBTree<BOX>* pt, BOX* pb, List<BOX*>& lres) {
	typedef BBTree<BOX> Tree;
	typedef BOX Box;
	typename Tree::Fun_flag fun = [&pb, lres](typename Tree::pNode pn, bool& flag) {
		if (!pn->box.are_overlapping(&pb)) {
			flag = false;
		} else {
			flag = true;
		}
		if (pn->is_leaf() && (flag != -1)) {
			lres.push_back(&(pn->box));
		}
	};

	//Int flag = 0;
	//Array<utPointer, 2> arr;
	//arr[0] = pb;
	//arr[1] = &lre
	//pt->PreOrder_flag(flag, cb_intersect_box_obj, &arr);
	//if (lres.empty()) {
	//	return false;
	//} else {
		// there are some boxes intersect with the obj in pb

	//}
}

}

#endif /* TS_INSTERSECT_H_ */
