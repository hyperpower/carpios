#ifndef STENCIL_HPP_
#define STENCIL_HPP_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "algebra/space.hpp"
#include "algebra/array_list.hpp"

#include <functional>
#include <math.h>

namespace carpio {

template<typename COO_VALUE, typename VALUE, St DIM, St DIMST>
class Stencil_ {
public:
	static const St Dim = DIMST;  //Dimension of stencil

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef const Node_<COO_VALUE, VALUE, DIM> const_Node;
	typedef const Node_<COO_VALUE, VALUE, DIM>* const_pNode;
	typedef Cell_<COO_VALUE, DIM> Cell;
	typedef Cell_<COO_VALUE, DIM> *pCell;
	typedef Data_<VALUE, DIM> Data;
	typedef Data *pData;
protected:
	/*
	 *  data
	 */
	template<typename CV, typename V, int D>
	struct sNode_ { //Stencil node
		Node_<CV, V, D>* pnode;
		int type;  //the type indicates that the node is new or already exist.
	};
	typedef sNode_<COO_VALUE, VALUE, DIM> sNode;
	typedef sNode_<COO_VALUE, VALUE, DIM>* psNode;

	SpaceT<sNode, Dim> _pnodes;
	ArrayListT<Axes> _axes;
	ArrayListT<St> _steps_f;
	ArrayListT<St> _steps_b;
protected:
	/*
	 * function
	 */
	pNode _neighbor_forward(pNode pn, Axes a) {
		ASSERT(pn != nullptr);
		pNode res = nullptr;
		if (a == _X_) {
			res = pn->get_neighbor_fast(_XP_);
		}
		if (a == _Y_) {
			res = pn->get_neighbor_fast(_YP_);
		}
		if (a == _Z_) {
			res = pn->get_neighbor_fast(_ZP_);
		}
		return res;
	}
	pNode _neighbor_backward(pNode pn, Axes a) {
		ASSERT(pn != nullptr);
		pNode res = nullptr;
		if (a == _X_) {
			res = pn->get_neighbor_fast(_XM_);
		}
		if (a == _Y_) {
			res = pn->get_neighbor_fast(_YM_);
		}
		if (a == _Z_) {
			res = pn->get_neighbor_fast(_ZM_);
		}
		return res;
	}
	St _get_idx_c(St sf, St sb) const {
		return sb;
	}
	void _construct_1d(pNode pnc, Axes a, St sf, St sb) {
		ASSERT(pnc != nullptr);
		ASSERT(Dim == 1);
		_axes[0] = a;
		_steps_f[0] = sf;
		_steps_b[0] = sb;
		_pnodes.reconstruct(sf + sb + 1);
		// set null
		for (St i = 0; i < _pnodes.size(); ++i) {
			_pnodes.at_1d(i).pnode = nullptr;
			_pnodes.at_1d(i).type = 0;   //which is not created by this class
		}
		_set_1d(pnc);
	}

	void _construct_2d_initial(Axes a1, St sf1, St sb1, Axes a2, St sf2,
			St sb2) {
		ASSERT(Dim == 2);
		_axes[0] = a1;
		_steps_f[0] = sf1;
		_steps_b[0] = sb1;
		_axes[1] = a2;
		_steps_f[1] = sf2;
		_steps_b[1] = sb2;
		_pnodes.reconstruct(sf1 + sb1 + 1, sf2 + sb2 + 1);
		// set null
		for (St i = 0; i < _pnodes.size(); ++i) {
			_pnodes.at_1d(i).pnode = nullptr;
			_pnodes.at_1d(i).type = 0;   //whish is not created by this class
		}
	}
	void _set_2d_1(pNode pn, St idx, Axes a, St sf, St sb) {
		ASSERT(Dim == 2);
		ASSERT(pn != nullptr);
		St ixc = (a == _axes[0]) ? sb : idx;
		St iyc = (a == _axes[0]) ? idx : sb;
		// set center node
		_pnodes(ixc, iyc).pnode = pn;
		// find neighbor forward;
		pNode pc = pn;
		for (St i = 0; i < sf; ++i) {
			pNode pnt = nullptr;
			pnt = _neighbor_forward(pc, a);
			if (pnt == nullptr) {
				break;
			} else {
				pNode& ref_pn = (
						(a == _X_) ?
								_pnodes(ixc + i + 1, iyc).pnode :
								_pnodes(ixc, iyc + i + 1).pnode);
				if (ref_pn != nullptr) {
					if (pnt->get_level() > ref_pn->get_level()) {
						ref_pn = pnt;
					}
				} else {
					ref_pn = pnt;
				}
				pc = pnt;
			}
		}
		// find neighbor backward;
		pc = pn;
		for (St i = 0; i < sb; ++i) {
			pNode pnt = nullptr;
			pnt = _neighbor_backward(pc, a);
			if (pnt == nullptr) {
				break;
			} else {
				pNode& ref_pn =
						(a == _X_) ?
								_pnodes(ixc - i - 1, iyc).pnode :
								_pnodes(ixc, iyc - i - 1).pnode;
				if (ref_pn != nullptr) {
					if (pnt->get_level() > ref_pn->get_level()) {
						ref_pn = pnt;
					}
				} else {
					ref_pn = pnt;
				}
				pc = pnt;
			}
		}
	}
	void _construct_2d(pNode pnc, Axes a1, St sf1, St sb1, Axes a2, St sf2,
			St sb2) {
		ASSERT(pnc != nullptr);
		_construct_2d_initial(a1, sf1, sb1, a2, sf2, sb2);
		// set center node
		_set_2d(pnc);
	}
public:
	/*
	 *  1d constructor
	 */
	Stencil_(pNode pnc, Axes a1, St sf1, St sb1) :
			_pnodes(), _axes(Dim), _steps_f(Dim), _steps_b(Dim) {
		ASSERT(Dim == 1);
		_construct_1d(pnc, a1, sf1, sb1);
	}
	Stencil_(pNode pnc, Axes a1, St sf1, St sb1, Axes a2, St sf2, St sb2) :
			_pnodes(), _axes(Dim), _steps_f(Dim), _steps_b(Dim) {
		ASSERT(Dim == 2);
		_construct_2d(pnc, a1, sf1, sb1, a2, sf2, sb2);
	}

protected:
	/*
	 * iterator
	 */

public:
	/*
	 * set
	 */
	void _set_1d(pNode pnc) {
		ASSERT(Dim == 1);
		// set center node
		_pnodes.at_1d(_steps_b[0]).pnode = pnc;
		//find neighbor forward
		pNode pc = pnc;
		for (St i = 0; i < _steps_f[0]; ++i) {
			pNode pnt = nullptr;
			pnt = _neighbor_forward(pc, _axes[0]);
			if (pnt == nullptr) {
				break;
			} else {
				_pnodes.at_1d(_steps_b[0] + i + 1).pnode = pnt;
				pc = pnt;
			}
		}
		//find neighbor backward
		pc = pnc;
		for (St i = 0; i < _steps_b[0]; ++i) {
			pNode pnt = nullptr;
			pnt = _neighbor_backward(pc, _axes[0]);
			if (pnt == nullptr) {
				break;
			} else {
				_pnodes.at_1d(_steps_b[0] - 1 - i).pnode = pnt;
				pc = pnt;
			}
		}
	}
	void _set_2d(pNode pnc) {
		// set center node
		Axes a1 = _axes[0];
		St sf1 = _steps_f[0];
		St sb1 = _steps_b[0];
		Axes a2 = _axes[1];
		St sf2 = _steps_f[1];
		St sb2 = _steps_b[1];

		_pnodes(sb1, sb2).pnode = pnc;
		_set_2d_1(pnc, sb2, a1, sf1, sb1);
		_set_2d_1(pnc, sb1, a2, sf2, sb2);
		// loop a1, find pnode on a2
		for (St i = 0; i < sb1; ++i) {
			pNode pc = _pnodes(i, sb2).pnode;
			if (pc != nullptr) {
				_set_2d_1(pc, i, a2, sf2, sb2);
			}
		}
		for (St i = sb1 + 1; i <= sb1 + sf1; ++i) {
			pNode pc = _pnodes(i, sb2).pnode;
			if (pc != nullptr) {
				_set_2d_1(pc, i, a2, sf2, sb2);
			}
		}
		// loop a2, find pnode on a1
		for (St i = 0; i < sb2; ++i) {
			pNode pc = _pnodes(sb1, i).pnode;
			if (pc != nullptr) {
				_set_2d_1(pc, i, a1, sf1, sb1);
			}
		}
		for (St i = sb2 + 1; i <= sb2 + sf2; ++i) {
			pNode pc = _pnodes(sb1, i).pnode;
			if (pc != nullptr) {
				_set_2d_1(pc, i, a1, sf1, sb1);
			}
		}
	}

	void set(pNode pnc) {
		ASSERT(pnc != nullptr);
		this->clear();
		if (Dim == 1) {
			_set_1d(pnc);
		} else if (Dim == 2) {
			_set_2d(pnc);
		} else {
			ASSERT_MSG(false, "unfinish");
		}
	}

	pNode operator()(St i, St j = 0, St k = 0) {
		return _pnodes(i, j, k).pnode;
	}
	const_pNode operator()(St i, St j = 0, St k = 0) const {
		return _pnodes(i, j, k).pnode;
	}

	pNode at_1d(St i) {
		return _pnodes.at_1d(i).pnode;
	}
	const_pNode at_1d(St i) const {
		return _pnodes.at_1d(i).pnode;
	}
	St size() const {
		return _pnodes.size();
	}
	St null_nodes() const {
		St res = 0;
		for (St i = 0; i < _pnodes.size(); ++i) {
			if (at_1d(i) == nullptr) {
				res++;
			}
		}
		return res;
	}
	St non_null_nodes() const {
		St res = 0;
		for (St i = 0; i < _pnodes.size(); ++i) {
			if (at_1d(i) != nullptr) {
				res++;
			}
		}
		return res;
	}
	St dim() const {
		return Dim;
	}
	/*
	 *  get
	 */
	pNode center_pnode() {
		if (Dim == 1) {
			return _pnodes(_steps_b[0], 0, 0).pnode;
		} else if (Dim == 2) {
			return _pnodes(_steps_b[0], _steps_b[1], 0).pnode;
		} else {
			return _pnodes(_steps_b[0], _steps_b[1], _steps_b[2]).pnode;
		}
	}
	const_pNode center_pnode() const {
		if (Dim == 1) {
			return _pnodes(_steps_b[0], 0, 0).pnode;
		} else if (Dim == 2) {
			return _pnodes(_steps_b[0], _steps_b[1], 0).pnode;
		} else {
			return _pnodes(_steps_b[0], _steps_b[1], _steps_b[2]).pnode;
		}
	}
	Axes axes(St i = 0) const {
		return _axes[i];
	}
	void clear() {
		for (St i = 0; i < _pnodes.size(); ++i) {
			_pnodes.at_1d(i).pnode = nullptr;
			_pnodes.at_1d(i).type = 0;   //whish is not created by this class
		}
	}
protected:
	bool _is_valid_axes(Axes a) {
		for (St i = 0; i < _axes.size(); ++i) {
			if (_axes[i] == a) {
				return true;
			}
		}
		return false;
	}
	St _to_arraylist_idx(Axes a) const {
		St i = 0;
		for (; i < _axes.size(); ++i) {
			if (_axes[i] == a) {
				return i;
			}
		}
		ASSERT_MSG(false, "Not valid axes");
		return i;
	}
public:
	pNode forward_pnode(St step, Axes a = _X_) {
		St ai = _to_arraylist_idx(a);
		if (!(step > 0 && step <= _steps_f[ai])) {
			return nullptr;
		}
		if (ai == 0) {
			return _pnodes(_steps_b[ai] + step, 0, 0).pnode;
		} else if (ai == 1) { //ai==1
			return _pnodes(0, _steps_b[ai] + step, 0).pnode;
		} else {
			return _pnodes(0, 0, _steps_b[ai] + step).pnode;
		}
	}
	const_pNode forward_pnode(St step, Axes a = _X_) const {
		St ai = _to_arraylist_idx(a);
		if (!(step > 0 && step <= _steps_f[ai])) {
			return nullptr;
		}
		if (ai == 0) {
			return _pnodes(_steps_b[ai] + step, 0, 0).pnode;
		} else if (ai == 1) { //ai==1
			return _pnodes(0, _steps_b[ai] + step, 0).pnode;
		} else {
			return _pnodes(0, 0, _steps_b[ai] + step).pnode;
		}
	}
	pNode backward_pnode(St step, Axes a = _X_) {
		St ai = _to_arraylist_idx(a);
		if (!(step > 0 && step <= _steps_b[ai])) {
			return nullptr;
		}
		if (ai == 0) {
			return _pnodes(_steps_b[ai] - step, 0, 0).pnode;
		} else if (ai == 1) { //ai==1
			return _pnodes(0, _steps_b[ai] - step, 0).pnode;
		} else {
			return _pnodes(0, 0, _steps_b[ai] - step).pnode;
		}
	}
	const_pNode backward_pnode(St step, Axes a = _X_) const {
		St ai = _to_arraylist_idx(a);
		if (!(step > 0 && step <= _steps_b[ai])) {
			return nullptr;
		}
		if (ai == 0) {
			return _pnodes(_steps_b[ai] - step, 0, 0).pnode;
		} else if (ai == 1) { //ai==1
			return _pnodes(0, _steps_b[ai] - step, 0).pnode;
		} else {
			return _pnodes(0, 0, _steps_b[ai] - step).pnode;
		}
	}
protected:
	pNode _get_pnode_1d(St step) {
		ASSERT(Dim == 1);
		return _pnodes(_steps_b[0] + step).pnode;
	}
	const_pNode _get_pnode_1d(St step) const {
		ASSERT(Dim == 1);
		return _pnodes(_steps_b[0] + step).pnode;
	}
	pNode _get_pnode_2d(St step1, St step2) {
		ASSERT(Dim == 2);
		return _pnodes(_steps_b[0] + step1, _steps_b[1] + step2).pnode;
	}
	const_pNode _get_pnode_2d(St step1, St step2) const {
		ASSERT(Dim == 2);
		return _pnodes(_steps_b[0] + step1, _steps_b[1] + step2).pnode;
	}
public:
	pNode get_pnode(St step1, St step2 = 0, St step3 = 0) {
		if (Dim == 1) {
			return _get_pnode_1d(step1);
		} else if (Dim == 2) {
			return _get_pnode_2d(step1, step2);
		} else {
			SHOULD_NOT_REACH;
			return nullptr; //unfinish
		}
	}

	const_pNode get_pnode(St step1, St step2 = 0, St step3 = 0) const {
		if (Dim == 1) {
			return _get_pnode_1d(step1);
		} else if (Dim == 2) {
			return _get_pnode_2d(step1, step2);
		} else {
			return nullptr; //unfinish
		}
	}

	/*
	 *  show
	 */
	void show() const {
		std::cout << " stencil dim = " << Dim << "\n";
		std::cout << " size        = " << size() << "\n";
		std::cout << " null        = " << null_nodes() << "\n";
		std::cout << " no null     = " << size() - null_nodes() << "\n";
		if (Dim == 1) {
			std::cout << " ---";
			for (St i = 0; i < _pnodes.size(); ++i) {
				if (i == _steps_b(0) && at_1d(i) != nullptr) {
					std::cout << "C";
				} else if (at_1d(i) != nullptr) {
					std::cout << "O";
				} else {
					std::cout << "X";
				}
				std::cout << "---";
			}
			std::cout << "--> " << ToString(_axes(0)) << "\n";
		}
	}

};

}

#endif
