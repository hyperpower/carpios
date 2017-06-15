#ifndef GRID_H_
#define GRID_H_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"

#include "algebra/space.hpp"

namespace carpio {

template<typename COO_VALUE, typename VALUE, St DIM>
class Grid_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef Grid_<COO_VALUE, VALUE, DIM> Self;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pSelf;
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pSelf;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM>& ref_Node;
	typedef const Node_<COO_VALUE, VALUE, DIM>& const_ref_Node;
	typedef const Node_<COO_VALUE, VALUE, DIM> const_Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef const Node_<COO_VALUE, VALUE, DIM> * const_pNode;

	typedef Path_<Dim> Path;
	typedef Face_<Node, pNode> Face;
	typedef Face_<Node, pNode> *pFace;
	typedef typename SpaceT<pNode, Dim>::reference reference;
	typedef typename SpaceT<pNode, Dim>::const_reference const_reference;
	typedef typename SpaceT<pNode, Dim>::size_type size_type;

	typedef void (*pfunction)(pNode, utPointer);

	typedef void (*pfunction_conditional)(arrayList &, pNode, utPointer);

	typedef typename SpaceT<pNode, Dim>::iterator iterator;
	typedef typename SpaceT<pNode, Dim>::const_iterator const_iterator;

	/*
	 *  data
	 */
	SpaceT<pNode, Dim> nodes;

	/*
	 *  constructor
	 */

	Grid_(St ni, cvt ox, cvt dx, //
			St nj = 0, cvt oy = 0, cvt dy = 0, //
			St nk = 0, cvt oz = 0, cvt dz = 0) :
			nodes(ni, nj, nk) {
		St i_1d = 0;
		for (St i = 0; i < ni; i++) {
			for (St j = 0; j < ((Dim >= 2) ? nj : 1); j++) {
				for (St k = 0; k < ((Dim == 3) ? nk : 1); k++) {
					// new
					Path p(1);
					nodes.at_1d(i_1d) =  //
							new Node(nullptr, // father
									0, // type
									0, //level
									i_1d, //root idx
									0, //child idx
									p, ox + (i + 0.5) * dx, 0.5 * dx, //x
									oy + (j + 0.5) * dy, 0.5 * dy, //y
									oz + (k + 0.5) * dz, 0.5 * dz); //z
					i_1d++;
				}
			}
		}
	}

protected:
	void _delete() {
		for (St i = 0; i < nodes.size(); i++) {
			if (nodes.at_1d(i) != nullptr) {
				delete nodes.at_1d(i);
			}
		}
	}

public:
	~Grid_() {
		_delete();
	}
	/*
	 *  index as 2 or 3 demension
	 */
	reference operator()(size_type i, size_type j = 0, size_type k = 0) {
		return nodes(i, j, k);
	}

	const_reference operator()(size_type i, size_type j = 0,
			size_type k = 0) const {
		return nodes(i, j, k);
	}
	/*
	 *  index as 1 demension
	 */
	reference at_1d(size_type i) {
		return nodes.at_1d(i);
	}

	const_reference at_1d(size_type i) const {
		return nodes.at_1d(i);
	}

	/*
	 *  size
	 */
	inline St size_i() const {
		return nodes.size_i();
	}

	inline St size_j() const {
		return nodes.size_j();
	}

	inline St size_k() const {
		return (Dim < 3) ? 0 : nodes.size_k();
	}

	inline bool empty() const {
		if (nodes.size() <= 0) {
			return true;
		} else {
			return false;
		}
	}

	St get_dim() const {
		return Dim;
	}

	St size() const {
		return nodes.size();
	}

	St get_num_root() const {
		St num = 0;
		for (auto iter = nodes.begin(); iter != nodes.end(); ++iter) {
			if ((*iter) != nullptr) {
				num++;
			}
		}
		return num;
	}

	const pNode get_last_root_pNode() const {
		for (double i = nodes.size() - 1; i >= 0; --i) {
			if (nodes.at_1d(St(i)) != nullptr) {
				return nodes.at_1d(St(i));
			}
		}
		return nullptr;
	}
	pNode get_last_root_pNode() {
		for (double i = nodes.size() - 1; i >= 0; --i) {
			if (nodes.at_1d(St(i)) != nullptr) {
				return nodes.at_1d(St(i));
			}
		}
		return nullptr;
	}
	const pNode get_first_root_pNode() const {
		for (St i = 0; i < nodes.size(); ++i) {
			if (nodes.at_1d(i) != nullptr) {
				return nodes.at_1d(i);
			}
		}
		return nullptr;
	}
	pNode get_first_root_pNode() {
		for (St i = 0; i < nodes.size(); ++i) {
			if (nodes.at_1d(i) != nullptr) {
				return nodes.at_1d(i);
			}
		}
		return nullptr;
	}

	St get_first_root_pNode_idx1d() const {
		St i = 0;
		for (; i < nodes.size(); ++i) {
			if (nodes.at_1d(i) != nullptr) {
				return i;
			}
		}
		return i;
	}

	St get_last_root_pNode_idx1d() const {
		size_type i = nodes.size() - 1;
		for (; i >= 0; --i) {
			if (nodes.at_1d(i) != nullptr) {
				return i;
			}
		}
		return i;
	}

	/*
	 *  iterator
	 */
	iterator begin() {
		return nodes.begin();
	}

	const_iterator begin() const {
		return nodes.begin();
	}

	iterator end() {
		return nodes.end();
	}

	const_iterator end() const {
		return nodes.end();
	}

	void connect_root() {
		for (St i = 0; i < nodes.size_i(); i++) {
			for (St j = 0; j < (Dim >= 2 ? nodes.size_j() : 1); j++) {
				for (St k = 0; k < (Dim == 3 ? nodes.size_k() : 1); k++) {
					pNode xp = nullptr, xm = nullptr, //x
							yp = nullptr, ym = nullptr, //y
							zp = nullptr, zm = nullptr; //z
					pNode cnode = nodes(i, j, k);
					if (cnode != nullptr) {
						// x m  and  p
						xm = nodes.check_idx_ijk(i - 1, j, k) ?
								nodes(i - 1, j, k) : nullptr;
						xp = nodes.check_idx_ijk(i + 1, j, k) ?
								nodes(i + 1, j, k) : nullptr;
						ym = nodes.check_idx_ijk(i, j - 1, k) ?
								nodes(i, j - 1, k) : nullptr;
						yp = nodes.check_idx_ijk(i, j + 1, k) ?
								nodes(i, j + 1, k) : nullptr;
						zm = nodes.check_idx_ijk(i, j, k - 1) ?
								nodes(i, j, k - 1) : nullptr;
						zp = nodes.check_idx_ijk(i, j, k + 1) ?
								nodes(i, j, k + 1) : nullptr;
						cnode->set_neighbor(xm, xp, ym, yp, zm, zp);
					}
				}
			}
		}
	}
	void connect_nodes() {
		for (St i = 0; i < nodes.size(); i++) {
			pNode pn = nodes.at_1d(i);
			if (pn != nullptr) {
				pn->connect_nodes();
			}
		}
	}
	/*
	 * set
	 */
	void set_data_index() {
		// index as iterator order
		int i = 0;
		for (iterator_leaf iter = this->begin_leaf(); iter != this->end_leaf();
				++iter) {
			Node& node = *iter;
			node.d_idx() = i;
			i++;
		}
	}
	/*
	 * new data
	 */
	void new_data_on_leaf(const St& nc, const St& nf, const St& nv,
			const St& nutp) {
		for (iterator_leaf iter = this->begin_leaf(); iter != this->end_leaf();
				++iter) {
			pNode pn = iter.get_pointer();
			if (pn != nullptr) {
				pn->new_data(nc, nf, nv, nutp);
			}
		}
	}
	void resize_data_on_leaf(const St& nc, const St& nf, const St& nv,
			const St& nutp) {
		for (iterator_leaf iter = this->begin_leaf(); iter != this->end_leaf();
				++iter) {
			pNode pn = iter.get_pointer();
			if (pn != nullptr) {
				pn->resize_data(nc, nf, nv, nutp);
			}
		}
	}
	/*
	 * get grid range
	 */
	cvt max(Axes a) const {
		cvt max;
		int f = 0;
		for (const_iterator_leaf iter = this->begin_leaf();
				iter != this->end_leaf(); ++iter) {
			const_pNode pn = iter.get_pointer();
			if (pn != nullptr) {
				cvt lm = pn->p(_P_, a);
				if (f == 0) {
					max = lm;
					f = 1;
				} else if (lm > max) {
					max = lm;
				}
			}
		}
		return max;
	}
	cvt min(Axes a) const {
		cvt min;
		int f = 0;
		for (const_iterator_leaf iter = this->begin_leaf();
				iter != this->end_leaf(); ++iter) {
			const_pNode pn = iter.get_pointer();
			if (pn != nullptr) {
				cvt lm = pn->p(_M_, a);
				if (f == 0) {
					min = lm;
					f = 1;
				} else if (lm < min) {
					min = lm;
				}
			}
		}
		return min;
	}
	/*
	 * for each
	 */
	typedef std::function<void(const_ref_Node)> Fun_const_ref_Node;
	typedef std::function<void(ref_Node)> Fun_ref_Node;
	typedef std::function<void(const_pNode)> Fun_const_pNode;
	typedef std::function<void(pNode)> Fun_pNode;

	void for_each_root(Fun_const_ref_Node fun) const {
		for (St i = 0; i < nodes.size(); i++) {
			const_pNode pn = nodes.at_1d(i);
			if (pn != nullptr) {
				fun((*pn));
			}
		}
	}

	void for_each_root(Fun_ref_Node fun) {
		for (St i = 0; i < nodes.size(); i++) {
			pNode pn = nodes.at_1d(i);
			if (pn != nullptr) {
				fun((*pn));
			}
		}
	}
	void for_each_root(Fun_const_pNode fun) const {
		for (St i = 0; i < nodes.size(); i++) {
			const_pNode pn = nodes.at_1d(i);
			if (pn != nullptr) {
				fun(pn);
			}
		}
	}

	void for_each_root(Fun_pNode fun) {
		for (St i = 0; i < nodes.size(); i++) {
			pNode pn = nodes.at_1d(i);
			if (pn != nullptr) {
				fun(pn);
			}
		}
	}

	void for_each_leaf(Fun_pNode fun) {
		// this function can be improved by using traversal
		//for (iterator_leaf iter = this->begin_leaf(); iter != this->end_leaf();
		//		++iter) {
		//	pNode pn = iter.get_pointer();
		//	if (pn != nullptr) {
		//		fun(pn);
		//	}
		//}
		// new function using traversal
		for (St i = 0; i < nodes.size(); i++) {
			pNode pn = nodes.at_1d(i);
			if (pn != nullptr) {
				TraversalLeaf(pn, fun);
			}
		}
	}

	void for_each_leaf(Fun_const_pNode fun) const {
		for (const_iterator_leaf iter = this->begin_leaf();
				iter != this->end_leaf(); ++iter) {
			const_pNode pn = iter.get_pointer();
			if (pn != nullptr) {
				fun(pn);
			}
		}
	}

	/*
	 *  iterator leaf node
	 */
protected:
	template<typename COV, typename V, St D, class _Ref, class _Ptr>
	class iterator_leaf_ {
	public:
		typedef COV cvt;
		typedef VALUE vt;
		typedef Node_<COV, V, D> Node;
		typedef Node_<COV, V, D> *pNode;
		typedef const Node_<COV, V, D> const_Node;
		typedef const Node_<COV, V, D>* const_pNode;

		typedef iterator_leaf_<COV, V, D, Node&, pNode> iterator;
		typedef iterator_leaf_<COV, V, D, const_Node&, const_pNode> const_iterator;
		typedef iterator_leaf_<COV, V, D, _Ref, _Ptr> Self;

		typedef Grid_<COV, V, D> Grid;
		typedef Grid_<COV, V, D>* pGrid;
		typedef const Grid_<COV, V, D>* const_pGrid;

		typedef Node value_type;
		typedef _Ptr pointer;
		typedef _Ref reference;

		const_pGrid _f;
		St _idx;
		_Ptr _ptr;

		iterator_leaf_() {
			_ptr = nullptr;
			_f = nullptr;
			_idx = 0;
		}
		iterator_leaf_(const_pGrid f, St idx, _Ptr ptr) :
				_f(f), _idx(idx), _ptr(ptr) {
		}
		iterator_leaf_(const iterator& _x) :
				_f(_x._f), _idx(_x._idx), _ptr(_x._ptr) {
		}
	protected:
		pNode _incr_root() {
			St count = 0;
			for (St ii = (_idx + 1) % _f->size(); count < _f->size();
					++ii, ++count) {
				pNode pt = _f->at_1d(ii);
				if (nullptr != pt) {
					_idx = ii;
					return pt;
				}
			}
			return nullptr;
		}
		void _incr() {
			pNode end = _f->get_last_root_pNode();
			if (_ptr == end) {
				return;  //will not increase
			}
			_Ptr s = GetpNodeSiblingPlus(_ptr);
			if (s->father != nullptr) {
				_ptr = GetFirstLeaf(s);
			} else {
				if (s == end) {
					_ptr = s;
				} else {
					pNode root_pt = _incr_root(); //this will change _idx
					_ptr = GetFirstLeaf(root_pt);
				}
			}
		}
	public:
		bool operator==(const iterator_leaf_& _x) const {
			return _ptr == _x._ptr;
		}
		bool operator!=(const iterator_leaf_& _x) const {
			return _ptr != _x._ptr;
		}

		reference operator*() const {
			return (*_ptr);
		}

		pointer operator->() const {
			return &(operator*());
		}

		Self & operator++() {
			this->_incr();
			return *this;
		}

		Self operator++(int) {
			Self __tmp = *this;
			this->_incr();
			return __tmp;
		}

		bool is_exist() {
			return _ptr != nullptr;
		}

		pointer get_pointer() {
			return _ptr;
		}

		const pointer get_pointer() const {
			return _ptr;
		}
	};

	template<typename COV, typename V, St D, class _Ref, class _Ptr>
	class iterator_face_ {
	public:
		typedef COV cvt;
		typedef VALUE vt;
		typedef Node_<COV, V, D> Node;
		typedef Node_<COV, V, D> *pNode;
		typedef const Node_<COV, V, D> const_Node;
		typedef const Node_<COV, V, D>* const_pNode;

		//typedef Face_<COV, V, D> Face;
		//typedef Face_<COV, V, D> *pFace;
		//typedef const Face_<COV, V, D> const_Face;
		//typedef const Face_<COV, V, D>* const_pFace;

		typedef iterator_face_<COV, V, D, Node&, pNode> iterator;
		typedef iterator_face_<COV, V, D, const_Node&, const_pNode> const_iterator;
		typedef iterator_face_<COV, V, D, _Ref, _Ptr> Self;

		typedef Grid_<COV, V, D> Grid;
		typedef Grid_<COV, V, D>* pGrid;
		typedef const Grid_<COV, V, D>* const_pGrid;

		typedef Node value_type;
		typedef _Ptr pointer;
		typedef _Ref reference;

		const_pGrid _pg;
		St _idx;

		Face_<Node, _Ptr> _face;

		iterator_face_() {
			_face = nullptr;
			_pg = nullptr;
			_idx = 0;
		}
		iterator_face_(const_pGrid f, St idx, const Face_<Node, _Ptr>& ptr) :
				_pg(f), _idx(idx), _face(ptr) {
		}
		iterator_face_(const iterator& _x) :
				_pg(_x._pg), _idx(_x._idx), _face(_x._face) {
		}
	protected:
		pNode _incr_root() {
			St count = 0;
			for (St ii = (_idx + 1) % _pg->size(); count < _pg->size();
					++ii, ++count) {
				pNode pt = _pg->at_1d(ii);
				if (nullptr != pt) {
					_idx = ii;
					return pt;
				}
			}
			return nullptr;
		}

		Direction face_order_dir(St idx) const {
			ASSERT(idx < 5);
			const Direction ARR[] = { _XM_, _XP_, _YM_, _YP_, _ZM_, _ZP_ };
			return ARR[idx];
		}

		St face_order(const Direction& dir) const {
			ASSERT(IsFaceDirection(dir));
			switch (dir) {
			case _XM_:
				return 0;
				break;
			case _XP_:
				return 1;
				break;
			case _YM_:
				return 2;
				break;
			case _YP_:
				return 3;
				break;
			case _ZM_:
				return 4;
				break;
			case _ZP_:
				return 5;
				break;
			default:
				ASSERT_MSG(false, " Error Dirction");
				return 0;
			}
		}
		void _incr() {
			pNode end = _pg->get_last_root_pNode();
			if ((_face.pori()) == end) {
				return;  //will not increase
			}
			//face order
			St fo = face_order(_face.dir());
			if (fo < NumFaces - 1) {
				_face.dir() = face_order_dir(fo + 1);
				_face.pnei() = (_face.pori())->get_neighbor(_face.dir());
				_face.ft() = GetFaceType(_face.pori(), _face.pnei());
				return;
			} else {
				pNode s = GetpNodeSiblingPlus(_face.pori());
				if (s->father != nullptr) {
					_face.pori() = GetFirstLeaf(s);
				} else {
					if (s == end) {
						_face.pori() = s;
					} else {
						pNode root_pt = _incr_root(); //this will change _idx
						_face.pori() = GetFirstLeaf(root_pt);
					}
				}
				if (_face.pori() != end) {
					_face.dir() = _XM_;
					_face.pnei() = _face.pori()->get_neighbor(_XM_);
					_face.ft() = GetFaceType(_face.pori(), _face.pnei());
				} else {            //the last face-----------------
					_face.dir() = _XM_;
					_face.pnei() = nullptr;
					_face.ft() = _Boundary_;
				}
			}

		}
	public:
		bool operator==(const iterator_face_& _x) const {
			return _face == _x._face;
		}
		bool operator!=(const iterator_face_& _x) const {
			return _face != _x._face;
		}

		const Face_<Node, _Ptr>& operator*() const {
			return _face;
		}

		const Face_<Node, _Ptr>* operator->() const {
			return &(operator*());
		}

		Face_<Node, _Ptr>& operator*() {
			return _face;
		}

		Face_<Node, _Ptr>* operator->() {
			return &(operator*());
		}

		Self & operator++() {
			this->_incr();
			return *this;
		}

		Self operator++(int) {
			Self __tmp = *this;
			this->_incr();
			return __tmp;
		}

		bool is_exist() {
			return _face != nullptr;
		}

		pointer get_pointer() {
			return _face;
		}

		const pointer get_pointer() const {
			return _face;
		}
	};

public:
	/*
	 *  typedef of iterator leaf
	 */
	typedef iterator_leaf_<cvt, vt, Dim, Node&, Node*> iterator_leaf;
	typedef iterator_leaf_<cvt, vt, Dim, const Node&, const Node*> const_iterator_leaf;

	typedef iterator_face_<cvt, vt, Dim, Node&, Node*> iterator_face;
	typedef iterator_face_<cvt, vt, Dim, const Node&, Node*> const_iterator_face;

	iterator_leaf begin_leaf() {
		size_type idx = this->get_first_root_pNode_idx1d();
		pNode pt = this->get_first_root_pNode();
		pNode pn = GetFirstLeaf(pt);
		return iterator_leaf(this, idx, pn);
	}
	const_iterator_leaf begin_leaf() const {
		size_type idx = this->get_first_root_pNode_idx1d();
		const Node* pt = this->get_first_root_pNode();
		const Node* pn = GetFirstLeaf(pt);
		const_pSelf pg = this;
		return const_iterator_leaf(pg, idx, pn);
	}
	iterator_leaf last_leaf() {
		size_type idx = this->get_last_root_pNode_idx1d();
		pNode pt = this->get_last_root_pNode();
		pNode pn = GetFirstLeaf(pt);
		return iterator_leaf(this, idx, pn);
	}
	const_iterator_leaf last_leaf() const {
		size_type idx = this->get_last_root_pNode_idx1d();
		pNode pt = this->get_first_root_pNode();
		pNode pn = GetFirstLeaf(pt);
		return const_iterator_leaf(this, idx, pn);
	}
	iterator_leaf end_leaf() {
		size_type idx = this->get_last_root_pNode_idx1d();
		pNode pt = this->get_last_root_pNode();
		return iterator_leaf(this, idx, pt);
	}
	const_iterator_leaf end_leaf() const {
		size_type idx = this->get_last_root_pNode_idx1d();
		const pNode pt = this->get_last_root_pNode();
		return const_iterator_leaf(this, idx, pt);
	}

	iterator_face begin_face() {
		size_type idx = this->get_first_root_pNode_idx1d();
		pNode pt = this->get_first_root_pNode();
		pNode pn = GetFirstLeaf(pt);
		pNode pnei = pn->get_neighbor(_XM_);
		Face f(pn, pnei, _XM_, GetFaceType(pn, pnei));
		return iterator_face(this, idx, f);
	}
	const_iterator_face begin_face() const {
		size_type idx = this->get_first_root_pNode_idx1d();
		const Node* pt = this->get_first_root_pNode();
		const Node* pn = GetFirstLeaf(pt);
		const_pSelf pg = this;
		return const_iterator_leaf(pg, idx, pn);
	}

	iterator_face end_face() {
		size_type idx = this->get_last_root_pNode_idx1d();
		pNode pt = this->get_last_root_pNode();
		//pNode pnei = pn->get_neighbor(_XM_);
		Face f(pt, nullptr, _XM_, _Boundary_);
		return iterator_face(this, idx, f);
	}
	const_iterator_face end_face() const {
		size_type idx = this->get_last_root_pNode_idx1d();
		const pNode pt = this->get_last_root_pNode();
		return const_iterator_leaf(this, idx, pt);
	}

	pNode get_pnode(const cvt& x, const cvt& y, const cvt& z = 0.0) {
		for (const_iterator iter = this->begin(); iter != this->end(); ++iter) {
			pNode pn = (*iter);
			if (pn != nullptr) {
				if (pn->cell->is_in_on(x, y, z)) {
					return GetpNodeAt(pn, x, y, z);
				}
			}
		}
		return nullptr;
	}

	/*
	 * get leaf on Axes
	 */
	std::list<pNode> get_leaf(Axes aix, const cvt& x) {
		std::list<pNode> ret;
		for (const_iterator iter = this->begin(); iter != this->end(); ++iter) {
			pNode pn = (*iter);
			std::list<pNode> lp;
			if (pn != nullptr) {
				if (pn->cell->is_in_on(aix, x, _co_)) {
					lp = Node::GetLeaf(pn, aix, x);
				}
			}
			ret.merge(lp);
		}
		return ret;
	}
	/*
	 * show
	 */
	St count_empty() const {
		St res = 0;
		for (int i = 0; i < nodes.size(); i++) {
			if (nodes.at_1d(i) == nullptr) {
				res++;
			}
		}
		return res;
	}
	St count_non_empty() const {
		St res = 0;
		for (St i = 0; i < nodes.size(); i++) {
			if (nodes.at_1d(i) != nullptr) {
				res++;
			}
		}
		return res;
	}
	St count_leaf() const {
		St res = 0;
		for (St i = 0; i < nodes.size(); i++) {
			if (nodes.at_1d(i) != nullptr) {
				res += nodes.at_1d(i)->count_leaf();
			}
		}
		return res;
	}

	void show_info() const {
		// copy format
		std::ios oldState(nullptr);
		oldState.copyfmt(std::cout);
		//
		std::cout << "=>Grid Info: <========\n";
		std::cout << "Dim          :" << Dim << "\n";
		std::cout << "Size         :" << this->size_i() << " x "
				<< this->size_j();
		if (Dim == 3) {
			std::cout << " x " << this->size_k() << "\n";
		} else {
			std::cout << "\n";
		}
		std::cout << "Num of Tree space    :" << this->size() << "\n";
		std::cout << "Num of non empty Tree:" << count_non_empty() << "\n";
		std::cout << "Num of leaf          :" << this->count_leaf()
				<< std::endl;

		//Leaf information---------------------------
		_IF_TRUE_RETRUN(this->empty());

		pNode pt = this->get_first_root_pNode();
		_IF_TRUE_RETRUN(pt == nullptr);
		int _maxlevel = pt->max_level();
		arrayList_int num_node(_maxlevel + 1);
		arrayList_int num_leaf(_maxlevel + 1);
		for (St i = 0; i < size(); i++) {
			pNode tree = this->nodes.at_1d(i);
			if (tree != nullptr) {
				for (int il = 0; il <= _maxlevel; il++) {
					num_node[il] += tree->count_level(il);
					num_leaf[il] += tree->count_leaf_at_level(il);
				}
			}
		}
		std::cout << "level   Num   leaf   ratio%\n";
		int totalnode = 0;
		int totalleaf = 0;
		for (int i = 0; i <= _maxlevel; i++) {
			int nn = num_node[i];
			std::cout.flags(std::ios::right);
			std::cout.width(4);
			std::cout << i;
			std::cout.width(7);
			std::cout << nn;
			std::cout.width(7);
			int nl = num_leaf[i];
			std::cout << nl;
			std::cout.width(7);
			std::cout.precision(1);
			std::cout.setf(std::ios::fixed, std::ios::floatfield);
			std::cout << Float(nl) / Float(nn) * 100 << std::endl;
			totalnode += nn;
			totalleaf += nl;
		}
		std::cout.flags(std::ios::right);
		std::cout.width(11);
		std::cout << totalnode;
		std::cout.width(7);
		std::cout << totalleaf;
		std::cout.width(7);
		std::cout.precision(1);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout << Float(totalleaf) / Float(totalnode) * 100 << std::endl;

		// restore format
		std::cout.copyfmt(oldState);
	}
};

/*
 *  Function out of class =================================================
 */

}
#endif
