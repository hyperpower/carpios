#ifndef NODE_H_
#define NODE_H_

//#include "../typedefine.hpp"
#include "domain_define.hpp"
#include "cell.hpp"
#include "data.hpp"
#include "path.hpp"
#include <functional>

#include <math.h>

namespace carpio {

enum NodeIdx {
	//=========================
	//   y
	//   |
	//   ---------------
	//   |  PM  |  PP  |
	//   |  2   |  3   |
	//   |  WN  |  NE  |
	//   ---------------
	//   |  SW  |  SE  |
	//   |  0   |  1   |
	//   |  MM  |  MP  |
	//   ------------------->x
	//
	//   ---------------  ---------------
	//   |  MPM |  MPP |  |  PPM |  PPP |
	//   |  2   |  3   |  |  6   |  7   |
	//   |  WNB |  NEB |  |  WNF |  NEF |
	//   ---------------  ---------------
	//   |  SWB |  SEB |  |  SWP |  SEP |
	//   |  0   |  1   |  |  4   |  5   |
	//   |  MMM |  MMP |  |  PMM |  PMP |
	//   ---------------  ---------------
	//=========================

	//2D
	_MM_ = 0,
	_MP_ = 1,
	_PM_ = 2,
	_PP_ = 3,
	//3D
	_MMM_ = 0,
	_MMP_ = 1,
	_MPM_ = 2,
	_MPP_ = 3,
	_PMM_ = 4,
	_PMP_ = 5,
	_PPM_ = 6,
	_PPP_ = 7,
};

inline bool is_x_p(St i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 6) == 7;
}

inline bool is_x_m(St i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 6) == 6;
}

inline bool is_y_p(St i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 5) == 7;
}

inline bool is_y_m(St i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 5) == 5;
}

inline bool is_z_p(St i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 3) == 7;
}

inline bool is_z_m(St i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 3) == 3;
}

inline bool is_on_direction(St i, const Direction& dir) {
	ASSERT(i >= 0 && i < 8);
	unsigned short hi = HI(dir);
	unsigned short lo = LO(dir);
	return (hi & i) == (hi & lo);
}

enum NodeType {
	_Normal_ = 1 << 0, _Ghost_ = 1 << 1, _Cut_ = 1 << 2,
};

#define _TEMPLATE_COOV_V_DIM_ template<typename COO_VALUE, typename VALUE, St DIM>
#define _COOV_V_DIM_ COO_VALUE, VALUE, DIM

_TEMPLATE_COOV_V_DIM_ class Node_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;
	static const St NumChildren = NumVertexes;

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef VALUE& ref_vt;
	typedef const VALUE& const_ref_vt;
	typedef Node_<_COOV_V_DIM_> Self;
	typedef Node_<_COOV_V_DIM_> *pSelf;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell_<COO_VALUE, Dim> *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef const Data* const_pData;
	typedef Self Node;
	typedef Self *pNode;
	typedef Path_<Dim> Path;
	typedef const Self const_Node;
	typedef const Self* const_pNode;

	typedef void (*pFun)(pNode, utPointer);

	typedef void (*pFun_Conditional)(arrayList &, pNode, utPointer);

protected:
	//
		int _node_type;
		St _level;
		St _root_idx;
		St _idx;
		Path _path;
	public:
		pNode father;
		pNode child[NumChildren];
		pNode neighbor[NumNeighbors];
		pCell cell;
		pData data;

	protected:
		int _height(const pNode Current) const {
			if (Current == nullptr) {
				return 0;
			}
			if (!Current->has_child()) {
				return 0;
			} else {
				ArrayListV<St> arrh(NumChildren);
				for (St i = 0; i < this->NumChildren; ++i) {
					arrh = _height(Current->child[i]);
				}
				return 1 + arrh.max();
			}
		}

		template<class Ret1, class Ret2, class Args1, class Args2>
		void _traversal_conditional(pNode pn,
				std::function<Ret1(bool[], pNode, Args1)> fun_con,
				Args1 argsc,
				std::function<Ret2(pNode, Args2)> fun,
				Args2 args) {
			if (pn == nullptr) {
				return;
			} else {
				fun(pn, args);
				if (pn->has_child()) {
					bool ist[NumChildren];
					for (St i = 0; i < NumChildren; i++) {
						ist[i] = false;
					}
					fun_con(ist, pn, argsc);
					for (St i = 0; i < NumChildren; i++) {
						pNode c = pn->child[i];
						if (c != nullptr && ist[i]) {
							_traversal_conditional(c, fun_con, argsc, fun, args);
						}
					}
				}
			}
		}

		void _traversal_conditional(pNode pn, pFun_Conditional pfun_con,
				pFun pfun, utPointer utp) {
			if (pn == nullptr) {
				return;
			} else {
				(*pfun)(pn, utp);
				_IF_TRUE_RETRUN(pn==nullptr);
				if (pn->has_child()) {
					arrayList avt(NumChildren);
					pfun_con(avt, pn, utp);
					for (int i = 0; i < NumChildren; i++) {
						pNode c = pn->child[i];
						if (c != nullptr && avt[i] == 1) {
							_traversal_conditional(c, pfun_con, pfun, utp);
						}
					}
				}
			}
		}

		void _traversal(pNode pn, pFun pfun, utPointer utp) {
			if (pn == nullptr) {
				return;
			} else {
				(*pfun)(pn, utp);
				_IF_TRUE_RETRUN(pn==nullptr);
				if (pn->has_child()) {
					for (St i = 0; i < NumChildren; i++) {
						pNode c = pn->child[i];
						if (c != nullptr) {
							_traversal(c, pfun, utp);
						}
					}
				}
			}
		}

		template<class Ret, class Args>
		void _traversal(pNode pn, std::function<Ret(pNode&, Args)> fun, Args &args) {
			if (pn == nullptr) {
				return;
			} else {
				fun(pn, args);
				_IF_TRUE_RETRUN(pn==nullptr);
				if (pn->has_child()) {
					for (St i = 0; i < NumChildren; i++) {
						pNode c = pn->child[i];
						if (c != nullptr) {
							_traversal(c, fun, args);
						}
					}
				}
			}
		}

		void _traversal(pNode& pn, std::function<void(pNode&)> fun) {
			if (pn == nullptr) {
				return;
			} else {
				fun(pn);
				_IF_TRUE_RETRUN(pn==nullptr);
				if (pn->has_child()) {
					for (St i = 0; i < NumChildren; i++) {
						pNode c = pn->child[i];
						if (c != nullptr) {
							_traversal(c, fun);
						}
					}
				}
			}
		}

		void _traversal(const_pNode& pn, std::function<void(const_pNode&)> fun) const{
			if (pn == nullptr) {
				return;
			} else {
				fun(pn);
				_IF_TRUE_RETRUN(pn==nullptr);
				if (pn->has_child()) {
					for (St i = 0; i < NumChildren; i++) {
						const_pNode c = pn->child[i];
						if (c != nullptr) {
							_traversal(c, fun);
						}
					}
				}
			}
		}

		template<class Ret, class Args>
		void _traversal(const_pNode pn, std::function<Ret(const_pNode, Args)> fun, Args &args) const {
			if (pn == nullptr) {
				return;
			} else {
				fun(pn, args);
				_IF_TRUE_RETRUN(pn==nullptr);
				if (pn->has_child()) {
					for (St i = 0; i < NumChildren; i++) {
						pNode c = pn->child[i];
						if (c != nullptr) {
							_traversal(c, fun, args);
						}
					}
				}
			}
		}

	public:
		/*
		 *  constructor
		 */
		Node_(pNode f, int nt, St level, St root_idx, St c_idx, Path _p,  //
				const vt &x, const vt &dhx,//
				const vt &y = 0.0, const vt &dhy = 0.0,//
				const vt &z = 0.0, const vt &dhz = 0.0) {
			_node_type = nt;
			_level = level;
			cell = new Cell(x, dhx, y, dhy, z, dhz);
			father = f;
			_root_idx = root_idx;
			_idx = c_idx;
			_path = _p;

			data = nullptr;
			for (St i = 0; i < this->NumChildren; i++) {
				child[i] = nullptr;
			}
			for (St i = 0; i < this->NumNeighbors; i++) {
				neighbor[i] = nullptr;
			}
		}
		Node_(pNode f, int nt, St level, St root_idx, St child_idx, Path _p, //
				const Cell& c) {
			_node_type = nt;
			_level = level;
			cell = new Cell(c);
			father = f;
			_root_idx = root_idx;
			_idx = child_idx;
			_path = _p;

			data = nullptr;
			for (St i = 0; i < this->NumChildren; i++) {
				child[i] = nullptr;
			}
			for (St i = 0; i < this->NumNeighbors; i++) {
				neighbor[i] = nullptr;
			}
		}
		/*
		 *  delete
		 */
	protected:
		/*
		 *  before using this function, making sure that this node is a leaf
		 */
		void _DeleteLeaf() {
			pNode f = this->father;
			if (f != nullptr) {
				f->child[_idx] = nullptr;
			}
			delete cell;
			if (data != nullptr) {
				delete data;
			}
		}

		void _Delete(pNode pn) {
			if (pn == nullptr) {
				return;
			}
			if (pn->has_child()) {
				for (St i = 0; i < NumChildren; i++) {
					pNode ch = pn->child[i];
					if (ch != nullptr) {
						_Delete(ch);
					}
				}
			} else { // is leaf
				pn->_DeleteLeaf();
			}
		}

	public:
		~Node_() {
			_Delete(this);
		}

		/*
		 * type
		 */
		inline int get_type() const {
			return _node_type;
		}

		inline void set_type(int type) {
			_node_type = type;
		}

		inline St get_level() const {
			return _level;
		}

		inline St get_idx() const {
			//return (_path >> int(pow(Dim, _level))) & (NumVertexes - 1);
			return _idx;
		}
		inline void set_idx(St i) {
			_idx = i;
		}

		inline Path get_path() const {
			return _path;
		}

		inline St get_root_idx() const {
			return _root_idx;
		}

		inline St height() const {
			return this->_height(this);
		}

		/*
		 *  child
		 */
		inline bool has_child() const {
			for (St i = 0; i < this->NumChildren; ++i) {
				if (this->child[i] != nullptr) {
					return true;
				}
			}
			return false;
		}
		inline bool has_child(St idx) const {
			ASSERT(idx < this->NumChildren);
			if (this->child[idx] != nullptr) {
				return true;
			}
			return false;
		}

		inline bool is_leaf() const {
			return !has_child();
		}

		inline bool is_root() const {
			if (this->father == nullptr) {
				return true;
			} else {
				return false;
			}
		}

		inline bool is_full_child() const {
			bool res = this->child[0] != nullptr;
			for (St i = 1; i < this->NumChildren; ++i) {
				res = res && (this->child[i] != nullptr);
			}
			return res;
		}

		inline bool is_alone() const {
			return this->is_leaf()&&this->is_root();
		}

		/*
		 *  count
		 */
		inline St count_children() const {
			St res = 0;
			for (St i = 0; i < this->NumChildren; ++i) {
				res += (this->child[i] != nullptr) ? 1 : 0;
			}
			return res;
		}

		inline St count_all() const {
			St res = 0;
			std::function<void(const_pNode, St)> fun = [&res](const_pNode pn, St dummy) {
				res++;
			};
			this->traversal(fun, res);
			return res;
		}

		inline St count_leaf() const {
			St res = 0;
			std::function<void(const_pNode, St)> fun = [&res](const_pNode pn, St dummy) {
				if(pn->is_leaf()) {
					res++;
				}
			};
			this->traversal(fun, res);
			return res;
		}

		inline St count_level(St le) const {
			St res = 0;
			std::function<void(const_pNode, St)> fun = [&res](const_pNode pn, St le) {
				if(pn->get_level()==le) {
					res++;
				}
			};
			this->traversal(fun, le);
			return res;
		}

		inline St count_leaf_at_level(St le) const {
			St res = 0;
			std::function<void(const_pNode, St)> fun = [&res](const_pNode pn, St le) {
				if(pn->get_level()==le && pn->is_leaf()) {
					res++;
				}
			};
			this->traversal(fun, le);
			return res;
		}
		/*
		 * Connect
		 */

		void connect_nodes() {
			std::function<void(pNode&, St)> fun = [&](pNode& pn, St le) {
				if(!pn->is_root()) {
					//set neighbor
					pn->neighbor[0] = pn->get_neighbor(_XM_);
					pn->neighbor[1] = pn->get_neighbor(_XP_);
					pn->neighbor[2] = pn->get_neighbor(_YM_);
					pn->neighbor[3] = pn->get_neighbor(_YP_);
					if(Dim == 3) {
						pn->neighbor[4] = pn->get_neighbor(_ZM_);
						pn->neighbor[5] = pn->get_neighbor(_ZP_);
					}
				}
			};
			St dummy =0;
			this->traversal(fun, dummy);
		}

		/*
		 *  new
		 */
	protected:
		Path _cal_this_path(St i) {
			Path res(Dim);
			res.clear();
			if(is_x_p(i)&& Dim>=1) {
				res.set(0);
			}
			if(is_y_p(i)&& Dim>=2) {
				res.set(1);
			}
			if(is_z_p(i)&& Dim>=3) {
				res.set(2);
			}
			return res;
		}
	public:
		void new_full_child() {
			if (!has_child()) {
				St ltmp = _level + 1;
				vt nhdx = this->cell->get_hd(_X_) * 0.5;
				vt nhdy = this->cell->get_hd(_Y_) * 0.5;
				vt nhdz = this->cell->get_hd(_Z_) * 0.5;
				vt cx = this->cell->get(_C_, _X_);
				vt cy = this->cell->get(_C_, _Y_);
				vt cz = this->cell->get(_C_, _Z_);
				for (St i = 0; i < this->NumChildren; ++i) {
					pNode f = this;
					int nt = 1;
					St l = ltmp;
					St ridx = _root_idx;
					St child_idx = i;
					Path p;
					if(!this->is_root()) {
						p = this->_path;
					}
					p.append(this->_cal_this_path(i));
					this->child[i] = new Node_(		//
							f, nt, l, ridx, child_idx, p,//
							cx + (is_x_p(i) ? nhdx : -nhdx), nhdx,//
							cy + (is_y_p(i) ? nhdx : -nhdx), nhdy,//
							cz + (is_z_p(i) ? nhdx : -nhdx), nhdz);
				}
			}
		}
		void new_child(St idx) {
			ASSERT(idx < this->NumChildren);
			if (!has_child()) {
				St ltmp = _level + 1;
				vt nhdx = this->cell->get_hd(_X_) * 0.5;
				vt nhdy = this->cell->get_hd(_Y_) * 0.5;
				vt nhdz = this->cell->get_hd(_Z_) * 0.5;
				vt cx = this->cell->get(_C_, _X_);
				vt cy = this->cell->get(_C_, _Y_);
				vt cz = this->cell->get(_C_, _Z_);
				pNode f = this;
				int nt = 1;
				St l = ltmp;
				St ridx = _root_idx;
				St npath = idx;
				Path p;
				if(!this->is_root()) {
					p = this->_path;
				}
				p.append(this->_cal_this_path(idx));
				this->child[idx] = new Node_( //
						f, nt, l, ridx, npath,//
						cx + (is_x_p(idx) ? nhdx : -nhdx), nhdx,//
						cy + (is_y_p(idx) ? nhdx : -nhdx), nhdy,//
						cz + (is_z_p(idx) ? nhdx : -nhdx), nhdz);
			}
		}

		/*
		 *  make sure point is in this node
		 */
		St which_child(const COO_VALUE &x,
				const COO_VALUE &y = 0,
				const COO_VALUE &z = 0) {
			St idx = 0;
			if (Dim == 3) {
				idx = IsInRange(this->cell->get(_M_, _Z_), z, this->cell->get(_C_, _Z_), _cc_) ? 0 : 1;
			}
			idx = idx << 1;
			if (Dim >= 2) {
				idx = idx + (IsInRange(this->cell->get(_M_, _Y_), y, this->cell->get(_C_, _Y_), _cc_) ? 0 : 1);
			}
			idx = idx << 1;
			idx = idx + (IsInRange(this->cell->get(_M_, _X_), x, this->cell->get(_C_, _X_), _cc_) ? 0 : 1);
			return idx;
		}

		/*
		 *  neighbor find
		 */
		void set_neighbor(
				pNode xm, pNode xp, //x
				pNode ym, pNode yp,//y
				pNode zm, pNode zp) { //z
			//       yp 3
			//      ______
			//     |      |
			//xm 0 |      | xp 1
			//     |______|
			//       ym 2
			neighbor[0] = xm;
			neighbor[1] = xp;
			if (Dim >= 2) {
				neighbor[2] = ym;
				neighbor[3] = yp;
			}
			if (Dim == 3) {
				neighbor[4] = zm;
				neighbor[5] = zp;
			}
		}
		void set_neighbor(pNode pn, Direction d) {
			switch (d) {
				case _XM_: {
					this->neighbor[0] = pn;
					break;
				}
				case _XP_: {
					this->neighbor[1] = pn;
					break;
				}
				case _YM_: {
					this->neighbor[2] = pn;
					break;
				}
				case _YP_: {
					this->neighbor[3] = pn;
					break;
				}
				case _ZM_: {
					ASSERT(Dim==3);
					this->neighbor[4] = pn;
					break;
				}
				case _ZP_: {
					ASSERT(Dim==3);
					this->neighbor[5] = pn;
					break;
				}
			}
		}

		inline bool is_adjacent(const Direction &d) const {
			// Direction on x y or z
			unsigned short hi = d >> 3;
			return ((hi & _idx) ^ (hi & d)) == 0;
		}

		inline St reflect(const Direction &d) const {
			// Direction on x y or z
			return _idx ^ (d >> 3);
		}

		inline bool has_diagonal_sibling(const Direction &d) const {
			unsigned short hi = d >> 3;
			return ((_idx ^ hi) & hi) == (LO(d) & hi);
		}

		inline bool is_out_corner(const Direction &d) const {
			return (_idx & HI(d)) == (HI(d) & LO(d));
		}

		inline Direction out_common_direction(const Direction &d) const {
			// return direction on x y or z
			static const unsigned short COMMON_AXES[4][4] = { {0, 2, 1, 0},
				{	2, 0, 0, 1},
				{	1, 0, 0, 2},
				{	0, 1, 2, 0}};
			std::function<unsigned short(unsigned short)> to2bit = [](unsigned short num) {
				unsigned short z = (num & 4) >> 2;
				return ((num & 1) << 1) + z;
			};
			std::function<unsigned short(unsigned short)> to3bit = [](unsigned short num) {
				unsigned short z = (num & 1) << 2;
				return z + ((num & 2) >> 1);
			};
			unsigned short hi = d >> 3;
			unsigned short lo = LO(d) & hi;
			unsigned short id = _idx & hi;
			unsigned short com;
			switch (hi) {
				case 3: //xy
				com = COMMON_AXES[id][lo];
				break;
				case 6://yz
				com = COMMON_AXES[id >> 1][lo >> 1];
				com = com << 1;
				break;
				case 5: { //zx
					com = COMMON_AXES[to2bit(id)][to2bit(lo)];
					com = to3bit(com);
					break;
				}
				default:
				ASSERT(false);
			}
			return (com << 3) + (lo & com);
		}

		inline St diagonal_idx(Direction d) const {
			return _idx ^ HI(d);
		}

	protected:
		pNode _get_root_neighbor_step(const Direction &d1, const Direction &d2) const {
			pNode pstep1 = this->get_root_neighbor(d1);
			if (pstep1 == nullptr) {
				return nullptr;
			} else {
				return pstep1->get_root_neighbor(d2);
			}
		}

		pNode _get_root_neighbor_step(const Direction &d1,
				const Direction &d2,
				const Direction &d3) const {
			pNode pstep = this->get_root_neighbor(d1);
			if (pstep == nullptr) {
				return nullptr;
			} else {
				pstep = pstep->get_root_neighbor(d2);
				if (pstep == nullptr) {
					return nullptr;
				} else {
					return pstep->get_root_neighbor(d3);
				}
			}
		}

		pNode _get_root_neighbor_path(const Direction &d1,
				const Direction &d2) const {
			pNode pt = this->_get_root_neighbor_step(d1, d2);
			if (pt != nullptr) {
				return pt;
			} else {
				return this->_get_root_neighbor_step(d2, d1);
			}
		}

		pNode _get_root_neighbor_path(const Direction &d1,
				const Direction &d2,
				const Direction &d3) const {
			pNode pt = this->_get_root_neighbor_step(d1, d2, d3);
			if (pt != nullptr) {
				return pt;
			}
			pt = this->_get_root_neighbor_step(d1, d3, d2);
			if (pt != nullptr) {
				return pt;
			}
			pt = _get_root_neighbor_step(d2, d1, d3);
			if (pt != nullptr) {
				return pt;
			}
			pt = _get_root_neighbor_step(d2, d3, d1);
			if (pt != nullptr) {
				return pt;
			}
			pt = _get_root_neighbor_step(d3, d1, d2);
			if (pt != nullptr) {
				return pt;
			}
			pt = _get_root_neighbor_step(d3, d2, d1);
			if (pt != nullptr) {
				return pt;
			}
			return nullptr;
		}

	public:
		pNode get_root_neighbor(const Direction &d) const {
			unsigned short hi = HI(d);
			ASSERT(hi != 0);
			switch (HI(d)) {
				case 1:   //x
				return this->neighbor[hi & LO(d)];
				case 2://y
				return this->neighbor[((hi & LO(d)) >> 1) + hi];
				case 3://xy
				return _get_root_neighbor_path(8 + LO(d), 16 + LO(d));
				case 4://z
				return this->neighbor[((hi & LO(d)) >> 2) + hi];
				case 5:
				return _get_root_neighbor_path(32 + LO(d), 8 + LO(d));
				case 6:
				return _get_root_neighbor_path(16 + LO(d), 32 + LO(d));
				case 7:
				return _get_root_neighbor_path(8 + LO(d), 16 + LO(d), 32 + LO(d));
				default:
				return nullptr;
			}
		}

		pNode get_adj_neighbor(const pNode Current, Direction d) {
			// face direction
			//std::cout << "--------" << "\n";
			//std::cout << "d n     " << d << "\n";
			//std::cout << "idx adj " << Current->get_idx() << "\n";
			pNode ca = nullptr;//common ancestor
			if (Current->father != nullptr
					&& Current->is_adjacent(d)) {
				ca = get_adj_neighbor(Current->father, d);
			} else {
				ca = Current->father;
			}
			pNode pt = nullptr;
			if (ca != nullptr && ca->has_child()) {
				pt = ca->child[Current->reflect(d)];
			} else if (ca == nullptr) {
				pt = Current->get_root_neighbor(d);
			} else {
				pt = ca;
			}
			return pt;
		}

		pNode get_cor_neighbor(const pNode Current, Direction d) {
			//std::cout << "--------" << "\n";
			//std::cout << "d n     " << d << "\n";
			//std::cout << "idx cor " << Current->get_idx() << "\n";
			pNode ca = nullptr;//common ancestor
			int flag = 0;
			if (Current->father != nullptr &&
					!Current->has_diagonal_sibling(d)) {
				//Find a common ancestor
				if (Current->is_out_corner(d)) {
					//	std::cout << "cor " << Current->get_idx() << "\n";
					ca = get_cor_neighbor(Current->father, d);
				} else {
					//	std::cout << "adj " << Current->get_idx() << "\n";
					ca = get_adj_neighbor(Current->father,
							Current->out_common_direction(d));
				}
			} else {
				std::cout << "dia " << Current->get_idx() << "\n";
				flag = 1;
				ca = Current->father;
			}
			//Follow opposite path to locate the neighbor
			pNode pt = nullptr;
			if (ca != nullptr && ca->has_child()) {
				pt = ca->child[Current->diagonal_idx(d)];
			} else if (ca == nullptr && flag == 1) {
				pt = Current->get_root_neighbor(d);
			} else {
				pt = ca;
			}
			return pt;
		}

		pNode get_neighbor_adj_cor(const pNode Current, //
				Direction d) {
			if (IsFaceDirection(d)) {
				return get_adj_neighbor(Current, d);
			}
			if (IsCornerDirection(d)) {
				return get_cor_neighbor(Current, d);
			}
			return nullptr;
		}

		pNode get_cor_neighbor_xyz(const pNode Current, Direction d) {
			//std::cout << "xyz " << Current->get_idx() << "\n";
			pNode ca = nullptr;//common ancestor
			if (Current->father != nullptr &&
					Current->is_out_corner(d)) {
				ca = get_cor_neighbor_xyz(Current->father, d);
			} else {
				Direction idx = Current->get_idx();
				if (idx==((~LO(d))&7)) {
					ca = Current->father;
				} else {
					Direction nd = (((~(idx ^ (d & 7))) & 7) << 3) + idx;
					ca = get_neighbor_adj_cor(Current->father, nd);
				}
			}
			pNode pt = nullptr;
			if (ca != nullptr && ca->has_child()) {
				pt = ca->child[Current->diagonal_idx(d)];
			} else if (ca == nullptr) {
				pt = Current->get_root_neighbor(d);
			} else {
				pt = ca;
			}
			return pt;
		}

		pNode get_neighbor(const pNode Current, Direction d) {
			ASSERT(d > 7);
			Direction nd = d & 63;
			if ( IsXYZDirection(nd)) {
				return get_cor_neighbor_xyz(Current, nd);
			} else {
				return get_neighbor_adj_cor(Current, nd);
			}
		}

		pNode get_neighbor(Direction d) {
			return get_neighbor(this, d);
		}
		pNode get_neighbor_adj_fast(const pNode Current, Direction d) {
			pNode res=nullptr;
			switch(d) {
				case _XM_: {res = Current->neighbor[0];break;}
				case _XP_: {res = Current->neighbor[1];break;}
				case _YM_: {ASSERT(Dim>=2);res = Current->neighbor[2];break;}
				case _YP_: {ASSERT(Dim>=2);res = Current->neighbor[3];break;}
				case _ZM_: {ASSERT(Dim>=3);res = Current->neighbor[4];break;}
				case _ZP_: {ASSERT(Dim>=3);res = Current->neighbor[5];break;}
				default:
				res = nullptr;break;
			}
			return res;
		}
		pNode get_neighbor_adj_cor_fast(const pNode Current, Direction d) {
			if (IsFaceDirection(d)) {
				return get_neighbor_adj_fast(Current, d);
			}
			if (IsCornerDirection(d)) {
				return get_cor_neighbor(Current, d);
			}
			return nullptr;
		}
		pNode get_neighbor_fast(const pNode Current, Direction d) {
			ASSERT(d > 7);
			Direction nd = d & 63;
			if ( IsXYZDirection(nd)) {
				return get_cor_neighbor_xyz(Current, nd);
			} else {
				return get_neighbor_adj_cor_fast(Current, nd);
			}
		}
		pNode get_neighbor_fast(Direction d) {
			return get_neighbor_fast(this, d);
		}
		/*
		 *  Traverse
		 */
		void traversal(pFun pfun, utPointer utp) {
			this->_traversal(this, pfun, utp);
		}

		template<class Args>
		void traversal(std::function<void(pNode&, Args)> fun, Args &args) {
			this->_traversal(this, fun, args);
		}

		template<class Args>
		void traversal(std::function<void(const_pNode, Args)> fun, Args &args) const {
			this->_traversal(this, fun, args);
		}

		template<class Ret1, class Ret2, class Args1, class Args2>
		void traversal_conditional(  //
				std::function<Ret1(bool[], pNode, Args1)> fun_con,//
				Args1 argsc,//
				std::function<Ret2(pNode, Args2)> fun,//
				Args2 args) {
			this->_traversal_conditional(this, fun_con, argsc, fun, args);
		}

		St max_level() const {
			St res = 0;
			std::function<void(const_pNode, St)> fun =
			[&res](const_pNode pn, St dummy) {
				if (pn->get_level()>res) {
					res = pn->get_level();
				}
			};
			St dummy;
			this->_traversal(this, fun, dummy);
			return res;
		}

		/*
		 *  overload the function of cell
		 */
		bool is_in_on(const vt x,
				const vt y = 0,
				const vt z = 0) const {
			bool res = this->cell->is_in_on(x,y,z);
			if(this->is_leaf()||res == false) {
				return res;
			}
			std::function<void(const_pNode, int)> fun =
			[&res,&x,&y,&z](const_pNode pn, int dummy) {
				if (pn->is_leaf() && res==false) {
					res = pn->cell->is_in_on(x,y,z);
				}
			};
			int dummy = 1;
			this->_traversal(this, fun, dummy);
			return res;
		}

		cvt cp(Axes axes) const {  //center point
			if( Dim <= 2 && axes == _Z_) {
				return 0.0;
			}
			return this->cell->get(_C_, axes);
		}
		cvt d(Axes axes) const {  //d
			return this->cell->get_d(axes);
		}
		cvt hd(Axes axes) const {  //d
			return this->cell->get_hd(axes);
		}
		cvt p(Orientation ori, Axes axes) const {  //point
			if( Dim <= 2 && axes == _Z_) {
				return 0.0;
			}
			return this->cell->get(ori, axes);
		}
		cvt p(Direction dir, Axes axes) const {
			if( Dim <= 2 && axes == _Z_) {
				return 0.0;
			}
			return this->cell->get(ToOrientation(dir, axes), axes);
		}
		cvt face_area(Direction dir) const {
			ASSERT(IsFaceDirection(dir));
			Axes a = FaceDirectionToAxes(dir);
			cvt a_f = 1.0;
			if (Dim == 2) {
				Axes dp = VerticalAxes2D(a);
				a_f = this->d(dp);
			} else {  //3d
				Axes dp1 = VerticalAxes1(a);
				Axes dp2 = VerticalAxes2(a);
				a_f = this->d(dp1) * this->d(dp2);
			}
			return a_f;
		}
		cvt volume() const {
			return this->cell->volume();
		}
		/*
		 *  data
		 *  overload data functions
		 */
		void new_data(const St& nc, const St& nf, const St& nv, const St& nutp) {
			if(this->data ==nullptr) {
				this->data = new Data(nc,nf,nv,nutp);
			} else {
				this->data->reconstruct(nc,nf,nv,nutp);
			}
		}
		void resize_data(const St& nc, const St& nf, const St& nv, const St& nutp) {
			if(this->data ==nullptr) {
				this->data = new Data(nc,nf,nv,nutp);
			} else {
				this->data->resize(nc,nf,nv,nutp);
			}
		}
		pData pdata() {
			ASSERT(this->data != nullptr);
			return this->data;
		}
		const_pData pdata() const {
			ASSERT(this->data != nullptr);
			return this->data->idx();
		}
		const int& d_idx() const {
			ASSERT(this->data != nullptr);
			return this->data->idx();
		}
		int& d_idx() {
			ASSERT(this->data != nullptr);
			return this->data->idx();
		}
		utPointer& utp(St i) {
			ASSERT(this->data != nullptr);
			return this->data->utp(i);
		}
		const_utPointer& utp(St i) const {
			ASSERT(this->data != nullptr);
			return this->data->utp(i);
		}
		St size_cd() const {
			ASSERT(this->data != nullptr);
			return this->data->size();
		}
		ref_vt cd(St i) { //center data
			ASSERT(this->data != nullptr);
			return this->data->center(i);
		}
		const_ref_vt cd(St i) const { //center data
			ASSERT(this->data != nullptr);
			return this->data->center(i);
		}
		vt cda(St i) const {
			if(this->data != nullptr) {
				return this->data->center(i);
			} else {
				// averange data from children
				bool nullflag = true;
				St count = 0;
				vt sum = 0.0;
				std::function<void(const_pNode, St)> fun = [&sum, &nullflag, &count](const_pNode pn, St i) {
					if(pn->is_leaf()) {
						if(pn->data != nullptr) {
							if(pn->data->has_center(i)) {
								nullflag = false;
								count ++;
								sum= sum+ pn->cd(i);
							}
						}
					}
				};
				this->traversal(fun,i);
				ASSERT_MSG(nullflag==false, " >! All children null");
				return sum/count;
			}
		}
		vt cdva(St i) const { //cell volume weight average
			if(this->data != nullptr) {
				return this->data->center(i);
			} else {
				// averange data from children
				bool nullflag = true;
				vt sum = 0.0;
				std::function<void(const_pNode, St)> fun = [&sum, &nullflag](const_pNode pn, St i) {
					if(pn->is_leaf()) {
						if(pn->data != nullptr) {
							if(pn->data->has_center(i)) {
								nullflag = false;
								sum+= pn->cd(i) * pn->cell->volume();
							}
						}
					}
				};
				this->traversal(fun,i);
				ASSERT_MSG(nullflag==false, " >! All children null");
				return sum/this->cell->volume();
			}
		}
		/*
		 *  show
		 */
		void show(const std::string& name = "" ) const {
			std::cout << "Node --- "<<name<<"\n";
			std::cout << "Dimension = " << Dim << "D\n";
			std::cout << "level      :" << this->_level << "\n";
			std::cout << "node type  :" << this->_node_type << "\n";
			std::cout << "idx        :" << this->_idx << "\n";
			std::cout << "path       :";
			this->_path.show(1);
			std::cout << "\n";
			std::cout << "CELL  show =========\n";
			std::cout << std::scientific;
			std::cout << "       min    " <<"     max    "<<"     d    "<< "     c    \n";
			std::cout << "x:" << p(_M_, _X_)<<" "<< p(_P_, _X_) <<" "<<d(_X_)<<" "<<cp(_X_)<< "\n";
			std::cout << "y:" << p(_M_, _Y_)<<" "<< p(_P_, _Y_) <<" "<<d(_Y_)<<" "<<cp(_Y_)<< "\n";
			if (Dim == 3) {
				std::cout << "z:" << p(_M_, _Z_)<<" "<< p(_P_, _Z_) <<" "<<d(_Z_)<<" "<<cp(_Z_)<< "\n";
			}
			std::cout << "DATA  show =========\n";
			if (nullptr == this->data) {
				std::cout << "No data \n";
			} else {
				//this->data->show_info();
			}
		}
		/*
		 * static function
		 */

		static pNode GetChild(pNode p, const Direction& dir, int assert_flag=1) {
			if(assert_flag == 1) {
				ASSERT(p!=nullptr);
				if(Dim == 2) {
					ASSERT(IsCornerDirection(dir));
					ASSERT(IsDirectionOn(dir, _X_) &&IsDirectionOn(dir, _Y_));
				} else { //dim =1
					ASSERT(IsVertexDirection(dir));
				}
			}
			unsigned short lo = LO(dir);
			return p->child[lo];
		}
		static pNode GetChild_CornerLeaf(pNode p, const Direction& dir) {
			ASSERT(p!=nullptr);
			if(Dim == 2) {
				ASSERT(IsCornerDirection(dir));
				ASSERT(IsDirectionOn(dir, _X_) &&IsDirectionOn(dir, _Y_));
			} else { //dim =1
				ASSERT(IsVertexDirection(dir));
			}
			pNode pc = GetChild(p, dir);
			while(!(pc==nullptr || pc->is_leaf())) {
				pc = GetChild(pc, dir, 0);   //without assert
			}
			return pc;
		}
		static std::list<pNode> GetLeaf(pNode pn, Axes a, cvt cor) {
			ASSERT(pn != nullptr);
			if (!pn->cell->is_in_on(a, cor, _co_)) {
				return std::list<pNode>();
			}
			std::function<void(bool[], pNode, cvt)> condition =
			[a](bool arr[], pNode pn, cvt cor ) {
				for(St i = 0; i<NumChildren;++i) {
					arr[i] = pn->child[i]->cell->is_in_on(a, cor,_co_);
				}
			};
			std::list<pNode> ret;
			std::function<void(pNode, int)> fun =
			[&ret](pNode pn, int dummy) {
				if (pn->is_leaf()) {
					ret.push_back(pn);
				}
			};
			//------------------
			int dummy = 0;
			pn->traversal_conditional(condition, cor , fun, dummy);
			return ret;
		}
	}
	;

	/*
	 *  functions out of class ================================================
	 */
	_TEMPLATE_COOV_V_DIM_ int GetDataIdx(const Node_<COO_VALUE, VALUE, DIM> *pn) {
		ASSERT(pn != nullptr);
		return pn->data->get_idx();
	}

	_TEMPLATE_COOV_V_DIM_
	Node_<COO_VALUE, VALUE
, DIM> *
GetpNodeAt(Node_<COO_VALUE, VALUE, DIM> *pn,
		const COO_VALUE &x,
		const COO_VALUE &y = 0,
		const COO_VALUE &z = 0) {
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	ASSERT(pn != nullptr);
	if (!pn->cell->is_in_on(x, y, z)) {
		return nullptr;
	}
	std::function<void(bool[], pNode, Float[])> condition =
	[](bool arr[], pNode pn, Float point[]) {
		St idx = pn->which_child(point[0],
				(pn->Dim >= 2) ? point[1] : 0,
				(pn->Dim == 3) ? point[2] : 0);
		arr[idx] = true;
	};
	pNode ret = nullptr;
	std::function<void(pNode, int)> fun =
	[&ret](pNode pn, int dummy) {
		if (pn->is_leaf()) {
			ret = pn;
		}
	};
	Float point[pn->Dim];
	point[0] = x;
	if (pn->Dim >= 2) {point[1] = y;}
	if (pn->Dim == 3) {point[2] = z;}
	int dummy = 0;
	pn->traversal_conditional(condition, point, fun, dummy);
	return ret;
}

_TEMPLATE_COOV_V_DIM_ const Node_<COO_VALUE, VALUE, DIM> *
GetpNodeSiblingPlus(const Node_<COO_VALUE, VALUE, DIM> *p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* f = p->father;
	if (f == nullptr) {
		return p;
	}
	for (St i = p->get_idx() + 1; i < Node::NumChildren; ++i) {
		Node* c = f->child[i];
		if (c != nullptr) {
			return c;
		}
	}
	return GetpNodeSiblingPlus(f);
}

_TEMPLATE_COOV_V_DIM_ Node_<COO_VALUE, VALUE, DIM> *
GetpNodeSiblingPlus(Node_<COO_VALUE, VALUE, DIM> *p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* f = p->father;
	if (f == nullptr) {
		return p;
	}
	for (St i = p->get_idx() + 1; i < Node::NumChildren; ++i) {
		Node* c = f->child[i];
		if (c != nullptr) {
			return c;
		}
	}
	return GetpNodeSiblingPlus(f);
}

_TEMPLATE_COOV_V_DIM_ Node_<COO_VALUE, VALUE, DIM> *
GetFirstChild(Node_<COO_VALUE, VALUE, DIM> * p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* c = nullptr;
	if(p==nullptr) {
		return c;
	}
	for (St i = 0; i < Node::NumChildren; ++i) {
		c = p->child[i];
		if (c != nullptr) {
			return c;
		}
	}
	return c;
}

_TEMPLATE_COOV_V_DIM_ const Node_<COO_VALUE, VALUE, DIM> *
GetFirstChild(const Node_<COO_VALUE, VALUE, DIM> * p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* c = nullptr;
	if(p==nullptr) {
		return c;
	}
	for (St i = 0; i < Node::NumChildren; ++i) {
		c = p->child[i];
		if (c != nullptr) {
			return c;
		}
	}
	return c;
}

_TEMPLATE_COOV_V_DIM_ Node_<COO_VALUE, VALUE, DIM> *
GetFirstLeaf(Node_<COO_VALUE, VALUE, DIM> * p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* c = GetFirstChild(p);
	if (c == nullptr) {
		return p;
	} else {
		Node* resc = c;
		while (c != nullptr) {
			resc = c;
			c = GetFirstChild(c);
		}
		return resc;
	}
}

_TEMPLATE_COOV_V_DIM_ const Node_<COO_VALUE, VALUE, DIM> *
GetFirstLeaf(const Node_<COO_VALUE, VALUE, DIM> * p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	const Node* c = GetFirstChild(p);
	if (c == nullptr) {
		return p;
	} else {
		const Node* resc = c;
		while (c != nullptr) {
			resc = c;
			c = GetFirstChild(c);
		}
		return resc;
	}
}

/*
 * face
 */
enum FaceType {
	_Null_ = -1, _Boundary_ = 0, _Equal_ = 1, _FineCoarse_ = 2, _CoarseFine_ = 3,
};

inline std::string ParseFaceType(const FaceType& t) {
	switch (t) {
	case _Boundary_:
		return "Boundary";
		break;
	case _Equal_:
		return "Equal";
		break;
	case _FineCoarse_:
		return "FineCorase";
		break;
	case _CoarseFine_:
		return "CoraseFine";
		break;
	default:
		break;
	}
	return "Null";
}

template<typename COO_VALUE, typename VALUE, St DIM>
FaceType GetFaceType(  //
		const Node_<_COOV_V_DIM_>* p,     //main node
const Node_<_COOV_V_DIM_>* pn) {
	if (p == nullptr) {
		return _Null_;
	}
	if (pn == nullptr) {
		return _Boundary_;
	}
	if (pn->get_type() == _Ghost_) {
		return _Boundary_;
	}
	if (p->get_level() > pn->get_level())
	return _FineCoarse_;
	if (p->get_level() == pn->get_level() && !pn->has_child())
	return _Equal_;
	return _CoarseFine_;
}

template<typename NODE, typename PNODE>
class Face_ {
public:
	typedef NODE Node;
	typedef PNODE pNode;

	typedef typename Node::vt vt;

	static const St Dim = Node::Dim;
	static const St NumFaces = Node::NumFaces;
	static const St NumVertexes = Node::NumVertexes;
	static const St NumNeighbors = Node::NumNeighbors;
	static const St NumChildren = Node::NumChildren;

	typedef Face_<Node, pNode> Face;
	typedef Face_<Node, pNode>* pFace;
	typedef Face_<Node, pNode>& ref_Face;
	typedef const Face_<Node, pNode>& const_ref_Face;
	typedef const Face_<Node, pNode>* const_pFace;
public:
	pNode pnode;
	pNode pneighbor;
	FaceType face_type;
	Direction direction;  //pnode ---> pneighbor

	Face_() :
			pnode(nullptr), pneighbor(nullptr), face_type(_Null_), direction(
					_XM_) {

	}
	Face_(pNode pn, pNode pnei, Direction d, FaceType ft) :
			pnode(pn), pneighbor(pnei), face_type(ft), direction(d) {

	}
	//NFace_(pNode pn, Direction d) {
	//	pnode = pn;
	//	pneighbor = pnode->get_neighbor(d);
	//	face_type = GetFaceFype(pnode, pneighbor);
	//	direction = d;
	//}
	Face_(const Face& a) :
			pnode(a.pnode), pneighbor(a.pneighbor), face_type(a.face_type), direction(
					a.direction) {
	}
	//
	//operator ==================================
	ref_Face& operator=(const Face& a) {
		pnode = a.pnode;
		pneighbor = a.pneighbor;
		face_type = a.face_type;
		direction = a.direction;
		return *this;
	}
	bool operator==(const Face& a) const {
		return (pnode == a.pnode) && (pneighbor == a.pneighbor)
				&& (face_type == a.face_type) && (direction == a.direction);
	}
	bool operator!=(const Face& a) const {
		return !((pnode == a.pnode) && (pneighbor == a.pneighbor)
				&& (face_type == a.face_type) && (direction == a.direction));
	}
	/*
	 * get
	 */
	pNode& pori() { //pNode origin
		return pnode;
	}
	pNode& pnei() { //pNode center
		return pneighbor;
	}
	FaceType& ft() {
		return face_type;
	}
	const FaceType& ft() const {
		return face_type;
	}
	Direction& dir() {
		return direction;
	}
	const Direction& dir() const {
		return direction;
	}
	vt area() const {
		ASSERT(pnode!=nullptr);
		ASSERT(pneighbor!=nullptr);
		vt res = pnode->face_area(direction);
		if (this->face_type == _CoarseFine_) {
			res = pneighbor->face_area(Opposite(this->direction));
		}
		return res;
	}

	void set(pNode pn, pNode pnei, const Direction& d, const FaceType& ft) {
		pnode = pn;
		pneighbor = pnei;
		face_type = ft;
		direction = d;
	}
	//show=======================================
	void show() const {
		std::cout << "Face show ===============\n";
		std::cout << "Dim        : " << Dim << '\n';
		std::cout << "Direction  : " << direction << '\n';
		std::cout << "Face type  : " << ParseFaceType(face_type) << '\n';
		std::cout << "Center p x : " << pnode->p(_C_, _X_) << '\n';
		std::cout << "       p y : " << pnode->p(_C_, _Y_) << '\n';
		if (Dim == 3) {
			std::cout << "       p z : " << pnode->p(_C_, _Z_) << '\n';
		}
	}

};

/*
 * Function out of calss ==================================
 */
template<class Node, class Ret, class Args>
void Traversal(Node*& pn, std::function<Ret(Node*&, Args)> fun, Args &args) {
	if (pn == nullptr) {
		return;
	} else {
		fun(pn, args);
		_IF_TRUE_RETRUN(pn == nullptr);
		if (pn->has_child()) {
			for (St i = 0; i < Node::NumChildren; i++) {
				Node*& c = pn->child[i];
				if (c != nullptr) {
					Traversal(c, fun, args);
				}
			}
		}
	}
}

template<class Node>
void Traversal(Node*& pn, std::function<void(Node*&)> fun) {
	if (pn == nullptr) {
		return;
	} else {
		fun(pn);
		_IF_TRUE_RETRUN(pn == nullptr);
		if (pn->has_child()) {
			for (St i = 0; i < Node::NumChildren; i++) {
				Node*& c = pn->child[i];
				if (c != nullptr) {
					Traversal(c, fun);
				}
			}
		}
	}
}

template<class Node>
void TraversalLeaf(Node*& pn, std::function<void(Node*)> fun) {
	if (pn == nullptr) {
		return;
	} else {
		if (pn->is_leaf()){
			fun(pn);
		}
		_IF_TRUE_RETRUN(pn == nullptr);
		if (pn->has_child()) {
			for (St i = 0; i < Node::NumChildren; i++) {
				Node*& c = pn->child[i];
				if (c != nullptr) {
					TraversalLeaf(c, fun);
				}
			}
		}
	}
}


}
//

#endif /* NODE_H_ */
