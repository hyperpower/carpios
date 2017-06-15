#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "path.hpp"
#include "shape.hpp"
#include "geometry/_line.hpp"
#include <functional>
#include <map>

namespace carpio {
/*
 * struct GhostID
 * This struct used to identify the ghost node
 */

template<typename COO_VALUE, typename VALUE, St DIM>
struct GhostID_ {
	//typedef int (*pFun_set_bc)(Node_<COO_VALUE, VALUE, DIM>*,
	//		GhostID_<COO_VALUE, VALUE, DIM>&, utPointer);

	St root_idx;     //the root idx of the origin node
	//st idx;        //the local idx of the origin node
	Path_<DIM> path; //the path of the origin node
	St step;         //the steps of ghost node, we can choose multiple ghost Node,
				     // usually step = 0
	Direction direction; //The direction only on x, y or z
//--------------------------------------------------------------------
	//int bc_type;
	//pFun_set_bc pfun_bc;
	St shape_idx;
	St seg_idx;
};

template<typename COO_VALUE, typename VALUE, St DIM>
struct GhostID_compare_ {
	typedef GhostID_<COO_VALUE, VALUE, DIM> Gid;
	bool operator()(const Gid& lhs, const Gid& rhs) const {
		if (lhs.root_idx < rhs.root_idx) {
			return true;
		} else if (lhs.root_idx == rhs.root_idx) {
			if (lhs.path < rhs.path) {
				return true;
			} else if (lhs.path == rhs.path) {
				if (int(lhs.direction) < int(rhs.direction)) {
					return true;
				} else if (int(lhs.direction) == int(rhs.direction)) {
					return lhs.step < rhs.step;
				}
			}
		}
		return false;
	}
};
/*
 * class Ghost
 *
 * It saves all the ghost nodes on a map structure.
 */

template<typename COO_VALUE, typename VALUE, int DIM>
class Ghost_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef Ghost_<COO_VALUE, VALUE, DIM> Self;
	typedef const Ghost_<COO_VALUE, VALUE, DIM> const_Self;
	typedef Ghost_<COO_VALUE, VALUE, DIM>* pSelf;
	typedef const Ghost_<COO_VALUE, VALUE, DIM>* const_pSelf;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell* pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data* pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM>* pNode;
	typedef const Node_<COO_VALUE, VALUE, DIM>* const_pNode;
	typedef void (*pfunction)(pNode, utPointer);
	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef Shape_<COO_VALUE, DIM> Shape;
	typedef Shape_<COO_VALUE, DIM> *pShape;
	typedef void (*pfunction_conditional)(arrayList&, pNode, utPointer);
	typedef GhostID_<COO_VALUE, VALUE, DIM> GhostID;
	typedef GhostID_compare_<COO_VALUE, VALUE, DIM> GhostID_compare;

	struct GhostVal {
		pNode pghost;
		St shape_idx;
		St seg_idx;
	};

	typedef std::pair<const GhostID, GhostVal> GhostNode;

	typedef std::map<GhostID, GhostVal, GhostID_compare> GhostMap;
	typedef typename GhostMap::iterator iterator;
	typedef typename GhostMap::const_iterator const_iterator;

protected:
	GhostMap _ghostmap;
	pGrid _pgrid;
	//
	/*
	 *  new ghost node
	 */
	pNode new_ghost_node(pNode pn, Direction dir) {
		int ghost_node_type = _Ghost_;
		//direction 4 5 6 7
		Cell c(*(pn->cell));
		switch (dir) {
		case _XM_:
			c.transfer(-c.get_d(_X_), 0.0, 0.0);
			break;
		case _YM_:
			c.transfer(0.0, -c.get_d(_Y_), 0.0);
			break;
		case _XP_:
			c.transfer(c.get_d(_X_), 0.0, 0.0);
			break;
		case _YP_:
			c.transfer(0.0, c.get_d(_Y_), 0.0);
			break;
		case _ZM_:
			c.transfer(0.0, 0.0, -c.get_d(_Z_));
			break;
		case _ZP_:
			c.transfer(0.0, 0.0, c.get_d(_Z_));
			break;
		default:
			ASSERT(false);
		}
		//pNode f, int nt, st level, st root_idx, st path,const Cell& c
		pNode ghostnode = new Node(  //
				nullptr,//
				ghost_node_type,  //
				pn->get_level(),  //
				pn->get_root_idx(),  //
				pn->get_idx(),       //
				pn->get_path(),      //
				c);
		ghostnode->father = pn;
		if (pn->data != nullptr) {
			ghostnode->data = new Data(*(pn->data));
		}
		return ghostnode;
	}
	/*
	 *  new ghost nodes outside of the boundary of Grid
	 */
	void new_ghost_nodes(pGrid pg) {
		for (typename Grid::iterator_face iter = pg->begin_face();
				iter != pg->end_face(); ++iter) {
			if (iter->ft() == _Boundary_) {
				pNode po = (*iter).pori();
				GhostID gid;
				gid.root_idx = po->get_root_idx(); //the root idx of the origin node
				gid.path = po->get_path(); //the local idx of the origin node
				gid.step = 0; //the steps of ghost node, we can choose multiple ghost Node,
							  // usually step =
				gid.direction = iter->dir(); //The direction only on x, y or z
				// --------------------------------
				//gid.bc_type = 0;
				//gid.pfun_bc = nullptr;
				//gid.shape_idx = -1;
				//gid.seg_idx = -1;
				pNode pghost = new_ghost_node(po, iter->dir());
				GhostVal gval = { pghost, 0, 0 };
				//change the index id ------------------------
				_ghostmap.insert(GhostNode(gid, gval));
			}
		}
	}

public:
	Ghost_(pGrid pg) :
			_pgrid(pg) {
	}
	~Ghost_() {
		for (typename GhostMap::iterator it = _ghostmap.begin();
				it != _ghostmap.end(); ++it) {
			if (it->second.pghost != nullptr) {
				delete it->second.pghost;
				it->second.pghost = nullptr;
			}
		}
	}
	void build() {
		ASSERT(_ghostmap.empty());
		new_ghost_nodes(_pgrid);
	}
protected:
	void _connect(const GhostNode& gn) {
		pNode pg = gn.second.pghost;
		pNode po = pg->father;
		Direction dir_o = gn.first.direction;
		po->set_neighbor(pg, dir_o);
	}
public:
	void connect() {
		// when the ghost is built, the ghost node are already connect to original node.
		// the father of ghost node is the original one,
		// so the find_neighbor_fast() can not be used on ghost node.
		// here, we connect original node to the ghost node
		// the find_neighbor_fast() can be used on boundary node
		for (typename GhostMap::iterator it = _ghostmap.begin();
				it != _ghostmap.end(); ++it) {
			this->_connect(*it);
		}
	}

	iterator begin() {
		return _ghostmap.begin();
	}
	const_iterator begin() const {
		return _ghostmap.begin();
	}
	iterator end() {
		return _ghostmap.end();
	}
	const_iterator end() const {
		return _ghostmap.end();
	}
	void set_boundary_index(iterator& iter, St si, St segi) {
		GhostVal& gval = iter->second;
		gval.shape_idx = si;
		gval.seg_idx = segi;
	}
	int get_shape_index(iterator& iter) const {
		return iter->first.shape_index;
	}
	int get_shape_index(const_iterator& iter) const {
		return iter->first.shape_index;
	}
	int get_seg_index(iterator& iter) const {
		return iter->first.shape_index;
	}
	int get_seg_index(const_iterator& iter) const {
		return iter->first.shape_index;
	}

	/*
	 * set boundary index
	 */
	void set_boundary_index(int shape_idx, const Shape& shape) {
		for (iterator iter = this->begin(); iter != this->end(); ++iter) {
			set_seg_index(shape_idx, shape, iter);
		}
	}

	St _choose_a_seg(std::list<St>& ls, const Shape& shape, iterator& iter,
			int flag) {
		ASSERT(!ls.empty());
		if (ls.size() == 1) {
			return *(ls.begin());
		} else {
			pNode pg = iter->second.pghost;
			for (auto i = ls.begin(); i != ls.end(); ++i) {
				typename Shape::Seg2D seg = shape.seg(*i);
				typename Shape::Poi2D poi(pg->cp(_X_), pg->cp(_Y_));
				int side = OnWhichSide3(seg, poi);
				if (side != flag) {
					i = ls.erase(i);
				}
			}
			if (ls.size() > 1) {
				// choose the closest one
				const GhostID& gid = iter->first;
				Orientation ori;
				Axes axi;
				Direction dir = gid.direction;
				FaceDirectionToOrientationAndAxes(dir, ori, axi);
				Axes vaxi = VerticalAxes2D(axi);
				auto iterres = ls.begin();
				cvt min_intercept;
				for (auto i = ls.begin(); i != ls.end(); ++i) {
					typename Shape::Seg2D seg = shape.seg(*i);
					typename Shape::Poi2D poi(pg->cp(_X_), pg->cp(_Y_));
					int side = OnWhichSide3(seg, poi);
					if (side != flag) {
						i = ls.erase(i);
					} else {
						Line_<cvt> l(seg.ps(), seg.pe());
						cvt intercept = Abs(
								l.cal(vaxi, pg->cp(vaxi)) - pg->cp(vaxi));
						if (i == ls.begin() || intercept < min_intercept) {
							min_intercept = intercept;
							iterres = i;
						}
					}
				}
				return *iterres;
			}
			return *(ls.begin());
		}
	}

	void set_seg_index(int shape_idx, const Shape& shape, iterator& iter) {
		ASSERT(Dim == 2);
		// get a line of ghost node
		Orientation ori;
		Axes axi;
		Direction dir = iter->first.direction;
		FaceDirectionToOrientationAndAxes(dir, ori, axi);
		Axes va = VerticalAxes2D(axi);
		//Orientation oori = Opposite(ori);
		std::list<St> ls;
		GhostVal& gval = iter->second;
		shape.find_seg_across(ls, va, gval.pghost->cp(va));
		if (!ls.empty()) {
			St seg_idx = _choose_a_seg(ls, shape, iter, -1);
			set_boundary_index(iter, shape_idx, seg_idx);
		} else {
			// first, find the closest vertex on the shape. The two segments
			// connect to the vertex will be chosen.
			St ver_idx = shape.find_closest_vertex(gval.pghost->cp(_X_),
					gval.pghost->cp(_Y_));
			shape.find_seg_connect_to_vertex(ls, ver_idx);
			// second, find
			St seg_idx = _choose_a_seg(ls, shape, iter, -1);
			set_boundary_index(iter, shape_idx, seg_idx);
		}

	}

	void for_each_node(std::function<void(GhostNode&)> fun) {
		for (iterator iter = this->begin(); iter != this->end(); ++iter) {
			fun(*iter);
		}
	}
	void for_each_node(std::function<void(const GhostNode&)> fun) const {
		for (const_iterator iter = this->begin(); iter != this->end(); ++iter) {
			fun((*iter));
		}
	}
	/*
	 * the data index of ghost node simply equal to negative of original node
	 */
	void set_data_idx() {
		std::function<void(GhostNode&)> fun = [](GhostNode& node) {
			GhostID gid = node.first;
			pNode pg = node.second.pghost;
			pNode po = pg->father;
			int id = po->d_idx() * 100 + gid.direction;
			// the date index of the ghost node is (
				pg->d_idx() = -(id);//negative added here;
			};
		for_each_node(fun);
	}

	void new_data(const St& nc, const St& nf, const St& nv, const St& nutp
			) {
		std::function<void(GhostNode&)> fun = [&](GhostNode& node) {
			pNode pg = node.second.pghost;
			pg->new_data(nc, nf, nv, nutp);
		};
		for_each_node(fun);
	}
	void resize_data(const St& nc, const St& nf, const St& nv, const St& nutp
			) {
		std::function<void(GhostNode&)> fun = [&](GhostNode& node) {
			pNode pg = node.second.pghost;
			pg->resize_data(nc, nf, nv, nutp);
		};
		for_each_node(fun);
	}

	iterator find(const GhostID& gid) {
		return _ghostmap.find(gid);
	}

	/*
	 * the static function
	 *
	 */
	static GhostID ToGhostID(const_pNode pn) {
		ASSERT(pn != nullptr);
		ASSERT(pn->get_type() == _Ghost_);

		//st root_idx;   //the root idx of the origin node
		//Path_<DIM> path; //the path of the origin node
		//st step;       //the steps of ghost node, we can choose multiple ghost Node,
		// usually step = 0
		//Direction direction; //The direction only on x, y or z
		GhostID res;
		res.step = 0;
		const_pNode pc = pn;  //the pnode close to original one
		while (pc->father->get_type() == _Ghost_) {
			pc = pc->father;
			res.step++;
		}
		const_pNode po = pc->father;
		res.root_idx = po->get_root_idx();
		res.path = po->get_path();
		for (St i = 0; i < NumFaces; i++) {
			if (po->neighbor[i] == pc) {
				res.direction = FaceDirectionInOrder(i);
			}
		}
		return res;
	}
}
;

template<typename COO_VALUE, typename VALUE>
class BoundaryCondition_ {
public:
	typedef VALUE vt;
	typedef COO_VALUE cvt;
	typedef std::function<vt(cvt, cvt, cvt)> Fun;
	typedef BoundaryCondition_<cvt, vt> Self;
	static const int _BC1_ = 1;
	static const int _BC2_ = 2;
protected:
	// data
	// 1 Bondary conditon type
	// 2 function
	int _type;
	Fun _function;
	Fun _function2;
	Fun _function3;
public:
	// Constructor
	BoundaryCondition_() {
		// default boundary condition is symmetric boundary condition
		_type = _BC2_;
		Fun df = [](cvt x, cvt y, cvt z) {return 0;};
		_function = df;
		_function2 = df;
		_function3 = df;
	}
	/*
	 * this constructor should not used to BC2
	 */
	BoundaryCondition_(int type, Fun fun) :
			_type(type), _function(fun), _function2(fun), _function3(fun) {
	}
	BoundaryCondition_(int type, Fun fun_x, Fun fun_y, Fun fun_z) :
			_type(type), _function(fun_x), _function2(fun_y), _function3(fun_z) {
	}
	BoundaryCondition_(const Self& self) :
			_type(self._type), _function(self._function), _function2(
					self._function2), _function3(self._function3) {
	}
	// get
	int get_type() const {
		return _type;
	}
	vt get_val(cvt x, cvt y, cvt z) const {
		return _function(x, y, z);
	}
	/*
	 * for the vector type boundary condition
	 * (x,y,z) is the location
	 * Axes    is axes
	 */
	vt get_val(cvt x, cvt y, cvt z, Axes a) const {
		switch (a) {
		case _X_: {
			return _function(x, y, z);
			break;
		}
		case _Y_: {
			return _function(x, y, z);
			break;
		}
		case _Z_: {
			return _function(x, y, z);
			break;
		}
		}
		SHOULD_NOT_REACH;
		return 0;
	}
	// set
	void set_function(Fun fun, Axes a = _X_) {
		switch (a) {
		case _X_: {
			_function = fun;
			break;
		}
		case _Y_: {
			_function2 = fun;
			break;
		}
		case _Z_: {
			_function3 = fun;
			break;
		}
		}
	}
	void set_function(Fun fun_x, Fun fun_y, Fun fun_z) {
		_function = fun_x;
		_function2 = fun_y;
		_function3 = fun_z;
	}
	void set_default_1_bc(const vt& val) {
		_type = _BC1_;
		Fun f = [val](cvt x, cvt y, cvt z) {return val;};
		_function = f;
		_function2 = f;
		_function3 = f;
	}
	void set_default_1_bc(Fun fun) {
		_type = _BC1_;
		Fun f = fun;
		_function = f;
		_function2 = f;
		_function3 = f;
	}
	void set_default_2_bc(const vt& val) {
		_type = _BC2_;
		Fun f = [val](cvt x, cvt y, cvt z) {return val;};
		_function = f;
		_function2 = f;
		_function3 = f;
	}
	void set_default_2_bc(Fun fun) {
		_type = _BC2_;
		Fun f = fun;
		_function = f;
		_function2 = f;
		_function3 = f;
	}

};

template<typename COO_VALUE, typename VALUE>
class BoundaryIndex_ {
protected:
	struct BCID {
		St shape_idx;
		St seg_idx;
		St val_idx;
	};

	struct BCID_compare {
		typedef BCID BCid;
		bool operator()(const BCid& lhs, const BCid& rhs) const {
			if (lhs.shape_idx < rhs.shape_idx) {
				return true;
			} else if (lhs.shape_idx == rhs.shape_idx) {
				if (lhs.seg_idx < rhs.seg_idx){
					return true;
				}else if(lhs.seg_idx == rhs.seg_idx){
					return lhs.val_idx < rhs.val_idx;
				}
			}
			return false;
		}

	};
public:
	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef BoundaryCondition_<cvt, vt> BoundaryCondition;
	typedef BoundaryCondition_<cvt, vt>* pBoundaryCondition;
	typedef const BoundaryCondition_<cvt, vt>* const_pBoundaryCondition;
	//typedef BCID_compare_<cvt,vt> BCID_compare;
	typedef std::pair<const BCID, pBoundaryCondition> BCNode;
	typedef std::map<BCID, pBoundaryCondition, BCID_compare> BCMap;
protected:
	BCMap _BCmap;
	const_pBoundaryCondition pdefault_BC;

	typedef typename BCMap::iterator iterator;
	typedef typename BCMap::const_iterator const_iterator;

public:
	//constructor
	BoundaryIndex_() :
			_BCmap() {
		pdefault_BC = new BoundaryCondition();
	}
	~BoundaryIndex_() {
		delete pdefault_BC;
	}
	//
	void insert(St si, St segi, St vi, pBoundaryCondition pbc) {
		//
		BCID key = { si, segi, vi };
		//key.seg_idx = segi;
		//key.shape_idx = si;
		//key.val_idx = vi;
		BCNode bcn(key, pbc);
		_BCmap.insert(bcn);
	}

	const_pBoundaryCondition find(St si, St segi, St vali) {
		BCID key;
		key.seg_idx = segi;
		key.shape_idx = si;
		key.val_idx = vali;
		iterator it = _BCmap.find(key);
		if (it != _BCmap.end()) {
			// found
			return (it->second);
		} else {
			// not found
			return pdefault_BC;
		}
	}
};

}

#endif
