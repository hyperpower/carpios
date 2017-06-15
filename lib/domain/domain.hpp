#ifndef DOMAIN_HPP_
#define DOMAIN_HPP_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "grid.hpp"
#include "adaptive.hpp"
#include "stencil.hpp"
#include "data.hpp"
#include "shape.hpp"
#include "boundary.hpp"

#include <cmath>

namespace carpio {
typedef Float Cvt;
typedef Float Vt;

typedef Cell_<Cvt, 2> Cell_2D;
typedef Cell_<Cvt, 3> Cell_3D;

typedef Data_<Float, 2> Data_2D;
typedef Data_<Float, 3> Data_3D;

typedef PData_<Float, Float, 2> PData_2D;
typedef PData_<Float, Float, 3> PData_3D;

typedef Node_<Float, Float, 2> Node_2D;
typedef Node_<Float, Float, 3> Node_3D;
typedef Node_2D* pNode_2D;
typedef Node_3D* pNode_3D;
typedef const Node_<Float, Float, 2>* const_pNode_2D;
typedef const Node_<Float, Float, 3>* const_pNode_3D;

typedef Grid_<Float, Float, 2> Grid_2D;
typedef Grid_<Float, Float, 3> Grid_3D;

typedef Ghost_<Float, Float, 2> Ghost_2D;
typedef Ghost_<Float, Float, 3> Ghost_3D;

typedef Stencil_<Float, Float, 2, 1> Stencil_2D1;
typedef Stencil_<Float, Float, 2, 2> Stencil_2D2;

template<typename COO_VALUE, typename VALUE, St DIM>
class Domain_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef Domain_<cvt, vt, Dim> Self;
	typedef Domain_<cvt, vt, Dim>* pSelf;
	typedef Grid_<cvt, vt, Dim> Grid;
	typedef Grid_<cvt, vt, Dim> *pGrid;
	typedef const Grid_<cvt, vt, Dim> * const_pGrid;
	typedef Grid_<cvt, vt, Dim>& ref_Grid;
	typedef const Grid_<cvt, vt, Dim>& const_ref_Grid;
	typedef Cell_<cvt, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<vt, Dim> Data;
	typedef Data *pData;
	typedef Node_<cvt, vt, Dim> Node;
	typedef Node_<cvt, vt, Dim>& ref_Node;
	typedef const Node_<cvt, vt, Dim>& const_ref_Node;
	typedef Node_<cvt, vt, Dim> *pNode;
	typedef const Node_<cvt, vt, Dim>* const_pNode;
	typedef Ghost_<cvt, vt, Dim> Ghost;
	typedef Ghost_<cvt, vt, Dim> *pGhost;
	typedef const Ghost_<cvt, vt, Dim> * const_pGhost;
	typedef Ghost_<cvt, vt, Dim>& ref_Ghost;
	typedef const Ghost_<cvt, vt, Dim>& const_ref_Ghost;

	typedef Adaptive_<cvt, vt, Dim> Adaptive;
	typedef Adaptive_<cvt, vt, Dim>& ref_Adaptive;
	typedef const Adaptive_<cvt, vt, Dim>& const_ref_Adaptive;
	typedef Adaptive_<cvt, vt, Dim> *pAdaptive;
	typedef BoundaryIndex_<cvt, vt> BoundaryIndex;
	typedef BoundaryIndex_<cvt, vt>* pBoundaryIndex;
	typedef const BoundaryIndex_<cvt, vt>* const_pBoundaryIndex;
	typedef BoundaryIndex_<cvt, vt>& ref_BoundaryIndex;
	typedef const BoundaryIndex_<cvt, vt>& const_ref_BoundaryIndex;

	typedef BoundaryCondition_<cvt, vt> BoundaryCondition;
	typedef BoundaryCondition_<cvt, vt>* pBoundaryCondition;
	typedef const BoundaryCondition_<cvt, vt>* const_pBoundaryCondition;

	typedef Face_<Node, pNode> Face;
	typedef Face_<Node, pNode> *pFace;
	typedef Shape_<cvt, Dim> Shape;
	typedef Shape_<cvt, Dim>* pShape;

	typedef std::function<vt(cvt, cvt, cvt)> Function;

	typedef std::function<void(const_ref_Node)> Fun_const_ref_Node;
	typedef std::function<void(ref_Node)> Fun_ref_Node;
	typedef std::function<void(const_pNode)> Fun_const_pNode;
	typedef std::function<void(pNode)> Fun_pNode;

public:
	// data
	// 2D : the domain is bounded by a shape
	//      the shape doesn't have holes.
	pShape _pshape_bound;
	// Innner solid
	std::list<pShape> _l_pshape_inner;
	// Grid: the calculation region
	pGrid _pgrid;
	// Adaptive
	pAdaptive _padaptive;
	// BC index
	pBoundaryIndex _pbindex;
	// Ghost nodes
	pGhost _pghost;

protected:

	void _new_grid(cvt UL) {
		// build grid ------------------
		Float max_x = _pshape_bound->max_x();
		Float max_y = _pshape_bound->max_y();
		Float min_x = _pshape_bound->min_x();
		Float min_y = _pshape_bound->min_y();
		St n_x = std::ceil((max_x - min_x) / UL);
		St n_y = std::ceil((max_y - min_y) / UL);
		_pgrid = new Grid(  //
				n_x, min_x, UL, //
				n_y, min_y, UL);
	}
	void _new_adaptive(St lmin = 1, St lmax = 1) {
		_padaptive = new Adaptive(_pgrid, lmin, lmax);
	}
	void _new_boundary_index() {
		_pbindex = new BoundaryIndex();
	}
	void _new_ghost() {
		_pghost = new Ghost(_pgrid);
	}

public:
	Domain_(pShape bound, cvt unit_length) :
			_pshape_bound(bound) {
		_new_grid(unit_length);
		_new_adaptive();
		_new_boundary_index();
		_new_ghost();
	}

	Domain_(pShape bound, cvt unit_length, St minl, St maxl) :
			_pshape_bound(bound) {
		_new_grid(unit_length);
		_new_adaptive(minl, maxl);
		_new_boundary_index();
		_new_ghost();
		//_padaptive->adapt_full();
	}

	Domain_(pShape bound, const std::list<pShape>& l_pshape_inner,
			cvt unit_length, St minl, St maxl) :
			_pshape_bound(bound), _l_pshape_inner(l_pshape_inner) {
		_new_grid(unit_length);
		_new_adaptive(minl, maxl);
		_new_boundary_index();
		_new_ghost();
		//_padaptive->adapt_full();
	}

	~Domain_() {
		if (_pghost != nullptr) {
			delete _pghost;
			_pghost = nullptr;
		}
		if (_pbindex != nullptr) {
			delete _pbindex;
			_pbindex = nullptr;
		}
		if (_padaptive != nullptr) {
			delete _padaptive;
			_padaptive = nullptr;
		}
		if (_pgrid != nullptr) {
			delete _pgrid;
			_pgrid = nullptr;
		}
	}

	void build() {
		// adapt the solid
		// bug ----> if the solids contact, the grid may not build well.
		_padaptive->adapt_bound_solid(*_pshape_bound, 1.0);
		for (typename std::list<pShape>::iterator iter =
				_l_pshape_inner.begin(); iter != _l_pshape_inner.end();
				++iter) {
			if ((*iter) != nullptr) {
				_padaptive->adapt_inner_solid(*(*iter));
			}
		}
		// connect nodes
		_pgrid->connect_root();
		_pgrid->connect_nodes();
		// new date on leaf, here intial set all the arrays on data
		// is empty, reset the data in actual calculation classes;
		_pgrid->new_data_on_leaf(0, 0, 0, 0);
		_pgrid->set_data_index();
		// build ghost
		_pghost->build();
		_pghost->connect();
		_pghost->set_data_idx();
		// set shape_idx and seg_idx in ghost node
		/*
		 * Method: calculates the ghost node belongs to which boundary segment.
		 * find the segment in shape which intersect with the line between ghost
		 * node and the boundary node
		 * choose nearest segment
		 * The line is axi align
		 * function : vt Intersect(Segment seg, Axis axis, vt ori)
		 */
		_pghost->set_boundary_index(0, *_pshape_bound);
		int i = 1;
		for (auto iter = _l_pshape_inner.begin(); iter != _l_pshape_inner.end();
				++iter) {
			pShape ps = (*iter);
			_pghost->set_boundary_index(i, *ps);
			i++;
		}
	}

public:

	pGrid pgrid() {
		return _pgrid;
	}
	const_pGrid pgrid() const {
		return _pgrid;
	}
	ref_Grid grid() {
		return *_pgrid;
	}
	const_ref_Grid grid() const {
		return *_pgrid;
	}
	ref_Adaptive adaptive() {
		return *_padaptive;
	}
	const_ref_Adaptive adaptive() const {
		return *_padaptive;
	}
	pGhost p_ghost() {
		return _pghost;
	}
	const_pGhost p_ghost() const {
		return _pghost;
	}
	ref_Ghost ghost() {
		return *_pghost;
	}
	const_ref_Ghost ghost() const {
		return *_pghost;
	}
	pBoundaryIndex p_b_index() {
		return _pbindex;
	}
	const_pBoundaryIndex p_b_index() const {
		return _pbindex;
	}
	ref_BoundaryIndex ref_b_index() {
		return (*_pbindex);
	}
	const_ref_BoundaryIndex ref_b_index() const {
		return (*_pbindex);
	}

	const_pBoundaryCondition find_bc(St si, St segi, St vali) const {
		return this->_pbindex->find(si, segi, vali);
	}
	/*
	 * set
	 */
	void set_val_grid(St idx, Function fun) {
		pGrid pg = this->_pgrid;
		for (auto it = pg->begin_leaf(); it != pg->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			pn->cd(idx) = fun(it->p(_C_, _X_), it->p(_C_, _Y_),
					it->p(_C_, _Z_));
		}
	}
	void set_val_ghost(St idx, Function fun) {
		typedef typename Ghost::GhostNode GhostNode;
		std::function<void(GhostNode&)> _fun = [&idx, &fun](GhostNode& node) {
			pNode pg = node.second.pghost;
			if (pg != nullptr) {
				pg->cd(idx) = fun(pg->p(_C_, _X_), pg->p(_C_, _Y_),
						pg->p(_C_, _Z_));
			}
		};
		this->p_ghost()->for_each_node(_fun);
	}
	void set_val(St idx, Function fun) {
		this->set_val_grid(idx, fun);
		this->set_val_ghost(idx, fun);
	}

	void set_val_ghost_by_bc(St idx) {  //bug bc 2
		typedef typename Ghost::GhostNode GhostNode;
		std::function<void(GhostNode&)> _fun = [&idx, this](GhostNode& node) {
			pNode pg = node.second.pghost;
			if (pg != nullptr) {
				// find in BoundaryIndex
				St si = node.second.shape_idx;
				St segi = node.second.seg_idx;
				const_pBoundaryCondition pbc = this->find_bc(si, segi, idx);
				if(pbc->get_type() == BoundaryCondition::_BC1_) {
					pg->cd(idx) = pbc->get_val(pg->p(_C_, _X_), pg->p(_C_, _Y_),
							pg->p(_C_, _Z_));
				} else { //bc 2
					typename Ghost::GhostID gid = Ghost::ToGhostID(pg);
					ASSERT(gid.step == 0);
					const_pNode po = pg->father;
					Axes a = FaceDirectionToAxes(gid.direction);
					vt val = pbc->get_val(pg->cp(_X_), pg->cp(_Y_),
							pg->cp(_Z_), a);
					cvt dl = (po->cp(a) - pg->cp(a));
					// pg = po - val*dl
					pg->cd(idx) = po->cd(idx) - val*dl;
				}
			}
		};
		this->p_ghost()->for_each_node(_fun);
	}
	/*
	 * new data
	 */
	void new_data(const St& nc, const St& nf, const St& nv, const St& nutp) {
		this->_pgrid->new_data_on_leaf(nc, nf, nv, nutp);
		this->_pghost->new_data(nc, nf, nv, nutp);
	}
	void resize_data(const St& nc, const St& nf, const St& nv, const St& nutp) {
		this->_pgrid->resize_data_on_leaf(nc, nf, nv, nutp);
		this->_pghost->resize_data(nc, nf, nv, nutp);
	}
};



}

#endif
