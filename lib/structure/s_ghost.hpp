/*
 * s_ghost.hpp
 *
 *  Created on: Jul 8, 2017
 *      Author: zhou
 */

#ifndef _S_GHOST_HPP_
#define _S_GHOST_HPP_


#include "s_define.hpp"
#include "s_grid.hpp"
#include "s_boundary.hpp"

namespace structure{

template<St DIM>
class Ghost_ {
public:
	static const St Dim = DIM;
	typedef BoundaryCondition BC;
	typedef std::shared_ptr<BoundaryCondition> spBC;
	typedef std::shared_ptr<const BoundaryCondition> spcBC;
	typedef BoundaryIndex BI;
	typedef std::shared_ptr<BI> spBI;
	typedef Grid_<Dim> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Index_<Dim> Index;
	typedef Stencil_<Dim> Stencil;
	typedef typename Stencil::Expression Expression;
	typedef typename Stencil::spExpression spExpression;
	typedef std::shared_ptr<Stencil> spStencil;

	typedef std::function<void(const Index&)> Fun_index;
protected:
	std::string _name;
	spBI _bi;
	spGrid _grid;
public:
	Ghost_(spGrid spg, spBI bi) :
			_grid(spg), _bi(bi) {
		_name = "Ghost";
	}

	const std::string& name() const {
		return _name;
	}

	virtual St which_boudary_seg(const Index& index) const {
		SHOULD_NOT_REACH;
		return 0;
	}
	virtual St which_boudary_shape(const Index& index) const {
		SHOULD_NOT_REACH;
		return 0;
	}

	virtual bool is_ghost(const Index& index) const {
		std::cout << "Boundary : is ghost\n";
		SHOULD_NOT_REACH;
		return false;
	}

	virtual bool IS_GHOST(const Index& INDEX) const {
		std::cout << "Boundary : is ghost\n";
		SHOULD_NOT_REACH;
		return false;
	}
	virtual bool is_normal(const Index& index) const {
		std::cout << "Boundary : is ghost\n";
		SHOULD_NOT_REACH;
		return false;
	}

	virtual bool IS_NORAML(const Index& INDEX) const {
		std::cout << "Boundary : is ghost\n";
		SHOULD_NOT_REACH;
		return false;
	}

	virtual void for_each_ghost(Fun_index fun) {
		std::cout << "Boundary : for each ghost\n";
		SHOULD_NOT_REACH;
	}

	//virtual void add_bc(St shape_idx, St seg_idx, const std::string&name,
	//		spBC spbc) {
	//	std::cout << "Boundary : add bc\n";
	//	SHOULD_NOT_REACH;
	//}

	//virtual spcBC find(const Index& index, const std::string&name) const {
	//	std::cout << "Boundary : find\n";
	//	SHOULD_NOT_REACH;
	//	return nullptr;
	//}

	virtual Index ori_index(const Index& index) const {
		SHOULD_NOT_REACH;
		return Index();
	}
	virtual St ori_axes(const Index& index) const {
		SHOULD_NOT_REACH;
		return 0;
	}

	virtual Orientation orientation(const Index& index) const {
		SHOULD_NOT_REACH;
		return _C_;
	}
	virtual spExpression substitute(spExpression spc,
			const std::string& vname) const {
		SHOULD_NOT_REACH;
		return nullptr;
	}

	virtual spcBC find(const Index& index, const std::string& vname) const {
		return nullptr;
	}

	virtual ~Ghost_() {

	}
protected:
	spExpression _substitute_ghost(spExpression spc, spcBC pbc,
			const Index& ghoidx, const Vt& coe,
			const std::string& vname) const {
		// find boundary condition
		//          bc
		//   ghost  |   node
		// ----o----|----o----
		//     0    |    1
		//    c0   cf   c1
		//    |-----dl---|
		//    |-----|----|
		//      d0    d1
		Index oriidx = this->ori_index(ghoidx);
		St a = this->ori_axes(ghoidx);
		Orientation ori = this->orientation(ghoidx);
		Vt c0 = this->_grid->c_(a, ghoidx[a]);
		Vt c1 = this->_grid->c_(a, oriidx[a]);
		Vt cf = this->_grid->f_(a, ori, oriidx[a]);
		Vt dl = std::abs(c1 - c0);
		Vt d0 = std::abs(cf - c0);
		Vt d1 = std::abs(c1 - cf);
		Vt val = pbc->get_val(               //
				0.0,//time ==========================
				this->_grid->f_(_X_, ori, oriidx[_X_]), //
				this->_grid->f_(_Y_, ori, oriidx[_Y_]), //
				this->_grid->f_(_Z_, ori, oriidx[_Z_]));
		if (pbc->get_type() == BC::_BC1_) {
			//		// Boundary condition 1
			//spc->insert(coe * val, ghoidx, 0);
			// -------
			spc->insert(-coe * d0 / d1, oriidx, 1);
			spc->insert(coe * dl * val / d1, oriidx, 0);
			// -------
		} else { // bc2
			// pg = po - val * dl;
			spc->insert(coe, oriidx, 1);
			spc->insert(-val * dl * coe, oriidx, 0);	//constant
		}
		return spc;
	}
}
;

template<St DIM>
class GhostRegular_: public Ghost_<DIM> {
public:
	static const St Dim = DIM;
	typedef BoundaryCondition BC;
	typedef std::shared_ptr<BoundaryCondition> spBC;
	typedef std::shared_ptr<const BoundaryCondition> spcBC;
	typedef BoundaryIndex BI;
	typedef std::shared_ptr<BI> spBI;
	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Index_<DIM> Index;
	typedef Stencil_<Dim> Stencil;
	typedef typename Stencil::Expression Expression;
	typedef typename Stencil::spExpression spExpression;

	typedef std::function<void(const Index&)> Fun_index;
protected:
	//BI _bi;
	//spGrid _grid;
	static const St _ShapeID = 0;
public:
	GhostRegular_(spGrid g, spBI bi) :
			Ghost_<DIM>(g, bi) {
		this->_name = "GhostRegular";
	}

	// y   ____3___
	// ^  |     /4 |
	// | 0|    /   |1
	// |  |___/____|
	// |     5 2
	// O-----> x
	///
	//z
	St which_boudary_seg(const Index& index) const {
		// get seg idx in BCID
		St ABI[3][2] = { { 0, 1 }, { 2, 3 }, { 4, 5 } };
		Index n = this->_grid->n();
		for (St d = 0; d < Dim; ++d) {
			Idx res = index.value(d);
			if (res < 0) {
				return ABI[d][0];
			} else if (res >= n.value(d)) {
				return ABI[d][1];
			}
		}
		SHOULD_NOT_REACH;
		return 0;
	}

	St which_boudary_shape(const Index& index) const {
		return _ShapeID;
	}

	spcBC find(const Index& index, const std::string& vname) const {
		St seg_idx = this->which_boudary_seg(index);
		St shape_idx = this->which_boudary_shape(index);
		return this->_bi->find(shape_idx, seg_idx, vname);
	}

	bool is_ghost(const Index& index) const {
		for (St d = 0; d < Dim; ++d) {
			Idx res = index.value(d);
			if (res < 0) {
				return true;
			} else if (res >= this->_grid->n().value(d)) {
				return true;
			}
		}
		return false;
	}

	bool IS_GHOST(const Index& INDEX) const {
		for (St d = 0; d < Dim; ++d) {
			Idx res = this->_grid->_idx(INDEX.value(d));
			if (res < 0) {
				return true;
			} else if (res >= this->_grid->n().value(d)) {
				return true;
			}
		}
		return false;
	}

	bool is_normal(const Index& index) const {
		return !is_ghost(index);
	}

	bool IS_NORMAL(const Index& INDEX) const {
		return !IS_GHOST(INDEX);
	}

	Index ori_index(const Index& index) const {
		// index is a ghost index
		ASSERT(is_ghost(index));
		for (St d = 0; d < Dim; ++d) {
			Idx res = index[d];
			if (res < 0) {
				Index res(index);
				res[d] = 0;
				return res;
			} else if (res >= this->_grid->n().value(d)) {
				Index res(index);
				res[d] = this->_grid->n().value(d) - 1;
				return res;
			}
		}
		SHOULD_NOT_REACH;
		return Index();
	}

	St ori_axes(const Index& index) const {
		// index is a ghost index
		ASSERT(is_ghost(index));
		for (St d = 0; d < Dim; ++d) {
			Idx res = index.value(d);
			if (res < 0) {
				return d;
			} else if (res >= this->_grid->n().value(d)) {
				return d;
			}
		}
		SHOULD_NOT_REACH;
		return 0;
	}

	Orientation orientation(const Index& index) const {
		// index is a ghost index
		ASSERT(is_ghost(index));
		for (St d = 0; d < Dim; ++d) {
			Idx res = index.value(d);
			if (res < 0) {
				return _M_;
			} else if (res >= this->_grid->n().value(d)) {
				return _P_;
			}
		}
		SHOULD_NOT_REACH;
		return _C_;
	}

	void for_each_ghost(Fun_index fun) {
		for (typename Grid::Ijk idx = this->_grid->begin_ijk(); !idx.is_end();
				++idx) {
			if (is_ghost(idx.current())) {
				fun(idx.current());
			}
		}
	}

	spExpression substitute(spExpression spc, const std::string& vname) const {
		for (auto iter = spc->begin(); iter != spc->end();) {
			const Vt& coe = iter->second;
			const std::pair<Index, short>& key = iter->first;
			if (this->is_ghost(key.first)) {
				auto itere = iter;
				++iter;
				if (coe != 0) {
					spcBC pbc = this->find(key.first, vname);
					// spExpression spc,
					// spcBC pbc,
					// const Index& index,
					// const Vt& coe,
					// const std::string& vname
					spc = this->_substitute_ghost(spc, pbc, key.first, coe,
							vname);
				}
				spc->erase(itere);
			} else {
				++iter;
			}
		}
		return spc;
	}

//spcBC find(const Index& index, const std::string& vname) const {
//	St seg_idx = which_boudary_seg(index);
//	return this->_bi.find(_ShapeID, seg_idx, vname);
//}

//void add_bc(St shape_idx, St seg_idx, const std::string&name, spBC spbc) {
//	this->_bi.insert(shape_idx, seg_idx, name, spbc);
//}

}
;


template<St DIM>
class Order_ {
public:
	static const St Dim = DIM;

	typedef Index_<DIM> Index;

	typedef Scalar_<DIM> CenterScalar;
	typedef std::shared_ptr<CenterScalar> spCenterScalar;

	typedef BoundaryCondition BC;
	typedef std::shared_ptr<BoundaryCondition> spBC;
	typedef std::shared_ptr<const BoundaryCondition> spcBC;

	typedef BoundaryIndex BI;
	typedef std::shared_ptr<BI> spBI;
	typedef std::shared_ptr<const BI> spcBI;
	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Ghost_<DIM> Ghost;
	typedef std::shared_ptr<Ghost> spGhost;
	typedef Stencil_<DIM> Stencil;
	typedef std::shared_ptr<Stencil> spStencil;
	typedef typename Stencil::Expression Expression;
	typedef std::shared_ptr<Expression> spExpression;
	typedef Data_<spExpression, DIM> CenterExp;
	typedef std::shared_ptr<CenterExp> spCenterExp;
protected:
	spGrid _grid;
	spGhost _ghost;
public:
	Order_(spGrid grid, spGhost ghost) :
			_grid(grid), _ghost(ghost) {
	}

	virtual St get_order(const Index& index) {
		SHOULD_NOT_REACH;
		return 0;
	}

	virtual St size() const {
		SHOULD_NOT_REACH;
		return 0;
	}

	virtual Index begin() const {
		SHOULD_NOT_REACH;
		return Index();
	}

	virtual Index end() const {
		SHOULD_NOT_REACH;
		return Index();
	}
	virtual Index next(const Index& index) const {
		SHOULD_NOT_REACH;
		return Index();
	}
	virtual ~Order_() {

	}
};

template<St DIM>
class OrderRegular_: public Order_<DIM> {
public:
	static const St Dim = DIM;

	typedef Index_<DIM> Index;
	typedef Scalar_<DIM> CenterScalar;
	typedef std::shared_ptr<CenterScalar> spCenterScalar;

	typedef BoundaryCondition BC;
	typedef std::shared_ptr<BoundaryCondition> spBC;
	typedef std::shared_ptr<const BoundaryCondition> spcBC;

	typedef BoundaryIndex BI;
	typedef std::shared_ptr<BI> spBI;
	typedef std::shared_ptr<const BI> spcBI;
	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Ghost_<DIM> Ghost;
	typedef std::shared_ptr<Ghost> spGhost;
	typedef Stencil_<DIM> Stencil;
	typedef std::shared_ptr<Stencil> spStencil;
	typedef typename Stencil::Expression Expression;
	typedef std::shared_ptr<Expression> spExpression;
	typedef Data_<spExpression, DIM> CenterExp;
	typedef std::shared_ptr<CenterExp> spCenterExp;
	typedef Ijk_<DIM> Ijk;
public:
	OrderRegular_(spGrid grid, spGhost ghost) :
			Order_<DIM>(grid, ghost) {
		ASSERT(ghost->name() == "GhostRegular");
	}

	St get_order(const Index& index) {
		Index n = this->_grid->n();
		return index.i() + index.j() * n.i() + index.k() * n.i() * n.j();
	}

	St size() const {
		St res = 1;
		Index n = this->_grid->n();
		for (St d = 0; d < Dim; ++d) {
			res *= n[d];
		}
		return res;
	}

	Index begin() const {
		return this->_grid->begin_ijk().current();
	}

	Index end() const {
		return this->_grid->end_ijk().current();
	}
	Index next(const Index& index) const {
		Ijk ijk(index, this->_grid->n());
		++ijk;
		return ijk.current();
	}
};


}



#endif /* LIB_STRUCTURE_S_GHOST_HPP_ */
