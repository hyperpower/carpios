/*
 * timestep.hpp
 *
 *  Created on: Oct 12, 2016
 *      Author: zhou
 */

#ifndef _TIMESTEP_HPP_
#define _TIMESTEP_HPP_

#include "calculation_define.hpp"
#include "expression.hpp"
#include "event.hpp"
#include "utility/clock.h"
#include "utility/format.h"
#include "io/mmio.h"

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM>
class TimeStep_Explicit_;

template<typename COO_VALUE, typename VALUE, int DIM>
class TimeStep_Implicit_;

template<typename COO_VALUE, typename VALUE, int DIM>
class TimeStep_CrankNicolson_;

template<typename COO_VALUE, typename VALUE, int DIM>
class TimeStep_ {
public:
	static const St Dim = DIM;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef TimeStep_<cvt, vt, Dim> Self;
	typedef TimeStep_<cvt, vt, Dim>& ref_Self;
	typedef TimeStep_<cvt, vt, Dim>* pSelf;
	typedef std::shared_ptr<Self> spTimeStep;

	typedef Domain_<cvt, vt, Dim> Domain;
	typedef Domain_<cvt, vt, Dim>& ref_Domain;
	typedef const Domain_<cvt, vt, Dim>& const_ref_Domain;
	typedef Domain_<cvt, vt, Dim>* pDomain;
	typedef const Domain_<cvt, vt, Dim>* const_pDomain;
	typedef BoundaryCondition_<cvt, vt> BoundaryCondition;
	typedef BoundaryCondition_<cvt, vt>* pBoundaryCondition;
	typedef const BoundaryCondition_<cvt, vt>* const_pBoundaryCondition;

	typedef Grid_<cvt, vt, Dim> Grid;
	typedef Grid_<cvt, vt, Dim> *pGrid;
	typedef const Grid_<cvt, vt, Dim> * const_pGrid;
	typedef Ghost_<cvt, vt, Dim> Ghost;
	typedef const Ghost_<cvt, vt, Dim> const_Ghost;
	typedef Ghost_<cvt, vt, Dim>* pGhost;
	typedef const Ghost_<cvt, vt, Dim>* const_pGhost;
	typedef typename Ghost::GhostNode GhostNode;
	typedef Cell_<cvt, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<vt, Dim> Data;
	typedef Data *pData;
	typedef Node_<cvt, vt, Dim> Node;
	typedef Node_<cvt, vt, Dim> *pNode;
	typedef const Node_<cvt, vt, Dim> *const_pNode;

	typedef Face_<Node, pNode> Face;
	typedef Face_<Node, pNode> *pFace;
	typedef const Face_<Node, pNode> * const_pFace;
	typedef Shape_<cvt, Dim> Shape;
	typedef Shape_<cvt, Dim>* pShape;

	typedef Expression_<cvt, vt, Dim> Exp;
	typedef Exp* pExp;
	typedef std::shared_ptr<Exp> spExp;
	typedef std::shared_ptr<Face> spFace;
	typedef std::function<spExp(pNode)> Fun_get_spExp;
	typedef std::function<void(spExp, pNode)> Fun_set_spExp;

protected:
	St _max_step;
	vt _dt; // dt

	ArrayListV<St> _arr_idx;
	St _inner_step;

	vt _tau; // coefficient

	//
	Fun_get_spExp _get_ex;  // get explicit expression
	Fun_get_spExp _get_im;  // get implicit expression

	Fun_set_spExp _set_rhs; // right hand side is constant
	Fun_set_spExp _set_lhs; // left hand side is expression

public:
	/*
	 * constructor
	 */
	TimeStep_() :
			_inner_step(0), _tau(1) {
		_arr_idx.resize(2);
	}

	virtual ~TimeStep_() {
	}

	vt dt() const {
		return _dt;
	}

	St max_step() const {
		return this->_max_step;
	}

	void set_dt(vt dt) {
		this->_dt = dt;
	}

	void set_max_step(St ms) {
		this->_max_step = ms;
	}

	vt tau() const {
		return _tau;
	}

	bool is_end() const {
		return _inner_step == _arr_idx.size() - 1;
	}

	St idx_n() const {
		return _arr_idx[0];
	}

	St idx_np() const {
		ASSERT(1 < _arr_idx.size());
		return _arr_idx[1];
	}

	St idx_nn() const {
		return _arr_idx[_arr_idx.size() - 1];
	}

	St idx_inner_step(St i) const {
		ASSERT(i < _arr_idx.size());
		return _arr_idx[i];
	}

	St idx_begin() const {
		return _arr_idx[0];
	}

	St idx_end() const {
		return _arr_idx[_arr_idx.size() - 1];
	}

	virtual St idx_old() const {
		SHOULD_NOT_REACH;
		return 0;
	}
	virtual St idx_new() const {
		SHOULD_NOT_REACH;
		return 0;
	}

	void inner_advance() {
		_inner_step += 1;
	}

	void next_step() {
		_inner_step = 0;
	}

	virtual spExp exp_rhs() {
		return spExp(nullptr);
	}

	virtual spExp exp_lhs() {
		return spExp(nullptr);
	}

	virtual bool do_solve() const {
		return false;
	}

	virtual spExp new_exp(pNode pn, spExp ex_spexp, spExp im_spexp = nullptr,
			vt ex = 1, vt im = 1) {

	}

	void set_unknow_idx(St idx) {
		this->_arr_idx[0] = idx;
	}
	St get_unknow_idx() const {
		return this->_arr_idx[0];
	}
	St size_inner_step() const {
		return this->_arr_idx.size();
	}
	void set_idx(St o_max_idx) {
		//the array begins from 1
		//the arr[0] is idx unknown
		for (St i = 1; i < _arr_idx.size(); i++) {
			this->_arr_idx[i] = o_max_idx + i;
		}
	}

	void set_tau(vt tau) {
		this->_tau = tau;
	}

	static spTimeStep Creator(const std::string& name = "explicit") {
		std::map<std::string, std::function<spTimeStep()> > map;
		typedef std::pair<std::string, std::function<spTimeStep()> > Term;
		map.insert(
				Term("explicit",
						[]() {return spTimeStep(new TimeStep_Explicit_<cvt,vt,Dim>());}));
		map.insert(
				Term("implicit",
						[]() {return spTimeStep(new TimeStep_Implicit_<cvt,vt,Dim>());}));
		map.insert(
						Term("CrankNicolson",
								[]() {return spTimeStep(new TimeStep_CrankNicolson_<cvt,vt,Dim>());}));

		return map.at(name)();
	}

};

template<typename COO_VALUE, typename VALUE, int DIM>
class TimeStep_Explicit_: public TimeStep_<COO_VALUE, VALUE, DIM> {
public:
public:
	typedef TimeStep_<COO_VALUE, VALUE, DIM> Base;
	typedef typename Base::vt vt;
	typedef typename Base::cvt cvt;

	typedef typename Base::Exp Exp;
	typedef typename Base::pExp pExp;
	typedef typename Base::spExp spExp;
	typedef typename Base::spFace spFace;
	typedef typename Base::Fun_get_spExp Fun_get_spExp;

	typedef typename Base::Node Node;
	typedef typename Base::pNode pNode;
	typedef typename Base::const_pNode const_pNode;

	TimeStep_Explicit_() :
			Base() {
	}

	bool do_solve() const {
		return false;
	}

	virtual St idx_old() const {
		return this->_arr_idx[0];
	}
	virtual St idx_new() const {
		return this->_arr_idx[1];
	}

	virtual spExp new_exp(pNode pn, spExp ex_spexp, spExp im_spexp = nullptr,
			vt ex = 1, vt im = 1) {
		// the method implict will treat explicit term as implict
		vt coe = this->_tau / this->_dt;
		St idx_n = this->idx_n();
		St v_n = pn->cda(idx_n);

		vt vexp = ex_spexp->substitute(idx_n) - coe * v_n;

		spExp res = spExp(new Exp());

		res->plus(coe, pn, 1);
		res->plus(vexp, pn, 0);
		return res;

	}

	virtual ~TimeStep_Explicit_() {
	}
	;

};

template<typename COO_VALUE, typename VALUE, int DIM>
class TimeStep_Implicit_: public TimeStep_<COO_VALUE, VALUE, DIM> {
public:
	typedef TimeStep_<COO_VALUE, VALUE, DIM> Base;
	typedef typename Base::vt vt;
	typedef typename Base::cvt cvt;

	typedef typename Base::Exp Exp;
	typedef typename Base::pExp pExp;
	typedef typename Base::spExp spExp;
	typedef typename Base::spFace spFace;
	typedef typename Base::Fun_get_spExp Fun_get_spExp;

	typedef typename Base::Node Node;
	typedef typename Base::pNode pNode;
	typedef typename Base::const_pNode const_pNode;

	TimeStep_Implicit_() :
			Base() {
	}

	bool do_solve() const {
		return true;
	}

	virtual St idx_old() const {
		return this->_arr_idx[0];
	}
	virtual St idx_new() const {
		return this->_arr_idx[1];
	}

	virtual spExp new_exp(pNode pn, spExp ex_spexp, spExp im_spexp = nullptr,
			vt ex = 1, vt im = 1) {
		ASSERT(im_spexp ==nullptr);
		// the method implict will treat explicit term as implict
		spExp res = spExp(new Exp(*ex_spexp));
		vt coe = this->_tau * pn->volume() / this->_dt;
		St idx_old = this->idx_old();
		vt v_n = pn->cd(idx_old);

		res->times(-1.0);
		res->plus(-coe * v_n, pn, 0); //negative add
		res->plus(coe, pn, 1);

		//if (pn->is_in_on(__X__, __Y__)) {
		//	std::cout << "---------------------------\n";
		//	res->show();
		//	std::cout<<" v n" << v_n<<std::endl;
		//}

		return res;

	}

	virtual ~TimeStep_Implicit_() {
	}
	;

};

template<typename COO_VALUE, typename VALUE, int DIM>
class TimeStep_CrankNicolson_: public TimeStep_<COO_VALUE, VALUE, DIM> {
public:
	typedef TimeStep_<COO_VALUE, VALUE, DIM> Base;
	typedef typename Base::vt vt;
	typedef typename Base::cvt cvt;

	typedef typename Base::Exp Exp;
	typedef typename Base::pExp pExp;
	typedef typename Base::spExp spExp;
	typedef typename Base::spFace spFace;
	typedef typename Base::Fun_get_spExp Fun_get_spExp;

	typedef typename Base::Node Node;
	typedef typename Base::pNode pNode;
	typedef typename Base::const_pNode const_pNode;

	TimeStep_CrankNicolson_() :
			Base() {
	}

	bool do_solve() const {
		return true;
	}

	virtual St idx_old() const {
		return this->_arr_idx[0];
	}
	virtual St idx_new() const {
		return this->_arr_idx[1];
	}

	virtual spExp new_exp(pNode pn, spExp ex_spexp, spExp im_spexp = nullptr,
			vt ex = 1, vt im = 1) {
		ASSERT(im_spexp ==nullptr);
		// the method implict will treat explicit term as implict
		spExp res = spExp(new Exp(*ex_spexp));
		vt coe = this->_tau * pn->volume() / this->_dt;
		St idx_old = this->idx_old();
		vt v_n = pn->cd(idx_old);

		res->times(0.5);    //coe1

		vt vex_exp = ex_spexp->substitute(idx_old);

		res->plus(0.5 * vex_exp, pn, 0);

		res->times(-1.0);
		res->plus(-coe * v_n, pn, 0); //negative add
		res->plus(coe, pn, 1);

		return res;

	}

	virtual ~TimeStep_CrankNicolson_() {
	}
	;

};

}

#endif /* CALCULATION_TIMESTEP_HPP_ */
