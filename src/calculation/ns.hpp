#ifndef _NS_H_
#define _NS_H_

#include "calculation_define.hpp"
#include "expression.hpp"
#include "utility/any.hpp"
#include "equation.hpp"
#include "io/mmio.h"

#include <vector>
#include <memory>
#include "utility/clock.h"

//#define __Debug__


namespace carpio {


//
//   ∇⋅U = 0
// without surface tension
//  d U      1
// ----- = ----- ( -∇P + mu*(∇U+∇UT) ) + Source(U)
//  d t     rho
// with surface tension
//  d U      1
// ----- = ----- ( -∇P + (mu (∇U+∇UT)) + delta * kappa ( δs n)) + Source(U)
//  d t     rho

template<typename COO_VALUE, typename VALUE, int DIM>
class NS_: public Equation_<COO_VALUE, VALUE, DIM> {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef Equation_<COO_VALUE, VALUE, DIM> Base;
	typedef Equation_<COO_VALUE, VALUE, DIM>& ref_Base;

	typedef NS_<COO_VALUE, VALUE, DIM> Self;
	typedef NS_<COO_VALUE, VALUE, DIM>& ref_Self;
	typedef NS_<COO_VALUE, VALUE, DIM>* pSelf;
	typedef const NS_<COO_VALUE, VALUE, DIM>* const_pSelf;

	typedef typename Base::vt vt;
	typedef typename Base::cvt cvt;

	typedef typename Base::Variables Variables;

	typedef typename Base::pDomain pDomain;
	typedef typename Base::spEvent spEvent;

	typedef typename Base::Grid Grid;
	typedef typename Base::pGrid pGrid;
	typedef typename Base::const_pGrid const_pGrid;

	typedef typename Base::Node Node;
	typedef typename Base::pNode pNode;
	typedef typename Base::const_pNode const_pNode;

	typedef typename Base::Exp Exp;
	typedef typename Base::pExp pExp;
	typedef typename Base::spExp spExp;
	typedef typename Base::spFace spFace;

	typedef typename Base::Face Face;
	typedef typename Base::pFace pFace;
	typedef typename Base::const_pFace const_pFace;

	typedef typename Base::Flag Flag;

	typedef ArrayListT<spExp> Arr_spExp;
	typedef ArrayListT<Any> Arr_Any;
	typedef typename Base::Map_FA Map_FA;
	typedef typename Base::pMap_FA pMap_FA;
	typedef typename Base::Iter_FA Iter_FA;
	typedef typename Base::const_Iter_FA const_Iter_FA;
	typedef typename Base::Pair_FA Pair_FA;
	typedef typename Base::Ret_FA Ret_FA;

	typedef typename Base::Mat Mat;
	typedef typename Base::Arr Arr;

	typedef typename Base::BoundaryCondition BoundaryCondition;
	typedef typename Base::pBoundaryCondition pBoundaryCondition;
	typedef typename Base::const_pBoundaryCondition const_pBoundaryCondition;

	typedef typename Base::Function Function;
	typedef typename Base::Fun_pNode Fun_pNode;
	typedef typename Base::Fun_pNode_spExp Fun_pNode_spExp;

public:
	/**
	 * @brief constructor
	 *
	 * @param  [pDomain] a point of domain you want to use.
	 *
	 * @param  [dt]      the difference of the time, default value is -1 which means
	 *                   the equation doesn't have time term.
	 * @param  [maxstep] the max time step. default value is 0.
	 *                   the equation doesn't have time term.
	 *
	 **/
	NS_(pDomain pd) :
			Base(pd) {
		//stand alone
		this->set_uniform_rho();
		this->set_uniform_mu();
		this->_construct_c_idx();
	}

	int _initial_ns() {
		ASSERT(this->is_stand_alone());
		// set time step
		ASSERT(this->has_time_term());

		this->_pdomain->resize_data( //
				this->max_idx(this->_vars_c) + 1, //
				0, //
				0, //
				this->max_idx(this->_vars_ut) + 1);

		// initial all the variables
		for (typename Variables::iterator it = this->_vars_c.begin();
				it != this->_vars_c.end(); ++it) {
			std::string key = it->first;
			if (key == "rho" || key == "mu") {
				init(key, 1.0);
			} else {
				init(key);  //defualt is 0.0
			}
		}

		this->_build_fun_get();

		return -1;
	}

	//virtual int initial() {
	//	std::cout << "  ns: initial \n";
	//	SHOULD_NOT_REACH;
	//	return -1;
	//}

	int finalize() {
		this->_delete_FA_on_leaf();
		this->_delete_NE_on_leaf();
		return -1;
	}

	/**
	 *  @brief set uniform rho
	 */
	void set_uniform_rho(vt rho = 1) {     //default
		spEvent pse(new Flag("uniform_rho", 1));
		this->_events["uniform_rho"] = pse;
		this->_values["uniform_rho"] = rho;
	}
	void unset_uniform_rho() {
		this->_events.erase("uniform_rho");
		this->_values.erase("uniform_rho");
	}
	/**
	 *  @brief set uniform mu
	 */
	void set_uniform_mu(vt v = 1) {     //default
		spEvent pse(new Flag("uniform_mu", 1));
		this->_events["uniform_mu"] = pse;
		this->_values["uniform_mu"] = v;
	}
	void unset_uniform_mu() {
		this->_events.erase("uniform_mu");
		this->_values.erase("uniform_mu");
	}


	/**
	 * @brief set value on center of the node
	 */
	void set_mu(Function fun) {
		this->unset_uniform_mu();
		//
		this->_functions["set_mu"] = fun;

		//this->set_val(fun, this->_vars_c["beta"]);
	}

	/**
	 * @brief set value on center of the node
	 */
	void set_rho(Function fun) {
		this->unset_uniform_rho();
		//
		this->_functions["set_rho"] = fun;

		//this->set_val(fun, this->_vars_c["beta"]);
	}

	void set_veo(Axes a, Function fun) {
		ASSERT(a < Dim);
		std::string _arrn[3] = { "u", "v", "w" };
		std::string key = "set_" + _arrn[int(a)];
		this->_functions[key] = fun;
		//this->set_val(fun, this->_vars_c["f"]);
	}

	void set_source(Axes a, Function fun) {
		ASSERT(a < Dim);
		std::string _arrn[3] = { "x", "y", "z" };
		std::string key = "set_s" + _arrn[int(a)];
		this->_functions[key] = fun;
		//this->set_val(fun, this->_vars_c["f"]);
	}

	void set_p(Function fun) {
		std::string key = "set_p";
		this->_functions[key] = fun;
	}

	void init(const std::string& vn, vt initv = 0.0) {
		if (this->has_event("uniform_" + vn)) {
			return;
		}
		if (this->has_function("set_" + vn)) {
			Function f = this->_functions["set_" + vn];
			this->set_val(f, this->_vars_c[vn]);
		} else {
			Function fun = [&initv](Float, Float, Float) {return initv;};
			std::string key = "set_" + vn;
			this->_functions[key] = fun;
			this->set_val(fun, this->_vars_c[vn]);
		}
	}

	St u_idx() const{
		return this->get_var_center_idx("u");
	}
	St v_idx() const{
		ASSERT(Dim>=2);
		return this->get_var_center_idx("v");
	}
	St w_idx() const{
		ASSERT(Dim>=3);
		return this->get_var_center_idx("w");
	}

protected:
	/*
	 * stand alone constructor
	 * constructor a NS class just used for solver only one group of NS equation
	 * with out surface tension
	 */
	void _construct_c_idx() {
		// velocity and source term
		switch (Dim) {
		case 1: {
			this->_vars_c["u"] = 0;
			this->_vars_c["sx"] = 1;
			this->_vars_c["p"] = 2;
			if (!this->has_event("uniform_rho")) {
				this->_vars_c["rho"] = 3;
			}
			if (!this->has_event("uniform_mu")) {
				this->_vars_c["mu"] = 4;
			}
		}
			break;
		case 2: {
			this->_vars_c["u"] = 0;
			this->_vars_c["v"] = 1;
			this->_vars_c["sx"] = 2;
			this->_vars_c["sy"] = 3;
			this->_vars_c["p"] = 4;
			if (!this->has_event("uniform_rho")) {
				this->_vars_c["rho"] = 5;
			}
			if (!this->has_event("uniform_mu")) {
				this->_vars_c["mu"] = 6;
			}
		}
			break;
		case 3: {
			this->_vars_c["u"] = 0;
			this->_vars_c["v"] = 1;
			this->_vars_c["w"] = 2;
			this->_vars_c["sx"] = 3;
			this->_vars_c["sy"] = 4;
			this->_vars_c["sz"] = 5;
			this->_vars_c["p"] = 6;
			if (!this->has_event("uniform_rho")) {
				this->_vars_c["rho"] = 7;
			}
			if (!this->has_event("uniform_mu")) {
				this->_vars_c["mu"] = 8;
			}
		}
			break;
		default:
			SHOULD_NOT_REACH;
			break;
		}
		// untype variables on the node
		this->_vars_ut["FA"] = 0;
		this->_vars_ut["NE"] = 1;
	}

	void _build_fun_get() {
		this->_get_node_spexp = [this](pNode pn) {
			ASSERT(pn!=nullptr);
			utPointer utp = pn->utp(this->_vars_ut["NE"]);
			Arr_spExp& arr_spexp = CAST_REF(Arr_spExp*, utp);
			return arr_spexp[0];
		};
	}

	void _update_face_adv_dif(St idx_val, St idx_gra, St idx_adv, St idx_dif) {
		Grid& grid = this->_pdomain->grid();
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			utPointer utp = pn->utp(this->_vars_ut["FA"]);
			ASSERT(utp != nullptr);
			Map_FA& mfe = CAST_REF(pMap_FA, utp);
			for (typename Map_FA::iterator it = mfe.begin(); it != mfe.end();
					it++) {
				Arr_Any& arr_any = it->second;
				spExp exp_val = any_cast<spExp>(arr_any[idx_val]);
				spExp exp_gra = any_cast<spExp>(arr_any[idx_gra]);
				arr_any[idx_adv] = _face_adv((*(it->first)), *exp_val);
				arr_any[idx_dif] = _face_dif((*(it->first)), *exp_val,
						*exp_gra);
			}
		}
	}

	void _update_face_adv(St idx_val, St idx_adv) {
		Grid& grid = this->_pdomain->grid();
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			utPointer utp = pn->utp(this->_vars_ut["FA"]);
			ASSERT(utp != nullptr);
			Map_FA& mfe = CAST_REF(pMap_FA, utp);
			for (typename Map_FA::iterator it = mfe.begin(); it != mfe.end();
					it++) {
				Arr_Any& arr_any = it->second;
				spExp exp_val = any_cast<spExp>(arr_any[idx_val]);
				arr_any[idx_adv] = _face_adv((*(it->first)), *exp_val);
			}
		}
	}

	void _update_face_dif(St idx_val, St idx_gra, St idx_dif) {
		Grid& grid = this->_pdomain->grid();
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			utPointer utp = pn->utp(this->_vars_ut["FA"]);
			ASSERT(utp != nullptr);
			Map_FA& mfe = CAST_REF(pMap_FA, utp);
			for (typename Map_FA::iterator it = mfe.begin(); it != mfe.end();
					it++) {
				Arr_Any& arr_any = it->second;
				spExp exp_val = any_cast<spExp>(arr_any[idx_val]);
				spExp exp_gra = any_cast<spExp>(arr_any[idx_gra]);
				arr_any[idx_dif] = _face_dif((*(it->first)), *exp_val,
						*exp_gra);
			}
		}
	}

	/*
	 * add up exp on FA
	 * pn node
	 * st idx on FA
	 */
	spExp _sum_face_exp(pNode pn, St idx) {
		//assert(pn->d());
		spExp exp(new Exp());
		utPointer& utp = pn->utp(this->_vars_ut["FA"]);
		ASSERT(utp != nullptr);
		Map_FA& mfe = CAST_REF(pMap_FA, utp);
		for (typename Map_FA::iterator it = mfe.begin(); it != mfe.end();
				it++) {
			Arr_Any& arr_any = it->second;
			spExp exp_f = any_cast<spExp>(arr_any[idx]);
			exp->plus(*exp_f);
		}
		return exp;
	}

	spExp _face_adv(Face& f, const Exp& exp_val) {
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(f.dir(), o, a);

		// get velocity on face
		St v_idx = this->veo_c_idx(a);
		vt veo_f = exp_val.substitute(v_idx);

		spExp sp_exp(new Exp());

		// times area
		int sign = IsP(o) ? 1 : -1;
		vt area = f.area();
		// first order up wind -----------------------------
		pNode pC = _find_C(&f, veo_f);
		sp_exp->insert(sign * area * veo_f, pC, 1.0);
		// center difference -------------------------------
		//exp_adv.plus(exp_val);
		//exp_adv.times(sign * area * veo_f);
#ifdef __Debug__
		pNode pn = f.pori();
		if (pn->is_in_on(__x__, __y__)) {
			std::cout << "_face_adv_  ========================\n";
			std::cout << " Direction   " << ToString(dir) << "\n";
			std::cout << " sign        " << sign << "\n";
			std::cout << " area        " << area << "\n";
			std::cout << " coe         " << sign * area * veo_f << "\n";
		}
#endif
		return sp_exp;
	}

	/*
	 * Input: 1  face
	 *        2  expression value on face
	 *        3  expression gradient on face
	 * Out  : 4  expression of diffusion on face
	 */
	spExp _face_dif(Face& f, const Exp& exp_val, const Exp& exp_gra) {
		// mu * ▽(phi(x,y))
		vt mu_f = 1.0;
		if (this->has_event("uniform_mu")) {
			mu_f = this->_values["uniform_mu"];
		} else {
			mu_f = exp_val.substitute(this->_vars_c["mu"]);
		}

		int sign = IsFacePDirection(f.dir()) ? 1 : -1;
		spExp sp_exp(new Exp(exp_gra));
		vt area = f.area();
		sp_exp->times(sign * area * mu_f);
#ifdef __Debug__
		pNode pn = f.pori();
		if (pn->is_in_on(__x__, __y__)) {
			std::cout << "_face_dif_  ========================\n";
			std::cout << "Direction   " << ToString(f.dir()) << "\n";
			std::cout << " sign        " << sign << "\n";
			std::cout << " mu_f        " << mu_f << "\n";
			std::cout << " coe         " << sign * mu_f << "\n";
			exp_gra.show();
		}
#endif
		return sp_exp;
	}

	vt _CFL_number(pNode pn, vt dt) {
		vt veo[Dim];
		vt cfl[Dim];
		for (St i = 0; i < Dim; i++) {
			veo[i] = pn->cdva(this->veo_c_idx(i));
			cfl[i] = veo[i] * dt / pn->d(ToAxes(i));
		}
		return Max(cfl, Dim);
	}

	pNode _find_C(pFace pface, Float veo_f) {
		pNode pori = pface->pori();
		pNode pnei = pface->pnei();
		Direction dir = pface->dir();
		if (IsFacePDirection(dir)) {
			if (veo_f > 0) {
				return pori;  //o
			} else {
				return pnei;
			}
		} else {
			if (veo_f > 0) {
				return pnei;
			} else {
				return pori;
			}
		}
		SHOULD_NOT_REACH;
		return nullptr;
	}

public:
	/*
	 * \brief get idx functions
	 */
	St c_idx(const std::string& name) const {
		return this->_vars_c[name];
	}

	St veo_c_idx(const St& dim) const {
		ASSERT(dim < Dim);
		std::string _arrn[3] = { "u", "v", "w" };
		auto it = this->_vars_c.find(_arrn[dim]);
		if (it != this->_vars_c.end()) {
			return it->second;
		} else {
			SHOULD_NOT_REACH;
			return 0;
		}
	}
	St veo_c_idx(Axes a) const {
		return veo_c_idx(St(a));
	}
	St pressure_idx() const {
		return this->_vars_c["p"];
	}

}
;

}

#endif
