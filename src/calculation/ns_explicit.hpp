#ifndef _NS_EXPLICIT_HPP_
#define _NS_EXPLICIT_HPP_

#include "ns.hpp"
#include <glog/logging.h>

#define __Debug__

#define __x__  -0.4
#define __y__  0.45

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM>
class NS_explicit_: public NS_<COO_VALUE, VALUE, DIM> {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef NS_<COO_VALUE, VALUE, DIM> Base;
	typedef NS_<COO_VALUE, VALUE, DIM>& ref_Base;

	typedef NS_explicit_<COO_VALUE, VALUE, DIM> Self;
	typedef NS_explicit_<COO_VALUE, VALUE, DIM>& ref_Self;
	typedef NS_explicit_<COO_VALUE, VALUE, DIM>* pSelf;
	typedef const NS_explicit_<COO_VALUE, VALUE, DIM>* const_pSelf;

	typedef typename Base::vt vt;
	typedef typename Base::cvt cvt;

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
protected:
	static const St _NE_veo_star = 0;
	static const St _FA_val = 0;
	static const St _FA_gra = 1;
	static const St _FA_adv = 2;
	static const St _FA_dif = 3;
public:
	NS_explicit_(pDomain pd) :
			Base(pd) {
		this->_construct_c_idx_tmp();
	}

	void set_solve_para(vt tol = 1e-6, int max_iter = 1000) {
		// solve para for pressure equation
		this->_values["tolerance"] = tol;
		this->_values["max_iter"] = max_iter;
	}

	int initial() {
		this->_initial_ns();

		this->_new_FA(4);
		this->_update_face_val_gra(_FA_val, _FA_gra);
		this->_update_face_adv_dif(_FA_val, _FA_gra, _FA_adv,  _FA_dif);
		this->_build_NE_on_leaf(2);
		return -1;

	}

	int run_one_step(St step) {
		std::cout << "ns one step " << step << "\n";
		_update_veo_star_grid(_NE_veo_star);
		// veo star on center
		for (St i = 0; i < Dim; ++i) {
			_solve_veo_star_c(ToAxes(i));
			_ghost_veo_star_c(ToAxes(i));
		}
		// solve pressure
		_solve_pressure();
		_ghost_pressure();
		// correct veo on center
		for (St i = 0; i < Dim; ++i) {
			_correct_veo_c(ToAxes(i));
		}
		return -1;
	}

	void _update_veo_star_grid(St idx) {
		pGrid pgrid = this->_pdomain->pgrid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			utPointer utp = pn->utp(this->_vars_ut["NE"]);
			Arr_spExp& arr_spexp = CAST_REF(Arr_spExp*, utp);
			arr_spexp[_NE_veo_star] = _update_exp_veo_star_node(pn);
		}
	}

	void _solve_veo_star_c(Axes a) {

	}

	void _ghost_veo_star_c(Axes a) {

	}

	void _solve_pressure() {

	}

	void _ghost_pressure() {

	}

	void _correct_veo_c(Axes a) {

	}

protected:
	void _construct_c_idx_tmp() {
		St midx = this->max_idx(this->_vars_c);
		// velocity and source term
		switch (Dim) {
		case 1: {
			this->_vars_c["us"] = midx + 1;
		}
			break;
		case 2: {
			this->_vars_c["us"] = midx + 1;
			this->_vars_c["vs"] = midx + 2;
		}
			break;
		case 3: {
			this->_vars_c["us"] = midx + 1;
			this->_vars_c["vs"] = midx + 2;
			this->_vars_c["ws"] = midx + 3;
		}
			break;
		default:
			SHOULD_NOT_REACH;
			break;
		}
	}

	/*
	 * this is the main function of the predict step
	 */
	spExp _update_exp_veo_star_node(pNode pn) {
		spExp exp(new Exp());
		// Traversal all the faces
		// cfl number check
		vt cfl = this->_CFL_number(pn, this->_timestep->dt());
		ASSERT(cfl < 1.0);
		// calculate the term separately
		// 1 ----- advection term
		spExp adv = this->_sum_face_exp(pn, _FA_adv);
		// divide volume and negative;
		adv->times(-1.0 / pn->volume());

		// 2 ----- diffusion term
		spExp dif = this->_sum_face_exp(pn, _FA_dif);
#ifdef __Debug__
		if (pn->is_in_on(__x__, __y__)) {
			LOG(INFO) << "u star mid =================\n";
			dif->show();
		}
#endif
		//divide volume
		dif->times(1.0 / pn->volume());

		// 3 ----- source term
		spExp src(new Exp());
		//_node_exp_src_term(pn, ai, exp_src);
		//   add up diffusion and source
		dif->plus(*src);
		if(this->has_event("uniform_rho")){
			dif->times(1.0 / this->_values["uniform_rho"]);
		}else{
			dif->times(1.0 / pn->cdva(this->_vars_c["rho"]));
		}
		dif->plus(*adv);   // add advection term
		dif->times(this->_timestep->dt());      // time
		// 4 u term
		exp->insert(1.0, pn, 1.0);
		exp->plus(*dif);
		return exp;
	}


}
;

}

#endif
