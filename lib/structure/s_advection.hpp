/*
 * s_advection.hpp
 *
 *  Created on: Jun 9, 2017
 *      Author: zhou
 */

#ifndef _S_ADVECTION_HPP_
#define _S_ADVECTION_HPP_

#include <structure/s_data.hpp>
#include "s_define.hpp"
#include "s_stencil.hpp"
#include "s_equation.hpp"
#include "s_solver.hpp"

namespace structure {

/*
 * the Poisson class
 */
template<St DIM>
class Advection_: public Equation_<DIM> {
public:
	static const St Dim = DIM;
	typedef Equation_<DIM> Base;
	typedef Equation_<DIM> Equation;
	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Index_<Dim> Index;

	typedef Scalar_<DIM> Scalar;
	typedef std::shared_ptr<Scalar> spScalar;
	typedef VectorCenter_<DIM> VectorCenter;
	typedef std::shared_ptr<VectorCenter> spVectorCenter;
	typedef VectorFace_<DIM> VectorFace;
	typedef std::shared_ptr<VectorFace> spVectorFace;

	typedef Stencil_<DIM> Stencil;
	typedef std::shared_ptr<Stencil> spStencil;
	typedef typename Stencil::Expression Expression;
	typedef typename Stencil::spExpression spExpression;

	// function define ------------------------------
	typedef std::function<Vt(Vt, Vt, Vt, Vt)> Function;
	typedef std::unordered_map<std::string, Function> Functions;
	// Vaules ---------------------------------------
	typedef std::unordered_map<std::string, Vt> Values;

	// Event ----------------------------------------
	typedef Event_<DIM> Event;
	typedef std::shared_ptr<Event> spEvent;
	typedef std::unordered_map<std::string, spEvent> Events;
	typedef EventFlag_<Dim> Flag;

	typedef Ghost_<Dim> Ghost;
	typedef std::shared_ptr<Ghost> spGhost;
	typedef typename Ghost::Fun_index Fun_index;

	// time -----------------------------------------
	typedef Time_<DIM> Time;
	typedef std::shared_ptr<Time> spTime;
	// solver ---------------------------------------
	typedef Solver_<DIM> Solver;
	typedef std::shared_ptr<Solver_<DIM> > spSolver;
	typedef std::map<std::string, spSolver> Solvers;
	// Any
	typedef carpio::Any Any;
	typedef std::map<std::string, Any> AFlag;

	typedef Operation_<DIM> Operation;
protected:
	typedef Vt (*VdotNabla)(const VectorFace&, const Scalar&, Index);
	VdotNabla fun_v_dot_nabla;

	spVectorCenter _veoc;
	spVectorFace _veof;

public:
	Advection_(spGrid spg) :
			Base(spg) {
		_default_setup();
	}

protected:
	void _default_setup() {
		spEvent pse(new Flag("stand_alone", 1));
		this->_events["stand_alone"] = pse;
		// default velocity is 0
		Function fun = [](Vt, Vt, Vt, Vt) {return 0;};
		std::string veo_name[] = { "u", "v", "w" };
		std::string vf_name[] = { "uf", "vf", "wf" };
		FOR_EACH_DIM
		{
			this->_functions[veo_name[d]] = fun;
			this->_functions[vf_name[d]] = fun;
		}
		//
		this->_functions["phi"] = fun;
		fun_v_dot_nabla = &Operation::VdotNabla_upwind1;
	}

	int initial() {
		ASSERT(this->has_event("time_term"));

		// if stand alone
		if (this->has_event("stand_alone")) {   // if stand alone
			this->_ghost = spGhost(
					new GhostRegular_<Dim>(this->_grid, this->_bi));

			spGrid spg = this->_grid;
			this->_css["phi"] = spScalar(new Scalar(spg));
			this->set_CS("phi", this->_functions["phi"]);
			std::string veo_name[] = { "u", "v", "w" };
			std::string veof_name[] = { "uf", "vf", "wf" };
			spScalar vc[] = { nullptr, nullptr, nullptr };
			spScalar vf[] = { nullptr, nullptr, nullptr };
			FOR_EACH_DIM
			{
				vc[d] = spScalar(new Scalar(spg));
				this->_css[veo_name[d]] = vc[d];
				this->set_CS(veo_name[d], this->_functions[veo_name[d]]);
				this->apply_bc(veo_name[d]);
			}
			FOR_EACH_DIM
			{
				vf[d] = spScalar(new Scalar(spg));
				this->_css[veof_name[d]] = vf[d];
				// this->set_CS(veof_name[d], this->_functions[veof_name[d]]);
			}
			this->_veoc = spVectorCenter(new VectorCenter(vc[0], vc[1], vc[2]));
			this->_veof = spVectorFace(new VectorFace(vf[0], vf[1], vf[2]));

			Operation::InterpolateC2F(this->_veoc, this->_veof);

		} else {

		}
		return -1;
	}

public:
	int set_advection_scheme(const std::string& name) {
		this->_aflags["AdvectionScheme_name"] = name;
		if (name == "upwind1") {
			fun_v_dot_nabla = &Operation::VdotNabla_upwind1;
			return 1;
		}
		if (name == "center") {
			fun_v_dot_nabla = &Operation::VdotNabla_center;
			return 1;
		}
		if (name == "center4") {
			fun_v_dot_nabla = &Operation::VdotNabla_center4;
			return 1;
		}
		if (name == "upwind2") {
			fun_v_dot_nabla = &Operation::VdotNabla_upwind2;
			return 1;
		}
		if (name == "VanLeer") {
			fun_v_dot_nabla = &Operation::VdotNabla_TVD_VanLeer;
			return 1;
		}
		if (name == "superbee") {
			fun_v_dot_nabla = &Operation::VdotNabla_TVD_superbee;
			return 1;
		}
		if (name == "WAHYD") {
			fun_v_dot_nabla = &Operation::VdotNabla_TVD_WAHYD;
			return 1;
		}
		if (name == "QUICK") {
			fun_v_dot_nabla = &Operation::VdotNabla_QUICK;
			return 1;
		}
		if (name == "CUI") {
			fun_v_dot_nabla = &Operation::VdotNabla_CUI;
			return 1;
		}
		SHOULD_NOT_REACH;
		return 0;
	}

	int run_one_step(St step) {
		// debug
		int stepo = 0;
		// restrict time

		// apply boundary condition : veo
		const std::string veo_name[] = { "u", "v", "w" };
		FOR_EACH_DIM
		{
			this->apply_bc(veo_name[d]);
		}
		this->apply_bc("phi");

		spVectorFace spvf = this->_veof;
		Vt dt = this->_time->dt();
		//std::cout << "dt = " << dt << std::endl;

		Grid& grid = *(this->_grid);
		Scalar& sphi = *(this->_css["phi"]);
		for (typename Grid::Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
			Index& index = ijk.current();
			//index.show();
			Vt adv  = fun_v_dot_nabla(*spvf, sphi, index);
			Vt phin = sphi(index);
			sphi(index) = phin - dt * adv;
		}
		return -1;
	}

}
;

}

#endif /* LIB_STRUCTURE_S_ADVECTION_HPP_ */
