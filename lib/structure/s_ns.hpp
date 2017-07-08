#ifndef _S_NS_HPP
#define _S_NS_HPP

#include <structure/s_data.hpp>
#include "s_define.hpp"
#include "s_stencil.hpp"
#include "s_equation.hpp"
#include "s_solver.hpp"
#include "s_vector.hpp"
#include "s_io_plotly.hpp"
#include "s_poisson.hpp"
#include "s_operation.hpp"

namespace structure {

template<St DIM>
class NS_: public Equation_<DIM> {
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

	// Variables ------------------------------------
	typedef std::map<std::string, spScalar> CenterScalars;
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
	typedef std::shared_ptr<Solver_<DIM>> spSolver;
	typedef std::map<std::string, spSolver> Solvers;

	typedef Poisson_<DIM> Poisson;
	typedef std::shared_ptr<Poisson_<DIM>> spPoisson;

	typedef carpio::Any Any;
	typedef Operation_<DIM> Operation;

protected:
	spVectorCenter _veoc;
	spVectorFace _veof;

	spPoisson _projection;

	typedef Vt (*VdotNabla)(const VectorFace&, const Scalar&, Index);
	VdotNabla fun_v_dot_nabla;

	Function _fun0;

	std::array<std::string, 3> _nv;
	std::array<std::string, 3> _nvf;

public:
	NS_(spGrid spg, spGhost pgh = nullptr) :
			Base(spg, pgh) {
		_default_setup();
	}

	/**
	 *  @brief set uniform beta
	 */
	void set_uniform_rho(Vt s = 1) {     //default
		spEvent pse(new Flag("uniform_rho", 1));
		this->_events["uniform_rho"] = pse;
		this->_values["uniform_rho"] = s;
	}
	void unset_uniform_rho() {
		this->_events.erase("uniform_rho");
		this->_values.erase("uniform_rho");
	}
	void set_rho(Function fun) {
		this->unset_uniform_beta();
		this->_functions["rho"] = fun;
	}
	void set_uniform_mu(Vt s = 1) {     //default
		spEvent pse(new Flag("uniform_mu", 1));
		this->_events["uniform_mu"] = pse;
		this->_values["uniform_mu"] = s;
	}
	void unset_uniform_mu() {
		this->_events.erase("uniform_mu");
		this->_values.erase("uniform_mu");
	}
	void set_mu(Function fun) {
		this->unset_uniform_beta();
		this->_functions["mu"] = fun;
	}

	void set_projection_solver(const std::string& name, const int& max_iter =
			1000,      //
			const Vt& tol = 1e-3,      //
			const Any& any = nullptr    //
			) {
		ASSERT(name == "Jacobi" //
		|| name == "IC_CGS" //
		|| name == "SOR");
		this->_aflags["SetProjectionSolver"] = name;
		this->_projection->set_solver(name, any, max_iter, tol);
	}

	void set_CFL_number(Vt val = 0.5) {
		spEvent pse(new Flag("CFL_number", 1));
		this->_events["CFL_number"] = pse;
		this->_values["CFL_number"] = val;
	}

	void set_diffusion_number(Vt val = 1e10, Vt coe = 1.0) {
		spEvent pse(new Flag("Diffusion_number", 1));
		this->_events["Diffusion_number"] = pse;
		this->_values["Diffusion_number"] = val;

		this->_values["Diffusion_coefficient"] = coe;
	}

	void restrict_time() {
		Vt rdt = Operation::TimeRestrict_CourantNumber(
				this->_values["CFL_number"], *(this->_veoc));
		//if (rdt < this->_time->dt()) {
		//	this->_time->set_dt(rdt);
		//}
		// if the diffusion number larger than 100
		// no restriction (default)
		Vt rdt2 = this->_time->dt();
		if (this->_values["Diffusion_number"] < 100) {  //
			rdt2 = Operation::TimeRestrict_DiffusionNumber(
					this->_values["Diffusion_number"], *(this->_grid),
					this->_values["Diffusion_coefficient"]); //should change
			//if (rdt2 < this->_time->dt()) {
			//	this->_time->set_dt(rdt2);
			//}
		}
		if (this->_time->dt() > std::min(rdt, rdt2)) {
			this->_time->set_dt(std::min(rdt, rdt2));
		}

		// std::cout<< "dt = " << rdt << std::endl;

	}

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
protected:
	void _default_setup() {
		this->_nv = {"u", "v", "w"};
		this->_nvf = {"uf", "vf", "wf"};
		// set time step
		set_CFL_number();
		set_diffusion_number();
		// mu and rho
		set_uniform_rho();//rho = 1
		set_uniform_mu();//mu = 1
		this->_projection = spPoisson(new Poisson(this->_grid));
		// default velocity is 0
		_fun0 = [](Vt, Vt, Vt, Vt) {return 0;};
		this->fun_v_dot_nabla = &Operation::VdotNabla_upwind1;
	}

	spVectorCenter _new_veoc(const std::string& nu, const std::string& nv,
			const std::string& nw) {
		std::string veo_name[] = {nu, nv, nw};
		spScalar vc[] = {nullptr, nullptr, nullptr};
		FOR_EACH_DIM
		{
			// new center velocity and set
			vc[d] = spScalar(new Scalar(this->_grid));
			this->_css[veo_name[d]] = vc[d];
			if (this->has_function(veo_name[d])) {
				this->set_CS(veo_name[d], this->_functions[veo_name[d]]);
			} else {
				this->set_CS(veo_name[d], _fun0);
			}
		}
		return spVectorCenter(new VectorCenter(vc[0], vc[1], vc[2]));
	}

	spVectorFace _new_veof(const std::string& nu, const std::string& nv,
			const std::string& nw) {
		std::string veo_name[] = {nu, nv, nw};
		spScalar vc[] = {nullptr, nullptr, nullptr};
		FOR_EACH_DIM
		{
			// new center velocity and set
			vc[d] = spScalar(new Scalar(this->_grid));
			this->_css[veo_name[d]] = vc[d];
			if (this->has_function(veo_name[d])) {
				this->set_CS(veo_name[d], this->_functions[veo_name[d]]);
			} else {
				this->set_CS(veo_name[d], _fun0);
			}
		}
		return spVectorFace(new VectorFace(vc[0], vc[1], vc[2]));
	}

	int _ns_initial() {
		// std::cout << "  NS: initial \n";
		// check time
		ASSERT(this->has_event("time_term"));

		// if stand alone
		spGrid spg = this->_grid;
		this->_css["p"] = spScalar(new Scalar(spg));
		this->_css["p"]->assign(0.0);
		if (!this->has_event("uniform_rho")) {
			this->_css["rho"] = spScalar(new Scalar(spg));
			this->set_CS("rho", this->_functions["rho"]);
		}
		if (!this->has_event("uniform_mu")) {
			this->_css["mu"] = spScalar(new Scalar(spg));
			this->set_CS("mu", this->_functions["mu"]);
		}

		this->_veoc = _new_veoc(_nv[0], _nv[1], _nv[2]);
		this->_veof = _new_veof(_nvf[0], _nvf[1], _nvf[2]);

		if(this->_ghost == nullptr) {
			this->_ghost = spGhost(new GhostRegular_<Dim>(this->_grid, this->_bi));
		}
		return -1;
	}
};

}

#endif
