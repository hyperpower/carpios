#ifndef _S_POISSON_HPP
#define _S_POISSON_HPP

#include <structure/s_data.hpp>
#include "s_define.hpp"
#include "s_stencil.hpp"
#include "s_equation.hpp"
#include "s_solver.hpp"

namespace structure {
//This file use to solve poisson equation
//
//   ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
//
//   alpha(x,y)•phi(x,y) + ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   alpha(x,y)•phi(x,y) + ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
//
//   tau (d phi(x,y) / dt) + alpha(x,y)•phi(x,y) + ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   tau (d phi(x,y) / dt) + alpha(x,y)•phi(x,y) + ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
template<St DIM>
class StencilPoisson_: public Stencil_<DIM> {
public:
	static const St Dim = DIM;
	typedef Index_<DIM> Index;
	typedef Stencil_<DIM> Stencil;

	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef std::function<Vt(const Index&, const short&, const short&)> FunGrid;
	typedef Scalar_<DIM> CenterScalar;
	typedef std::shared_ptr<CenterScalar> spCenterScalar;
	typedef typename Stencil::Expression Expression;
	typedef typename Stencil::spExpression spExpression;

protected:
	FunGrid _d; // distance to other center
	FunGrid _a; // face area
	FunGrid _s; // face area
	FunGrid _betaf;  //beta on face
	FunGrid _source; // source term

	spCenterScalar _cs_beta;   // scalar field of beta
	spCenterScalar _cs_source;   // scalar field of beta
	Vt _beta, _src;
public:
	StencilPoisson_(spGrid g, spCenterScalar csb, spCenterScalar source, Vt b =
			1, Vt s = 0) :
			Stencil(g) {
		_cs_beta = csb;
		_cs_source = source;
		_beta = b;
		_src = s;
		_build_default_a();
		_build_default_s();
		_build_default_beta();
		_build_default_d();
		_build_default_source();
	}
	virtual Vt coe(const Index& index, short dim, short ori) const {
		ASSERT(dim < DIM);
		if (ori == _C_) {
			Vt coec = 0;
			for (St d = 0; d < Dim; d++) {
				Vt coep = -_coe_other(index, d, _P_);
				Vt coem = -_coe_other(index, d, _M_);
				coec += coem;
				coec += coep;
			}
			return coec;
		} else {
			return _coe_other(index, dim, ori);
		}
	}
	virtual Vt src(const Index& index, short dim, short ori) const {
		ASSERT(dim < DIM);
		ASSERT(ori == _C_);
		Vt v = 1.0;
		for (St d = 0; d < Dim; d++) {
			v *= _s(index, dim, ori);
		}
		return v * _source(index, dim, ori);
	}

	virtual spExpression coe_row(const Index& index) const {
		spExpression res(new Expression());
		// all the variables are on left hand side
		for (St d = 0; d < Dim; d++) {
			Vt coep = coe(index, d, _P_);
			Vt coem = coe(index, d, _M_);
			Index idxp = index.p(d);
			Index idxm = index.m(d);
			res->insert(coep, idxp, 1);
			res->insert(coem, idxm, 1);
		}
		Vt coec = coe(index, 0, _C_);
		res->insert(coec, index, 1);
		// constant
		Vt s = src(index, 0, _C_);
		res->insert(s, index, 0);
		return res;
	}

	virtual ~StencilPoisson_() {

	}

protected:
	Vt _coe_other(const Index& index, short dim, short ori) const {
		Vt betaf = _betaf(index, dim, ori);
		Vt a = _a(index, dim, ori);
		Vt d = _d(index, dim, ori);
		return -betaf * a / d;
	}
	void _build_default_d() {
		FunGrid fun =
				[this](const Index& index, const short& dim, const short& ori) {
					Vt c = this->_grid->c_(dim, index.value(dim));
					int atoadd[] = {-1, 1, 0};
					Vt co = this->_grid->c_(dim, index.value(dim) + atoadd[ori]);
					return std::abs(c - co);
				};
		this->_d = fun;
	}
	void _build_default_a() {
		FunGrid fun =
				[this](const Index& index, const short& dim, const short& ori) {
					St d1[] = {1,2,0};
					St d2[] = {2,0,1};
					Vt s1 = this->_grid->s_(d1[dim], index.value(d1[dim]));
					Vt s2 = this->_grid->s_(d2[dim], index.value(d2[dim]));
					return s1 * s2;
				};
		this->_a = fun;
	}
	void _build_default_s() {
		FunGrid fun =
				[this](const Index& index, const short& dim, const short& ori) {
					Vt s = this->_grid->s_(dim, index.value(dim));
					return s;
				};
		this->_s = fun;
	}
	void _build_default_beta() {
		// get average value on face
		FunGrid fun;
		if (_cs_beta == nullptr) {
			fun =
					[this](const Index& index, const short& dim, const short& ori) {
						return this->_beta;
					};
		} else {
			fun =
					[this](const Index& index, const short& dim, const short& ori) {
						Vt v = this->_cs_beta->val(index);
						Vt df = this->_grid->df_(dim, index.value(dim));
						int atoadd[] = {-1, 1, 0};
						Vt vo = this->_cs_beta->val(index.p(dim));
						Vt dfo = this->_grid->df_(dim, index.value(dim) + atoadd[ori]);
						return (df* v + dfo * vo)/ (df + dfo);
					};
		}
		this->_betaf = fun;
	}
	void _build_default_source() {
		// get average value on face
		FunGrid fun;
		if (_cs_source == nullptr) {
			fun =
					[this](const Index& index, const short& dim, const short& ori) {
						return this->_src;
					};
		} else {
			fun =
					[this](const Index& index, const short& dim, const short& ori) {
						Vt v = this->_cs_source->val(index);
						return v;
					};
		}
		this->_source = fun;
	}
}
;

template<St DIM>
class StencilHelmholtz_: public StencilPoisson_<DIM> {
public:
	static const St Dim = DIM;
	typedef Index_<DIM> Index;
	typedef Stencil_<DIM> Stencil;
	typedef StencilPoisson_<DIM> StencilPoisson;
	typedef StencilPoisson Base;

	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef std::function<Vt(const Index&, const short&, const short&)> FunGrid;
	typedef Scalar_<DIM> CenterScalar;
	typedef std::shared_ptr<CenterScalar> spCenterScalar;
	typedef typename Stencil::Expression Expression;
	typedef typename Stencil::spExpression spExpression;

protected:
	spCenterScalar _cs_alpha;   // scalar field of beta
	Vt _alpha;
public:
	StencilHelmholtz_(spGrid g, spCenterScalar csalpha, spCenterScalar csb,
			spCenterScalar source, Vt a = 1, Vt b = 1, Vt s = 0) :
			Base(g, csb, source, b, s) {
		_cs_alpha = csalpha;
		_alpha = a;
	}
	Vt coe(const Index& index, short dim, short ori) const {
		ASSERT(dim < DIM);
		if (ori == _C_) {
			Vt coec = 0;
			for (St d = 0; d < Dim; d++) {
				Vt coep = -this->_coe_other(index, d, _P_);
				Vt coem = -this->_coe_other(index, d, _M_);
				coec += coem;
				coec += coep;
			}
			return coec;
		} else {
			return this->_coe_other(index, dim, ori);
		}
	}
	Vt src(const Index& index, short dim, short ori) const {
		ASSERT(dim < DIM);
		ASSERT(ori == _C_);
		Vt v = 1.0;
		for (St d = 0; d < Dim; d++) {
			v *= this->_s(index, d, ori);
		}
		return v * this->_source(index, dim, ori);
	}

	spExpression coe_row(const Index& index) const {
		spExpression res(new Expression());
		// all the variables are on left hand side
		for (St d = 0; d < Dim; d++) {
			Vt coep = coe(index, d, _P_);
			Vt coem = coe(index, d, _M_);
			Index idxp = index.p(d);
			Index idxm = index.m(d);
			res->insert(coep, idxp, 1);
			res->insert(coem, idxm, 1);
		}
		Vt coec = coe(index, 0, _C_);
		Vt v = 1.0;
		for (St d = 0; d < Dim; d++) {
			v *= this->_s(index, d, _P_);
		}
		res->insert(coec + v * _alpha, index, 1);   //
		// constant
		Vt s =  v * this->_source(index, 0, _C_);
		res->insert(s, index, 0);
		return res;
	}

};

/*
 * the Poisson class
 */
template<St DIM>
class Poisson_: public Equation_<DIM> {
public:
	static const St Dim = DIM;
	typedef Equation_<DIM> Base;
	typedef Equation_<DIM> Equation;
	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Index_<Dim> Index;

	typedef Scalar_<DIM> CenterScalar;
	typedef std::shared_ptr<CenterScalar> spCenterScalar;

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
	typedef std::map<std::string, spCenterScalar> CenterScalars;
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

public:
	Poisson_(spGrid spg, spGhost sph = nullptr) :
			Base(spg, sph) {
		_default_setup();
	}

	int initial() {
		//std::cout << "  Poisson: initial \n";

		if (this->has_event("stand_alone")) {   // if stand alone
			spGrid spg = this->_grid;
			this->_css["phi"] = spCenterScalar(new CenterScalar(spg));
			this->_css["phi"]->assign(0.0);
			if (!this->has_event("uniform_beta")) {
				this->_css["beta"] = spCenterScalar(new CenterScalar(spg));
				this->set_CS("beta", this->_functions["beta"]);
			}
			if (!this->has_event("uniform_source")) {
				this->_css["source"] = spCenterScalar(new CenterScalar(spg));
				this->set_CS("source", this->_functions["source"]);
			}
			if (this->has_event("time_term")) {
				this->_css["phis"] = spCenterScalar(new CenterScalar(spg));
				this->_css["phis"]->assign(0.0);
			}
		} else {

		}

		if (this->has_event("alpha_term")) {
			this->_stencils["poisson"] = _new_helmholtz_stencil();
		} else {
			this->_stencils["poisson"] = _new_poisson_stencil();
		}

		if ( this->_ghost == nullptr){
			this->_ghost = spGhost(new GhostRegular_<Dim>(this->_grid, this->_bi));
		}
		_init_solver();

		return -1;
	}

	void _init_solver() {
		spSolver spsolver;
		if (this->has_flag("SetSolver")) {
			std::string sn = carpio::any_cast<std::string>(
					this->_aflags["SetSolver"]);
			if (sn == "Jacobi") {
				spsolver = spSolver(
						new Solver_Jacobi_<Dim>(this->_grid, this->_ghost,
								this->_stencils["poisson"], "phi",
								this->_css["phi"]));
			} else if (sn == "IC_CGS") {
				spsolver = spSolver(
						new Solver_IC_CGS_<Dim>(this->_grid, this->_ghost,
								this->_stencils["poisson"], "phi",
								this->_css["phi"]));
			} else if (sn == "SOR") {
				ASSERT(this->has_flag("SOR_omega"));
				Vt sn = carpio::any_cast<Vt>(this->_aflags["SOR_omega"]);
				spsolver = spSolver(
						new Solver_SOR_<Dim>(this->_grid, this->_ghost,
								this->_stencils["poisson"], "phi",
								this->_css["phi"], sn));
			}
		} else {
			// default solver
			spsolver = spSolver(
					new Solver_Jacobi_<Dim>(this->_grid, this->_ghost,
							this->_stencils["poisson"], "phi",
							this->_css["phi"]));
		}
		// max iter and tolerance -----------------------------------
		if (this->has_flag("SetSolver_max_iter")) {
			int val = carpio::any_cast<int>(
					this->_aflags["SetSolver_max_iter"]);
			spsolver->set_max_iter(val);
		}
		if (this->has_flag("SetSolver_tolerance")) {
			Vt val = carpio::any_cast<Vt>(this->_aflags["SetSolver_tolerance"]);
			spsolver->set_tolerance(val);
		}

		if (this->has_flag("OutputSolverResidual")) {
			FILE* pf = carpio::any_cast<std::FILE*>(
					this->_aflags["OutputSolverResidual"]);
			spsolver->set_output_residual(pf);
		}
		this->_solvers["main"] = spsolver;
	}

	int finalize() {
		//std::cout << "  Poisson: finalize \n";
		return -1;
	}

	int solve() {
		// build coe matrix
		// std::cout << "  Poisson: solve \n";
		spSolver solver = this->_solvers["main"];
		//solver->set_output_residual();
		//solver->set_max_iter(1000);
		solver->solve();
		// build matrix ----
		return -1;
	}

	int run_one_step(St step) {
		this->apply_bc("phi");
		_run_one_step_explicit(step);
		return 1;
	}

	int _run_one_step_explicit(St step) {
		CenterScalar& phi = *(this->get_CS("phi"));
		CenterScalar& phis = *(this->get_CS("phis"));
		if (this->has_event("uniform_beta")) {
			Vt beta = this->_values["uniform_beta"];
			Operation::NablaMuNabla(beta, phi, phis);
		} else {
			CenterScalar& beta = *(this->get_CS("beta"));
			Operation::NablaMuNabla(beta, phi, phis);
		}

		if (this->has_event("uniform_source")) {
			Vt source = this->_values["uniform_source"];
			Operation::Add(phis, source);
		} else {
			CenterScalar& source = *(this->get_CS("source"));
			Operation::Add(phis, source);
		}
		Operation::Multiply(phis, (this->_time->dt()));
		Operation::Add(phi, phis);
		return 1;
	}

// set function ===
	/**
	 *  @brief set uniform beta
	 */
	void set_uniform_beta(Vt beta = 1) {     //default
		spEvent pse(new Flag("uniform_beta", 1));
		this->_events["uniform_beta"] = pse;
		this->_values["uniform_beta"] = beta;
	}

	void set_uniform_alpha(Vt a = 1) {
		spEvent psa(new Flag("alpha_term", 1));
		this->_events["alpha_term"] = psa;
		spEvent pse(new Flag("uniform_alpha", 1));
		this->_events["uniform_alpha"] = pse;
		this->_values["uniform_alpha"] = a;
	}
protected:
	void unset_uniform_beta() {
		this->_events.erase("uniform_beta");
		this->_values.erase("uniform_beta");
	}
	void unset_uniform_alpha() {
		this->_events.erase("uniform_alpha");
		this->_values.erase("uniform_alpha");
	}
public:
	void set_beta(Function fun) {
		this->unset_uniform_beta();
		this->_functions["beta"] = fun;
	}
	void set_alpha(Function fun) {
		spEvent pse(new Flag("alpha_term", 1));
		this->_events["alpha_term"] = pse;
		this->unset_uniform_beta();
		this->_functions["alpha"] = fun;
	}

	/**
	 *  @brief set uniform beta
	 */
	void set_uniform_source(Vt s = 0) {     //default
		spEvent pse(new Flag("uniform_source", 1));
		this->_events["uniform_source"] = pse;
		this->_values["uniform_source"] = s;
	}
	void unset_uniform_source() {
		this->_events.erase("uniform_source");
		this->_values.erase("uniform_source");
	}
	void set_source(Function fun) {
		this->unset_uniform_source();
		this->_functions["source"] = fun;
	}

	void set_time(const St& max_step, const Vt& dt, const Vt& tau = 1) {
		this->_time = spTime(new Time(max_step, dt));
		spEvent pse(new Flag("time_term", 1));
		this->_events["time_term"] = pse;
		this->_values["tau"] = tau;
	}

	void set_depend( //
			spCenterScalar spphi, //
			spCenterScalar spbeta,  //
			spCenterScalar spsource, //
			Vt beta = 1, Vt source = 0) {
		this->_events.erase("stand_alone");
		ASSERT(spphi != nullptr);
		this->_css["phi"] = spphi; // set phi
		if (spbeta != nullptr) {  //  set beta
			unset_uniform_beta();
			this->_css["beta"] = spbeta;
		} else {
			set_uniform_beta(beta);
		}
		if (spsource != nullptr) {  // set source
			unset_uniform_source();
			this->_css["source"] = spsource;
		} else {
			set_uniform_source(source);
		}
	}

	void set_solver( //
			const std::string& name,                 //
			const Any& any = nullptr,   //
			const int& max_iter = 1000,      //
			const Vt& tol = 1e-3       //
			) {
		// name should
		ASSERT(name == "Jacobi" //
		|| name == "IC_CGS" //
		|| name == "SOR");
		this->_aflags["SetSolver"] = name;
		this->_aflags["SetSolver_max_iter"] = max_iter;
		this->_aflags["SetSolver_tolerance"] = tol;
		if (name == "SOR") {
			this->_aflags["SOR_omega"] = any;
		}
	}

protected:
	void _default_setup() {
		spEvent pse(new Flag("stand_alone", 1));
		this->_events["stand_alone"] = pse;
		set_uniform_beta();  //bata = 1
		set_uniform_source(); //source = 0

	}

	spStencil _new_poisson_stencil() {
		return spStencil(
				new StencilPoisson_<DIM>(
						this->_grid, //
						this->has_event("uniform_beta") ?
								nullptr : this->_css["beta"], //
						this->has_event("uniform_source") ?
								nullptr : this->_css["source"],
						this->has_event("uniform_beta") ?
								this->_values["uniform_beta"] : 1,
						this->has_event("uniform_source") ?
								this->_values["uniform_source"] : 0));
	}

	spStencil _new_helmholtz_stencil() {
		return spStencil(
				new StencilHelmholtz_<DIM>(
						this->_grid, //
						this->has_event("uniform_alpha") ?
								nullptr : this->_css["alpha"], //
						this->has_event("uniform_beta") ?
								nullptr : this->_css["beta"], //
						this->has_event("uniform_source") ?
								nullptr : this->_css["source"], //
						this->has_event("uniform_alpha") ?
								this->_values["uniform_alpha"] : 1,  //
						this->has_event("uniform_beta") ?
								this->_values["uniform_beta"] : 1, //
						this->has_event("uniform_source") ?
								this->_values["uniform_source"] : 0));
	}

};

}

#endif
