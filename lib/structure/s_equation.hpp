#ifndef _S_EQUATION_HPP
#define _S_EQUATION_HPP

#include <structure/s_data.hpp>
#include "s_define.hpp"
#include "s_boundary.hpp"
#include "s_event.hpp"
#include "s_grid.hpp"
#include "s_stencil.hpp"
#include "s_time.hpp"
#include "s_solver.hpp"
#include "s_operation.hpp"

namespace structure {

template<St DIM>
class Equation_ {
public:
	static const St Dim = DIM;
	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;

	typedef Scalar_<DIM> Scalar;
	typedef Scalar_<DIM>& ref_Scalar;
	typedef const Scalar_<DIM>& ref_const_Scalar;

	typedef std::shared_ptr<Scalar> spScalar;
	typedef std::shared_ptr<const Scalar> spcScalar;

	typedef BoundaryCondition BC;
	typedef std::shared_ptr<BoundaryCondition> spBC;
	typedef std::shared_ptr<const BoundaryCondition> spcBC;

	typedef BoundaryIndex BI;
	typedef std::shared_ptr<BI> spBI;
	typedef std::shared_ptr<const BI> spcBI;

	// function define ------------------------------
	typedef std::function<Vt(Vt, Vt, Vt, Vt)> Function;
	typedef std::unordered_map<std::string, Function> Functions;
	// Vaules ---------------------------------------
	typedef std::unordered_map<std::string, Vt> Values;
	// Stencil --------------------------------------
	typedef Index_<DIM> Index;
	typedef Stencil_<DIM> Stencil;
	typedef std::shared_ptr<Stencil> spStencil;
	typedef std::map<std::string, spStencil> Stencils;

	// Variables ------------------------------------
	typedef std::map<std::string, spScalar> CenterScalars;
	// Event ----------------------------------------
	typedef Event_<DIM> Event;
	typedef std::shared_ptr<Event> spEvent;
	typedef std::unordered_map<std::string, spEvent> Events;
	typedef EventFlag_<Dim> Flag;
	typedef Ghost_<DIM> Ghost;
	typedef std::shared_ptr<Ghost> spGhost;
	typedef typename Stencil::Expression Expression;
	typedef std::shared_ptr<Expression> spExpression;

	// time -----------------------------------------
	typedef Time_<DIM> Time;
	typedef std::shared_ptr<Time> spTime;

	// solver ---------------------------------------
	typedef Solver_<DIM> Solver;
	typedef std::shared_ptr<Solver> spSolver;
	typedef std::map<std::string, spSolver> Solvers;
	// Any
	typedef carpio::Any Any;
	typedef std::map<std::string, Any> AFlag;

protected:
	spGrid _grid;

	spGhost _ghost;
	spBI _bi;
	// time relates variables
	spTime _time;

	// Flag
	Events _events;       // _events
	Functions _functions; // independent of grid
	Values _values;       // values for equation
	Stencils _stencils;   // stencil
	Solvers _solvers;
	AFlag _aflags;        // other types of data put in this map

	// Variables
	CenterScalars _css;  //!< variables on the center of node
	//Variables _vars_ut; //!< untype variables on the node

	//Fun_pNode_spExp _get_node_spexp;
public:
	/*
	 * default constructor
	 *   2D
	 *
	 * utp
	 * utp map          0
	 */
	Equation_(spGrid pf) :
			_grid(pf) {
		_bi = spBI(new BI());
		this->_time = nullptr;
		//_timestep = spTimeStep(nullptr);
		//this->set_stand_alone();
		//_get_node_spexp = nullptr;
	}

	virtual int run_one_step(St step) {
		std::cout << step << "  Equation: run one step \n";
		return -1;
	}

	virtual int initial() {
		std::cout << "  Equation: initial \n";
		SHOULD_NOT_REACH;
		return -1;
	}

	virtual int finalize() {
		std::cout << "  Equation: finalize \n";
		return -1;
	}

	virtual int solve() {
		std::cout << "  Equation: solve \n";
		return -1;
	}

	void run() {
		// the equation don't have time
		if (!this->has_event("time_term")) {
			initial();
			// start events
			run_events(0, 0.0, Event::START);

			solve();

			run_events(1, 0.0, Event::END);
			finalize();
		} else {
			//Vt t = this->_time->_current_time();
			//St step = this->_time->_current_step();
			// events before calculation
			initial();
			run_events(this->_time->current_step(), //
					this->_time->current_step(),    //
					Event::START);
			//Vt ms = this->_time->max_step();
			//Vt dt = this->_time->dt();
			// loop
			while (!this->_time->is_end() && (!this->has_event("_STOP_"))) {
				//
				// events before each step
				run_events(this->_time->current_step(),  //
						this->_time->current_time(),     //
						Event::BEFORE);

				// run one step =================
				run_one_step(this->_time->current_step());
				// ==============================

				// events after each step
				run_events(this->_time->current_step(),  //
						this->_time->current_time(),     //
						Event::AFTER);
				//
				this->_time->advance();
			}
			// events after calculation
			run_events(this->_time->current_step(),    //
					this->_time->current_time(),       //
					Event::END);
			finalize();
		}
	}

	void run_events(St step, Vt t, int fob) {
		for (auto& event : this->_events) {
			if (event.second->_do_execute(step, t, fob)) {
				event.second->execute(step, t, fob, this);
			}
		}
	}

	virtual spGhost get_ghost() {
		return this->_ghost;
	}

	virtual spGrid get_grid() {
		return this->_grid;
	}

	/**
	 * @brief this function check the events
	 *        if flags contain key and value == val return true
	 *        the default val == 1
	 */
	bool has_event(const std::string& key, const int& val = 1) const {
		auto it = this->_events.find(key);
		if (it != this->_events.end()) {
			if ((*it).second->flag() == val) {
				return true;
			}
		}
		return false;
	}

	void add_event_flag(const std::string& key, const int& val = 1) {
		spEvent pse(new Flag(key, 1));
		this->_events[key] = pse;
	}

	bool has_function(const std::string& key) const {
		auto it = this->_functions.find(key);
		if (it != this->_functions.end()) {
			return true;
		}
		return false;
	}

	bool has_value(const std::string& key) const {
		auto it = this->_values.find(key);
		if (it != this->_values.end()) {
			return true;
		}
		return false;
	}

	bool has_flag(const std::string& key) const {
		auto it = this->_aflags.find(key);
		if (it != this->_aflags.end()) {
			return true;
		}
		return false;
	}

	void set_function(const std::string& name, Function fun) {
		this->_functions[name] = fun;
	}

	spScalar get_CS(const std::string& vname) {
		auto it = this->_css.find(vname);
		if (it != this->_css.end()) {
			return this->_css.at(vname);
		} else {
			std::cerr << ">! Center scalar not found :" << vname << "\n";
			SHOULD_NOT_REACH;
		}
		return nullptr;
	}

	spcScalar get_CS(const std::string& vname) const {
		auto it = this->_css.find(vname);
		if (it != this->_css.end()) {
			return this->_css.at(vname);
		} else {
			std::cerr << ">! Center scalar not found :" << vname << "\n";
			SHOULD_NOT_REACH;
		}
		return nullptr;
	}

	void set_CS(const std::string& vname, Function fun) {
		spScalar spcs = this->get_CS(vname);
		Vt ct = 0;
		if (this->_time != nullptr) {
			ct = this->_time->current_time();
		}
		for (typename Grid::Ijk IJK = _grid->begin_IJK(); !IJK.is_end();
				++IJK) {
			typename Grid::Ijk ijk = _grid->to_ijk(IJK);
			Vt val = fun(ct, _grid->c_(_X_, ijk.i()), _grid->c_(_Y_, ijk.j()),
					_grid->c_(_Z_, ijk.k()));
			spcs->val(ijk) = val;
		}
	}

	void set_CS(const std::string& vname, Vt val) {
		spScalar spcs = this->get_CS(vname);
		for (typename Grid::Ijk IJK = _grid->begin_IJK(); !IJK.is_end();
				++IJK) {
			spcs->VAL(IJK) = val;
		}
	}

	void initial_CS(const std::string& vname, Vt val) {
		Function fun = [val](Vt, Vt, Vt, Vt) {return val;};
		this->_functions[vname] = fun;
	}

	void initial_CS(const std::string& vname, Function fun) {
		this->_functions[vname] = fun;
	}

	void add_CS(const std::string& vname, Function fun) {
		auto it = this->_css.find(vname);
		if (it != this->_css.end()) {
			return;
		} else {
			spGrid spg = this->_grid;
			this->_css[vname] = spScalar(new Scalar(spg));
			this->set_CS(vname, fun);
		}
	}

	void add_CS(const std::string& vname, Vt v) {
		auto it = this->_css.find(vname);
		if (it != this->_css.end()) {
			return;
		} else {
			spGrid spg = this->_grid;
			this->_css[vname] = spScalar(new Scalar(spg));
			this->set_CS(vname, v);
		}
	}

	void set_time(const St& max_step, const Vt& dt) {
		this->_time = spTime(new Time(max_step, dt));
		spEvent pse(new Flag("time_term", 1));
		this->_events["time_term"] = pse;
	}

	/*
	 * Boundary condition
	 */

	spcBC find(const Index& index, const std::string& vname) const {
		St seg_idx = _ghost->which_boudary_seg(index);
		St shape_idx = _ghost->which_boudary_shape(index);
		return this->_bi->find(shape_idx, seg_idx, vname);
	}

	void add_bc(St shape_idx, St seg_idx, const std::string&name, spBC spbc) {
		this->_bi->insert(shape_idx, seg_idx, name, spbc);
	}

	void apply_bc(const std::string&name) {
		spScalar sps = this->_css[name];
		Scalar& scalar = (*sps);
		for (typename Grid::Ijk IJK = _grid->begin_IJK(); !IJK.is_end();
				++IJK) {
			typename Grid::Ijk ijk = _grid->to_ijk(IJK);
			Index index = ijk.current();
			if (this->_ghost->is_ghost(index)) {
				spcBC spbc = find(index, name);
				Orientation ori = this->_ghost->orientation(index);
				Index oriidx = this->_ghost->ori_index(index);
				St a = this->_ghost->ori_axes(index);
				Vt val = spbc->get_val(               //
						this->_time->current_time(), //time ==========================
						this->_grid->f_(_X_, ori, oriidx[_X_]), //
						this->_grid->f_(_Y_, ori, oriidx[_Y_]), //
						this->_grid->f_(_Z_, ori, oriidx[_Z_]));
				Vt dl = std::abs(
						this->_grid->c_(a, index[a])
								- this->_grid->c_(a, oriidx[a]));
				if (spbc->get_type() == BC::_BC1_) {
					// Boundary condition 1
					Vt inter = carpio::Interp_<Vt, Vt>::Linear(
							this->_grid->c_(a, index[a]), // coordinate ghost
							this->_grid->f_(a, this->_ghost->orientation(index),
									oriidx[a]),  // coordinate face
							this->_grid->c_(a, oriidx[a]), // coordinate original center
							val,                            // face value
							scalar(oriidx));               // original value
					scalar(index) = inter;                  // ghost value
				} else {
					// Boundary condition 2
					Vt po = scalar(oriidx);
					Vt pg = po - val * dl;
					(*sps)(index) = pg;
				}
			}
		}
	}

	void set_output_time(int is = -1, int ie = -1, int istep = -1,
			int flag = -1, const std::string& fn = "") {
		// output to screen
		spEvent pse(new EventOutputTime_<Dim>(is, ie, istep, flag, fn));
		this->_events["OutputTime"] = pse;
	}

	void set_center_scalar(const std::string& vname, Function fun, int is = -1,
			int ie = -1, int istep = -1, int flag = -1) {
		// output to screen
		spEvent pse(
				new EventSetCenterScalar_<Dim>(vname, fun, is, ie, istep,
						flag));
		this->_events["SetCenterScalar"] = pse;
	}

// the vname and fn is the key
//
	void set_output_cs2file(const std::string& vname, const std::string& fn,
			int is = -1, int ie = -1, int istep = -1, int flag = -1) {
		spEvent pse(
				new EventOutputCenterScalar_<Dim>(vname, fn, is, ie, istep,
						flag));
		std::stringstream sst;
		sst << "OutputCenterScalar_" << vname << "_" << fn;
		this->_events[sst.str()] = pse;
	}

	void set_output_error(const std::string& vexa, const std::string& vval,
			int is = -1, int ie = -1, int istep = -1, int flag = -1,
			const std::string& fn = "") {
		spEvent pse(
				new EventOutputError_<Dim>(vexa, vval, is, ie, istep, flag,
						fn));
		std::stringstream sst;
		sst << "OutputError_" << vexa << "_" << vval;
		this->_events[sst.str()] = pse;
	}

	void set_output_solver_residual(const std::string& fn) {
		ASSERT(!this->has_flag("OutputSolverResidual"));
		FILE* fp = std::fopen(fn.c_str(), "w");
		if (!fp) {
			std::perror("File opening failed: file for output solver residual");
			return;
		}
		this->_aflags["OutputSolverResidual"] = fp;
	}

	void set_stop_cs(const std::string& vname, Vt e1 = -1, Vt e2 = -1, Vt e3 =
			-1, int is = -1, int ie = -1, int istep = -1, int flag = -1) {
		spEvent pse(
				new EventStopScalarError_<Dim>(is, ie, istep, flag, vname, e1,
						e2, e3));
		std::stringstream sst;
		sst << "StopError_" << vname;
		this->_events[sst.str()] = pse;
	}

	void show_events() const {
		std::cout << "Events:  (" << this->_events.size() << ")\n";
		for (auto iter = this->_events.begin(); iter != this->_events.end();
				++iter) {
			std::cout << "  " << iter->first << std::endl;
		}
	}

	void show_functions() const {
		std::cout << "Functions:  (" << this->_functions.size() << ")\n";
		for (auto iter = this->_functions.begin();
				iter != this->_functions.end(); ++iter) {
			std::cout << "  " << iter->first << std::endl;
		}
	}

	void show_values() const {
		std::cout << "Values:  (" << this->_values.size() << ")\n";
		for (auto iter = this->_values.begin(); iter != this->_values.end();
				++iter) {
			std::cout << "  " << std::setw(15) << iter->first << " : "
					<< iter->second << std::endl;
		}
	}

	virtual ~Equation_() {
		if (this->has_flag("OutputSolverResidual")) {
			FILE* pf = carpio::any_cast<std::FILE*>(
					this->_aflags["OutputSolverResidual"]);
			std::fclose(pf);
		}
	}

}
;

}

#endif
