#ifndef _EQUATION_H_
#define _EQUATION_H_

#include "calculation_define.hpp"
#include "expression.hpp"
#include "event.hpp"
#include "utility/clock.h"
#include "utility/format.h"
#include "io/mmio.h"
#include "timestep.hpp"

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

#include "utility/any.hpp"

//#define __Debug__

#define __x__  -0.4
#define __y__  0.45

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM>
class Equation_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Equation_<cvt, vt, Dim> Self;
	typedef Equation_<cvt, vt, Dim>& ref_Self;
	typedef Equation_<cvt, vt, Dim>* pSelf;
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

	typedef Stencil_<cvt, vt, Dim, Dim> Stencil;

	typedef Interpolate_<cvt, vt, Dim> Interpolate;

	// function define ------------------------------
	typedef std::function<vt(cvt, cvt, cvt)> Function;
	typedef std::unordered_map<std::string, Function> Functions;

	typedef std::unordered_map<std::string, vt> Values;

	typedef std::function<void(pNode)> Fun_pNode;

	// Events ----------------------------------------
	typedef Event_<cvt, vt, Dim> Event;
	typedef std::shared_ptr<Event> spEvent;
	typedef std::unordered_map<std::string, spEvent> Events;

	typedef EventFlag_<cvt, vt, Dim> Flag;

	// Variables ------------------------------------
	typedef std::map<std::string, St> Variables;

	typedef Expression_<cvt, vt, Dim> Exp;
	typedef Exp* pExp;
	typedef std::shared_ptr<Exp> spExp;
	typedef std::function<spExp(pNode)> Fun_pNode_spExp;

	typedef std::shared_ptr<Face> spFace;

	typedef typename Exp::Term Term;
	typedef typename Exp::iterator ExpIter;
	typedef typename Exp::const_iterator const_ExpIter;

	// Time step ------------------------------------
	typedef TimeStep_<cvt, vt, Dim> TimeStep;
	typedef std::shared_ptr<TimeStep> spTimeStep;

protected:
	typedef MatrixSCR_<vt> Mat;
	typedef ArrayListV<vt> Arr;

	typedef ArrayListT<spExp> Arr_spExp;
	typedef ArrayListT<spExp>* pArr_spExp;

	typedef ArrayListT<Any> Arr_Any;

	struct comp_spFace {
		bool operator()(const spFace& lhs, const spFace& rhs) {
			// just compare pface
			if (lhs->pori() < rhs->pori()) {
				return true;
			} else if (lhs->pori() == rhs->pori()) {
				if (lhs->pnei() < rhs->pnei()) {
					return true;
				} else if (lhs->pnei() == rhs->pnei()) {
					if (lhs->dir() < rhs->dir()) {
						return true;
					} else {
						return false;
					}
				} else {
					return false;
				}
			} else {
				return false;
			}
		}
	};

	typedef std::map<spFace, Arr_Any, comp_spFace> Map_FA;
	typedef std::map<spFace, Arr_Any, comp_spFace>* pMap_FA;
	typedef typename std::map<spFace, Arr_Any, comp_spFace>::iterator Iter_FA;
	typedef typename std::map<spFace, Arr_Any, comp_spFace>::const_iterator const_Iter_FA;
	typedef std::pair<spFace, Arr_Any> Pair_FA;
	typedef std::pair<Iter_FA, bool> Ret_FA;

	void _new_map_fa(pNode pn, St idx) {
		utPointer& utp = pn->utp(idx);
		if (utp == nullptr) {
			utp = new Map_FA();
		} else {
			Map_FA& mfe = CAST_REF(pMap_FA, utp);
			mfe.clear();
		}
	}

	// new node expression
	void _new_ne(pNode pn, St len, St idx) {
		utPointer& utp = pn->utp(idx);
		if (utp == nullptr) {
			utp = new Arr_spExp(len);
			Arr_spExp& arrne = CAST_REF(pArr_spExp, utp);
			for (St i = 0; i < len; i++) {
				spExp sp(new Exp());
				arrne[i] = sp;
			}
		} else {
			Arr_spExp& arrne = CAST_REF(pArr_spExp, utp);
			for (St i = 0; i < len; i++) {
				arrne[i]->clear();
			}
		}
	}

	void _delete_map_fa(pNode pn, St idx) {
		utPointer& utp = pn->utp(idx);
		if (utp != nullptr) {
			pMap_FA mfe = CAST(pMap_FA, utp);
			mfe->clear();
			delete mfe;
			utp = nullptr;
		} else {
			SHOULD_NOT_REACH;
		}
	}

	void _delete_ne(pNode pn, St idx) {
		utPointer& utp = pn->utp(idx);
		if (utp != nullptr) {
			pArr_spExp arrne = CAST(pArr_spExp, utp);
			for (St i = 0; i < arrne->size(); i++) {
				(*arrne)[i]->clear();
			}
			delete arrne;
			utp = nullptr;
		} else {
			SHOULD_NOT_REACH;
		}
	}

	void _face_val_gra(Face& pf, Exp& exp_val, Exp& exp_gra) {
		// face direction
		pNode pori = pf.pori();
		Direction dir = pf.dir();
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(dir, o, a);

		// interpolate on face
		exp_val = Interpolate::OnFace(pori, dir, 1);
		exp_gra = Interpolate::OnFace_Gradient(pori, dir, 2);
	}

	void _new_arr_any(pNode pn, St ut_idx, St len = 2) {
		utPointer utp = pn->utp(ut_idx);
		ASSERT(utp != nullptr);
		Map_FA& mfe = CAST_REF(pMap_FA, utp);
		for (St i = 0; i < NumFaces; i++) {
			Direction dir = FaceDirectionInOrder(i);
			pNode pnei = pn->get_neighbor_fast(dir);
			// The face need to be calculate
			// 1 the boundary face
			// 2 the equal face on P direction;
			// 3 the fine to coarse face
			FaceType ft = GetFaceType(pn, pnei);
			Arr_Any arr_any(len);
			for (St ii = 0; ii < arr_any.size(); ++ii) {
				spExp sp(new Exp());
				arr_any[ii] = sp;
			}
			spFace spface(new Face(pn, pnei, dir, ft));
			if ((ft == _Boundary_)				//1
			|| (ft == _Equal_ && (IsFacePDirection(dir)))				//2
					|| (ft == _FineCoarse_))				//3
					{
				// work on spface
				mfe.insert(Pair_FA(spface, arr_any));
			}
			// case 2 3
			if ((ft == _Equal_ && (IsFacePDirection(dir)))				//2
			|| (ft == _FineCoarse_))				//3
					{
				// work on pnei
				utPointer utp_nei = pnei->utp(ut_idx);
				ASSERT(utp_nei != nullptr);
				Map_FA& mfe_nei = CAST_REF(pMap_FA, utp_nei);
				FaceType oft = (ft == _FineCoarse_) ? _CoarseFine_ : _Equal_;
				spFace spface_nei(new Face(pnei, pn, Opposite(dir), oft));
				mfe_nei.insert(Pair_FA(spface_nei, arr_any));
			}
		}
	}

	void _update_face_val_gra(St idx_val = 0, St idx_gra = 1) {
		Grid& grid = _pdomain->grid();
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			utPointer utp = pn->utp(this->_vars_ut["FA"]);
			ASSERT(utp != nullptr);
			Map_FA& mfe = CAST_REF(pMap_FA, utp);
			for (typename Map_FA::iterator it = mfe.begin(); it != mfe.end();
					it++) {
				FaceType ft = it->first->ft();
				Direction dir = it->first->dir();
				Arr_Any& arr_any = it->second;
				if ((ft == _Boundary_)				//1
				|| (ft == _Equal_ && (IsFacePDirection(dir)))				//2
						|| (ft == _FineCoarse_))				//3
						{
					spExp exp_val = any_cast<spExp>(arr_any[idx_val]);
					spExp exp_gra = any_cast<spExp>(arr_any[idx_gra]);
					_face_val_gra((*(it->first)), *exp_val, *exp_gra);
				}
			}
		}
	}

	void _new_FA(St len = 2) {
		Grid& grid = _pdomain->grid();
		// new map FE on leaf
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_new_map_fa(pn, this->_vars_ut["FA"]);
		}
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_new_arr_any(pn, this->_vars_ut["FA"], len);
		}
	}

	void _build_NE_on_leaf(St len) {
		Grid& grid = _pdomain->grid();
		// new map FE on leaf
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_new_ne(pn, len, this->_vars_ut["NE"]);
		}
	}

	void _delete_FA_on_leaf() {
		Grid& grid = _pdomain->grid();
		// delete map FE on leaf
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_delete_map_fa(pn, this->_vars_ut["FA"]);
		}
	}

	void _delete_NE_on_leaf() {
		Grid& grid = _pdomain->grid();
		// delete map FE on leaf
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_delete_ne(pn, this->_vars_ut["NE"]);
		}
	}

protected:
	//matrx
	typedef ArrayListV<St> Arr_st;

	pDomain _pdomain;

	// time relates variables
	// vt _dt;
	// vt _max_step;
	spTimeStep _timestep;

	// Flag
	Events _events;       // _events
	Functions _functions; // independent of grid
	Values _values;       // values for equation

	// Variables
	Variables _vars_c;  //!< variables on the center of node
	Variables _vars_ut; //!< untype variables on the node

	Fun_pNode_spExp _get_node_spexp;

#ifdef __Debug__
	std::fstream debug_log;
#endif
public:
	/*
	 * default constructor
	 *   2D
	 *
	 * utp
	 * utp map          0
	 */
	Equation_(pDomain pf) :
			_pdomain(pf) {
		_timestep = spTimeStep(nullptr);
		this->set_stand_alone();
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
		if (!this->has_time_term()) {
			initial();
			// start events
			run_events(0, 0.0, Event::START);

			solve();

			run_events(1, 0.0, Event::END);
			finalize();
		} else {
			vt t = 0.0;
			St step = 0;
			// events before calculation
			initial();
			run_events(0, 0.0, Event::START);
			vt ms = this->_timestep->max_step();
			vt dt = this->_timestep->dt();
			// loop
			for (; step < ms; ++step) {
				//
				// events before each step
				run_events(step, t, Event::BEFORE);

				// run one step =================
				run_one_step(step);
				// ==============================

				// events after each step
				run_events(step, t, Event::AFTER);
				//
				t += dt;
			}
			// events after calculation
			run_events(ms, t, Event::END);
			finalize();
		}
	}

	virtual ~Equation_() {

	}

	void run_events(St step, vt t, int fob) {
		for (auto& event : this->_events) {
			if (event.second->_do_execute(step, t, fob)) {
				event.second->execute(step, t, fob, _pdomain);
			}
		}
	}

public:
	/**
	 * @brief set value on center of node
	 */
	void set_val_grid(Function pfun, St idx) {
		_pdomain->set_val_grid(idx, pfun);
	}

	void set_val_grid(Function pfun, const std::string& name) {
		St idx = this->get_var_center_idx(name);
		set_val_grid(pfun, idx);
	}

	void set_val_ghost(Function pfun, St idx) {
		_pdomain->set_val_ghost(idx, pfun);
	}

	void set_val_ghost(Function pfun, const std::string& name) {
		St idx = this->get_var_center_idx(name);
		set_val_ghost(pfun, idx);
	}

	void set_val(Function pfun, St idx) {
		this->set_val_grid(pfun, idx);
		this->set_val_ghost(pfun, idx);
	}

	void set_val(Function pfun, const std::string& name) {
		St idx = this->get_var_center_idx(name);
		set_val(pfun, idx);
	}

	/**
	 * @brief this function returns the max index of the valuables in
	 *        ns, the max index makes sure the array size is ok.
	 */
	St max_idx(const Variables& vars) {
		St max = 0;
		for (auto& item : vars) {
			if (item.second > max) {
				max = item.second;
			}
		}
		return max;
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

	void set_output_time(int is = -1, int ie = -1, int istep = -1,
			int flag = -1, std::FILE * f = nullptr) {
		// output to screen
		spEvent pse(new EventOutputTime_<cvt, vt, Dim>(is, ie, istep, flag, f));
		this->_events["OutputTime"] = pse;
	}

	bool is_stand_alone() const {
		return this->has_event("stand_alone");
	}
	/**
	 *  @brief add event
	 */
	void add_event(const std::string& name, spEvent pse) {
		if (!this->has_event(name)) {
			this->_events[name] = pse;
		} else {
			SHOULD_NOT_REACH;
		}
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

	void show_variables_c() const {
		std::cout << "Vars c:  (" << this->_vars_c.size() << ")\n";
		for (auto iter = this->_vars_c.begin(); iter != this->_vars_c.end();
				++iter) {
			std::cout << "  " << std::setw(5) << iter->first << " : "
					<< iter->second << std::endl;
		}
	}

	void show_variables_ut() const {
		std::cout << "Vars ut:  (" << this->_vars_ut.size() << ")\n";
		for (auto iter = this->_vars_ut.begin(); iter != this->_vars_ut.end();
				++iter) {
			std::cout << "  " << std::setw(5) << iter->first << " : "
					<< iter->second << std::endl;
		}
	}
	/**
	 *  @brief set time term
	 *         default --> stand_alone
	 */
	void set_stand_alone() {
		spEvent pse(new Flag("stand_alone", 1));
		this->_events["stand_alone"] = pse;
	}

	/*
	 * Does the equation has term of time
	 */
	bool has_time_term() const {
		return this->has_event("time_term");
	}

	/**
	 *  @brief set time term
	 *         default --> no time term
	 */
	void set_time_term(vt dt, St max_step, const std::string& scheme =
			"explicit") {
		spEvent pse(new EventFlag_<cvt, vt, Dim>("time_term", 1));
		this->_events["time_term"] = pse;
		//
		//spEvent pdt(new FlagV("dt", 1, dt));
		//this->_events["dt"] = pdt;
		//this->_dt = dt;
		//this->_max_step = max_step;
		//spEvent pms(new FlagV("max_step", 1, max_step));
		//this->_events["max_step"] = pms;

		this->_timestep = TimeStep::Creator(scheme);

		this->_timestep->set_dt(dt);
		this->_timestep->set_max_step(max_step);
	}
	void unset_time_term() {
		this->_events.erase("time_term");
		this->_timestep = nullptr;
	}

	/*
	 * Boundary condition
	 */
	void set_boundary_condition(St si, St segi, St vi, pBoundaryCondition pbc) {
		_pdomain->_pbindex->insert(si, segi, vi, pbc);
	}

	/*
	 * Show valuable table
	 */
	void show_var_center() {
		for (auto& term : this->_vars_c) {
			fmt::print("({:<20} : {})\n", term.first, term.second);
		}
	}

	void show_var_ut() {
		for (auto& term : this->_vars_ut) {
			fmt::print("({:<20} : {})\n", term.first, term.second);
		}
	}

	/*
	 * add var center
	 */
	void add_var_center(const std::string& name) {
		St max = this->max_idx(this->_vars_c);
		this->_vars_c[name] = max + 1;
		this->_pdomain->resize_data(this->max_idx(this->_vars_c) + 1, 0, 0,
				this->max_idx(this->_vars_ut) + 1);
	}

	St get_var_center_idx(const std::string& name) const {
		Variables::const_iterator got = this->_vars_c.find(name);
		ASSERT_MSG(got != this->_vars_c.end(),
				("No var named : " + name).c_str());
		return got->second;
	}

protected:

	/*
	 *
	 */
	const_pBoundaryCondition _find_bc(const_pNode pg, St vali) {
		//input a ghost node
		BoundaryCondition res;
		typename Ghost::GhostID gid = Ghost::ToGhostID(pg);
		// find in map ghost
		typename Ghost::iterator iter = this->_pdomain->p_ghost()->find(gid);
		ASSERT(iter != this->_pdomain->p_ghost()->end());
		// find in BoundaryIndex
		St si = iter->second.shape_idx;
		St segi = iter->second.seg_idx;
		return this->_pdomain->find_bc(si, segi, vali);
	}

	int _exp_substitute_ghost(pNode pn, Exp& exp, St vali) {
		for (ExpIter iter = exp.begin(); iter != exp.end();) {
			if (Exp::IsGhostNode(iter)) {
				ExpIter itere = iter;
				++iter;
				// find boundary condition
				const_pNode pg = Exp::GetpNode(itere);
				if (!Exp::IsZeroCoe(itere)) { // if coe = 0 the term will be deleted
					const_pBoundaryCondition pbc = _find_bc(pg, vali);
					if (pbc->get_type() == BoundaryCondition::_BC1_) {
						// Boundary condition 1
						vt val = pbc->get_val(pg->cp(_X_), pg->cp(_Y_),
								pg->cp(_Z_));
						//vt vol = pg->volume();
						vt coe = Exp::GetCoe(itere);
						exp.insert(coe * val, pn, 0);
					} else { // bc2
						//
						//SHOULD_NOT_REACH;
						typename Ghost::GhostID gid = Ghost::ToGhostID(pg);
						ASSERT(gid.step == 0);
						const_pNode po = pg->father;
						Axes a = FaceDirectionToAxes(gid.direction);
						vt val = pbc->get_val(pg->cp(_X_), pg->cp(_Y_),
								pg->cp(_Z_), a);
						cvt dl = (pg->cp(a) - po->cp(a));
						vt coe = Exp::GetCoe(itere);
						// pg = po - val*dl
						exp.insert(coe, po, 1);
						exp.insert(-val * dl * coe, po, 0);			//constant
						//SHOULD_NOT_REACH;
						//vt val = pbc->get_val(pg->cp(_X_), pg->cp(_Y_),
						//		pg->cp(_Z_));
						//exp.insert(Term(val, pn, 0));
					}
				}
				exp.erase(itere);
			} else {
				++iter;
			}
		}
		exp.concise();
		return 1;
	}

	vt _CFL_number(pNode pn, vt dt, St vidx) {
		vt veo[Dim];
		vt cfl[Dim];
		for (St i = 0; i < Dim; i++) {
			veo[i] = pn->cdva(vidx);
			cfl[i] = veo[i] * dt / pn->d(ToAxes(i));
		}
		return Max(cfl, Dim);
	}

	void _arr_to_grid(const Arr& x, St idx) {
		pGrid pgrid = this->_pdomain->pgrid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			it->cd(idx) = x[it->d_idx()];
		}
	}
	void _grid_to_arr(Arr& x, St idx) {
		pGrid pgrid = this->_pdomain->pgrid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			x[it->d_idx()] = it->cd(idx);
		}
	}

	void _build_node_exp(std::function<void(pNode, Exp&)> fun_node_exp,
			St idx) {
		Grid& grid = _pdomain->grid();
		// new map FE on leaf
		for (typename Grid::iterator_leaf it = grid.begin_leaf();
				it != grid.end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			utPointer& utp = pn->utp(this->_vars_ut["NE"]);
			if (utp != nullptr) {
				Arr_spExp& arrne = CAST_REF(pArr_spExp, utp);
				fun_node_exp(pn, *(arrne[idx]));
			} else {
				SHOULD_NOT_REACH;
			}
		}
	}

	void _build_matrix(Mat& mat, Arr& b, //
			Fun_pNode_spExp fun) {
		//2 Treverse face
		//  to calculate the expression for solver
		typedef std::list<St> ListST;
		typedef std::list<vt> ListVT;
		ListST l_rowptr;
		l_rowptr.push_back(0);
		ListST l_colid;
		ListVT l_val;
		ListVT l_b;
		int countnz = 0;

		pGrid pgrid = _pdomain->pgrid();

		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			spExp spexp = fun(pn);
			Exp& exp = (*spexp);
			exp.concise();

			int fconst = 0;
			for (typename Exp::iterator ite = exp.begin(); ite != exp.end();
					++ite) {
				if (!Exp::IsConstant(ite)) {
					const_pNode pn = Exp::GetpNode(ite);
					vt val = Exp::GetCoe(ite);
					l_colid.push_back(pn->d_idx());
					l_val.push_back(val);
					countnz++;
				} else {
					fconst = 1;
					vt val = Exp::GetCoe(ite);
					l_b.push_back(-val);    //!!!!! negative added here
				}
			}
			if (fconst == 0) {
				l_b.push_back(0);
			}
			l_rowptr.push_back(countnz);
			//exp.show();
		}
		//copy list to array  ======================
		St nr = l_rowptr.size() - 1;
		St nz = l_val.size();
		ASSERT(nz == l_colid.size());
		ASSERT(nr <= nz);
		mat.newsize(nr, nr, nz);
		b.reconstruct(l_b.size());
		int i = 0;
		for (typename ListST::iterator it = l_colid.begin();
				it != l_colid.end(); ++it) {
			mat.col_ind(i) = (*it);
			i++;
		}
		i = 0;
		for (typename ListST::iterator it = l_rowptr.begin();
				it != l_rowptr.end(); ++it) {
			mat.row_ptr(i) = (*it);
			i++;
		}
		i = 0;
		for (typename ListVT::iterator it = l_val.begin(); it != l_val.end();
				++it) {
			mat.val(i) = (*it);
			i++;
		}
		i = 0;
		for (typename ListVT::iterator it = l_b.begin(); it != l_b.end();
				++it) {
			b[i] = (*it);
			i++;
		}
	}

}
;

}

#endif
