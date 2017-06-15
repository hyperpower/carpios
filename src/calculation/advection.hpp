#ifndef _ADVECTION_H_
#define _ADVECTION_H_
//#define __Debug__
#include "calculation_define.hpp"
#include "expression.hpp"

#include <vector>
#include <iomanip>

namespace carpio {
//This file use to solve advection equation
// 1 non-conservation from
//
//       d(phi)              d( phi)       2D --version
//    u ------------ +  v -------------  = 0
//          dx                 dy
//
//      d(phi)       d(phi)       d(phi)       3D --version
//    u -------- + v -------  + w -------  = 0
//        dx           dy          dz
//
// 2 conservation from
//
//      d( u phi)            d( v phi)       2D --version
//     --------------- +   --------------  = 0
//          dx                 dy
//
//      d(u phi)     d( v phi)       d(w phi)       3D --version
//     ----------- + ---------  +  ----------  = 0
//            dx        dy             dz
//
// 2 with time
//
//   d(phi)    d(u phi)    d(v phi)       2D --version
//   ------ +  ------- +  -------  = 0
//     dt         dx         dy
//
//   d(phi)     d(u phi)  d(v phi)   d(w phi)       3D --version
//   ------ +  ------- +  -------  + -------  = 0
//     dt         dx         dy          dz

/*
 * the Advection class
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class Advection_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Advection_<COO_VALUE, VALUE, DIM> Self;
	typedef Advection_<COO_VALUE, VALUE, DIM>& ref_Self;
	typedef Domain_<COO_VALUE, VALUE, DIM> Domain;
	typedef Domain_<COO_VALUE, VALUE, DIM>& ref_Domain;
	typedef const Domain_<COO_VALUE, VALUE, DIM>& const_ref_Domain;
	typedef Domain_<COO_VALUE, VALUE, DIM>* pDomain;
	typedef const Domain_<COO_VALUE, VALUE, DIM>* const_pDomain;
	typedef typename Domain::BoundaryIndex::BoundaryCondition BoundaryCondition;
	typedef typename Domain::BoundaryIndex::pBoundaryCondition pBoundaryCondition;
	typedef typename Domain::BoundaryIndex::const_pBoundaryCondition const_pBoundaryCondition;

	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pGrid;
	typedef Ghost_<COO_VALUE, VALUE, DIM> Ghost;
	typedef const Ghost_<COO_VALUE, VALUE, DIM> const_Ghost;
	typedef Ghost_<COO_VALUE, VALUE, DIM>* pGhost;
	typedef const Ghost_<COO_VALUE, VALUE, DIM>* const_pGhost;
	typedef typename Ghost::GhostNode GhostNode;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef const Node_<COO_VALUE, VALUE, DIM> *const_pNode;

	typedef typename Domain::Face Face;
	typedef typename Domain::pFace pFace;
	typedef Shape_<COO_VALUE, DIM> Shape;
	typedef Shape_<COO_VALUE, DIM>* pShape;

	typedef Stencil_<COO_VALUE, VALUE, DIM, DIM> Stencil;

	typedef Interpolate_<COO_VALUE, VALUE, DIM> Interpolate;

	typedef std::function<vt(cvt, cvt, cvt)> Function;

	typedef Expression_<COO_VALUE, VALUE, DIM> Exp;
	typedef Exp* pExp;
	typedef typename Exp::Term Term;
	typedef typename Exp::iterator ExpIter;
	typedef typename Exp::const_iterator const_ExpIter;

	typedef std::function<vt(vt)> Limiter;

protected:
	//matrix
	typedef MatrixSCR_<vt> Mat;
	typedef ArrayListV<vt> Arr;
	typedef ArrayListV<St> Arr_st;

	pDomain _pdomain;
	// time ---------------------------
	vt _dt;
	St _max_step;

	St _scheme;
	std::vector<Limiter> _v_limiter;

	Arr_st _phi_idx;
	St _v_idx[Dim];

	St _utp_idx;

	// flag ---------------------------
	int _flag_conserve; //if 0 --->  unconserve
						//   1 --->    conserve   (default)
	int _flag_time;     //if 0 --->  no time term (default)
						//   1 --->  with time

public:
	/*
	 * constructor
	 */
	Advection_(pDomain pf, St scheme = 0, int flag_time = 0) :
			_pdomain(pf), _scheme(scheme) {
		_push_limiter();
		_utp_idx = 0;
		for (St i = 0; i < Dim; ++i) {
			_v_idx[i] = i;
		}
		if (flag_time == 1) {
			_phi_idx.reconstruct(3);
			_phi_idx[0] = Dim;
			_phi_idx[1] = Dim+1;
			_phi_idx[2] = Dim+2;
			_pdomain->new_data(this->max_idx() + 1, 0, 0, 1);
		} else {
			_phi_idx.reconstruct(1);
			_phi_idx[0] = Dim;
			_pdomain->new_data(this->max_idx() + 1, 0, 0, 1);
		}
		_flag_conserve = 1;
		_flag_time = flag_time;
	}

	/*
	 * set
	 */
	void set_nonconserve() {
		_flag_conserve = 0;
	}

	void set_v_grid(Function pfun, Axes a) {
		ASSERT(a < Dim);
		_pdomain->set_val_grid(_v_idx[a], pfun);
	}

	void set_v_ghost(Function pfun, Axes a) {
		_pdomain->set_val_ghost(_v_idx[a], pfun);
	}
	void set_v(Function pfun, Axes a) {
		this->set_v_grid(pfun, a);
		this->set_v_ghost(pfun, a);
	}
	void set_phi_ghost() {
		// the boundary must be set
		_pdomain->set_val_ghost_by_bc(_phi_idx[0]);
	}
	// can be used to set initial boundary condition
	void set_phi(Function pfun) {
		_pdomain->set_val(_phi_idx[0], pfun);
	}

	void set_time(vt dt, St max_step) {
		_dt = dt;
		_max_step = max_step;
	}
	/*
	 * get pDomain
	 */
	pDomain get_pDomain(){
		return _pdomain;
	}

	const_pDomain get_pDomain() const{
		return _pdomain;
	}


	/*
	 * this function return the max index of valuable used in
	 * poisson, the max index make sure the array size is ok.
	 */
	St max_idx() {
		St max = 0;
		St m1 = Max(_v_idx, Dim);
		St m2 = _phi_idx.max();
		max = Max(max, m1);
		max = Max(max, m2);
		return max;
	}
	/*
	 * Boundary condition
	 */
	void set_boundary_condition(St si, St segi, St vi, pBoundaryCondition pbc) {
		_pdomain->_pbindex->insert(si, segi, vi, pbc);
	}
	/*
	 * this function can be deleted
	 */
	int _face_scheme_fou(pFace pface, Exp& exp) {
		//face type
		switch (pface->ft()) {
		//_Null_ = -1, _Boundary_ = 0, _Equal_ = 1, _FineCoarse_ = 2, _CoarseFine_ = 3,
		case _Null_:
			SHOULD_NOT_REACH;
			break;
		case _Boundary_: {
			this->_face_scheme_equal_fou(pface, exp);
			break;
		}
		case _Equal_: {
			this->_face_scheme_equal_fou(pface, exp);
			break;
		}
		case _FineCoarse_: {
			this->_face_scheme_equal_fou(pface, exp);
			break;
		}
		case _CoarseFine_: {
			SHOULD_NOT_REACH;
			std::cout << "Corase Fine function unfinish\n";
			break;
		}
		default:
			return -1;
		}
		return pface->ft();
	}

	/*
	 * this function can be deleted
	 */
	int _face_scheme_tvd(pFace pface, Exp& exp, St iphi) {
		//face type
		switch (pface->ft()) {
		//_Null_ = -1, _Boundary_ = 0, _Equal_ = 1, _FineCoarse_ = 2, _CoarseFine_ = 3,
		case _Null_:
			SHOULD_NOT_REACH;
			break;
		case _Boundary_: {
			this->_face_scheme_equal_tvd(pface, exp, iphi);
			break;
		}
		case _Equal_: {
			this->_face_scheme_equal_tvd(pface, exp, iphi);
			break;
		}
		case _FineCoarse_: {
			this->_face_scheme_equal_tvd(pface, exp, iphi);
			break;
		}
		case _CoarseFine_: {
			SHOULD_NOT_REACH;
			std::cout << "Corase Fine function unfinish\n";
			break;
		}
		default:
			return -1;
		}
		return pface->ft();
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

	int _face_scheme_equal_fou(pFace pface, Exp& exp) {
		// face direction
		pNode pori = pface->pori();
		Direction dir = pface->dir();
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(dir, o, a);
		// interpolate veo_f
		Exp expv = Interpolate::OnFace(pori, dir, 1);
		Float veo_f = expv.subsitute(_v_idx[a]);

		//get U C D ------------------------------
		pNode pC = _find_C(pface, veo_f);
		//
		// exp should be empty before insert
		if (_flag_conserve == 1) {
			exp.insert(veo_f, pC, 1.0);
		} else {
			exp.insert(1.0, pC, 1.0);
		}
		return 1;
	}

	int _face_scheme_equal_tvd(pFace pface, Exp& exp, St iphi) {
		// face direction
		pNode pori = pface->pori();
		pNode pnei = pface->pnei();
		Direction dir = pface->dir();
		Orientation o;
		Axes a;
		FaceDirectionToOrientationAndAxes(dir, o, a);
		// interpolate veo_f
		Exp expv = Interpolate::OnFace(pori, dir, 1);
		Float veo_f = expv.subsitute(_v_idx[a]);

		//get U C D ------------------------------
		pNode pU = nullptr;
		pNode pC = nullptr;
		pNode pD = nullptr;
		Direction dircu = dir;
		if (IsP(o)) {
			if (veo_f > 0) {
				pC = pori;
				pD = pnei;
				dircu = Opposite(dir);
			} else {
				pC = pnei;
				pD = pori;
				dircu = dir;
			}
		} else {
			if (veo_f > 0) {
				pC = pnei;
				pD = pori;
				dircu = dir;
			} else {
				pC = pori;
				pD = pnei;
				dircu = Opposite(dir);
			}
		}
		pU = pC->get_neighbor_fast(dircu);
		vt cor = 0;
		if (pU != nullptr) {
			// cal limiter
			//st iphi = _phi_idx[0];
			vt vU = pU->cda(iphi);
			vt vC = pC->cda(iphi);
			vt vD = pD->cda(iphi);
			// cal \Psi(r)
			vt r = (vC - vU) / (vD - vC + SMALL);
			vt psi = limiter(r, _scheme);
			cor = 0.5 * psi * (vD - vC);
		}
		if (_flag_conserve == 1) {
			exp.insert(cor * veo_f, pC, 0.0);
			exp.insert(1.0 * veo_f, pC, 1.0);
		} else {
			exp.insert(cor, pC, 0.0);
			exp.insert(1.0, pC, 1.0);
		}
		return 1;
	}
	/*
	 * build face exp to utp
	 */
protected:
	typedef std::list<std::pair<pFace, pExp> > ListFE;
	typedef std::pair<pFace, pExp> PairFE;
public:
	void _face_exp_fou(pNode& pn) {
		_new_list_face_exp(pn, _utp_idx);
		utPointer& utp = pn->utp(_utp_idx);
		for (St i = 0; i < NumFaces; i++) {
			Direction dir = FaceDirectionInOrder(i);
			pNode pnei = pn->get_neighbor_fast(dir);
			// The face need to be calculate
			// 1 the boundary face
			// 2 the equal face on P direction;
			// 3 the fine to coarse face
			FaceType ft = GetFaceType(pn, pnei);
			pExp pexp = new Exp();
			bool df = true;
			if ((ft == _Boundary_)				//1
			|| (ft == _Equal_ && (IsFacePDirection(dir)))				//2
					|| (ft == _FineCoarse_))				//3
					{
				// work on pn
				pFace pf = new Face(pn, pnei, dir, ft);
				_face_scheme_fou(pf, *pexp);
				ListFE& lpexp = CAST_REF(ListFE*, utp);
				PairFE pfe(pf, pexp);
				lpexp.push_back(pfe);
				df = false;
			}
			// case 2 3
			if ((ft == _Equal_ && (IsFacePDirection(dir)))				//2
			|| (ft == _FineCoarse_))				//3
					{
				// work on pnei
				_new_list_face_exp(pnei, _utp_idx);
				utPointer& utpn = pnei->utp(_utp_idx);
				FaceType oft = (ft == _FineCoarse_) ? _CoarseFine_ : _Equal_;
				pFace fn = new Face(pnei, pn, Opposite(dir), oft);
				pExp pexpn = new Exp(*pexp);
				//pexpn->times(-1.0);
				ListFE& lpexp = (*CAST(ListFE*, utpn));
				PairFE pairn(fn, pexpn);
				lpexp.push_back(pairn);
			}
			if (df == true) {
				delete pexp;
			}
		}
	}
	void _face_exp_tvd(pNode& pn, St iphi) {
		_new_list_face_exp(pn, _utp_idx);
		utPointer& utp = pn->utp(_utp_idx);
		for (St i = 0; i < NumFaces; i++) {
			Direction dir = FaceDirectionInOrder(i);
			pNode pnei = pn->get_neighbor_fast(dir);
			// The face need to be calculate
			// 1 the boundary face
			// 2 the equal face on P direction;
			// 3 the fine to coarse face
			FaceType ft = GetFaceType(pn, pnei);
			pExp pexp = new Exp();
			bool df = true;
			if ((ft == _Boundary_)				//1
			|| (ft == _Equal_ && (IsFacePDirection(dir)))				//2
					|| (ft == _FineCoarse_))				//3
					{
				// work on pn
				pFace pf = new Face(pn, pnei, dir, ft);
				_face_scheme_tvd(pf, *pexp, iphi);
				ListFE& lpexp = CAST_REF(ListFE*, utp);
				PairFE pfe(pf, pexp);
				lpexp.push_back(pfe);
				df = false;
			}
			// case 2 3
			if ((ft == _Equal_ && (IsFacePDirection(dir)))				//2
			|| (ft == _FineCoarse_))				//3
					{
				// work on pnei
				_new_list_face_exp(pnei, _utp_idx);
				utPointer& utpn = pnei->utp(_utp_idx);
				FaceType oft = (ft == _FineCoarse_) ? _CoarseFine_ : _Equal_;
				pFace fn = new Face(pnei, pn, Opposite(dir), oft);
				pExp pexpn = new Exp(*pexp);
				//pexpn->times(-1.0);
				ListFE& lpexp = (*CAST(ListFE*, utpn));
				PairFE pairn(fn, pexpn);
				lpexp.push_back(pairn);
			}
			if (df == true) {
				delete pexp;
			}
		}
	}
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
	int _exp_substitute_ghost(pNode pn, Exp& exp, St phii = 0) {
		for (ExpIter iter = exp.begin(); iter != exp.end();) {
			if (Exp::IsGhostNode(iter)) {
				ExpIter itere = iter;
				++iter;
				// find boundary condition
				const_pNode pg = Exp::GetpNode(itere);
				if (!Exp::IsZeroCoe(itere)) { // if coe = 0 the term will be deleted
					const_pBoundaryCondition pbc = _find_bc(pg, _phi_idx[phii]);
					if (pbc->get_type() == BoundaryCondition::_BC1_) {
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

	void _node_exp_conserve(pNode pn, Exp& exp) {
		//assert(pn->d());
		ListFE& lpexp = CAST_REF(ListFE*, pn->utp(_utp_idx));
		Exp FF[NumFaces];
		Arr_st countCF(NumFaces);
		countCF.assign(0);
		for (typename ListFE::iterator iter = lpexp.begin();
				iter != lpexp.end(); ++iter) {
			Face* pface = iter->first;
			Exp* pexp = iter->second;

			ASSERT(pface->pori() == pn); //
			Direction dir = pface->dir();
			// times area
			int sign = IsFacePDirection(dir) ? 1 : -1;
			vt area = pface->area();
			if (pface->ft() == _CoarseFine_) {
				area = pface->pnei()->face_area(Opposite(dir));
			}
			pexp->times(sign * area);
			FF[FaceDirectionInOrder(dir)].plus(*pexp);
			countCF[FaceDirectionInOrder(dir)]++;
		}
		for (St i = 0; i < NumFaces; i++) {
			exp.plus(FF[i]);
#ifdef __Debug__
			if (pn->is_in_on(0.2, 0.4)) {
				std::cout<< "face i = "<<i<<"\n";
				_show_exp(exp, *_pdomain);
			}
#endif
		}
		//exp.times(-1); //negative;
#ifdef __Debug__
		if (pn->is_in_on(0.2, 0.4)) {
			Exp expc(exp);
			expc.divide(pn->volume());
			std::cout << "expc ------- " << "\n";
			_show_exp(expc, *_pdomain);
		}
#endif
		//
		//
		_delete_list_face_exp(pn, _utp_idx);
	}

	vt _CFL_number(pNode pn, vt dt) {
		vt veo[Dim];
		vt cfl[Dim];
		for (St i = 0; i < Dim; i++) {
			veo[i] = pn->cdva(_v_idx[i]);
			cfl[i] = veo[i] * dt / pn->d(ToAxes(i));
		}
		return Max(cfl, Dim);
	}

	void _node_exp(pNode pn, Exp& exp) {
		//assert(pn->d());
		ListFE& lpexp = CAST_REF(ListFE*, pn->utp(_utp_idx));
		Exp FF[NumFaces];
		Arr_st countCF(NumFaces);
		countCF.assign(0);
		for (typename ListFE::iterator iter = lpexp.begin();
				iter != lpexp.end(); ++iter) {
			Face* pface = iter->first;
			Exp* pexp = iter->second;
			ASSERT(pface->pori() == pn); //
			Direction dir = pface->dir();
			FF[FaceDirectionInOrder(dir)].plus(*pexp);
			countCF[FaceDirectionInOrder(dir)]++;
		}
		for (St i = 0; i < NumFaces; i++) {
#ifdef __Debug__
			if (pn->is_in_on(0.2, 0.4)) {
				std::cout << "ff -------- " << i << "\n";
				_show_exp(FF[i], *_pdomain);
			}
#endif
			if (countCF[i] > 1) {
				// SHOULD_NOT_REACH;
				// get average face velocity;
				FF[i].times(1.0 / Float(countCF[i]));
			}
		}
		//
		// x direction
		Axes a = _X_;
		St ip = FaceDirectionInOrder(a, _P_);
		St im = FaceDirectionInOrder(a, _M_);
		Exp exp_x(FF[ip]);
		exp_x.minus(FF[im]);
		Float l = pn->d(a);
		Float u = pn->cdva(_v_idx[a]);
		Float CFL_x = u / l;
		//ASSERT_MSG(CFL_x < 1, "CFL x > 0.5");
		exp_x.times(-CFL_x);  //negative

		// y direction
		a = _Y_;
		ip = FaceDirectionInOrder(a, _P_);
		im = FaceDirectionInOrder(a, _M_);
		Exp exp_y(FF[ip]);
		exp_y.minus(FF[im]);
		l = pn->d(a);
		u = pn->cdva(_v_idx[a]);
		Float CFL_y = u / l;
		//ASSERT_MSG(CFL_y < 0.5, " CFL y > 0.5");
		exp_y.times(-CFL_y); //negative
		exp.plus(exp_x);
		exp.plus(exp_y);
		//
		_delete_list_face_exp(pn, _utp_idx);
	}

	void _show_exp(Exp& exp, Domain& domain) {
		exp.show();
		std::list<Gnuplot_actor> lga;
		Gnuplot_actor ga;
		GnuplotActor_LeafNodes(ga, domain.grid());
		lga.push_back(ga);
		GnuplotActor_Expression(ga, exp);
		lga.push_back(ga);
		Gnuplot gp;
		gp.set_equal_ratio();
		//gp.set_xrange(2.0,3.0);
		//gp.set_yrange(1.5,2.5);
		//gp.set_cbrange(-2.0, 2.0);
		//gp.plot(lga);
		//delete shape
	}

	int _bulid_matrix_fou(Mat& mat, Arr& b) { //
		//1 Traverse face
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_face_exp_fou(pn);
		}

		typedef std::list<St> ListST;
		typedef std::list<vt> ListVT;
		ListST l_rowptr;
		l_rowptr.push_back(0);
		ListST l_colid;
		ListVT l_val;
		ListVT l_b;
		int countnz = 0;
		//2 build matrix
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			Exp exp;
			pNode pn = it.get_pointer();
			if (_flag_conserve == 1) {
				_node_exp_conserve(pn, exp);
			} else {
				_node_exp(pn, exp);
			}

			_exp_substitute_ghost(pn, exp);
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
		return 1;
	}

	int _bulid_matrix_tvd(Mat& mat, Arr& b) { //
		//1 Traverse face
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			pNode pn = it.get_pointer();
			_face_exp_tvd(pn, _phi_idx[0]);
		}

		typedef std::list<St> ListST;
		typedef std::list<vt> ListVT;
		ListST l_rowptr;
		l_rowptr.push_back(0);
		ListST l_colid;
		ListVT l_val;
		ListVT l_b;
		int countnz = 0;
		//2 build matrix
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			Exp exp;
			pNode pn = it.get_pointer();
#ifdef __Debug__
			if (pn->is_in_on(0.2, 0.4)) {
				std::cout << "stop\n";
			}
#endif
			if (_flag_conserve == 1) {
				_node_exp_conserve(pn, exp);
			} else {
				_node_exp(pn, exp);
			}
			// debug

			_exp_substitute_ghost(pn, exp);
			exp.concise();
#ifdef __Debug__
			if (pn->is_in_on(0.2, 0.4)) {
				//pn->show();
				_show_exp(exp, *_pdomain);
			}
#endif
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
		return 1;
	}

	int solve_fou(vt tol = 1e-6, int max_iter = 1000, int info = 0) {
		std::list<vt> lr;
		int rcode = this->solve_fou(tol, max_iter, lr, info);
		//std::cout<< "max iter "<<max_iter<< " tol "<< tol <<"\n";
		return rcode;
	}

	int solve_tvd(vt tol = 1e-6, int max_iter = 1000, int info = 0) {
		std::list<vt> lr;
		int rcode = this->solve_tvd(tol, max_iter, lr, info);
		//std::cout<< "max iter "<<max_iter<< " tol "<< tol <<"\n";
		return rcode;
	}

	int solve_fou(vt& tol, int& max_iter, std::list<vt>& lr, int info) {
		Mat mat;
		Arr b;
		this->_bulid_matrix_fou(mat, b);
		Arr x(b.size());
		// initial x
		//_grid_to_arr(x);
		//x.show();
		//defalt set up ========
		//std::list<Float> lr;	//list residual
		//solver =======================
		int sf = Dia_BiCGSTAB(mat, x, b, max_iter, tol, lr);
		if (info == 1) {
			std::cout << "iter n " << lr.size() << " residual " << tol << "\n";
		}
		if (sf != 0) {
			std::cerr << " >! Poisson solve failed \n";
			return -1;
		}
		//gnuplot_show_ylog(lr);
		//put the value back
		_arr_to_grid(x);
		return 1;
	}

	int solve_tvd(vt& tol, int& max_iter, std::list<vt>& lr, int info = 0) {
		if (info == 1) {
			std::cout << "Solve tvd scheme " << this->_scheme << "\n";
		}
		Mat mat;
		Arr b;
		this->_bulid_matrix_fou(mat, b);
		Arr x(b.size());
		// initial x
		//_grid_to_arr(x);
		//x.show();
		//defalt set up ========
		std::list<vt> lr1;	 //list residual
		//solver =======================
		vt tol1 = tol * 1000;
		int max_iter1 = 1000;
		int sf = Dia_BiCGSTAB(mat, x, b, max_iter1, tol1, lr1);
		if (info == 1) {
			std::cout << "Step 1 --> solve first order upwind\n";
			std::cout << "SP mat  nx = " << mat.size1() << " ny = "
					<< mat.size2() << " nonzero " << mat.NumNonzeros() << "\n";
			std::cout.precision(5);
			std::cout << std::scientific << "iter num " << lr1.size()
					<< " residual " << tol1 << "\n";
		}
		if (sf != 0) {
			std::cerr << " >! Poisson solve failed \n";
			return -1;
		}
		lr.push_back(tol1);  //
		//gnuplot_show_ylog(lr);
		//put the value back
		_arr_to_grid(x);
		//std::cout << "finish fou =====\n";
		//
		int ic = 0;
		if (info == 1) {
			std::cout << "Step 2 --> solve tvd\n";
		}
		vt stresid = tol + 0.1;
		while (ic < max_iter && stresid > tol) {   // revise check the residual
			Mat mat2;
			Arr b2;
			this->_bulid_matrix_tvd(mat2, b2);
			Arr x2(b2.size());
			lr.clear();
			// put the value to x
			_grid_to_arr(x2);
			//solver =======================
			std::list<vt> lr2;
			vt tol2 = tol/10000.0;
			int max_iter2 = 100 + ic * 10;
			if (info == 1) {
				stresid = Residual(mat2, x2, b2);
				//std::cout << "SP mat  nx = " << mat2.size1() << " ny = "
				//		<< mat2.size2() << " nonzero " << mat2.NumNonzeros()
				//		<< "\n";
				std::cout << ic << " r start " << stresid;
			}
			sf = Dia_BiCGSTAB(mat2, x2, b2, max_iter2, tol2, lr2);
			if (info == 1) {
				std::cout << " iter num " << lr2.size() << " r end " << tol2
						<< "\n";
			}
			//gnuplot_show_ylog(lr);
			if (sf != 0 && sf != 1) {
				std::cerr << " >! solver failed \n";
				return -1;
			}
			lr.push_back(tol2);
			// put the value back
			_arr_to_grid(x2);
			ic++;
		}
		max_iter = ic;
		return 0;
	}

	void advance(int info = 0) {
		St istep = 0;
		while (istep < _max_step) {
			if (info == 1) {
				std::cout << std::setw(5);
				std::cout << "i = " << istep;
			}
			//for (int i = 0; i < step; i++) {
			// Prediction step --------------------------------
			//1 Traverse face
			pGrid pgrid = this->_pdomain->p_grid();
			for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
					it != pgrid->end_leaf(); ++it) {
				pNode pn = it.get_pointer();
				_face_exp_fou(pn);
			}
			//2 add time
			for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
					it != pgrid->end_leaf(); ++it) {
				pNode pn = it.get_pointer();
				Exp exp;
				_node_exp_conserve(pn, exp); //advance half dt
				// cfl number check
				vt cfl = _CFL_number(pn, _dt * 0.5);
				ASSERT(cfl < 1.0);
				//
				exp.times(-0.5 * _dt / pn->volume()); //negative;
				exp.insert(1.0, pn, 1);   //phi n
				pn->cd(_phi_idx[1]) = exp.subsitute(_phi_idx[0]);
			}
			// Correction step add limiter: limite variables ------------
			for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
					it != pgrid->end_leaf(); ++it) {
				pNode pn = it.get_pointer();
				_face_exp_tvd(pn, _phi_idx[1]);
			}
			for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
					it != pgrid->end_leaf(); ++it) {
				pNode pn = it.get_pointer();
				Exp exp;
				_node_exp_conserve(pn, exp); //advance half dt
				// cfl number check
				vt cfl = _CFL_number(pn, _dt * 0.5);
				ASSERT(cfl < 1.0);
				//
				exp.times(-1.0 * _dt / pn->volume());  //negative;
				vt valn = exp.subsitute(_phi_idx[1]);
				pn->cd(_phi_idx[2]) = pn->cd(_phi_idx[0]) + valn;
				//
			}
			// Put back -------------------------------------------------
			for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
					it != pgrid->end_leaf(); ++it) {
				pNode pn = it.get_pointer();
				pn->cd(_phi_idx[0]) = pn->cd(_phi_idx[2]);
			}
			if (info == 1) {
				std::cout << "\n";
			}
			istep++;
		}
	}

protected:
	void _arr_to_grid(const Arr& x) {
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			it->cd(_phi_idx[0]) = x[it->d_idx()];
		}
	}
	void _grid_to_arr(Arr& x) {
		pGrid pgrid = this->_pdomain->p_grid();
		for (typename Grid::iterator_leaf it = pgrid->begin_leaf();
				it != pgrid->end_leaf(); ++it) {
			x[it->d_idx()] = it->cd(_phi_idx[0]);
		}
	}

	int _new_list_face_exp(pNode& pn, St utp_idx) {
		utPointer& utp = pn->utp(utp_idx);
		if (utp == nullptr) {
			utp = new ListFE();   //new !!!!!
			return 1;
		}
		return 0;
	}
	int _delete_list_face_exp(pNode pn, St utp_idx) {
		utPointer& utp = pn->utp(utp_idx);
		if (utp != nullptr) {
			ListFE& lpexp = CAST_REF(ListFE*, utp);
			for (typename ListFE::iterator iter = lpexp.begin();
					iter != lpexp.end(); ++iter) {
				delete iter->first;
				delete iter->second;
			}
			lpexp.clear();
			delete &lpexp;
			utp = nullptr;
			return 1;
		}
		return 0;
	}

	void _push_limiter() {
		Limiter fou = [](vt r) {return 0;};

		Limiter minmod = [](vt r) {
			return Max(0.0, Min(1.0, r));
		};

		Limiter superbee = [](vt r) {
			return Max(0.0, Min(2.0 * r, 1.0), Min(r, 2.0));
		};

		Limiter vanLeer = [](vt r) {
			return (r + Abs(r)) / (1 + r);
		};

		_v_limiter.push_back(fou);
		_v_limiter.push_back(minmod);
		_v_limiter.push_back(superbee);
		_v_limiter.push_back(vanLeer);
	}

	inline Float limiter(Float r, int i) {
		ASSERT(i < _v_limiter.size());
		return _v_limiter[i](r);
	}

}
;
}

#endif

