#ifndef _S_SOLVER_HPP
#define _S_SOLVER_HPP

#include <structure/s_data.hpp>
#include "s_define.hpp"
#include <map>
#include <utility>
#include <stdio.h>
#include "utility/any.hpp"
#include "s_ghost.hpp"

namespace structure {

template<St DIM>
class Solver_ {
public:
	static const St Dim = DIM;
	typedef Index_<DIM> Index;
	typedef Equation_<DIM> Equation;
	typedef std::shared_ptr<Equation> spEquation;

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
	typedef Order_<DIM> Order;
	typedef std::shared_ptr<Order_<DIM> > spOrder;

	typedef std::map<std::string, carpio::Any> Flags;
protected:
	Flags _flags;

	spGrid _grid;
	spGhost _ghost;
	spStencil _stencil;

	spCenterScalar _spcs; // the unknow variable
	spCenterExp _spc_exp;
	std::string _vname;
	spOrder _order;

	//St _max_iter;
	//Vt _residual;

	int _max_iter;
	Vt _tol;

public:
	Solver_(spGrid spe, spGhost spg, spStencil sps, const std::string& vn,
			spCenterScalar spphi) :
			_grid(spe), _ghost(spg), _stencil(sps), _vname(vn), _spcs(spphi) {
		_spc_exp = spCenterExp(new CenterExp(spe));
		//_spcs = spCenterScalar(new CenterScalar(spe));
		_build_order();
		_max_iter = 1000;
		_tol = 1e-4;
	}

	Arr _copy_to_x() const {
		St n = this->_order->size();
		Arr x(n);
		for (Index ijk = _order->begin(); ijk != _order->end();
				ijk = _order->next(ijk)) {
			x[_order->get_order(ijk)] = _spcs->val(ijk);
		}
		return std::move(x);
	}

	void _copy_to_cs(const Arr& x) {
		for (Index idx = _order->begin(); idx != _order->end();
				idx = _order->next(idx)) {
			_spcs->val(idx) = x[_order->get_order(idx)];
		}
	}

	void _build_matrix(Mat& mat, Arr&b) {
		//2 Treverse face
		//  to calculate the expression for solver
		typedef std::list<St> ListST;
		typedef std::list<Vt> ListVT;
		ListST l_rowptr;
		l_rowptr.push_back(0);
		ListST l_colid;
		ListVT l_val;
		ListVT l_b;
		int countnz = 0;

		for (Index ijk = _order->begin(); ijk != _order->end();
				ijk = _order->next(ijk)) {
			Expression& exp = *(_spc_exp->val(ijk));

			int fconst = 0;
			for (typename Expression::iterator ite = exp.begin();
					ite != exp.end(); ++ite) {
				if (!Expression::IsConstant(ite)) {
					Index index = Expression::GetTerm(ite);
					Vt val = Expression::GetCoe(ite);
					l_colid.push_back(_order->get_order(index));
					l_val.push_back(val);
					countnz++;
				} else {
					fconst = 1;
					Vt val = Expression::GetCoe(ite);
					l_b.push_back(-val);    //!!!!! negative added here
				}
			}
			if (fconst == 0) {
				l_b.push_back(0);
			}
			l_rowptr.push_back(countnz);
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

	void _build_exp() {
		for (Index ijk = _order->begin(); ijk != _order->end();
				ijk = _order->next(ijk)) {
			_spc_exp->val(ijk) = this->_stencil->coe_row(ijk);
			_spc_exp->val(ijk) = this->_ghost->substitute(_spc_exp->val(ijk),
					this->_vname);
			_spc_exp->val(ijk)->concise();
		}
	}

	void set_max_iter(int max) {
		this->_max_iter = max;
	}

	void set_tolerance(Vt tol) {
		this->_tol = tol;
	}

	void set_output_residual(std::FILE* pf) {
		this->_flags["OutputResidual"] = pf;
	}

	bool has_flag(const std::string& key) const {
		auto it = this->_flags.find(key);
		if (it != this->_flags.end()) {
			return true;
		}
		return false;
	}

	void output_residual(std::FILE* f, std::list<Vt>& lr) {
		if (f != nullptr) {
			std::fseek(f, 0, SEEK_END);
			int n = 0;
			for (auto v : lr) {
				fmt::print(f, "{0:10d},{1:10.7e}\n", n, v);
				n++;
			}
			fmt::print(f, "\n"); //append a empty line
		}
	}

	virtual void solve() {
		SHOULD_NOT_REACH;
	}

	virtual ~Solver_() {

	}
protected:
	void _build_order() {
		if (this->_ghost->name() == "GhostRegular") {
			this->_order = spOrder(
					new OrderRegular_<DIM>(this->_grid, this->_ghost));
		} else {
			SHOULD_NOT_REACH;
		}
	}

}
;

template<St DIM>
class Solver_Jacobi_: public Solver_<DIM> {
public:
	static const St Dim = DIM;
	typedef Equation_<DIM> Equation;
	typedef std::shared_ptr<Equation> spEquation;
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

	typedef std::map<std::string, carpio::Any> Flags;
public:
	Solver_Jacobi_(spGrid spe, spGhost spg, spStencil sps,
			const std::string& vn, spCenterScalar spcs) :
			Solver_<DIM>(spe, spg, sps, vn, spcs) {
	}

	void solve() {
		this->_build_exp();
		Mat mat;
		Arr b;
		this->_build_matrix(mat, b);
		Arr x = this->_copy_to_x();
		int maxi = this->_max_iter;
		Vt tol = this->_tol;
		std::list<Vt> lr;
		//
		typedef carpio::Solver_Jacobi_<Vt> Sol;
		Sol solver(maxi, tol);
		solver.solve(mat, x, b);
		//int res_code = carpio::Jacobi(mat,          // A  The matrix
		//		x,            // x
		//		b,            // b
		//		maxi,         // max iter
		//		tol,          // Tolerance
		//		lr);          // list residual
		//std::cout << "num iter = " << solver.num_iter();
		//std::cout << "   residual = " << solver.residual() << "\n";

		//std::cout << maxi << "  " << tol << "  " << res_code << "\n";
		this->_copy_to_cs(x);

		if (this->has_flag("OutputResidual")) {
			std::FILE* f = carpio::any_cast<std::FILE*>(
					this->_flags["OutputResidual"]);
			this->output_residual(f, lr);
		}
	}

};

template<St DIM>
class Solver_IC_CGS_: public Solver_<DIM> {
public:
	static const St Dim = DIM;
	typedef Equation_<DIM> Equation;
	typedef std::shared_ptr<Equation> spEquation;
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
public:
	Solver_IC_CGS_(spGrid spe, spGhost spg, spStencil sps,
			const std::string& vn, spCenterScalar spcs) :
			Solver_<DIM>(spe, spg, sps, vn, spcs) {
	}

	void solve() {
		this->_build_exp();
		Mat mat;
		Arr b;
		this->_build_matrix(mat, b);
		Arr x = this->_copy_to_x();
		int maxi = this->_max_iter;
		Vt tol = this->_tol;
		std::list<Vt> lr;
		int res_code = carpio::IC_CGS(mat,          // A  The matrix
				x,            // x
				b,            // b
				maxi,         // max iter
				tol,          // Tolerance
				lr);          // list residual
		//std::cout << maxi << "  " << tol << "  " << res_code << "\n";
		this->_copy_to_cs(x);

		if (this->has_flag("OutputResidual")) {
			std::FILE* f = carpio::any_cast<std::FILE*>(
					this->_flags["OutputResidual"]);
			this->output_residual(f, lr);
		}
	}
};

template<St DIM>
class Solver_SOR_: public Solver_<DIM> {
public:
	static const St Dim = DIM;
	typedef Equation_<DIM> Equation;
	typedef std::shared_ptr<Equation> spEquation;
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

	typedef std::map<std::string, carpio::Any> Flags;

protected:
	Vt _omega;
public:
	Solver_SOR_(spGrid spe, spGhost spg, spStencil sps,
			const std::string& vn, spCenterScalar spcs, Vt omega = 1.0) :
				Solver_<DIM>(spe, spg, sps, vn, spcs) {
		this->_omega = omega;
	}

	void solve() {
		this->_build_exp();
		Mat mat;
		Arr b;
		this->_build_matrix(mat, b);
		Arr x = this->_copy_to_x();
		int maxi = this->_max_iter;
		Vt tol = this->_tol;
		std::list<Vt> lr;
		//
		typedef carpio::Solver_SOR_<Vt> Sol;
		Sol solver(maxi, tol, this->_omega);
		solver.solve(mat, x, b);
		//int res_code = carpio::Jacobi(mat,          // A  The matrix
		//		x,            // x
		//		b,            // b
		//		maxi,         // max iter
		//		tol,          // Tolerance
		//		lr);          // list residual
		//std::cout << "num iter = " << solver.num_iter();
		//std::cout << "   residual = " << solver.residual() << "\n";

		//std::cout << maxi << "  " << tol << "  " << res_code << "\n";
		this->_copy_to_cs(x);

		if (this->has_flag("OutputResidual")) {
			std::FILE* f = carpio::any_cast<std::FILE*>(
					this->_flags["OutputResidual"]);
			this->output_residual(f, lr);
		}
	}
};

}

#endif
