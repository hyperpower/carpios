#ifndef _EXPRESSION_H_
#define _EXPRESSION_H_

#include "calculation_define.hpp"
#include <math.h>
#include <fstream>
namespace carpio {
/*
 * the Expression class
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class Expression_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Expression_<COO_VALUE, VALUE, DIM> Self;
	typedef Expression_<COO_VALUE, VALUE, DIM> pSelf;
	typedef const Expression_<COO_VALUE, VALUE, DIM> const_Self;
	typedef Expression_<COO_VALUE, VALUE, DIM>& ref_Self;
	typedef const Expression_<COO_VALUE, VALUE, DIM>& const_ref_Self;

	typedef Domain_<COO_VALUE, VALUE, DIM> Domain;
	typedef Domain_<COO_VALUE, VALUE, DIM>& ref_Domain;
	typedef const Domain_<COO_VALUE, VALUE, DIM>& const_ref_Domain;
	typedef Domain_<COO_VALUE, VALUE, DIM>* pDomain;
	typedef const Domain_<COO_VALUE, VALUE, DIM>* const_pDomain;

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

	struct pNode_less: public std::binary_function<const_pNode, const_pNode,
			bool> {
		bool operator()(const_pNode __x, const_pNode __y) const {
			return __x->d_idx() < __y->d_idx();
		}
	};
	typedef Polynomial_<vt, const_pNode, int, IsZero_<vt>, IsZero_<int>,
			pNode_less, std::less<int> > Exp;
	typedef Exp* pExp;
	typedef const Exp& const_ref_Exp;
	typedef typename Exp::Term Term;
	//typedef typename Exp::Value Value;
	//typedef typename Exp::Key Key;
	typedef typename Exp::iterator iterator;
	typedef typename Exp::const_iterator const_iterator;

protected:
	Exp _exp;
public:
	Expression_() :
			_exp() {
	}
	Expression_(const_ref_Self rhs) :
			_exp(rhs._exp) {
	}
	ref_Self operator=(const_ref_Self& rhs) {
		//std::cout<<"operator="<<std::endl;
		this->_exp = rhs._exp;
		return *this;
	}
	/*
	 *  iterator
	 */
	iterator begin() {
		return _exp.begin();
	}
	const_iterator begin() const {
		return _exp.begin();
	}
	iterator end() {
		return _exp.end();
	}
	const_iterator end() const {
		return _exp.end();
	}
	/*
	 * capacity
	 */
	bool empty() const {
		return _exp.empty();
	}
	St size() const {
		return _exp.size();
	}
	//Value& operator[](const Key& k) {
	//	return _exp[k];
	//}
	//const Value& operator[](const Key& k) const {
	//	return _exp[k];
	//}

	static bool IsConstant(const const_iterator& iter) {
		return Exp::IsConstant(iter);
	}
	static bool IsConstant(const iterator& iter) {
		return Exp::IsConstant(iter);
	}

	static bool IsZeroCoe(const iterator& iter) {
		return Exp::IsZeroCoe(iter);
	}

	static const_pNode GetpNode(const_iterator iter) {
		return iter->first.first;
	}

	static vt GetCoe(const iterator& iter) {
		return iter->second;
	}
	static vt GetCoe(const_iterator& iter) {
		return iter->second;
	}

	static int GetExp(const iterator& iter) {
		return iter->first.second;
	}

	static int GetExp(const_iterator& iter) {
		return iter->first.second;
	}

	static bool IsGhostNode(const iterator& iter) {
		const_pNode pn = iter->first.first;
		ASSERT(pn != nullptr);
		return (pn->get_type() == _Ghost_);
	}

protected:
	typedef std::pair<iterator, bool> _ret;
public:

	_ret insert(const Term &x) {
		return _exp.insert(x);
	}

	_ret insert(vt coe, const_pNode pn, int exp) {
		Term t(coe, pn, exp);
		return insert(t);
	}

	iterator find(const_pNode pn, int exp) {
		Term t(1.0, pn, exp);
		return _exp.find(t);
	}

	const_iterator find(const_pNode pn, int exp) const {
		Term t(1.0, pn, exp);
		return _exp.find(t);
	}

	void erase(iterator& iter) {
		_exp.erase(iter);
	}

	void clear() {
		_exp.clear();
	}

	void plus(const_ref_Self exp) {
		_exp.plus(exp._exp);
	}

	void plus(const vt& coe, const_pNode pn, const int& exp) {
		this->insert(coe, pn, exp);
	}

	void minus(const_ref_Self exp) {
		_exp.minus(exp._exp);
	}
	void times(const vt& rhs) {   //overload operator*
		_exp.times(rhs);
	}
	void divide(const vt& rhs) {  //overload operator/
		_exp.divide(rhs);
	}

	void trim_zero() {
		_exp.trim_zero();
	}

	void merge_const() {
		_exp.merge_const();
	}

	void concise() {
		_exp.concise();
	}

	void show() const {
		// all the valuable should overload <<
		std::cout << "  size = " << _exp.size() << "\n";
		for (const_iterator iter = _exp.begin(); iter != _exp.end(); ++iter) {
			std::cout.flags(std::ios::right);
			std::cout.width(9);
			std::cout << iter->second << " ";
			const_pNode pn = iter->first.first;
			std::cout.width(5);
			std::cout << pn->d_idx() << " ";
			std::cout.width(3);
			std::cout << iter->first.second << " ";
			if (Self::IsConstant(iter)) {
				std::cout << "CONST";
			}
			std::cout << "\n";
		}
		//_exp.show();
	}
	void show(std::fstream& fs) const {
		std::ios oldState(nullptr);
		oldState.copyfmt(fs);
		// all the valuable should overload <<
		fs << "  size = " << _exp.size() << "\n";
		for (const_iterator iter = _exp.begin(); iter != _exp.end(); ++iter) {
			fs.flags(std::ios::right);
			fs.width(9);
			fs << iter->second << " ";
			const_pNode pn = iter->first.first;
			fs.width(5);
			fs << pn->d_idx() << " ";
			fs.width(3);
			fs << iter->first.second << " ";
			if (Self::IsConstant(iter)) {
				fs << "CONST";
			}
			fs << "\n";
		}
		fs.copyfmt(oldState);
		//_exp.show();
	}
	void show_substitute(St idx) const {
		std::ios oldState(nullptr);
		oldState.copyfmt(std::cout);
		// all the valuable should overload <<
		std::cout << " size = " << _exp.size() << "\n";
		for (const_iterator iter = _exp.begin(); iter != _exp.end(); ++iter) {
			std::cout.flags(std::ios::right);
			std::cout.width(11);
	std::cout.precision(3);
			std::cout << std::scientific << iter->second << " ";
			const_pNode pn = iter->first.first;
			std::cout.width(7);
			std::cout << pn->d_idx() << " ";
			std::cout.width(5);
			std::cout << iter->first.second << " ";
			std::cout.width(7);
			std::cout << pn->cd(idx) << " ";
			if (Self::IsConstant(iter)) {
				std::cout << "CONST";
			}
			std::cout << "\n";
		}
		std::cout.copyfmt(oldState);
		//_exp.show();
	}

	vt substitute(St idx) const {
		// get the value of the expression
		// idx is the center data index on pnode
		vt res = 0;
		for (const_iterator iter = this->begin(); iter != this->end(); ++iter) {
			const_pNode pn = GetpNode(iter);
			vt coe = GetCoe(iter);
			int exp = GetExp(iter);
			if (exp == 0) {
				res += coe;
			} else if (exp == 1) {
				vt v = pn->cdva(idx);
				res += coe * v;
			} else {
				vt v = pn->cdva(idx);
				res += coe * pow(v, exp);
			}
		}
		return res;
	}
};



}

#endif
