#ifndef _S_EXPRESSION_HPP
#define _S_EXPRESSION_HPP

#include "s_define.hpp"
#include <map>
#include <utility>

namespace structure {

template<St DIM>
class Expression_: public carpio::Polynomial_<Vt, Index_<DIM>, short,
		carpio::IsZero_<Vt>,
		carpio::IsZero_<short>,
		Index_compare_<DIM>,
		std::less<short> > {
public:

	typedef carpio::Polynomial_<Vt, Index_<DIM>, short,
			carpio::IsZero_<Vt>,
			carpio::IsZero_<short>,
			Index_compare_<DIM>,
			std::less<short> > Base;
	typedef Index_<DIM> Index;
	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Index_compare_<DIM> Index_compare;
	//typedef std::unordered_map<Index, Vt, Index_compare> Coerow;
	//typedef std::shared_ptr<Coerow> spCoerow;

	typedef Expression_<DIM> Expression;
	typedef std::shared_ptr<Expression> spExpression;

	typedef typename Base::iterator iterator;
	typedef typename Base::const_iterator const_iterator;

public:
	Expression_(){

	}


};
}
#endif
