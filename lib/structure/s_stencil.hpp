#ifndef _S_STENCIL_HPP
#define _S_STENCIL_HPP

#include "s_define.hpp"
#include "s_expression.hpp"
#include <map>
#include <utility>

namespace structure {

template<St DIM>
class Stencil_ {
public:
	typedef Index_<DIM> Index;

	typedef Grid_<DIM> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Index_compare_<DIM> Index_compare;
	//typedef std::unordered_map<Index, Vt, Index_compare> Coerow;
	//typedef std::shared_ptr<Expression> spExpression;

	typedef Expression_<DIM> Expression;
	typedef std::shared_ptr<Expression> spExpression;
protected:
	spGrid _grid;
public:
	Stencil_(spGrid grid) :
			_grid(grid) {
	}

	virtual Vt coe(const Index& index, short dim, short ori) const {
		std::cout << "Stencil coefficient \n";
		return 0.0;
	}
	virtual Vt src(const Index& index, short dim, short ori) const {
		std::cout << "Stencil source \n";
		return 0.0;
	}

	virtual spExpression coe_row(const Index& index) const {
		std::cout << "Stencil ceo row \n";
		return nullptr;
	}

	virtual ~Stencil_() {
	}

	static void Show(const Expression& cr){
		for(auto iter = cr.begin(); iter!=cr.end(); ++iter){
			const Vt& v = iter->second;
			const std::pair<Index, short>& key = iter->first;
			std::cout<<"I:"<< key.first <<" E:" << key.second <<" V:"<<v<<"\n";
		}
	}

};

}

#endif
