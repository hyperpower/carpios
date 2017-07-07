/*
 * s_mask.hpp
 *
 *  Created on: Jul 7, 2017
 *      Author: zhou
 */

#ifndef _S_MASK_HPP_
#define _S_MASK_HPP_

#include "s_define.hpp"
#include "s_grid.hpp"

namespace structure {

template<St DIM, class ADATA>
class Mask_ {
public:
	static const St Dim = DIM;
	typedef Grid_<Dim> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Data_<int, DIM> Data;
	typedef ADATA Append;
	typedef std::shared_ptr<AData> spAppend;
	typedef std::shared_ptr<Data> spData;
	typedef typename Grid::Poi Poi;
	typedef typename Grid::Index Index;
	typedef typename Grid::Ijk Ijk;
	typedef carpio::SpaceT<Vt, DIM> Space;
public:
	spData   _spmask;
	spAppend _spapp;
public:
	Mask_() :
		 _spmask(nullptr), _spapp(nullptr){
	}







protected:

};
}



#endif /* LIB_STRUCTURE_S_MASK_HPP_ */
