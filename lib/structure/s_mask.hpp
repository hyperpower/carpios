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
#include "s_data.hpp"

#include <unordered_map>

namespace structure {

enum CellType {
	_NORMAL, _BOUNDARY, _GHOST, _MPI
};

template<St DIM, class ADATA>
class Mask_ {
public:
	static const St Dim = DIM;
	typedef Grid_<Dim> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef typename Grid::Poi Poi;
	typedef typename Grid::Index Index;
	typedef typename Grid::Ijk Ijk;
	typedef Data_<int, DIM> MaskValue;
	typedef ADATA AData;
	typedef std::unordered_map<Index, AData, Index_hash_<DIM> > Append;
	typedef std::shared_ptr<Append> spAppend;
	typedef std::shared_ptr<MaskValue> spMaskValue;

	typedef carpio::SpaceT<Vt, DIM> Space;
public:
	spMaskValue _spmask;
	spAppend _spapp;

	AData _defalt_app;
public:
	Mask_() :
			_spmask(nullptr), _spapp(nullptr) {
	}

	inline const int& v_mask(const Idx& i,   //
			const Idx& j = 0, //
			const Idx& k = 0) const {
		return _spmask->val(i, j, k);
	}

	inline int& v_mask(const Idx& i,   //
			const Idx& j = 0,  //
			const Idx& k = 0) {
		return _spmask->val(i, j, k);
	}

	inline const AData& v_append( //
			const Idx& i,   //
			const Idx& j = 0, //
			const Idx& k = 0) const {
		if (_spapp == nullptr) {
			return _defalt_app;
		} else {
			Index idx(i, j, k);
			auto it = _spapp->find(idx);
			if (it == _spapp->end()) {
				return _defalt_app;
			} else {
				return it->second;
			}
		}
	}
	inline AData& v_append( //
			const Idx& i,   //
			const Idx& j = 0, //
			const Idx& k = 0){
		if (_spapp == nullptr) {
			return _defalt_app;
		} else {
			Index idx(i, j, k);
			auto it = _spapp->find(idx);
			if (it == _spapp->end()) {
				return _defalt_app;
			} else {
				return it->second;
			}
		}
	}

protected:

};
}

#endif /* LIB_STRUCTURE_S_MASK_HPP_ */
