#ifndef _S_VECTOR_DATA_HPP_
#define _S_VECTOR_DATA_HPP_

#include <structure/s_data.hpp>
#include "s_define.hpp"
#include "s_grid.hpp"

namespace structure {

template<St DIM>
class Vector_ {
public:
	static const St Dim = DIM;
	typedef Grid_<Dim> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef std::shared_ptr<const Grid> spcGrid;

	typedef typename Grid::Poi Poi;
	typedef typename Grid::Index Index;
	typedef typename Grid::Ijk Ijk;
	typedef Scalar_<DIM> Data;
	typedef Scalar_<DIM>& ref_Data;
	typedef const Scalar_<DIM>& const_ref_Data;
typedef std::shared_ptr<Data> spData;
typedef std::shared_ptr<const Data> spcData;

protected:
	typedef std::array<spData, DIM> ArrspD;
	typedef typename ArrspD::size_type _st;
	ArrspD _arr;

public:
	Vector_(spData spd1, spData spd2 = nullptr, spData spd3 = nullptr) {
		spData a[] = { spd1, spd2, spd3 };
		FOR_EACH_DIM
		{
			ASSERT(a[d] != nullptr);
			_arr[d] = a[d];
		}
	}

	spGrid get_grid() {
		return _arr[0]->get_grid();
	}
	spcGrid get_grid() const {
		return _arr[0]->get_grid();
	}

	ref_Data operator[](const _st& dim) {
		return *(this->_arr[dim]);
	}

	const_ref_Data operator[](const _st& dim) const {
		return *(this->_arr[dim]);
	}

	spData get_spdata(const _st& dim) {
		return this->_arr[dim];
	}

	spcData get_spdata(const _st& dim) const {
		return this->_arr[dim];
	}

};

template<St DIM>
class VectorCenter_: public Vector_<DIM> {
public:
	static const St Dim = DIM;
	typedef Vector_<DIM> Base;
	typedef Grid_<Dim> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef typename Grid::Poi Poi;
	typedef typename Grid::Index Index;
	typedef typename Grid::Ijk Ijk;
	typedef Scalar_<DIM> Data;
	typedef Scalar_<DIM>& ref_Data;
	typedef const Scalar_<DIM>& const_ref_Data;
	typedef std::shared_ptr<Data> spData;

public:

	VectorCenter_(spData spd1, spData spd2 = nullptr, spData spd3 = nullptr) :
			Base(spd1, spd2, spd3) {

	}
	// max( sqrt(u*u + v*v))
	Vt max_magnitude() const {
		Vt max = 0;
		spGrid grid = this->_arr[0]->get_grid();
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			Vt sum = 0;
			for (St d = 0; d < DIM; ++d) {
				Index idx = ijk.current();
				sum += this->_arr[d]->val(idx) * this->_arr[d]->val(idx);
			}
			Vt sqrtsum = std::sqrt(sum);
			max = std::max(max, sqrtsum);
		}
		return max;

	}

};

template<St DIM>
class VectorFace_: public Vector_<DIM> {
public:
	static const St Dim = DIM;
	typedef Vector_<DIM> Base;
	typedef Grid_<Dim> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef typename Grid::Poi Poi;
	typedef typename Grid::Index Index;
	typedef typename Grid::Ijk Ijk;
	typedef Scalar_<DIM> Data;
	typedef Scalar_<DIM>& ref_Data;
	typedef const Scalar_<DIM>& const_ref_Data;
	typedef std::shared_ptr<Data> spData;

public:
	VectorFace_(spData spd1, spData spd2 = nullptr, spData spd3 = nullptr) :
			Base(spd1, spd2, spd3) {

	}

};

}

#endif
