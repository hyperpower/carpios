#ifndef _S_CENTER_DATA_HPP_
#define _S_CENTER_DATA_HPP_

#include "s_define.hpp"
#include "s_grid.hpp"

namespace structure {

template<class T, St DIM>
class Data_ {
public:
	static const St Dim = DIM;
	typedef Grid_<Dim> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef T Type;
	typedef T& ref_Type;
	typedef const T& const_ref_Type;
	typedef typename Grid::Poi Poi;
	typedef typename Grid::Index Index;
	typedef typename Grid::Ijk Ijk;
	typedef carpio::SpaceT<T, DIM> Space;
	typedef Data_<T, DIM> Data;
protected:
	Space _data;
	spGrid _grid;
public:
	Data_(spGrid p) {
		_grid = p;
		_init_space();
	}

	Data_(const Data& d) {
		_grid = d._grid;
		_data = d._data;
	}

	spGrid get_grid() const {
		return _grid;
	}

	St num_cells() const {
		return _grid->num_cells();
	}

	// work with grid
	Vt c_(const St& dim, const Index& idx) const {
		return _grid->c_(dim, idx);
	}
	Vt s_(const St& dim, const Index& idx) const {
		return _grid->s_(dim, idx);
	}
	Vt hs_(const St& dim, const Index& idx) const {
		return _grid->hs_(dim, idx);
	}

	const_ref_Type val(const Idx& i, const Idx& j = 0, const Idx& k = 0) const {
		return _data(_IDX(i), _IDX(j), _IDX(k));
	}

	ref_Type val(const Idx& i, const Idx& j = 0, const Idx& k = 0) {
		return _data(_IDX(i), _IDX(j), _IDX(k));
	}

	const_ref_Type VAL(const Idx& i, const Idx& j = 0, const Idx& k = 0) const {
		return _data(i, j, k);
	}

	ref_Type VAL(const Idx& i, const Idx& j = 0, const Idx& k = 0) {
		return _data(i, j, k);
	}

	const_ref_Type operator()(const Idx& i, const Idx& j = 0,
			const Idx& k = 0) const {
		return _data(_IDX(i), _IDX(j), _IDX(k));
	}

	ref_Type operator()(const Idx& i, const Idx& j = 0, const Idx& k = 0) {
		return _data(_IDX(i), _IDX(j), _IDX(k));
	}

	const_ref_Type operator()(const Index& index) const {
		return val(index.value(0), index.value(1), index.value(2));
	}

	ref_Type operator()(const Index& index) {
		return val(index.value(0), index.value(1), index.value(2));
	}

	const_ref_Type val(const Index& index) const {
		return val(index.value(0), index.value(1), index.value(2));
	}

	ref_Type val(const Index& index) {
		return val(index.value(0), index.value(1), index.value(2));
	}

	const_ref_Type VAL(const Index& INDEX) const {
		return VAL(INDEX.value(0), INDEX.value(1), INDEX.value(2));
	}

	ref_Type VAL(const Index& INDEX) {
		return VAL(INDEX.value(0), INDEX.value(1), INDEX.value(2));
	}

	const_ref_Type val(const Ijk& ijk) const {
		return val(ijk.i(), ijk.j(), ijk.k());
	}

	ref_Type val(const Ijk& ijk) {
		return val(ijk.i(), ijk.j(), ijk.k());
	}

	const_ref_Type VAL(const Ijk& IJK) const {
		Index INDEX = IJK.current();
		return VAL(INDEX.value(0), INDEX.value(1), INDEX.value(2));
	}

	ref_Type VAL(const Ijk& IJK) {
		Index INDEX = IJK.current();
		return VAL(INDEX.value(0), INDEX.value(1), INDEX.value(2));
	}

	void assign(const T& t) {
		_data.assign(t);
	}

	inline St _IDX(Idx i) const {
		return _grid->_IDX(i);
	}
	inline Idx _idx(St I) const {
		return _grid->_idx(I);
	}
protected:
	void _init_space() {
		Index ng = _grid->N();
		_data.reconstruct(ng.value(0), ng.value(1), ng.value(2));
	}
};
template<St DIM>
class Scalar_: public Data_<Vt, DIM> {
public:
	static const St Dim = DIM;
	typedef Grid_<Dim> Grid;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Data_<Vt, DIM> Data;
	typedef Scalar_<DIM> Scalar;
	typedef typename Grid::Poi Poi;
	typedef typename Grid::Index Index;
	typedef typename Grid::Ijk Ijk;
	typedef carpio::SpaceT<Vt, DIM> Space;

public:
	Scalar_(spGrid p) :
			Data(p) {
		this->assign(0);
	}



	Scalar_(const Scalar& s) :
			Data(s) {
	}

protected:

};

}

#endif
