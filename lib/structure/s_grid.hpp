#ifndef _S_GRID_HPP
#define _S_GRID_HPP

#include "s_define.hpp"
#include <limits>

namespace structure {

template<St DIM>
class Grid_ {
public:
	static const St Dim = DIM;
	typedef Poi_<Dim> Poi;
	typedef Index_<Dim> Index;
	typedef Ijk_<Dim> Ijk;
	static const St NumVertex = DIM == 1 ? 2 : (DIM == 2 ? 4 : 8);
	static const St NumFace = DIM == 1 ? 2 : (DIM == 2 ? 4 : 6);
protected:
	Poi _min, _max;
	Index _n, _ng;  // number of node and number of node and ghost
	Idx _gl;        // ghost layer
	Arr _cs[Dim];   // cell size
	Arr _c[Dim];    // coordinate center
public:
	Grid_(const Poi& min, const Poi& max, const Index& n, const Idx& gl) {
		_min = min;
		_max = max;
		_n = n;
		_gl = gl;
		for (St i = 0; i < Dim; ++i) {
			_ng[i] = _n[i] + 2 * _gl;
		}
		_set_uniform_grid();
	}
	Grid_(const Poi& min, const Idx& gl, const Arr& csx, //
			const Arr& csy = Arr(), //
			const Arr& csz = Arr()) {
		_min = min;
		Index n(csx.size(), csy.size(), csz.size());
		_n = n;
		_gl = gl;
		for (St i = 0; i < Dim; ++i) {
			_ng[i] = _n[i] + 2 * _gl;
		}
		_set_non_uniform_grid(csx, csy, csz);
	}
	// index ==================================
	const Index& n() const {
		return _n;
	}
	const Index& N() const {
		return _ng;
	}
	Idx n(St dim) const {
		return dim < Dim ? _n[dim] : 0;
	}

	Idx N(St dim) const {
		return dim < Dim ? _ng[dim] : 0;
	}

	// size ===================================
	inline Vt s_(St dim, Idx i) const {
		if (dim < Dim) {
			return _cs[dim][_IDX(i)];
		} else {
			return 1;
		}
	}
	inline Vt S_(St dim, Idx I) const {
		if (dim < Dim) {
			return _cs[dim][I];
		} else {
			return 1;
		}
	}
	inline Vt s_(St dim, Index idx) const {
		if (dim < Dim) {
			return _cs[dim][_IDX(idx[dim])];
		} else {
			return 1;
		}
	}
	// half size
	inline Vt hs_(St dim, Idx i) const {
		return (dim < Dim) ? _cs[dim][_IDX(i)] * 0.5 : 0.5;
	}
	inline Vt hs_(St dim, Index idx) const {
		return hs_(dim, idx[dim]);
	}
	inline Vt HS_(St dim, Idx I) const {
		return (dim < Dim) ? _cs[dim][I] * 0.5 : 0.5;
	}
	// volume ==================================
	inline Vt volume(Index idx) const {
		Vt res = 1.0;
		FOR_EACH_DIM
		{
			res *= this->s_(d, idx);
		}
		return res;
	}
	// center ==================================
	Poi c(Idx i, Idx j = 0, Idx k = 0) const {
		Poi res;
		Idx ai[] = { i, j, k };
		for (St d = 0; d < Dim; ++d) {
			res[d] = _c[d][_IDX(ai[d])];
		}
		return std::move(res);
	}
	Poi C(Idx i, Idx j = 0, Idx k = 0) const {
		Poi res;
		Idx ai[] = { i, j, k };
		for (St d = 0; d < Dim; ++d) {
			res[d] = _c[d][ai[d]];
		}
		return std::move(res);
	}
	Vt c_(const St& dim, const Idx& idx) const {
		return (dim < Dim) ? _c[dim][_IDX(idx)] : 0;
	}
	Vt C_(const St& dim, const Idx& IDX) const {
		return (dim < Dim) ? _c[dim][IDX] : 0;
	}
	Vt c_(const St& dim, const Index& idx) const {
		return (dim < Dim) ? _c[dim][_IDX(idx[dim])] : 0;
	}

	Poi c(Ijk ijk) const {
		return c(ijk.i(), ijk.j(), ijk.k());
	}
	Poi C(Ijk IJK) const {
		return C(IJK.i(), IJK.j(), IJK.k());
	}

	Poi c(Index index) const {
		return c(index.i(), index.j(), index.k());
	}
	Poi C(Index INDEX) const {
		return C(INDEX.i(), INDEX.j(), INDEX.k());
	}

	Vt cx(Idx i, Idx j = 0, Idx k = 0) const {
		return _c[0][_IDX(i)];
	}

	Vt cy(Idx i, Idx j = 0, Idx k = 0) const {
		ASSERT(Dim >= 2);
		return _c[1][_IDX(j)];
	}
	Vt cz(Idx i, Idx j = 0, Idx k = 0) const {
		ASSERT(Dim >= 3);
		return _c[2][_IDX(k)];
	}

	St num_cells() const {
		St res = 1;
		for (St d = 0; d < Dim; ++d) {
			res *= _n[d];
		}
		return res;
	}

	St num_vertex() const {
		return NumVertex;
	}
	St num_face() const {
		return NumFace;
	}
	Vt min_size() const {
		Vt mins = std::numeric_limits<Vt>::max();
		for (St d = 0; d < Dim; ++d) {
			mins = std::min(mins, this->_cs[d].min());
		}
		return mins;
	}

	Vt min_size(St dim) const {
		Vt mins = std::numeric_limits<Vt>::max();
		mins = std::min(mins, this->_cs[dim].min());
		return mins;
	}

	// face  ===================================
	Poi f(St dim, bool fb, Idx i, Idx j = 0, Idx k = 0) const {
		Poi pc = c(i, j, k);
		Idx ai[] = { i, j, k };
		Vt halfs = hs_(dim, ai[dim]);
		if (fb) { //right face
			pc[dim] = pc[dim] + halfs;
		} else {  //left face
			pc[dim] = pc[dim] - halfs;
		}
		return pc;
	}

	Vt f_(St dim, int fb, Idx idx) const {
		Vt halfs = hs_(dim, idx);
		Vt cen = c_(dim, idx);
		if (fb == _P_) { //right face
			return cen + halfs;
		} else if (fb == _M_) {  //left face
			return cen - halfs;
		} else {
			SHOULD_NOT_REACH;
			return 0;
		}
	}
	Vt F_(St dim, int fb, Idx IDX) const {
		Vt halfs = HS_(dim, IDX);
		Vt cen = C_(dim, IDX);
		if (fb == _P_) { //right face
			return cen + halfs;
		} else if (fb == _M_) {  //left face
			return cen - halfs;
		} else {
			SHOULD_NOT_REACH;
			return 0;
		}
	}

	// distance to face
	Vt df_(St dim, Idx idx) const {
		Vt halfs = hs_(dim, idx);
		Vt cen = c_(dim, idx);
		return std::abs(halfs - cen);
	}
	// vertex ================================
	Poi v(Idx order, Idx i, Idx j = 0, Idx k = 0) const {
		static const short VERTEX_IDX[][3] = { { _M_, _M_, _M_ }, //
				{ _P_, _M_, _M_ }, //
				{ _M_, _P_, _M_ }, //
				{ _P_, _P_, _M_ }, //
				{ _M_, _M_, _P_ }, //
				{ _P_, _M_, _P_ }, //
				{ _M_, _P_, _P_ }, //
				{ _P_, _P_, _P_ }  //
		};
		return v(i, VERTEX_IDX[order][0], j, VERTEX_IDX[order][1], k,
				VERTEX_IDX[order][2]);
	}
	Poi v(Idx order, Ijk ijk) const {
		return v(order, ijk.i(), ijk.j(), ijk.k());
	}
	Poi v(Idx i, short oi, Idx j = 0, short oj = 0, Idx k = 0,
			short ok = 0) const {
		Poi res;
		Idx ai[] = { i, j, k };
		short ao[] = { oi, oj, ok };
		for (St d = 0; d < Dim; ++d) {
			Vt cen = c_(d, ai[d]);
			Vt halfs = hs_(d, ai[d]);
			if (ao[d] == _P_) {
				res[d] = cen + halfs;
			} else if (ao[d] == _M_) {
				res[d] = cen - halfs;
			} else {
				res[d] = cen;
			}
		}
		return res;
	}

	Ijk begin_ijk() const {
		Idx a[] = { 0, 0, 0 };
		for (St d = 0; d < Dim; ++d) {
			a[d] = _n[d];
		}
		return Ijk(0, a[0], 0, a[1], 0, a[2]);
	}

	Ijk begin_IJK() const {
		Idx a[] = { 0, 0, 0 };
		for (St d = 0; d < Dim; ++d) {
			a[d] = _ng[d];
		}
		return Ijk(0, a[0], 0, a[1], 0, a[2]);
	}

	Ijk end_ijk() const {
		Ijk res = begin_ijk();
		res.to_end();
		return res;
	}

	Ijk end_IJK() const {
		Ijk res = begin_IJK();
		res.to_end();
		return res;
	}

	Ijk to_ijk(const Ijk& IJK) const {
		Index C = IJK.current();
		Index c(C.value(0) - _gl, C.value(1) - _gl, C.value(2) - _gl);
		return Ijk(c, this->_n);
	}

	inline St _IDX(Idx i) const {
		return St(i + _gl);
	}
	inline Idx _idx(St I) const {
		return Idx(I - _gl);
	}

	// Point is in the range
	inline bool is_in_on(Poi p) {
		bool res = true;
		FOR_EACH_DIM
		{
			Vt min = this->_min[d];
			Vt max = this->_max[d];
			Vt val = p[d];
			res = res && carpio::IsInRange(min, val, max, carpio::_cc_);
		}
		return res;
	}
	//  find closest index on the negative direction
	//  for example:
	//
	//      1.2   1.4   1.6   1.8   2.0
	//  -----+-----+-----+-----+-----+
	//       0     1     2     3     4
	//                      ^
	//               cor = 1.5
	//  return 2
	inline Idx find_close_idx_m(St dim, Vt cor) {
		ASSERT(dim < Dim);
		St res = 0;
		// this->_c is a order array
		for (St i = 0; i < this->_c[dim].size(); i++) {
			if (this->_c[dim][i] > cor) {
				res = i;
				break;
			}
		}
		res = (res == 0) ? 0 : res - 1;
		return _idx(res);
	}
	//  find closest index on the positive direction
	//  for example:
	//
	//      1.2   1.4   1.6   1.8   2.0
	//  -----+-----+-----+-----+-----+
	//       0     1     2     3     4
	//                      ^
	//               cor = 1.5
	//  return 3
	inline Idx find_close_idx_p(St dim, Vt cor) {
		ASSERT(dim < Dim);
		St res = 0;
		// this->_c is a order array
		for (St i = 0; i < this->_c[dim].size(); i++) {
			if (this->_c[dim][i] > cor) {
				res = i;
				break;
			}
		}
		return _idx(res);
	}

protected:
	void _set_uniform_grid() {
		for (St d = 0; d < Dim; ++d) {
			Vt l = _max[d] - _min[d];
			ASSERT(l > 0);
			Vt ds = l / _n[d];
			_cs[d].reconstruct(_ng[d]);
			_cs[d].assign(ds);
			_c[d].reconstruct(_ng[d]);
			Vt csum = 0;
			for (typename Arr::size_type ii = 0; ii < _c[d].size(); ++ii) {
				_c[d][ii] = csum + 0.5 * _cs[d][ii];
				csum += _cs[d][ii];
			}
			// shift
			Vt llg = 0;
			for (typename Arr::size_type ii = 0; ii < _gl; ++ii) {
				llg += _cs[d][ii];
			}
			Vt sh = _min[d] - llg;
			for (typename Arr::size_type i = 0; i < _c[d].size(); ++i) {
				_c[d][i] += sh;
			}
		}
	}

	void _set_non_uniform_grid(const Arr& csx, //
			const Arr& csy, //
			const Arr& csz) {
		const Arr* csn[3] = { &csx, &csy, &csz };
		for (St d = 0; d < Dim; ++d) {
			const Arr& acs = *(csn[d]);
			Arr::size_type acsn = acs.size();
			ASSERT(acsn > 1);
			_cs[d].reconstruct(_ng[d]);
			Vt len = 0;
			for (typename Arr::size_type ii = 0; ii < _cs[d].size(); ++ii) {
				if (ii < _gl) {
					_cs[d][ii] = acs[0];
				} else if (ii - _gl >= acsn) {
					_cs[d][ii] = acs[acsn - 1];
				} else {
					_cs[d][ii] = acs[ii - _gl];
					len += acs[ii - _gl];
				}
			}
			_c[d].reconstruct(_ng[d]);
			Vt csum = 0;
			for (typename Arr::size_type ii = 0; ii < _c[d].size(); ++ii) {
				_c[d][ii] = csum + 0.5 * _cs[d][ii];
				csum += _cs[d][ii];
			}
			// shift
			Vt llg = 0;
			for (typename Arr::size_type ii = 0; ii < _gl; ++ii) {
				llg += _cs[d][ii];
			}
			Vt sh = _min[d] - llg;
			for (typename Arr::size_type i = 0; i < _c[d].size(); ++i) {
				_c[d][i] += sh;
			}
			_max[d] = _min[d] + len;
		}
	}

};

}

#endif
