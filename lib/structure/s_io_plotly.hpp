#ifndef _S_IO_PLOTLY_HPP
#define _S_IO_PLOTLY_HPP

#include "s_grid.hpp"
#include "s_define.hpp"
#include "s_vector.hpp"

#include "io/plotly.h"

namespace structure {
typedef carpio::Plotly Plotly;
typedef std::shared_ptr<carpio::Plotly_actor> spPA;
typedef std::shared_ptr<carpio::Plotly_actor_scatter> spPA_scatter;
typedef std::shared_ptr<carpio::Plotly_actor_scatter3d> spPA_scatter3d;
typedef std::shared_ptr<carpio::Plotly_actor_mesh3d> spPA_mesh3d;
typedef std::shared_ptr<carpio::Plotly_actor_heatmap> spPA_heatmap;

typedef std::list<Vt> Listd;

template<St DIM>
spPA_scatter3d WireFrame(const Grid_<DIM>& grid) {
	typedef Grid_<DIM> Grid;
	Listd lx;
	Listd ly;
	Listd lz;
	short order[] = { 0, 1, 3, 2, 6, 4, 5, 7 };
	short order3[] = { 0, 1, 1, 3, 3, 2, 2, 0, 2, 6, 3, 7, 1, 5, 0, 4, 4, 5, 5,
			7, 7, 6, 6, 4 };
	for (typename Grid::Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
		if (DIM == 2) {
			for (short i = 0; i < grid.num_vertex(); ++i) {
				typename Grid::Poi p = grid.v(order[i], ijk);
				lx.push_back(p.value(0));
				ly.push_back(p.value(1));
				lz.push_back(p.value(2));
			}
			typename Grid::Poi p = grid.v(0, ijk);
			lx.push_back(p.value(0));
			ly.push_back(p.value(1));
			lz.push_back(p.value(2));
		}
		if (DIM == 3) {
			for (short i = 0; i < 23; i += 2) {
				typename Grid::Poi p = grid.v(order3[i], ijk);
				lx.push_back(p.value(0));
				ly.push_back(p.value(1));
				lz.push_back(p.value(2));
				typename Grid::Poi pp = grid.v(order3[i + 1], ijk);
				lx.push_back(pp.value(0));
				ly.push_back(pp.value(1));
				lz.push_back(pp.value(2));
			}
		}
	}
	spPA_scatter3d res;
	if (DIM == 2) {
		res = spPA_scatter3d(new carpio::Plotly_actor_scatter3d(lx, ly, lz, 5));
	} else if (DIM == 3) {
		res = spPA_scatter3d(new carpio::Plotly_actor_scatter3d(lx, ly, lz, 2));
	}
	res->set_mode("lines");
	return res;
}

template<St DIM>
spPA_scatter WireFrame2(const Grid_<DIM>& grid, St dim = 2, Idx idx = 0) {
	typedef Grid_<DIM> Grid;
	Listd lx;
	Listd ly;
	Listd lz;
	short order[] = { 0, 1, 3, 2, 6, 4, 5, 7 };
	short order3[] = { 0, 1, 1, 3, 3, 2, 2, 0, 2, 6, 3, 7, 1, 5, 0, 4, 4, 5, 5,
			7, 7, 6, 6, 4 };
	for (typename Grid::Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
		if (DIM == 3) {
			if (ijk.current().value(dim) == idx) {
				for (short i = 0; i < grid.num_vertex(); ++i) {
					typename Grid::Poi p = grid.v(order[i], ijk);
					lx.push_back(p.value(0));
					ly.push_back(p.value(1));
					lz.push_back(p.value(2));
				}
				typename Grid::Poi p = grid.v(0, ijk);
				lx.push_back(p.value(0));
				ly.push_back(p.value(1));
				lz.push_back(p.value(2));
			}
		} else {
			for (short i = 0; i < grid.num_vertex(); ++i) {
				typename Grid::Poi p = grid.v(order[i], ijk);
				lx.push_back(p.value(0));
				ly.push_back(p.value(1));
				lz.push_back(p.value(2));
			}
			typename Grid::Poi p = grid.v(0, ijk);
			lx.push_back(p.value(0));
			ly.push_back(p.value(1));
			lz.push_back(p.value(2));
		}
	}
	spPA_scatter res;
	res = spPA_scatter(new carpio::Plotly_actor_scatter(lx, ly, 5));
	res->set_mode("lines");
	return res;
}

template<St DIM>
spPA_scatter3d WireFrameGhost(const Grid_<DIM>& grid,
		const Ghost_<DIM>& ghost) {
	typedef Grid_<DIM> Grid;
	Listd lx;
	Listd ly;
	Listd lz;
	short order[] = { 0, 1, 3, 2, 6, 4, 5, 7 };
	short order3[] = { 0, 1, 1, 3, 3, 2, 2, 0, 2, 6, 3, 7, 1, 5, 0, 4, 4, 5, 5,
			7, 7, 6, 6, 4 };
	for (typename Grid::Ijk IJK = grid.begin_IJK(); !IJK.is_end(); ++IJK) {
		typename Grid::Ijk ijk = grid.to_ijk(IJK);
		if (!ghost.is_ghost(ijk.current())) {
			continue;
		}
		if (DIM == 2) {
			for (short i = 0; i < grid.num_vertex(); ++i) {
				typename Grid::Poi p = grid.v(order[i], ijk);
				lx.push_back(p.value(0));
				ly.push_back(p.value(1));
				lz.push_back(p.value(2));
			}
			typename Grid::Poi p = grid.v(0, ijk);
			lx.push_back(p.value(0));
			ly.push_back(p.value(1));
			lz.push_back(p.value(2));
		}
		if (DIM == 3) {
			for (short i = 0; i < 23; i += 2) {
				typename Grid::Poi p = grid.v(order3[i], ijk);
				lx.push_back(p.value(0));
				ly.push_back(p.value(1));
				lz.push_back(p.value(2));
				typename Grid::Poi pp = grid.v(order3[i + 1], ijk);
				lx.push_back(pp.value(0));
				ly.push_back(pp.value(1));
				lz.push_back(pp.value(2));
			}
		}
	}
	spPA_scatter3d res;
	if (DIM == 2) {
		res = spPA_scatter3d(new carpio::Plotly_actor_scatter3d(lx, ly, lz, 5));
	} else if (DIM == 3) {
		res = spPA_scatter3d(new carpio::Plotly_actor_scatter3d(lx, ly, lz, 2));
	}
	res->set_mode("lines");
	return res;
}

template<St DIM>
spPA_scatter3d CenterPoints(const Grid_<DIM>& grid) {
	typedef Grid_<DIM> Grid;
	Listd lx;
	Listd ly;
	Listd lz;
	for (typename Grid::Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
		typename Grid::Poi p = grid.c(ijk);
		lx.push_back(p.value(0));
		ly.push_back(p.value(1));
		lz.push_back(p.value(2));

	}
	spPA_scatter3d res;
	res = spPA_scatter3d(new carpio::Plotly_actor_scatter3d(lx, ly, lz, 1));
	res->set_mode("points");
	return res;
}

template<St DIM>
spPA_scatter ScalarCenter(const Scalar_<DIM>& v, St dim = 2, Idx i = 0) {
	typedef Grid_<DIM> Grid;
	typedef Scalar_<DIM> Scalar;
	Listd lx, ly, lz, lv;
	const Grid_<DIM>& g = *(v.get_grid());
	for (typename Grid::Ijk ijk = g.begin_ijk(); !ijk.is_end(); ++ijk) {
		typename Grid::Poi p = g.c(ijk);
		lx.push_back(p.value(0));
		ly.push_back(p.value(1));
		lz.push_back(p.value(2));
		lv.push_back(v.val(ijk));
	}
	spPA_scatter res;
	if (DIM == 3) {
		Listd* al[] = { &lx, &ly, &lz };
		St dim1, dim2;
		Normal(dim, dim1, dim2);
		res = spPA_scatter(
				new carpio::Plotly_actor_scatter(*(al[dim1]), *(al[dim2])));
		res->set_colorscale(lv);
		res->set_add_val(lv);
		res->set_mode("markers");
	} else if (DIM == 2) {
		res = spPA_scatter(new carpio::Plotly_actor_scatter(lx, ly));
		res->set_colorscale(lv);
		res->set_add_val(lv);
		res->set_mode("markers");
	} else { // DIM == 1
		res = spPA_scatter(new carpio::Plotly_actor_scatter(lx, lv));
		res->set_mode("lines + markers");
	}
	return res;

}

template<St DIM>
spPA_scatter SCALARCENTER(const Scalar_<DIM>& v, St dim = 2, Idx i = 0) {
	typedef Grid_<DIM> Grid;
	typedef Scalar_<DIM> Scalar;
	Listd lx;
	Listd ly;
	Listd lz;
	Listd lv;
	const Grid_<DIM>& g = *(v.get_grid());
	for (typename Grid::Ijk ijk = g.begin_IJK(); !ijk.is_end(); ++ijk) {
		typename Grid::Poi p = g.C(ijk);
		lx.push_back(p.value(0));
		ly.push_back(p.value(1));
		lz.push_back(p.value(2));
		lv.push_back(v.VAL(ijk.current()));
	}
	spPA_scatter res;
	if (DIM == 3) {
		Listd* al[] = { &lx, &ly, &lz };
		St dim1, dim2;
		Normal(dim, dim1, dim2);
		res = spPA_scatter(
				new carpio::Plotly_actor_scatter(*(al[dim1]), *(al[dim2])));
	} else if (DIM == 2) {
		res = spPA_scatter(new carpio::Plotly_actor_scatter(lx, ly));
	} else { // DIM == 1
		res = spPA_scatter(new carpio::Plotly_actor_scatter(lx, lv));
		res->set_mode("lines + markers");
		return res;
	}
	res->set_colorscale(lv);
	res->set_add_val(lv);
	res->set_mode("markers");
	return res;
}

template<St DIM>
spPA_scatter3d ScalarCenter3d(const Scalar_<DIM>& v, St dim = 2, Idx i = 0) {
	typedef Grid_<DIM> Grid;
	typedef Scalar_<DIM> Scalar;
	Listd lx;
	Listd ly;
	Listd lz;
	Listd lv;
	const Grid_<DIM>& g = *(v.get_grid());
	for (typename Grid::Ijk ijk = g.begin_ijk(); !ijk.is_end(); ++ijk) {
		typename Grid::Poi p = g.c(ijk);
		lx.push_back(p.value(0));
		ly.push_back(p.value(1));
		lz.push_back(p.value(2));
		lv.push_back(v.val(ijk));
	}
	spPA_scatter3d res;
	if (DIM == 3) {
		Listd* al[] = { &lx, &ly, &lz };
		St dim1, dim2;
		Normal(dim, dim1, dim2);
		res = spPA_scatter3d(
				new carpio::Plotly_actor_scatter3d(*(al[dim1]), *(al[dim2]),
						lv));
	} else {
		res = spPA_scatter3d(new carpio::Plotly_actor_scatter3d(lx, ly, lv));
	}
	res->set_colorscale(lv, 7);
	//res->set_add_val(lv);
	res->set_mode("markers");
	return res;
}

template<St DIM>
spPA_scatter VectorFace(const VectorFace_<DIM>& vf, St dim = 2, Idx i = 0) {
	typedef Grid_<DIM> Grid;
	typedef Scalar_<DIM> Scalar;
	Listd lx;
	Listd ly;
	Listd lz;
	Listd lv;
	const Grid_<DIM>& g = *(vf.get_grid());
	for (typename Grid::Ijk ijk = g.begin_ijk(); !ijk.is_end(); ++ijk) {
		for (St d = 0; d < DIM; ++d) {
			Vt hd[] = { 0, 0, 0 };
			hd[d] = g.hs_(d, ijk[d]);
			typename Grid::Poi p = g.c(ijk);
			lx.push_back(p.value(0) + hd[0]);
			ly.push_back(p.value(1) + hd[1]);
			lz.push_back(p.value(2) + hd[2]);
			lv.push_back(vf[d].val(ijk));
			//std::cout<<" v = "<<vf[d].val(ijk) << "\n";
			if (ijk[d] == 0) {
				lx.push_back(p.value(0) - hd[0]);
				ly.push_back(p.value(1) - hd[1]);
				lz.push_back(p.value(2) - hd[2]);
				lv.push_back(vf[d].val(ijk.current().m(d)));
			}
		}
	}
	spPA_scatter res;
	if (DIM == 3) {
		Listd* al[] = { &lx, &ly, &lz };
		St dim1, dim2;
		Normal(dim, dim1, dim2);
		res = spPA_scatter(
				new carpio::Plotly_actor_scatter(*(al[dim1]), *(al[dim2])));
	} else {
		res = spPA_scatter(new carpio::Plotly_actor_scatter(lx, ly));
	}
	//std::cout << " len = " << lv.size() << "\n";
	res->set_colorscale(lv, 11);
	res->set_add_val(lv);
	res->set_mode("markers");
	return res;
}

template<St DIM>
spPA_scatter VectorCenter(const VectorCenter_<DIM>& vf, St dim = 2, Idx i = 0) {
	typedef Grid_<DIM> Grid;
	typedef Scalar_<DIM> Scalar;
	Listd lx;
	Listd ly;
	Listd lz;
	Listd lv;
	const Grid_<DIM>& g = *(vf.get_grid());
	// find max velocity
	Vt veoabsmax = SMALL;
	for (typename Grid::Ijk ijk = g.begin_ijk(); !ijk.is_end(); ++ijk) {
		for (St d = 0; d < DIM; ++d) {
			if (std::abs(vf[d].val(ijk)) > veoabsmax) {
				veoabsmax = vf[d].val(ijk);
			}
		}
	}
	for (typename Grid::Ijk ijk = g.begin_ijk(); !ijk.is_end(); ++ijk) {
		for (St d = 0; d < DIM; ++d) {
			Vt hd[] = { 0, 0, 0 };
			hd[d] = g.hs_(d, ijk[d]) * 0.5 * (vf[d].val(ijk) / veoabsmax);
			typename Grid::Poi p = g.c(ijk);
			lx.push_back(p.value(0) + hd[0]);
			ly.push_back(p.value(1) + hd[1]);
			lz.push_back(p.value(2) + hd[2]);
			lv.push_back(vf[d].val(ijk));
			//std::cout<<" v = "<<vf[d].val(ijk) << "\n";
		}
	}
	spPA_scatter res;
	if (DIM == 3) {
		Listd* al[] = { &lx, &ly, &lz };
		St dim1, dim2;
		Normal(dim, dim1, dim2);
		res = spPA_scatter(
				new carpio::Plotly_actor_scatter(*(al[dim1]), *(al[dim2])));
	} else {
		res = spPA_scatter(new carpio::Plotly_actor_scatter(lx, ly));
	}
	//std::cout << " len = " << lv.size() << "\n";
	res->set_colorscale(lv, 11);
	res->set_add_val(lv);
	res->set_mode("markers");
	return res;
}

template<St DIM>
spPA_scatter VECTORCENTER(const VectorCenter_<DIM>& vf, St dim = 2, Idx i = 0) {
	typedef Grid_<DIM> Grid;
	typedef Scalar_<DIM> Scalar;
	Listd lx, ly, lz, lv;
	const Grid_<DIM>& g = *(vf.get_grid());
	for (typename Grid::Ijk IJK = g.begin_IJK(); !IJK.is_end(); ++IJK) {
		for (St d = 0; d < DIM; ++d) {
			Vt hd[] = { 0, 0, 0 };
			hd[d] = g.HS_(d, IJK[d]) * 0.5;
			typename Grid::Poi p = g.C(IJK);
			lx.push_back(p.value(0) + hd[0]);
			ly.push_back(p.value(1) + hd[1]);
			lz.push_back(p.value(2) + hd[2]);
			lv.push_back(vf[d].VAL(IJK.current()));
			//std::cout<<" v = "<<vf[d].val(IJK) << "\n";
		}
	}
	spPA_scatter res;
	if (DIM == 3) {
		Listd* al[] = { &lx, &ly, &lz };
		St dim1, dim2;
		Normal(dim, dim1, dim2);
		res = spPA_scatter(
				new carpio::Plotly_actor_scatter(*(al[dim1]), *(al[dim2])));
	} else {
		res = spPA_scatter(new carpio::Plotly_actor_scatter(lx, ly));
	}
	//std::cout << " len = " << lv.size() << "\n";
	res->set_colorscale(lv, 11);
	res->set_add_val(lv);
	res->set_mode("markers");
	return res;
}

template<St DIM>
spPA_heatmap Heatmap( //
		const Grid_<DIM>& grid, //
		const Scalar_<DIM>& csfield, //
		St dim1 = _X_, St dim2 = _Y_, Idx idx = 0) {
	typedef Grid_<DIM> Grid;
	typedef Index_<DIM> Index;
	Listd lx;
	Listd ly;
	Listd lz;
	lx.push_back(grid.f_(dim1, _M_, 0));
	for (Idx i = 0; i < grid.n(dim1); i++) {
		Vt f = grid.f_(dim1, _P_, i);
		lx.push_back(f);
	}
	ly.push_back(grid.f_(dim2, _M_, 0));
	for (Idx j = 0; j < grid.n(dim2); j++) {
		Vt f = grid.f_(dim2, _P_, j);
		ly.push_back(f);
	}

	for (Idx j = 0; j < grid.n(dim2); j++) {
		for (Idx i = 0; i < grid.n(dim1); i++) {
			Index index;
			index.set(dim1, i);
			index.set(dim2, j);
			index.set(Normal(dim1, dim2), idx);
			//
			lz.push_back(csfield.val(index));
		}
	}

	spPA_heatmap res;
	res = spPA_heatmap(new carpio::Plotly_actor_heatmap(lx, ly, lz));
	return res;
}

template<St DIM>
spPA_heatmap HEATMAP( //
		const Scalar_<DIM>& csfield, //
		St dim1 = _X_, St dim2 = _Y_, Idx idx = 0) {
	typedef Grid_<DIM> Grid;
	typedef Index_<DIM> Index;
	Listd lx;
	Listd ly;
	Listd lz;
	Grid_<DIM>& grid = *(csfield.get_grid());
	lx.push_back(grid.F_(dim1, _M_, 0));
	for (Idx i = 0; i < grid.N(dim1); i++) {

		Vt f = grid.F_(dim1, _P_, i);
		lx.push_back(f);
	}
	ly.push_back(grid.F_(dim2, _M_, 0));
	for (Idx j = 0; j < grid.N(dim2); j++) {

		Vt f = grid.F_(dim2, _P_, j);
		ly.push_back(f);
		//std::cout<< j <<" , "<< f <<std::endl;
	}
	std::cout << lx.size() << std::endl;
	std::cout << ly.size() << std::endl;

	for (Idx j = 0; j < grid.N(dim2); j++) {
		for (Idx i = 0; i < grid.N(dim1); i++) {
			Index index;

			index.set(dim1, i);
			index.set(dim2, j);
			index.set(Normal(dim1, dim2), idx);
			//

			lz.push_back(csfield.VAL(index));

		}
	}

	std::cout << lz.size() << std::endl;
	spPA_heatmap res;
	res = spPA_heatmap(new carpio::Plotly_actor_heatmap(lx, ly, lz));
	return res;
}

}

#endif
