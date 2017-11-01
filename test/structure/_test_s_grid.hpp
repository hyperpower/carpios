#ifndef _TEST_S_GRID_HPP
#define _TEST_S_GRID_HPP

#include "gtest/gtest.h"
#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_io_plotly.hpp"
#include "utility/random.h"
#include <math.h>
#include <memory>

namespace structure {
const short DIM = 2;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef Index_<DIM> Index;
typedef Poisson_<DIM> Poisson;

TEST(sgrid, index) {
	Index_<2> i1;
	Index_<2> i2(1, 2);
	i1 = i2;
	Index_<2> i3(i1);
	i3.show();

	Ijk_<2> cur(2, 3, 3, 4);
	cur.show();
	++cur;
	++cur;
	cur.show();
	--cur;
	cur.show();
	std::cout << "\n";
	Ijk_<2> cur2(0, 3, 0, 3);
	cur2.show();
	--cur2;
	--cur2;
	cur2.show();
	std::cout << (cur2 == cur) << "\n";
}

TEST(sgrid, construct) {
	Poi pmin(2.0, 2.0, 0.0);
	Poi pmax(6.0, 10.0, 3.0);
	Index mn(4, 5, 3);
	Grid grid(pmin, pmax, mn, 2);
	std::cout << "Number of Face     :" << Grid::NumFace << "\n";
	std::cout << "Number of Vertices :" << Grid::NumVertex << "\n";

	for (Grid::Ijk idx = grid.begin_ijk(); !idx.is_end(); ++idx) {
		//	idx.show();
		Poi pc = grid.c(idx);
		//	pc.show();
		//	std::cout<<"\n";
	}
	Poi ver = grid.v(3, 0, 0);
	ver.show();

	Plotly plt;
	plt.add(CenterPoints(grid));
	plt.add(WireFrame(grid));
	plt.title("Structure");
	plt.size(1200, 800);
	plt.plot();
}

Arr random_arr(Vt len, Arr::size_type n) {
	Arr arr(n);
	Vt dx = len / double(n);
	for (Arr::size_type i = 0; i < arr.size(); i++) {
		arr[i] = dx;
	}
	for (Arr::size_type i = 0; i < arr.size() - 1; i++) {
		Vt d2 = arr[i] + arr[i + 1];
		Vt dx1 = carpio::Random::nextFloat(0.7 * d2 / 2.0, 1.2 * d2 / 2.0);
		//Vt dx1 = d2 / 2. - 0.5;
		arr[i] = dx1;
		arr[i + 1] = d2 - dx1;
		ASSERT(arr[i] + arr[i + 1] - d2 < 1e-6);
	}
	Vt nlen =0;
	for (Arr::size_type i = 0; i < arr.size(); i++) {
		nlen += arr[i];
	}
	ASSERT(nlen - len < 1e-6);
	return arr;
}

TEST(sgrid, construct_non) {
	Poi pmin(2.0, 2.0, 0.0);
	Poi pmax(6.0, 10.0, 3.0);
	Arr acx = random_arr(10, 4);
	Arr acy = random_arr(7, 6);
	//for (Arr::size_type i = 0; i < acx.size(); i++) {
	//	acx[i] = carpio::Random::nextFloat(0.3, 0.5);
	//}
	//for (Arr::size_type i = 0; i < acy.size(); i++) {
	//	acy[i] = carpio::Random::nextFloat(0.3, 0.8);
	//}

	Grid grid(pmin, 2, acx, acy);
	std::cout << "Number of Face :" << Grid::NumFace << "\n";
	for (Grid::Ijk idx = grid.begin_ijk(); !idx.is_end(); ++idx) {
		//	idx.show();
		Poi pc = grid.c(idx);
		//	pc.show();
		//	std::cout<<"\n";
	}

	Plotly plt;
	plt.add(CenterPoints(grid));
	plt.add(WireFrame(grid));
	//plt.add(WireFrameGhost(gri));
	plt.title("Structure");
	plt.size(1200, 800);
//	plt.plot();
}

}

#endif
