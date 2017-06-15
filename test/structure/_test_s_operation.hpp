#ifndef _TEST_S_OPERATION_HPP
#define _TEST_S_OPERATION_HPP

#include "gtest/gtest.h"
#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_io_plotly.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <memory>

namespace structure {
const short DIM = 1;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid_<DIM> > spGrid;
typedef Scalar_<DIM> Scalar;
typedef std::shared_ptr<Scalar_<DIM> > spScalar;
typedef Index_<DIM> Index;
typedef Operation_<DIM> Op;

TEST(operation, set) {
	Poi pmin(1.0, 1.0, 0.0);
	Poi pmax(5.0, 5.0, 3.0);
	Index mn(10, 10, 3);
	spGrid spg(new Grid(pmin, pmax, mn, 3));
	spScalar sps(new Scalar(spg));

	typename Op::Function fun = [](Vt t, Vt x, Vt y, Vt z) {
		Vt cx = 3.0;
		Vt cy = 3.0;
		Vt r  = 0.5;
		return (x - cx) * (x - cx) + (y - cy)*(y - cy) - r*r;
	};

	Op::Set(sps, fun);
	Vt x = 4.3;
	Vt y = 4.1;
	Vt z = 4.6;

	Poi p(x, y, z);

	Vt val = Op::InterpolateCoordinate(sps, p[0], p[1], p[2]);
	std::cout << " val = " << val <<std::endl;
	Vt ext = fun(0.0, p[0], p[1], p[2]);
	std::cout << " ext = " << ext <<std::endl;
	Plotly plt;
	//plt.add(Heatmap(grid, *(poisson.get_CS("phi"))));
	//	auto actor1 = ScalarCenter(*(poisson.get_CS("phi")));
	auto actor1 = ScalarCenter(*sps);
	actor1->set_colorscale_range(0, 1);
	actor1->set_colorscale_name("RdBu");
	plt.add(actor1);
	//auto actor2 = VectorFace(*vfv);
	//actor2->set_colorscale_range(0, 1);
	//actor2->set_colorscale_name("RdBu");
	//plt.add(actor2);
	//plt.add(WireFrame2(grid));
	//plt.add(VectorFace(*spvf));
	//plt.add(VectorCenter(*spvc));
	//plt.add(CenterPoints(grid));
	//plt.add(WireFrame(grid));
	//plt.add(WireFrameGhost(grid, *(poisson.get_ghost())));
	//plt.title("Structure");
	plt.size(800, 800);
	//plt.set_output_file("fig_phi");
	//plt.plot();

}

}
#endif
