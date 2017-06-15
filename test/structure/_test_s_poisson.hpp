#ifndef _TEST_S_POISSON_HPP
#define _TEST_S_POISSON_HPP

#include "gtest/gtest.h"
#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_io_plotly.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <memory>

namespace structure {
const short DIM = 2;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid> spGrid;
typedef Index_<DIM> Index;
typedef Poisson_<DIM> Poisson;

TEST(DISABLED_sgrid, poisson) {
	Poi pmin(1.0, 1.0, 0.0);
	Poi pmax(5.0, 5.0, 3.0);
	Index mn(10, 10, 3);
	Grid grid(pmin, pmax, mn, 3);
	std::cout << "Number of Face :" << Grid::NumFace << "\n";
	Poi ver = grid.v(3, 0, 0);
	Poisson poisson(std::make_shared<Grid>(grid));

	Poisson::spBC bc1(new Poisson::BC()), bc2(new Poisson::BC());
	bc1->set_default_1_bc(0);
	bc2->set_default_1_bc(1);
	poisson.add_bc(0, 0, "phi", bc1);
	poisson.add_bc(0, 1, "phi", bc1);
	poisson.add_bc(0, 2, "phi", bc2);
	poisson.add_bc(0, 3, "phi", bc1);

	typename Poisson::Function fun = [](Vt t, Vt x, Vt y, Vt z) {
		return x < 2.0 ? 3: 10;
	};
	poisson.set_beta(fun);
	poisson.set_solver("IC_CGS");

	poisson.run();

	typedef Operation_<DIM> Op;
	auto spvf = Op::Grad(poisson.get_CS("phi"));
	auto spvc = Op::GradCenter(poisson.get_CS("phi"));
	auto vfv = Op::InterpolateC2F(poisson.get_CS("phi"));
	auto vfdiv = Op::Div(vfv);
	auto v = Op::InterpolateCoordinate(poisson.get_CS("phi"), 0, 0, 0);

	Plotly plt;
	//plt.add(Heatmap(grid, *(poisson.get_CS("phi"))));
//	auto actor1 = ScalarCenter(*(poisson.get_CS("phi")));
	auto actor1 = ScalarCenter(*vfdiv);
	actor1->set_colorscale_range(0, 1);
	actor1->set_colorscale_name("RdBu");
	plt.add(actor1);
	auto actor2 = VectorFace(*vfv);
	//actor2->set_colorscale_range(0, 1);
	//actor2->set_colorscale_name("RdBu");
	plt.add(actor2);
	plt.add(WireFrame2(grid));
	//plt.add(VectorFace(*spvf));
	//plt.add(VectorCenter(*spvc));
	//plt.add(CenterPoints(grid));
	//plt.add(WireFrame(grid));
	//plt.add(WireFrameGhost(grid, *(poisson.get_ghost())));
	//plt.title("Structure");
	plt.size(800, 800);
	//plt.set_output_file("fig_phi");
	plt.plot();

	// test output
	//Output("a.txt", *(poisson.get_CS("phi")));
}

TEST(DISABLED_sgrid, h) {
	Poi pmin(0.0, 0.0);
	Poi pmax(1.0, 1.0);
	Index mn(20, 20);
	spGrid grid = spGrid(new Grid(pmin, pmax, mn, 3));

	Poisson poisson(grid);
	// set boundary condition
	Vt alpha = 10;

	Poisson::BC::Fun exact = [](Vt t, Vt x, Vt y, Vt z) {
		return 2 * x * x + y * y;
	};

	Poisson::Function source = [alpha](Vt t, Vt x, Vt y, Vt z) {
		return 6 - alpha *( 2 * x * x + y * y);
	};
	poisson.set_source(source);
	poisson.set_uniform_alpha(-alpha);

	Poisson::spBC bc0(new Poisson::BC());
	Poisson::BC::Fun fbc0 = [](Vt t, Vt x, Vt y, Vt z) {
		return 2 * x * x + y * y;
	};
	bc0->set_default_1_bc(fbc0);
	poisson.add_bc(0, 0, "phi", bc0);
	poisson.add_bc(0, 1, "phi", bc0);
	poisson.add_bc(0, 2, "phi", bc0);
	poisson.add_bc(0, 3, "phi", bc0);

	poisson.add_CS("exact", exact);

	//poisson.set_output_solver_residual("residual");
	poisson.set_solver("IC_CGS", 0.0, 1000, 1e-9);

	poisson.run();

	auto spphi = poisson.get_CS("phi");
	int res = Output("center_phi", *spphi);

	auto spexa = poisson.get_CS("exact");

	typedef Operation_<DIM> Op;
	auto error = Op::new_Minus(spexa, spphi);

	//output
	Plotly plt;
	//plt.add(Heatmap(grid, *(poisson.get_CS("phi"))));
	//auto actor1 = ScalarCenter(*(poisson.get_CS("exact")));
	//auto actor1 = ScalarCenter(*(error));
	auto actor1 = ScalarCenter(*(spphi));
	//auto actor1 = ScalarCenter(*vfdiv);
	actor1->set_colorscale_range(0, 1);
	actor1->set_colorscale_name("RdBu");
	plt.add(actor1);
	//auto actor2 = VectorFace(*vfv);
	//actor2->set_colorscale_range(0, 1);
	//actor2->set_colorscale_name("RdBu");
	//plt.add(actor2);
	plt.add(WireFrame2(*grid));
	//plt.add(VectorFace(*spvf));
	//plt.add(VectorCenter(*spvc));
	//plt.add(CenterPoints(grid));
	//plt.add(WireFrame(grid));
	//plt.add(WireFrameGhost(grid, *(poisson.get_ghost())));
	//plt.title("Structure");
	plt.size(800, 800);
	//plt.set_output_file("fig_phi");
	plt.plot();
}

TEST(sgrid, t) {
	Poi pmin(0.0, 0.0);
	Poi pmax(1.0, 1.0);
	Index mn(20, 20);
	spGrid grid = spGrid(new Grid(pmin, pmax, mn, 3));

	Poisson poisson(grid);
	// set boundary condition
	// poisson.set_source(source);
	// poisson.set_uniform_alpha(-alpha);

	Poisson::spBC bc0(new Poisson::BC());
	Poisson::spBC bc1(new Poisson::BC());
	bc0->set_default_1_bc(0);
	bc1->set_default_1_bc(1);
	poisson.add_bc(0, 0, "phi", bc1);
	poisson.add_bc(0, 1, "phi", bc0);
	poisson.add_bc(0, 2, "phi", bc0);
	poisson.add_bc(0, 3, "phi", bc0);

	poisson.set_time(100, 0.000125);
	// source
	typename Poisson::Function fun = [](Vt t, Vt x, Vt y, Vt z) {
			return (x > 0.1 && x < 0.3 && y > 0.1 && y < 0.3) ? 10: 0;
		};
	//poisson.set_source(fun);

	//poisson.set_output_solver_residual("residual");
	//poisson.set_solver("IC_CGS", 0.0, 1000, 1e-9);
	poisson.set_output_time(-1, -1, 10,
			Poisson::Event::START | Poisson::Event::AFTER);

	poisson.run();

	auto spphi = poisson.get_CS("phi");
	int res = Output("center_phi", *spphi);

	//auto spexa = poisson.get_CS("exact");

	//typedef Operation_<DIM> Op;
	//auto error = Op::new_Minus(spexa, spphi);

	//output
	Plotly plt;
	//plt.add(Heatmap(grid, *(poisson.get_CS("phi"))));
	//auto actor1 = ScalarCenter(*(poisson.get_CS("exact")));
	//auto actor1 = ScalarCenter(*(error));
	auto actor1 = ScalarCenter(*(spphi));
	//auto actor1 = ScalarCenter(*vfdiv);
	actor1->set_colorscale_range(0, 1);
	actor1->set_colorscale_name("RdBu");
	plt.add(actor1);
	//auto actor2 = VectorFace(*vfv);
	//actor2->set_colorscale_range(0, 1);
	//actor2->set_colorscale_name("RdBu");
	//plt.add(actor2);
	plt.add(WireFrame2(*grid));
	//plt.add(VectorFace(*spvf));
	//plt.add(VectorCenter(*spvc));
	//plt.add(CenterPoints(grid));
	//plt.add(WireFrame(grid));
	//plt.add(WireFrameGhost(grid, *(poisson.get_ghost())));
	//plt.title("Structure");
	plt.size(800, 800);
	//plt.set_output_file("fig_phi");
	plt.plot();
}

}
#endif
