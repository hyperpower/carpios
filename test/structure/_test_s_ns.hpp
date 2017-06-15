#ifndef _TEST_S_NS_HPP
#define _TEST_S_NS_HPP

#include "gtest/gtest.h"
#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_ns_explicit.hpp"
#include "structure/s_ns_kim.hpp"
#include "structure/s_io_plotly.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <memory>
#include <vtkVersion.h>

namespace structure {

const short DIM = 2;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid_<DIM> > spGrid;
typedef Index_<DIM> Index;

typedef NS_kim_<DIM> NS;

TEST(sgrid, ns) {
	Poi pmin(-0.5, -0.5, 0.0);
	Poi pmax(0.5, 0.5, 0.0);
	Index mn(10, 10, 2);
	Grid grid(pmin, pmax, mn, 1);
	spGrid spgrid(new Grid(pmin, pmax, mn, 1));
	Poi ver = grid.v(3, 0, 0);
	NS ns(spgrid);

	// set time
	ns.set_time(10, 1e-1);

	NS::spBC bc0(new NS::BC()), bc1(new NS::BC());
	bc0->set_default_1_bc(0);
	bc1->set_default_1_bc(1);
	ns.add_bc(0, 0, "u", bc0);
	ns.add_bc(0, 1, "u", bc0);
	ns.add_bc(0, 2, "u", bc0);
	ns.add_bc(0, 3, "u", bc1);
	ns.add_bc(0, 0, "v", bc0);
	ns.add_bc(0, 1, "v", bc0);
	ns.add_bc(0, 2, "v", bc0);
	ns.add_bc(0, 3, "v", bc0);

	//typename NS::Function fun = [](Vt t, Vt x, Vt y, Vt z) {
	//	return  x*x + y*y - 4;
	//};
	Vt Re = 1000;
	ns.set_uniform_rho(Re);

	ns.set_output_time(3, -1, 1,
			NS::Event::AFTER | NS::Event::END);
	//ns.set_output_cs2file("u", "u", 0, -1, 1,
	//			NS::Event::START |NS::Event::AFTER);

	//ns.show_events();
	ns.set_projection_solver("IC_CGS", 100, 1e-6);
	//ns.set_stop_cs("u", 1e-2,1e-2,1e-2, 0, -1, 1, NS::Event::AFTER);

	ns.run();
	//ns.set_CS("v", fun);
	//ns.apply_bc("v");

	Plotly plt;
	//plt.add(SCALARCENTER(*(ns.get_CS("us"))));
	plt.add(VectorFace(*(ns._veo_f())));
	//plt.add(SCALARCENTER(*(ns.get_CS("u"))));
	//plt.add(SCALARCENTER(*(ns.get_CS("v"))));
	plt.add(WireFrame2(*(ns.get_grid())));
	//	auto actor1 = ScalarCenter(*(ns.get_CS("phi")));
	plt.size(800, 800);
	plt.plot();
}

}

#endif
