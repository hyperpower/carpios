#ifndef _TEST_S_NS_HPP
#define _TEST_S_NS_HPP

#include "gtest/gtest.h"
#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_advection.hpp"
#include "structure/s_io_plotly.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <memory>

namespace structure {

const short DIM = 2;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid_<DIM> > spGrid;
typedef Index_<DIM> Index;
typedef Advection_<DIM> Eq;

TEST(sgrid, advection) {
	Vt dt = 0.1;
	Vt v  = 1;
	int step = 1;
	Vt xs  = 50;
	Vt xe  = xs + step * dt * v;
	std::cout<< "x end = " << xe << std::endl;


	Poi pmin(0.0, 0.0, 0.0);
	Poi pmax(200, 200, 0.0);
	Index mn(30, 30, 2);
	spGrid spgrid(new Grid(pmin, pmax, mn, 2));
	Eq eq(spgrid);

	// set time
	eq.set_time(step, dt);

	eq.initial_CS("u", 1);
	eq.initial_CS("v", 1);
	Eq::Function fun_init_phi = [xs](Vt t , Vt x, Vt y, Vt z) {
		return (x >= xs - 10 && x <= xs + 10) ? 1 : 0;
	};
	eq.initial_CS("phi", fun_init_phi);

	Eq::spBC bc0(new Eq::BC()), bc1(new Eq::BC());
	bc0->set_default_1_bc(0);
	bc1->set_default_1_bc(1);
	eq.add_bc(0, 0, "phi", bc0);
	eq.add_bc(0, 1, "phi", bc0);
	//eq.add_bc(0, 2, "phi", bc0);
	//eq.add_bc(0, 3, "phi", bc0);

	Eq::Function fun_exact = [xs,v](Vt t , Vt x, Vt y, Vt z) {
		Vt xe  = xs + t * v;
		return (x >= xe - 10 && x <= xe + 10) ? 1 : 0;
	};

	eq.add_CS("exact", fun_exact);

	//typename Eq::Function fun = [](Vt t, Vt x, Vt y, Vt z) {
	//	return  x*x + y*y - 4;
	//};


	eq.set_output_time(0, -1, 1, Eq::Event::START | Eq::Event::AFTER);
	//eq.set_output_error("exact", "phi", 0, -1, 1, Eq::Event::START | Eq::Event::AFTER);
	//eq.set_center_scalar("exact", fun_exact, 0, -1, 1, Eq::Event::START | Eq::Event::BEFORE);
	//eq.set_output_cs2file("u", "u", 0, -1, 1,
	//			Eq::Event::START |Eq::Event::AFTER);

	//eq.show_events();
	//eq.set_advection_scheme("center");
	eq.run();
	//eq.set_CS("v", fun);
	//eq.apply_bc("v");

	Plotly plt;
	plt.add(SCALARCENTER(*(eq.get_CS("phi"))));
    //plt.add(SCALARCENTER(*(eq.get_CS("exact"))));
	//plt.add(VectorFace(*(eq._veo_f())));
	//plt.add(SCALARCENTER(*(eq.get_CS("u"))));
	//plt.add(SCALARCENTER(*(eq.get_CS("v"))));
	//plt.add(WireFrame2(*(eq.get_grid())));
	//	auto actor1 = ScalarCenter(*(eq.get_CS("phi")));
	plt.size(800, 800);
	plt.plot();
}

}

#endif
