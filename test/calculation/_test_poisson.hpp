#ifndef __TEST_POISSON_H_
#define __TEST_POISSON_H_


#include "domain/domain.hpp"
#include "../../src/calculation/poisson.hpp"
#include "io/gnuplot_actor.h"
#include "../../src/calculation/timestep.hpp"
#include <math.h>
#include "gtest/gtest.h"

namespace carpio {
Float coe_set_b(Float x, Float y, Float z) {
	return 1;
}
Float coe_set_f(Float x, Float y, Float z) {
	if (IsInRange(-0.25, x, 0.0, _oo_)) {
		return 1;
	}
	return 0;
}

TEST(DISABLED_Poisson, unigrid) {
	const St dim = 2;
	typedef Poisson_<Float, Float, dim> Poisson;
	// new shape--------------------
	Shape2D shape;
	Float x1 = -0.5, y1 = -0.5, x2 = 0.5, y2 = 0.5;
	CreatCube(shape, x1, y1, x2, y2);
	// CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 5, 10);
	domain.build();

	Poisson poisson(&domain);
	//poisson.set_time_term(0.1,10,"explicit");
	poisson.set_output_time(0, 1, 1,
			Poisson::Event::START | Poisson::Event::END);
	Poisson_<Float, Float, dim>::Function fun_beta =
			[](Float, Float, Float) {return 1;};
	poisson.set_beta(fun_beta);
	//poisson.set_alpha_term(fun_beta);

	Poisson_<Float, Float, dim>::BoundaryCondition bc1;
	bc1.set_default_1_bc(1);
	Poisson_<Float, Float, dim>::BoundaryCondition bc2;
	bc2.set_default_1_bc(5);
	//bc.set_default_1_bc(exact_fun_2);
	poisson.set_bc_phi(0, 0, &bc1);
	poisson.set_bc_phi(0, 1, &bc1);
	poisson.set_bc_phi(0, 2, &bc1);
	poisson.set_bc_phi(0, 3, &bc2);
	//poisson.set_output_time(0,10, 2);
	poisson.run();
	// output
	// show ================================
	//GnuplotActor::list_spActor lga;
	//lga.push_back(
	//		GnuplotActor::LeafNodesContour((*domain.pgrid()),
	//				poisson.phi_idx()));
	Gnuplot gp;
	gp.add(GnuplotActor::LeafNodesContour((*domain.pgrid()),
					poisson.phi_idx()));
	gp.set_equal_ratio();
	//gp.set_xrange(2.0,3.0);
	//gp.set_yrange(1.5,2.5);
	//gp.set_xlogscale(10);
	//gp.set_ylogscale(10);
	gp.plot();

	std::cout << "fin --------------\n";

}

TEST(Poisson, adpgrid_source) {
	const St dim = 2;
	typedef Poisson_<Float, Float, dim> Poisson;
	// new shape--------------------
	Shape2D shape;
	Float x1 = -0.5, y1 = -0.5, x2 = 0.5, y2 = 0.5;
	CreatCube(shape, x1, y1, x2, y2);
	Shape2D cube;
	CreatCube(cube, -0.3, -0.3, -0.1, -0.1);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 5, 6);
	domain.adaptive().adapt_full();
	domain.adaptive().adapt_shape_inner(cube);
	domain.build();

	Poisson poisson(&domain);
	//poisson.set_time_term(0.1,10,"explicit");
	poisson.set_output_time(0, 1, 1,
			Poisson::Event::START | Poisson::Event::END);
	Poisson_<Float, Float, dim>::Function fun_beta =
			[](Float, Float, Float) {return 1;};
	poisson.set_beta(fun_beta);
	Poisson_<Float, Float, dim>::Function fun_source =
			[](Float x, Float y, Float z) {
				bool res = IsInRange(-0.3,x,-0.1,_cc_);
				res = res && IsInRange(-0.3,y,-0.1,_cc_);
				return res ? -20: 0;
			};
	poisson.set_f(fun_source);
	//poisson.set_alpha_term(fun_beta);

	Poisson_<Float, Float, dim>::BoundaryCondition bc1;
	bc1.set_default_1_bc(0);
	Poisson_<Float, Float, dim>::BoundaryCondition bc2;
	bc2.set_default_1_bc(0.25);
	//bc.set_default_1_bc(exact_fun_2);
	poisson.set_bc_phi(0, 0, &bc1);
	poisson.set_bc_phi(0, 1, &bc2);
	poisson.set_bc_phi(0, 2, &bc2);
	poisson.set_bc_phi(0, 3, &bc1);
	//poisson.set_output_time(0,10, 2);
	poisson.run();
	// output
	// show ================================
	GnuplotActor::list_spActor lga;
	lga.push_back(
			GnuplotActor::LeafNodesContour((*domain.pgrid()),
					poisson.phi_idx()));
	//lga.push_back(
	//		GnuplotActor::LeafNodes((*domain.pgrid())));
	Gnuplot gp;
	gp.set_equal_ratio();
	//gp.set_xrange(2.0,3.0);
	//gp.set_yrange(1.5,2.5);
	//gp.set_xlogscale(10);
	//gp.set_ylogscale(10);
	gp.plot(lga);

	std::cout << "fin --------------\n";

}

TEST(Poisson, adp_time_term) {
	const St dim = 2;
	St timestep = 300;
	typedef Poisson_<Float, Float, dim> Poisson;
	// new shape--------------------
	Shape2D shape;
	Float x1 = -0.5, y1 = -0.5, x2 = 0.5, y2 = 0.5;
	CreatCube(shape, x1, y1, x2, y2);
	Shape2D cube;
	CreatCube(cube, -0.3, -0.3, -0.1, -0.1);
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 5, 6);
	domain.adaptive().adapt_full();
	domain.adaptive().adapt_shape_inner(cube);
	domain.build();

	Poisson poisson(&domain);
	//poisson.set_time_term(0.00001, timestep, "explicit");
	poisson.set_time_term(0.00001, timestep, "CrankNicolson");
	poisson.set_output_time(0, 10, 1,
			Poisson::Event::START | Poisson::Event::END
					| Poisson::Event::AFTER);
	Poisson_<Float, Float, dim>::Function fun_beta =
			[](Float, Float, Float) {return 2;};
	poisson.set_beta(fun_beta);
	Poisson_<Float, Float, dim>::Function fun_source =
			[](Float x, Float y, Float z) {
				bool res = IsInRange(-0.3,x,-0.1,_cc_);
				res = res && IsInRange(-0.3,y,-0.1,_cc_);
				return res ? -20.0: 0;
			};

	Poisson_<Float, Float, dim>::Function fun_phi =
			[](Float, Float, Float) {return 0;};
	poisson.set_phi(fun_phi);
	poisson.set_f(fun_source);

	//poisson.set_alpha_term(fun_beta);

	Poisson_<Float, Float, dim>::BoundaryCondition bc1;
	bc1.set_default_1_bc(0);
	Poisson_<Float, Float, dim>::BoundaryCondition bc2;
	bc2.set_default_1_bc(1);
	//bc.set_default_1_bc(exact_fun_2);
	poisson.set_bc_phi(0, 0, &bc1);
	poisson.set_bc_phi(0, 1, &bc2);
	poisson.set_bc_phi(0, 2, &bc2);
	poisson.set_bc_phi(0, 3, &bc1);
	//poisson.set_output_time(0,10, 2);

	typedef EventTraceCenterValue_<Float, Float, dim> EventTraceCenterValue;
	EventTraceCenterValue eventtcv(0, timestep, 1, Poisson::Event::AFTER);
	Poisson::spEvent spe = std::make_shared<EventTraceCenterValue>(eventtcv);
	poisson.add_event("trace center", spe);

	//---------------------------
	poisson.run();
	// output
	// show ================================
	GnuplotActor::list_spActor lga;
	lga.push_back(
			GnuplotActor::LeafNodesContour((*domain.pgrid()),
					poisson.phi_idx()));
	//lga.push_back(
	//		GnuplotActor::LeafNodes((*domain.pgrid())));
	Gnuplot gp;
	//gp.set_terminal_jpeg("phi.jpeg");
	gp.set_equal_ratio();
	//gp.set_xrange(2.0,3.0);
	//gp.set_yrange(1.5,2.5);
	//gp.set_xlogscale(10);
	//gp.set_ylogscale(10);
	gp.plot(lga);

	std::cout << "fin --------------\n";

}

}
#endif
