#ifndef _TEST_NS_EXPLICIT_HPP_
#define _TEST_NS_EXPLICIT_HPP_

#include "gtest/gtest.h"
#include "domain/domain.hpp"
#include "io/gnuplot_actor.h"
#include "../../src/calculation/poisson.hpp"
#include "../../src/calculation/ns_explicit.hpp"
#include "../../src/calculation/timestep.hpp"
#include <math.h>

namespace carpio {

TEST(ns, unigrid) {
	const St dim = 2;
	typedef NS_explicit_<Float, Float, dim> NSe;
	typedef NSe::BoundaryCondition BC;
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

	// boundary condition
	BC bc0, bc1;
	bc0.set_default_1_bc(0.0);
	bc1.set_default_1_bc(1.0);

	NSe nse(&domain);
	nse.show_events();
	nse.show_values();
	nse.show_functions();
	nse.show_variables_c();
	nse.show_variables_ut();
	nse.set_boundary_condition(0, 0, nse.get_var_center_idx("u"), &bc0);
	nse.set_boundary_condition(0, 1, nse.get_var_center_idx("u"), &bc0);
	nse.set_boundary_condition(0, 2, nse.get_var_center_idx("u"), &bc1);
	nse.set_boundary_condition(0, 3, nse.get_var_center_idx("u"), &bc0);
	nse.set_boundary_condition(0, 0, nse.get_var_center_idx("v"), &bc0);
	nse.set_boundary_condition(0, 1, nse.get_var_center_idx("v"), &bc0);
	nse.set_boundary_condition(0, 2, nse.get_var_center_idx("v"), &bc0);
	nse.set_boundary_condition(0, 3, nse.get_var_center_idx("v"), &bc0);

	nse.run();
	// show ================================
	Gnuplot gp;
	gp.add(
			GnuplotActor::LeafNodesContour((*domain.pgrid()),
					nse.get_var_center_idx("u")));
	//lga.push_back(
	//		GnuplotActor::LeafNodes((*domain.pgrid())));
	gp.set_equal_ratio();
	//gp.set_xrange(2.0,3.0);
	//gp.set_yrange(1.5,2.5);
	//gp.set_xlogscale(10);
	//gp.set_ylogscale(10);
	gp.plot();

	std::cout << "fin --------------\n";

}

}

#endif
