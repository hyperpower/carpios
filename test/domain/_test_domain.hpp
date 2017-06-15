#ifndef _TEST_DOMAIN_HPP_
#define _TEST_DOMAIN_HPP_

#include "../src/carpio_define.hpp"
//#include "../geometry/geometry.hpp"
//#include "../utility/clipper.hpp"
#include "domain/domain.hpp"
#include "algebra/matrix_SparCompRow.hpp"
#include "algebra/matrix_SparCompCol.hpp"
#include "io/gnuplot.h"
#include "io/gnuplot_actor.h"
#include <iostream>
#include <cmath>

#include "gtest/gtest.h"

namespace carpio {

TEST(Domain, try1) {
	//
	const St dim = 2;
	// new shape--------------------
	Shape2D shape;
	Float x1 = 0, y1 = 0, r = 1.5;
	CreatCircle(shape, x1, y1, r, 359);
	//CreatCube(shape, 1.5, 1.5, 2.5, 2.5);
	// shape is out bound
	//
	// define unit length
	Float UL = 0.499;
	// build grid ------------------
	Float max_x = shape.max_x();
	Float max_y = shape.max_y();
	Float min_x = shape.min_x();
	Float min_y = shape.min_y();
	St n_x = std::ceil((max_x - min_x) / UL);
	St n_y = std::ceil((max_y - min_y) / UL);
	std::cout << n_x << " " << n_y << "\n";
	Grid_<Float, Float, dim> g(n_x, min_x, UL, //
			n_y, min_y, UL);
	// g.show_info();
	// build adaptive
	Adaptive_<Float, Float, dim> adp(&g, 2, 2);
	std::cout << " here0 ----\n";
	adp.adapt_bound_solid(shape);
	std::cout << " here3 ----\n";
	// show ================================
	Gnuplot gp;
	gp.add(GnuplotActor::RootNodes(g));
	gp.add(GnuplotActor::LeafNodes(g));
	//lga.push_back(ga);
	gp.set_equal_ratio();
	gp.plot();
	//delete shape
}

TEST(Domain, try2) {
	//
	const St dim = 2;
	// new shape--------------------
	Shape2D shape;
	//CreatCircle(shape, 0.0, 0.0, 1.5, 359);
	Float x1 = 1.5, y1 = 1.5, x2 = 3.5, y2 = 3.5;
	CreatCube(shape, x1, y1, x2, y2);
	//CreatCube(shape, 1.5, 1.5, 3.5, 3.5);
	// shape is out bound
	//
	// define unit length
	Float UL = 0.25;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 2, 3);
	domain.build();
	std::cout<<"build -----\n";
	// show ================================
	GnuplotActor::list_spActor lga;

	lga.push_back(GnuplotActor::LeafNodes(domain.grid()));
	//lga.push_back(ga);
	lga.push_back(GnuplotActor::GhostNodes(domain.ghost()));
	//lga.push_back(GnuplotActor::GhostNodesContour_BoundaryIndex(domain.ghost()));
	lga.push_back(GnuplotActor::Shape( shape, 0));

	Gnuplot gp;
	gp.set_equal_ratio();
	//gp.plot(lga);
	//delete shape
}

}

#endif
