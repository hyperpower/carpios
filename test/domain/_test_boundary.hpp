#ifndef __TEST_BOUNDARY_H_
#define __TEST_BOUNDARY_H_

#include "../io/gnuplot.h"

#include "../domain/domain.hpp"
#include "../calculation/advection.hpp"
#include "gtest/gtest.h"
#include <math.h>

namespace carpio {

TEST(DISABLED_Boundary, square) {
	//
	const St dim = 2;
	// new shape--------------------
	Shape2D shape;
	CreatCircle(shape, 0.001, 0.001, 1.501, 10);
	//Float x1 = 1.5, y1 = 1.5, x2 = 3.5, y2 = 3.5;
	//CreatCube(shape, x1, y1, x2, y2);
	// shape is out bound
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 2, 4);
	domain.build();
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_LeafNodes(ga, domain.grid());
	lga.push_back(ga);
	GnuplotActor_GhostNodes(ga, domain.ghost());
	lga.push_back(ga);
	GnuplotActor_GhostNodesContour_BoundaryIndex(ga, domain.ghost());
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);

	Gnuplot gp;
	gp.set_equal_ratio();
	gp.plot(lga);
	//delete shape
}

TEST(DISABLED_Boundary, inner_soild) {
	//
	const St dim = 2;
	// new shape--------------------
	Shape2D shape;
	CreatCircle(shape, 0.001, 0.001, 1.501, 10);

	Shape2D inner_shape;
	CreatCircle(inner_shape, 0.4, 0.6, 0.4, 3);
	std::list<Shape2D*> lis;
	lis.push_back(&inner_shape);
	Shape2D inner_shape2;
	CreatCircle(inner_shape2, .5, -0.1, 0.3, 6);
	lis.push_back(&inner_shape2);
	//Float x1 = 1.5, y1 = 1.5, x2 = 3.5, y2 = 3.5;
	//CreatCube(shape, x1, y1, x2, y2);
	// shape is out bound
	// define unit length
	Float UL = 1.0;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, lis, UL, 2, 4);
	domain.build();
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_LeafNodes(ga, domain.grid());
	lga.push_back(ga);
	GnuplotActor_GhostNodes(ga, domain.ghost());
	lga.push_back(ga);
	GnuplotActor_GhostNodesContour_BoundaryIndex(ga, domain.ghost());
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);

	Gnuplot gp;
	gp.set_equal_ratio();
	gp.plot(lga);
	//delete shape
}

TEST(Boundary, index_boundary) {
	//
	const St dim = 2;
	// new shape--------------------
	Shape2D shape;
	CreatCircle(shape, 0.001, 0.001, 1.503, 10);
	//Float x1 = 1.5, y1 = 1.5, x2 = 3.5, y2 = 3.5;
	//CreatCube(shape, x1, y1, x2, y2);
	// shape is out bound
	// define unit length
	Float UL = 1.5;
	// build grid ------------------
	Domain_<Float, Float, dim> domain(&shape, UL, 3, 4);
	domain.build();
	// show ================================
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_LeafNodes(ga, domain.grid());
	lga.push_back(ga);
	GnuplotActor_GhostNodes(ga, domain.ghost());
	lga.push_back(ga);
	GnuplotActor_GhostNodesContour_BoundaryIndex(ga, domain.ghost());
	lga.push_back(ga);
	GnuplotActor_Shape2D(ga, shape, 0);
	lga.push_back(ga);

	Gnuplot gp;
	gp.set_equal_ratio();
	gp.plot(lga);
	//delete shape
}

}
#endif
