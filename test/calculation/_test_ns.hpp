#ifndef _TEST_NS_HPP_
#define _TEST_NS_HPP_



#include "gtest/gtest.h"
#include "domain/domain.hpp"
#include "io/gnuplot_actor.h"
#include "../../src/calculation/poisson.hpp"
#include "../../src/calculation/ns.hpp"
#include "../../src/calculation/timestep.hpp"
#include <math.h>

namespace carpio{

TEST(ns, unigrid) {
	const St dim = 2;
	typedef NS_<Float, Float, dim> NS;
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

	NS ns(&domain);


}


}


#endif
