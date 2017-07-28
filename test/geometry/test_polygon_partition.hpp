/*
 * test_polygon.hpp
 *
 *  Created on: Jun 12, 2017
 *      Author: zhou
 */

#ifndef _TEST_POLYGON_PARTITION_HPP_
#define _TEST_POLYGON_PARTITION_HPP_
#include "gtest/gtest.h"
#include "geometry/geometry.hpp"
#include "utility/random.h"
#include "utility/clock.h"

#include <functional>
#include <io/plotly.h>
#include <io/plotly_actor.h>

#include <math.h>
#include <string>
#include <memory>

#include <string>
#include <limits.h>
#include <unistd.h>

namespace carpio {

typedef Polygon_<double> Polygon;
typedef Point_<double, 2> Point;
typedef PolygonPartition_<double, 2> PP;

typedef GnuplotActor_<double, 2> GpActor;
using namespace std;

TEST(Polygon_p, isconvex) {
	// test is convex
	Point p1(0, 0);
	Point p2(0, 1);
	Point p3(0.5, 0.5);
	bool res = PP::IsConvex(p1, p2, p3);
	std::cout << "P1 = " << p1 << "\n";
	std::cout << "P2 = " << p2 << "\n";
	std::cout << "P3 = " << p3 << "\n";
	std::cout << "Is convex = " << res << "\n";



	Gnuplot gp;
	gp.add(GpActor::Lines(res.contour(0)));
	//gp.plot();
}

}

#endif /* TEST_GEOMETRY_TEST_POLYGON_HPP_ */
