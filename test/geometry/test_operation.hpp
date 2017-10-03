#ifndef _TEST_OPERATION_HPP_
#define _TEST_OPERATION_HPP_
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

const St dim = 3;
typedef Polygon_<double> Polygon;
typedef Contour_<double> Contour;
typedef Point_<double, dim> Point;
typedef Point_<double, 3> Point3;
typedef Point_<double, 2> Point2;
typedef PointChain_<double, 2> PointChain;
typedef PolygonPartition_<double, dim> PP;
typedef Operation_<double, 3> Op3;
typedef Operation_<double, 2> Op2;

typedef GPA_Geometry_<double, 2> GpActor;

TEST(operation, tripleproduct){
	Point3 x(1, 0, 0);
	Point3 y(0, 1, 0);
	Point3 z(0, 0, 1);

	std::cout<< "x      = " << x <<std::endl;
	std::cout<< "y      = " << y <<std::endl;
	std::cout<< "z      = " << z <<std::endl;
	std::cout<< "Volume = x, y, z " << Op3::TripleScalar(x, y, z)<< std::endl;
	std::cout<< "change order " << std::endl;
	std::cout<< "         y, x, z " << Op3::TripleScalar(y, x, z)<< std::endl;
	std::cout<< "         z, x, y " << Op3::TripleScalar(z, x, y)<< std::endl;

}

TEST(operation, CCW){
	Point2 x(0, 0 );
	Point2 y(1, 0);
	Point2 z(0, 1);

	std::cout<< "x      = " << x <<std::endl;
	std::cout<< "y      = " << y <<std::endl;
	std::cout<< "z      = " << z <<std::endl;
	std::cout<< "Volume = x, y, z " << Op2::IsCCW(x, y, z)<< std::endl;
	std::cout<< "change order " << std::endl;
	std::cout<< "         y, x, z " << Op2::IsCCW(y, x, z)<< std::endl;
	std::cout<< "         z, x, y " << Op2::IsCCW(z, x, y)<< std::endl;

}

}

#endif
