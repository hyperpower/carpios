/*
 * test_polygon.hpp
 *
 *  Created on: Jun 12, 2017
 *      Author: zhou
 */

#ifndef _TEST_contour_HPP_
#define _TEST_contour_HPP_
#include "gtest/gtest.h"
#include "geometry/geometry.hpp"
#include "utility/random.h"
#include "utility/clock.h"
#include "utility/format.h"

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

typedef Polygon_<double> Polygons;
using namespace std;


//inline std::string getexepath()
//{
//  char result[ PATH_MAX ];
//  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
//  return std::string( result, (count > 0) ? count : 0 );
//}

TEST(contour, read) {
	fmt::print("{}", "-------test contour-----------");

}
}

#endif /* TEST_GEOMETRY_TEST_POLYGON_HPP_ */
