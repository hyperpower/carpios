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

typedef Polygon_<double> Polygon;
typedef Contour_<double> Contour;
typedef Point_<double, 2> Point;
typedef PointChain_<double, 2> PointChain;
typedef PolygonPartition_<double, 2> PP;
typedef Operation_<double, 2> Op;

typedef GPA_Geometry_<double, 2> GpActor;



}

#endif
