#ifndef _GEOMETRY_HPP_
#define _GEOMETRY_HPP_

#include <array>
#include "geometry_define.hpp"
#include "_point.hpp"
#include "_line.hpp"
#include "_segment.hpp"
#include "_plane.hpp"
#include "_polygon.hpp"
#include "_contour.hpp"
#include "_bbox.hpp"
#include "_actor_gnuplot.hpp"
#include "_sweep.hpp"
//#include "_polygon_boolean.hpp"
#include "_operation.hpp"
#include "_creation.hpp"

namespace carpio {

typedef Point_<Float, 2> Point_2D;
typedef Point_<Float, 3> Point_3D;

typedef Segment_<Float, 2> Segment_2D;
typedef Segment_<Float, 3> Segment_3D;

typedef Polygon_<Float> Polygon;

}

#endif
