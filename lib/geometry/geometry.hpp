#ifndef _GEOMETRY_HPP_
#define _GEOMETRY_HPP_

#include <array>
#include "geometry_define.hpp"
#include "objects/_objects.hpp"
#include "io/_actor_gnuplot.hpp"
#include "io/_io_file.hpp"
#include "operations/_sweep.hpp"
#include "operations/_polygon_boolean.hpp"
#include "operations/_polygon_partition.hpp"
#include "operations/_operation.hpp"
#include "operations/_intersection.hpp"
#include "operations/_creation.hpp"
#include "operations/_tri_tri_intersect.hpp"

namespace carpio {

typedef Point_<Float, 2> Point_2D;
typedef Point_<Float, 3> Point_3D;

typedef Segment_<Float, 2> Segment_2D;
typedef Segment_<Float, 3> Segment_3D;

//typedef Polygon_<Float> Polygon;

}

#endif
