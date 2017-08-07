
#ifndef _OBJECTS_HPP_
#define _OBJECTS_HPP_

//#include <array>
#include <geometry/objects/_box.hpp>
#include "../geometry_define.hpp"
#include "_point.hpp"
#include "_box.hpp"
#include "_line.hpp"
#include "_segment.hpp"
#include "_plane.hpp"
#include "_polygon.hpp"
#include "_contour.hpp"
#include "_point_chain.hpp"
#include "_triangle.hpp"
#include "_edge.hpp"
#include "_triface.hpp"
#include "_trisurface.hpp"
#include "_vertex.hpp"
#include "_bbtree.hpp"

namespace carpio {

/// Geometry objects
///  Point
///  Line
///  Plane
///  Segment
///  Contour
///  Triangle
///  PointChain
///  Polygon

template <typename T>
struct IsGeoObj {
  static const bool value = false;
};

template<typename TYPE, St DIM>
struct IsGeoObj<Point_<TYPE, DIM> > {
  static const bool value = true;
};
template<typename TYPE>
struct IsGeoObj<Line_<TYPE> > {
  static const bool value = true;
};
template<typename TYPE>
struct IsGeoObj<Plane_<TYPE> > {
  static const bool value = true;
};
template<typename TYPE>
struct IsGeoObj<Contour_<TYPE> > {
  static const bool value = true;
};
template<typename TYPE>
struct IsGeoObj<Polygon_<TYPE> > {
  static const bool value = true;
};
template<typename TYPE, St DIM>
struct IsGeoObj<Segment_<TYPE, DIM> > {
  static const bool value = true;
};
template<typename TYPE, St DIM>
struct IsGeoObj<PointChain_<TYPE, DIM> > {
  static const bool value = true;
};
template<typename TYPE, St DIM>
struct IsGeoObj<Triangle_<TYPE, DIM> > {
  static const bool value = true;
};

}

#endif
