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

TEST(Operation, Intersectiontype) {

	Point p1(0, 0);
	Point p2(1, 0);

	Point p3(1, 0);
	Point p4(2, 0);

	int type = Op::IntersectType(p1, p2, p3, p4);
	std::cout << "Type = " << ToString_SegmentIntersectType(type) << "\n";
	bool b = Op::IsSegmentIntersect(p1, p2, p3, p4);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_NORMAL);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT);
	ASSERT_EQ(b, true);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT_2);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT_2);
	ASSERT_EQ(b, false);
}

TEST(Operation, Intersectiontype2) {

	Point p1(0, 0);
	Point p2(1.5, 0);

	Point p3(1, 0);
	Point p4(2, 0);

	int type = Op::IntersectType(p1, p2, p3, p4);
	std::cout << "Type = " << ToString_SegmentIntersectType(type) << "\n";
	bool b = Op::IsSegmentIntersect(p1, p2, p3, p4);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_NORMAL);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT_2);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT_2);
	ASSERT_EQ(b, true);
}

TEST(Operation, Intersectiontype3) {

	Point p1(0, 0);
	Point p2(1.5, 0);

	Point p3(0.5, 0);
	Point p4(2, 1);

	int type = Op::IntersectType(p1, p2, p3, p4);
	std::cout << "Type = " << ToString_SegmentIntersectType(type) << "\n";
	bool b = Op::IsSegmentIntersect(p1, p2, p3, p4);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_NORMAL);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT);
	ASSERT_EQ(b, true);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT_2);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT_2);
	ASSERT_EQ(b, false);
}

TEST(Operation, Intersectiontype4) {

	// Intersection normally
	Point p1(0, 0);
	Point p2(1.5, 0);

	Point p3(0.5, -0.1);
	Point p4(0.5, 1);

	int type = Op::IntersectType(p1, p2, p3, p4);
	std::cout << "Type = " << ToString_SegmentIntersectType(type) << "\n";
	bool b = Op::IsSegmentIntersect(p1, p2, p3, p4);
	ASSERT_EQ(b, true);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_NORMAL);
	ASSERT_EQ(b, true);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT_2);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT_2);
	ASSERT_EQ(b, false);
}

TEST(Operation, Intersectiontype5) {

	// point to point intersection
	// only one
	Point p1(0, 0);
	Point p2(1.5, 0);

	Point p3(1.5, 0);
	Point p4(0.5, 1);

	int type = Op::IntersectType(p1, p2, p3, p4);
	std::cout << "Type = " << ToString_SegmentIntersectType(type) << "\n";
	bool b = Op::IsSegmentIntersect(p1, p2, p3, p4);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_NORMAL);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT);
	ASSERT_EQ(b, true);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT_2);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT_2);
	ASSERT_EQ(b, false);
}

TEST(Operation, Intersectiontype6) {

	// segment is same
	Point p1(0, 0);
	Point p2(1.5, 0);

	Point p3(1.5, 0);
	Point p4(0, 0);

	int type = Op::IntersectType(p1, p2, p3, p4);
	std::cout << "Type = " << ToString_SegmentIntersectType(type) << "\n";
	bool b = Op::IsSegmentIntersect(p1, p2, p3, p4);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_NORMAL);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT);
	ASSERT_EQ(b, false);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_POINT_2);
	ASSERT_EQ(b, true);
	b = Op::IsSegmentIntersect(p1, p2, p3, p4, INTERSECT_POINT_SEGMENT_2);
	ASSERT_EQ(b, false);
}

}

#endif
