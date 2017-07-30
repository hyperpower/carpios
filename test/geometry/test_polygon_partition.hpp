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
typedef Contour_<double> Contour;
typedef Point_<double, 2> Point;
typedef PointChain_<double, 2> PointChain;
typedef PolygonPartition_<double, 2> PP;
typedef Operation_<double, 2> Op;
typedef Creation_<double, 2> Cr;

typedef GPA_Geometry_<double, 2> GpActor;
typedef IOFile_Geometry_<double, 2> IOF;
using namespace std;

//template<class C, class CT>
//void top_(C&, CT&);
//
//template<class C>
//void top(C& c){
//	top_(c, C::Tag());
//	std::cout<<"top\n";
//}
//
//template<class C>
//void top_(C& c, Tag_point){
//	std::cout<<"top point\n";
//}

TEST(Polygon_p, DISABLED_isconvex) {
	// test is convex
	// this test comfired that
	// IsConvex == IsCCW
	for (int i = 1; i < 10000; i++) {
		Point p1(Random::nextDouble(-1000, 1000),
				Random::nextDouble(-1000, 1000));
		Point p2(Random::nextDouble(-1000, 1000),
				Random::nextDouble(-1000, 1000));
		Point p3(Random::nextDouble(-1000, 1000),
				Random::nextDouble(-1000, 1000));
		bool res = PP::IsConvex(p1, p2, p3);

		bool res2 = Op::IsCCW(p1, p2, p3);
		bool resr = Op::IsCCW(p3, p2, p1);
		ASSERT_EQ(res, res2);
//		std::cout << "P1 = " << p1 << "\n";
//		std::cout << "P2 = " << p2 << "\n";
//		std::cout << "P3 = " << p3 << "\n";
//		std::cout << "Is convex   = " << res << "\n";
//		std::cout << "Is CCW      = " << res2 << "\n";

//		PointChain lp;
//		lp.push_back(p1);
//		lp.push_back(p2);
//		lp.push_back(p3);
//		lp.set_close();
//		Gnuplot gp;
//		gp.add(GpActor::Lines(lp));
//      gp.plot();
	}
}

TEST(Polygon, DISABLED_inside) {
	Point p1(0, 0);
	Point p2(1, 0);
	Point p3(0.5, 0.5);
	Point p(0.6, 0.6);
	bool res = Op::IsInOn(p1, p2, p3, p);
	std::cout << "Is inside   = " << res << "\n";

	PointChain lp;
	lp.push_back(p1);
	lp.push_back(p2);
	lp.push_back(p3);
	lp.set_close();
	Gnuplot gp;
	gp.add(GpActor::Lines(lp));
	gp.add(GpActor::Points(p));
	//gp.plot();
}

TEST(Polygon, DISABLED_intersect) {
	Point P1(-3.57639, -1.8857);
	Point P2(-9.8111, -5.23389);
	Point P3(-8.1276, -7.13455);
	Point P4(-3.90168, 3.49251);
	bool resb = Op::IsSegmentBoxIntersect(P1, P2, P3, P4);
	std::cout << "Is box ntersects   = " << resb << "\n";
	bool res = Op::IsSegmentIntersect(P1, P2, P3, P4);
	std::cout << "Is Intersects   = " << res << "\n";

	PointChain lp;
	lp.push_back(P1);
	lp.push_back(P2);
	PointChain lp2;
	lp2.push_back(P3);
	lp2.push_back(P4);

	Gnuplot gp;
	gp.add(GpActor::Lines(lp, 1));
	gp.add(GpActor::Lines(lp2));
	gp.plot();
}

TEST(Polygon, DISABLED_ec1) {
	Point p1(0, 0);
	Point p2(1, 0);
	Point p3(0.5, 1);
	Point p4(0.25, 1);
	Point p5(0.12, 0.3);
	PointChain lp;
	lp.push_back(p1);
	lp.push_back(p2);
	lp.push_back(p3);
	lp.push_back(p4);
	lp.push_back(p5);
	lp.set_close();

	bool res = lp.is_simple();

	std::list<PointChain> lres;
	int r = PP::EerClipping(lp, lres);
	std::cout << "return code = " << r << std::endl;
	ASSERT_EQ(r, 1);
	Gnuplot gp;
	gp.add(GpActor::Lines(lp, -1, 0));
	int i = 0;
	for (auto& pc : lres) {
		gp.add(GpActor::Lines(pc, -1, i));
		i++;
	}
	//gp.plot();
}

TEST(Polygon, is_simple) {
	PointChain lp;
	lp.push_back(Point(-9.8111, -5.23389));
	lp.push_back(Point(-8.1276, -7.13455));
	lp.push_back(Point(-3.90168, 3.49251));
	lp.push_back(Point(-9.91016, -5.12197));
	lp.push_back(Point(-5.67673, 8.65985));
	lp.push_back(Point(0.222038, -2.21205));
	lp.push_back(Point(-3.57639, -1.8857));
	lp.set_close();

	bool res = Op::IsSimple(lp.begin(), lp.end(), true);

	std::cout << "Is Simple : " << res << std::endl;
	std::cout << "Perimeter : " << lp.perimeter() << std::endl;
	Gnuplot gp;
	gp.add(GpActor::Lines(lp, -1, 0));
	//gp.plot();
}

TEST(Polygon, random) {
	PointChain lp;
	lp.set_close();
	Random::randomSeed();
	Cr::RandomSimplePointChain(lp, 8);
	bool res = lp.is_simple();
	std::cout << "Is Simple : " << res << std::endl;
	Gnuplot gp;
	gp.add(GpActor::Lines(lp));
	//gp.plot();
}

TEST(Polygon, DISABLED_random_ec) {
	PointChain lp;
	lp.set_close();
	Random::randomSeed();
	//Random::seed();
	Cr::RandomSimplePointChain(lp, 8);
	lp.show();
	bool res = lp.is_simple();
	std::cout << "Is Simple : " << res << std::endl;

	std::list<PointChain> lres;
	int r = PP::EerClipping(lp, lres);
	std::cout << "return code = " << r << std::endl;
	//ASSERT_EQ(r, 1);
	Gnuplot gp;
	gp.add(GpActor::Lines(lp, -1, 0));
//	int i = 0;
//	for (auto& pc : lres) {
//		gp.add(GpActor::Lines(pc, -1, i));
//		i++;
//	}

	gp.plot();
}

TEST(Polygon, EC2) {
	PointChain lp;
	lp.set_close();
	lp.push_back(Point(1.98886, -7.84437));
	lp.push_back(Point(-5.55701, -1.26642));
	lp.push_back(Point(-4.90507, -4.6697));
	lp.push_back(Point(1.37726, -9.60244));
	lp.push_back(Point(3.30717, -9.23728));
	lp.push_back(Point(-8.35034, 7.40719));
	lp.push_back(Point(-4.61044, -1.70564));
	lp.push_back(Point(0.989662, -6.21943));

	std::list<PointChain> lres;
	int r = PP::EerClipping(lp, lres);
	std::cout << "return code = " << r << std::endl;
	//ASSERT_EQ(r, 1);
	Gnuplot gp;
	int i = 0;
	for (auto& pc : lres) {
		gp.add(GpActor::Lines(pc, -1, i));
		i++;
	}
	gp.add(GpActor::Lines(lp, -1, 0));

	//gp.plot();
}

TEST(Polygon, read) {
	PointChain pc;
	IOF::ReadPointChain("./man", pc);
	pc.set_close();
	std::cout << "size = " << pc.size() << std::endl;
	std::list<PointChain> lres;
	int r = PP::EerClipping(pc, lres);

	std::cout << "return code = " << r << std::endl;
	//ASSERT_EQ(r, 1);
	Gnuplot gp;
	int i = 0;
	for (auto& t : lres) {
		gp.add(GpActor::Lines(t, -1, i));
		i++;
	}
	gp.add(GpActor::Lines(pc, -1, 0));

	IOF::WritePointChain("./man2", pc);

	gp.plot();
}

}

#endif /* TEST_GEOMETRY_TEST_POLYGON_HPP_ */
