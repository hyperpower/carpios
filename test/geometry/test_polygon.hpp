/*
 * test_polygon.hpp
 *
 *  Created on: Jun 12, 2017
 *      Author: zhou
 */

#ifndef _TEST_POLYGON_HPP_
#define _TEST_POLYGON_HPP_
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

typedef GnuplotActor_<double, 2> GpActor;
using namespace std;

std::string getexepath() {
	char result[ PATH_MAX];
	ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
	return std::string(result, (count > 0) ? count : 0);
}

TEST(Polygon, signedarea) {
	std::cout << "test Signed Area\n";
	typedef Point_<double, 2> Point;
	typedef Operation_<double, 2> Op;
	Point p1(0.0, 0.0);
	std::cout << "P1 = " << p1 << "\n";
	Point p2(1.0, 0.0);
	std::cout << "P2 = " << p2 << "\n";
	Point p3(0.0, 1.0);
	std::cout << "P3 = " << p3 << "\n";

	double res = Op::SignedArea(p1, p2, p3);
	std::cout << "Order      = " << "1 2 3 \n";
	std::cout << "SignedArea = " << res << std::endl;
	ASSERT_EQ(res, 1);

	res = Op::SignedArea(p2, p1, p3);
	std::cout << "Order      = " << "2 1 3 \n";
	std::cout << "SignedArea = " << res << std::endl;
	ASSERT_EQ(res, -1);
}

TEST(Polygon, sweepevent) {
	typedef double Vt;
	const St Dim = 2;
	typedef SweepEvent_<Vt, Dim> SweepEvent;
	typedef SweepEventComp_<Vt, Dim> SweepEventComp;
	typedef Segment_<Vt, Dim> Segment;
	typedef Operation_<Vt, Dim> Op;
	typedef Point_<double, 2> Point;
	std::priority_queue<SweepEvent*, std::vector<SweepEvent*>, SweepEventComp> eq;

	Point p1(0.0, 0.0);
	Point p2(1.0, 0.0);
	Point p3(0.0, 1.0);
	Point p4(0.5, 1.0);

	Segment s1(p1, p2);
	Segment s2(p3, p4);
	s1.show();
	s2.show();
	//1 const Point& pp,
	//2 bool         b,
	//3 PolygonType  apl,
	//4 SweepEvent*  o,
	//5 EdgeType     t= NORMAL
	SweepEvent* e1 = new SweepEvent(s1.ps(), true, SUBJECT, nullptr);
	SweepEvent* e2 = new SweepEvent(s1.pe(), true, SUBJECT, e1);
	e1->other = e2;

	if (e1->p.x() < e2->p.x()) {
		e2->left = false;
	} else if (e1->p.x() > e2->p.x()) {
		e1->left = false;
	} else if (e1->p.y() < e2->p.y()) { // the line segment is vertical. The bottom endpoint is the left endpoint
		e2->left = false;
	} else {
		e1->left = false;
	}

	std::cout << "Queue :" << eq.size() << "\n";

	e1->show();
	e2->show();

	eq.push(e1);
	eq.push(e2);

	std::cout << "Popping out elements...\n";
	while (!eq.empty()) {
		SweepEvent* et = eq.top();
		et->show();
		//std::cout << '\n';
		eq.pop();
	}
	std::cout << '\n';
}

TEST(Polygon, trivial1) {
	typedef double Vt;
	const St Dim = 2;
	typedef Segment_<Vt, Dim> Segment;
	typedef Operation_<Vt, Dim> Op;
	typedef Creation_<Vt, Dim> Cr;
	typedef Point_<double, 2> Point;
	typedef Clip_<Vt> Clip;
	Polygon poly;
	Cr::Cube(poly, 0.0, 0.0, 1.0, 1.0);
	Polygon poly2;
	// poly2 is empty
	//Cr::Cube(poly2, 0.5, 0.5, 1.5, 1.5);

	Clip clip(poly, poly2);
	Polygon res;
	clip.compute(DIFFERENCE, res);

	Gnuplot gp;
	gp.add(GpActor::Lines(res.contour(0)));
	//gp.plot();
}

TEST(Polygon, read) {
	//std::cout<< " path : " << argc << std::endl;
	string workdir =
			"/home/zhou/git-workspace/carpios/test/input_files/polygons/samples";
	//string workdir = "./";
	std::string fn = "/polygonwithholes";
	Polygon p(workdir + fn);
	cout << "Number of contours: " << p.ncontours() << '\n';
	cout << "Number of vertices: " << p.nvertices() << '\n';
	//p.computeHoles();
	// show information
	for (int i = 0; i < p.ncontours(); i++) {
		cout << "--- new contour ---\n";
		cout << "Identifier : " << i << "  ";
		cout << (p.contour(i).external() ? "External" : "Internal")
				<< " contour\n";
		cout << "Orientation: "
				<< (p.contour(i).clockwise() ? "clockwise" : "counterclockwise")
				<< '\n';
		cout << "Holes identifiers: ";
		for (int j = 0; j < p.contour(i).nholes(); j++)
			cout << p.contour(i).hole(j) << "  ";
		cout << '\n';
		cout << "Vertices: \n";
		for (int j = 0; j < p.contour(i).nvertices(); j++)
			cout << "   " << p.contour(i).vertex(j) << "\n";
		cout << '\n';
	}

	Gnuplot gp;
	gp.add(GpActor::Lines(p.contour(0)));
	gp.add(GpActor::Lines(p.contour(1)));
	//gp.plot();
}

TEST(Polygon, normal) {
	typedef double Vt;
	const St Dim = 2;
	typedef Segment_<Vt, Dim> Segment;
	typedef Operation_<Vt, Dim> Op;
	typedef Creation_<Vt, Dim> Cr;
	typedef Point_<double, 2> Point;
	typedef Clip_<Vt> Clip;
	Polygon poly;
	Cr::Cube(poly, 0.0, 0.0, 1.0, 1.0);
	Polygon poly2;
	Cr::Cube(poly2, 0.5, 0.5, 1.5, 1.5);

	Clip clip(poly, poly2);
	Polygon res;
	clip.compute(INTERSECTION, res);

	Gnuplot gp;
	gp.add(GpActor::Lines(poly.contour(0), 0, 0));
	gp.add(GpActor::Lines(poly2.contour(0), 0, 1));
	gp.add(GpActor::Lines(res.contour(0), 0, 2));
	//gp.plot();
}

TEST(Polygon, line_touch) {
	typedef double Vt;
	const St Dim = 2;
	typedef Segment_<Vt, Dim> Segment;
	typedef Operation_<Vt, Dim> Op;
	typedef Creation_<Vt, Dim> Cr;
	typedef Point_<double, 2> Point;
	typedef Clip_<Vt> Clip;
	Polygon poly;
	Cr::Cube(poly, 0.0, 0.0, 1.0, 1.0);
	Polygon poly2;
	Cr::Cube(poly2, 1.0, 0.5, 1.5, 1.5);

	Clip clip(poly, poly2);
	Polygon res;
	clip.compute(INTERSECTION, res);

	Gnuplot gp;
	gp.add(GpActor::Lines(poly.contour(0), 0, 0));
	gp.add(GpActor::Lines(poly2.contour(0), 0, 1));
	if (res.ncontours() > 0) {
		gp.add(GpActor::Lines(res.contour(0), 0, 2));
	}
	//gp.plot();
}

TEST(Polygon, point_line_touch) {
	typedef double Vt;
	const St Dim = 2;
	typedef Segment_<Vt, Dim> Segment;
	typedef Operation_<Vt, Dim> Op;
	typedef Creation_<Vt, Dim> Cr;
	typedef Point_<double, 2> Point;
	typedef Clip_<Vt> Clip;
	Polygon poly;
	Cr::Cube(poly, 0.0, 0.0, 1.0, 1.0);
	Polygon poly2;
	Cr::Triangle(poly2, Point(1.0, 0.5), Point(1.5,0.5), Point(1.3, 1.0));

	Clip clip(poly, poly2);
	Polygon res;
	clip.compute(INTERSECTION, res);

	Gnuplot gp;
	gp.add(GpActor::Lines(poly.contour(0), 0, 0));
	gp.add(GpActor::Lines(poly2.contour(0), 0, 1));
	if (res.ncontours() > 0) {
		gp.add(GpActor::Lines(res.contour(0), 0, 2));
	}
	//gp.plot();
}



}

#endif /* TEST_GEOMETRY_TEST_POLYGON_HPP_ */
