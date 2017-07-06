/*
 * test_polygon.hpp
 *
 *  Created on: Jun 12, 2017
 *      Author: zhou
 */

#ifndef _TEST_POLYGON_HPP_
#define _TEST_POLYGON_HPP_
#include "gtest/gtest.h"
#include "geometry/_polygon.hpp"
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

typedef Polygons_<double> Polygons;
using namespace std;


std::string getexepath()
{
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  return std::string( result, (count > 0) ? count : 0 );
}

TEST(Polygon, read) {
	std::cout<< " path : " <<getexepath() << std::endl;
	string workdir = "./test/input_files/polygons/samples";
	std::string fn = "/polygonwithholes";
	Polygons p(workdir + fn);
	cout << "Number of contours: " << p.ncontours() << '\n';
	cout << "Number of vertices: " << p.nvertices() << '\n';
	p.computeHoles();
	// show information
	for (int i = 0; i < p.ncontours(); i++) {
		cout << "--- new contour ---\n";
		cout << "Identifier: " << i << "  ";
		cout << (p.contour(i).external() ? "External" : "Internal")
				<< " contour\n";
		cout << "Orientation: "
				<< (p.contour(i).clockwise() ? "clockwise" : "counterclockwise")
				<< '\n';
		cout << "Holes identifiers: ";
		for (int j = 0; j < p.contour(i).nholes(); j++)
			cout << p.contour(i).hole(j) << "  ";
		cout << '\n';
		cout << "Vertices: ";
		for (int j = 0; j < p.contour(i).nvertices(); j++)
			cout << p.contour(i).vertex(j) << "  ";
		cout << '\n';
	}

}
}

#endif /* TEST_GEOMETRY_TEST_POLYGON_HPP_ */
