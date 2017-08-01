#ifndef _CREATION_HPP_
#define _CREATION_HPP_

#include "../geometry_define.hpp"
#include <array>
#include "../objects/_objects.hpp"
#include "../../utility/random.h"
#include <cmath>

namespace carpio {

template<typename TYPE, St DIM>
class Creation_ {
public:
	static const St Dim = DIM;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Point;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Segment_<TYPE, DIM>& ref_Segment;
	typedef const Segment_<TYPE, DIM>& const_ref_Segment;
	typedef PointChain_<TYPE, DIM> PointChain;

	typedef Polygon_<TYPE> Polygon;
	typedef Contour_<TYPE> Contour;
	public:
	//static void FromFile(Polygon& res, const std::string& filename) {
	//	Polygon(filename);
	//}

	static void Cube(Polygon& res, const Vt minx, const Vt miny, const Vt maxx,
			const Vt maxy) {
		ASSERT(Dim == 2);
		ASSERT(minx < maxx);
		ASSERT(miny < maxy);
		res.clear();
		std::vector<Point> vers;
		std::vector<St> holes;
		vers.push_back(Point(minx, miny));
		vers.push_back(Point(maxx, miny));
		vers.push_back(Point(maxx, maxy));
		vers.push_back(Point(minx, maxy));

		Contour con(vers, holes, true, true, true);
		res.push_back(con);
	}

	static void Triangle(Polygon& res, const Point& p1, const Point& p2,
			const Point& p3) {
		ASSERT(Dim == 2);
		res.clear();
		std::vector<Point> vers;
		std::vector<St> holes;
		vers.push_back(p1);
		vers.push_back(p2);
		vers.push_back(p3);
		Contour con(vers, holes, true, false, false);
		res.push_back(con);
	}

	static void RandomSimplePointChain(PointChain& pc, int num,
			const Point& min = Point(-10, -10, -10),
			const Point& max = Point(10, 10, 10)) {
		ASSERT(num < 15);
		pc.clear();
		for (int i = 0; i < num; i++) {
			Point p(Random::nextDouble(min[_X_], max[_X_]),
					Random::nextDouble(min[_Y_], max[_Y_]),
					Random::nextDouble(min[_Z_], max[_Z_]));
			pc.push_back(p);
			if(pc.size()>3){
				// check simple
				if(!(pc.is_simple())){
					pc.pop_back();
					--i;
				}
			}
			//std::cout<< "i = " << i << "p = "<< p <<std::endl;
		}

	}

};
}

#endif


