#ifndef _CREATION_HPP_
#define _CREATION_HPP_

#include "geometry_define.hpp"
#include <array>
#include "_point.hpp"
#include "_segment.hpp"
#include "_polygon.hpp"
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


	typedef Polygon_<TYPE> Polygon;
public:
	static void FromFile(Polygon& res, const std::string& filename){
		return Polygon(filename);
	}




};
}


#endif









}

#endif
