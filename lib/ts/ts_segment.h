/*
 * ts_segment.h
 *
 *  Created on: May 31, 2015
 *      Author: zhou
 */

#ifndef TS_SEGMENT_H_
#define TS_SEGMENT_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"

namespace TS {

template<class TYPE, st DIM> class Segment;
template<class TYPE, st DIM> class Surface;
template<class TYPE, st DIM> class Edge;

template<class TYPE, st DIM>
class Segment {
public:
	typedef Segment<TYPE, DIM> self_class;
	typedef st size_type;
	typedef Point<TYPE, DIM> Poi;
	typedef std::shared_ptr<Poi> pPoi;
	typedef Vertex<TYPE, DIM> Ver;
	typedef std::shared_ptr<Ver> pVer;
	typedef Segment<TYPE, DIM> Seg;
	typedef std::shared_ptr<Seg> pSeg;
	typedef Edge<TYPE, DIM> Edg;
	typedef std::shared_ptr<Edg> pEdg;
	typedef Triangle<TYPE, DIM> Tri;
	typedef std::shared_ptr<Tri> pTri;
	typedef Surface<TYPE, DIM> Sur;
	typedef std::shared_ptr<Sur> pSur;
	typedef List<pSeg> list_pSeg;
	typedef List<pVer> list_pVer;
	typedef List<pTri> list_pTri;
	typedef List<pSur> list_pSur;
	typedef List<pSeg> list_pEdg;
public:
	pVer v1;
	pVer v2;
public:
	Segment() {
		v1 = nullptr;
		v2 = nullptr;
	}
	Segment(const pVer& a, const pVer& b) {
		assert(a != nullptr);
		assert(b != nullptr);
		v1 = a;
		v2 = b;
		//v1->segments.push_back(pe);
		//v2->segments.push_back(pe);
	}
	Segment(const pVer& a, const pVer& b, pEdg pe) {
		assert(a != nullptr);
		assert(b != nullptr);
		v1 = a;
		v2 = b;
		v1->segments.push_back(pe);
		v2->segments.push_back(pe);
	}

};

//function out of class
/**
 * gts_segments_are_intersecting:
 * @s1: a #GtsSegment.
 * @s2: a #GtsSegment.
 *
 * Returns: %GTS_IN if @s1 and @s2 are intersecting, %GTS_ON if one of the
 * endpoints of @s1 (resp. @s2) lies on @s2 (resp. @s1), %GTS_OUT otherwise.
 */
template<class TYPE, st DIM>
Intersect segments_are_intersecting( //
		const Segment<TYPE, DIM>& s1, //
		const Segment<TYPE, DIM>& s2) {
	typename Segment<TYPE, DIM>::Poi * p1, *p2, *p3, *p4;
	Float d1, d2, d3, d4;

	p1 = s1.v1;
	p2 = s1.v2;
	p3 = s2.v1;
	p4 = s2.v2;
	d1 = point_orientation(p1, p2, p3);
	d2 = point_orientation(p1, p2, p4);
	if ((d1 > 0.0 && d2 > 0.0) || (d1 < 0.0 && d2 < 0.0))
		return TS_OUT;
	d3 = point_orientation(p3, p4, p1);
	d4 = point_orientation(p3, p4, p2);
	if ((d3 > 0.0 && d4 > 0.0) || (d3 < 0.0 && d4 < 0.0))
		return TS_OUT;
	if (d1 == 0.0 || d2 == 0.0 || d3 == 0.0 || d4 == 0.0)
		return TS_ON;
	return TS_IN;
}

template<class TYPE, st DIM>
inline bool segment_connect(  //
		Segment<TYPE, DIM> s, //
		Vertex<TYPE, DIM> e1, //
		Vertex<TYPE, DIM> e2) { //
	return ((*(s.v1)) == e1 && (*(s.v2)) == e2)
			|| ((*(s.v1)) == e2 && (*(s.v2)) == e1);
}
template<class TYPE, st DIM>
inline bool segments_are_identical(     //
		Segment<TYPE, DIM>* s1,   //
		Segment<TYPE, DIM>* s2) { //
	return ((*(s1->v1)) == (*(s2->v1)) && (*(s1->v2)) == (*(s2->v2)))
			|| ((*(s1->v1)) == (*(s2->v2)) && (*(s1->v2)) == (*(s2->v1)));
}
template<class TYPE, st DIM>
inline bool segments_touch(  //
		Segment<TYPE, DIM>* s1, //
		Segment<TYPE, DIM>* s2) {
	return (*(s1->v1) == *(s2->v1) || *(s1->v1) == *(s2->v2)
			|| *(s1->v2) == *(s2->v1) || *(s1->v2) == *(s2->v2));
}
}

#endif /* TS_TS_SEGMENT_H_ */
