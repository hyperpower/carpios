/*
 * ts_operation.h
 *
 *  Created on: Jul 19, 2017
 *      Author: zhou
 */

#ifndef _TS_OPERATION_H_
#define _TS_OPERATION_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_edge.h"
#include "ts_face.h"
#include <fstream>
#include <sstream>
#include <math.h>

namespace TS {

template<class TYPE, st DIM>
class Operation {
	static const st Dim = DIM;
	typedef TYPE vt;
	typedef Point<TYPE, Dim> Poi;
	typedef std::shared_ptr<Poi> spPoi;
	typedef Vertex<TYPE, Dim> Ver;
	typedef std::shared_ptr<Ver> spVer;
	typedef Segment<TYPE, Dim> Seg;
	typedef std::shared_ptr<Seg> spSeg;
	typedef Edge<TYPE, Dim> Edg;
	typedef std::shared_ptr<Edg> spEdg;
	typedef Triangle<TYPE, Dim> Tri;
	typedef std::shared_ptr<Tri> spTri;
	typedef Face<TYPE, Dim> Fac;
	typedef std::shared_ptr<Fac> spFac;
	typedef Surface<TYPE, Dim> Sur;
	typedef std::shared_ptr<Sur> spSur;
public:

	vt Distance2(const Poi& p1, const Poi& p2) {
		if (Dim == 2) {
			return (p1[0] - p2[0]) * (p1[0] - p2[0])
					+ (p1[1] - p2[1]) * (p1[1] - p2[1]);
		} else {
			return (p1[0] - p2[0]) * (p1[0] - p2[0])
					+ (p1[1] - p2[1]) * (p1[1] - p2[1])
					+ (p1[2] - p2[2]) * (p1[2] - p2[2]);
		}
	}

	//function out of class
	vt Distance(const Poi& p1, const Poi& p2) {
		if (Dim == 2) {
			return std::sqrt(
					(p1[0] - p2[0]) * (p1[0] - p2[0])
							+ (p1[1] - p2[1]) * (p1[1] - p2[1]));
		} else {
			return std::sqrt(
					(p1[0] - p2[0]) * (p1[0] - p2[0])
							+ (p1[1] - p2[1]) * (p1[1] - p2[1])
							+ (p1[2] - p2[2]) * (p1[2] - p2[2]));
		}
	}

	/**
	 * gts_point_segment_distance2:
	 * @p: a #GtsPoint.
	 * @s: a #GtsSegment.
	 *
	 * Returns: the square of the minimun Euclidean distance between @p and @s.
	 */
	TYPE Distance2( //
			const Poi& p, //
			const Seg& s) {
		TYPE t,ns2,x,y,z;
		const Poi* p1,p2;

		p1 = s->v1;
		p2 = s->v2;
		ns2 = Distance2(*p1, *p2);
		if (ns2 == 0.0)
			return Distance2(*p, *p1);
		t = ((p2->x() - p1->x()) * (p->x() - p1->x())
				+ (p2->y() - p1->y()) * (p->y() - p1->y())
				+ (p2->z() - p1->z()) * (p->z() - p1->z())) / ns2;
		if (t > 1.0)
			return Distance2(*p, *p2);
		if (t < 0.0)
			return Distance2(*p, *p1);
		x = (1. - t) * p1->x() + t * p2->x() - p->x();
		y = (1. - t) * p1->y() + t * p2->y() - p->y();
		z = (1. - t) * p1->z() + t * p2->z() - p->z();
		return x * x + y * y + z * z;
	}

	/**
	 * point_segment_distance:
	 * @p: a #GtsPoint.
	 * @s: a #GtsSegment.
	 *
	 * Returns: the minimun Euclidean distance between @p and @s.
	 */
	vt Distance( //
			const Poi& p, //
			const Seg& s) {
		return sqrt(Distance2(p, s));
	}

	bool IsCCW(const Poi& p1, const Poi& p2, const Poi& p3) {
		vt tmp;
		tmp = (p3.y() - p1.y()) * (p2.x() - p1.x())
				- (p3.x() - p1.x()) * (p2.y() - p1.y());
		if (tmp > 0)
			return true;
		else
			return false;
	}

	bool IsCW(const Poi& p1, const Poi& p2, const Poi& p3) {
		vt tmp;
		tmp = (p3.y() - p1.y()) * (p2.x() - p1.x())
				- (p3.x() - p1.x()) * (p2.y() - p1.y());
		if (tmp > 0)
			return true;
		else
			return false;
	}

};
}

#endif /* _TS_OPERATION_H_ */
