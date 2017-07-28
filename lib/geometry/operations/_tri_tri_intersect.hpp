#ifndef _TRI_TRI_INTERSECT_H_
#define _TRI_TRI_INTERSECT_H_

#include "../geometry_define.hpp"

#include "_tri_tri_intersect_97.hpp"
#include "_tri_tri_intersect_02.hpp"
#include "../objects/_triangle.hpp"

namespace carpio {

// * parameters: vertices of triangle 1: V0,V1,V2
// *             vertices of triangle 2: U0,U1,U2
// * result    : returns 1 if the triangles intersect, otherwise 0

template<class VALUE>
int TriTriIsect_Raw_97( //
		VALUE V0[3], VALUE V1[3], VALUE V2[3], //
		VALUE U0[3], VALUE U1[3], VALUE U2[3]) //
		{
	return TriTriIsect97::NoDivTriTriIsect(V0, V1, V2, U0, U1, U2);
}

template<class VALUE>
int TriTriIsect_Raw_02( //
		VALUE V0[3], VALUE V1[3], VALUE V2[3], //
		VALUE U0[3], VALUE U1[3], VALUE U2[3]) //
		{
	return TriTriIsect02::tri_tri_overlap_test_3d(V0, V1, V2, U0, U1, U2);
}

template<class VALUE>
int TriTriIsect_97(const Triangle_<VALUE, 3>& tri1,
		const Triangle_<VALUE, 3>& tri2) {
	if (tri1 == tri2) {
		return 0;  // equal is not intersect
	}
	VALUE V0[3], V1[3], V2[3];
	VALUE U0[3], U1[3], U2[3];
	for (typename Triangle_<VALUE, 3>::size_type i = 0; i < 3; i++) {
		V0[i] = tri1.p(0, ToAxes(i));
		V1[i] = tri1.p(1, ToAxes(i));
		V2[i] = tri1.p(2, ToAxes(i));

		U0[i] = tri2.p(0, ToAxes(i));
		U1[i] = tri2.p(1, ToAxes(i));
		U2[i] = tri2.p(2, ToAxes(i));
	}
	return TriTriIsect_Raw_97(V0, V1, V2, U0, U1, U2);
}

template<class VALUE>
int TriTriIsect_02(const Triangle_<VALUE, 3>& tri1,
		const Triangle_<VALUE, 3>& tri2) {
	VALUE V0[3], V1[3], V2[3];
	VALUE U0[3], U1[3], U2[3];
	if (tri1 == tri2) {
		return 0;  // equal is not intersect
	}
	for (typename Triangle_<VALUE, 3>::size_type i = 0; i < 3; i++) {
		V0[i] = tri1.p(0, ToAxes(i));
		V1[i] = tri1.p(1, ToAxes(i));
		V2[i] = tri1.p(2, ToAxes(i));

		U0[i] = tri2.p(0, ToAxes(i));
		U1[i] = tri2.p(1, ToAxes(i));
		U2[i] = tri2.p(2, ToAxes(i));
	}
	return TriTriIsect_Raw_02(V0, V1, V2, U0, U1, U2);
}

template<class VALUE>
int TriTriIsect_Segment(const Triangle_<VALUE, 3>& tri1,
		const Triangle_<VALUE, 3>& tri2, bool& coplane, Point_<VALUE, 3>& ps,
		Point_<VALUE, 3>& pe) {
	VALUE V0[3], V1[3], V2[3];
	VALUE U0[3], U1[3], U2[3];
	if (tri1 == tri2) {
		return 0;  // equal is not intersect
	}
	for (typename Triangle_<VALUE, 3>::size_type i = 0; i < 3; i++) {
		V0[i] = tri1.p(0, ToAxes(i));
		V1[i] = tri1.p(1, ToAxes(i));
		V2[i] = tri1.p(2, ToAxes(i));

		U0[i] = tri2.p(0, ToAxes(i));
		U1[i] = tri2.p(1, ToAxes(i));
		U2[i] = tri2.p(2, ToAxes(i));
	}
	VALUE start[3], end[3];
	int res = TriTriIsect02::tri_tri_intersection_test_3d(V0, V1, V2, U0, U1,
			U2, coplane, start, end);
	for (int i = 0; i < 3; i++) {
		ps[i] = start[i];
		pe[i] = end[i];
	}
	return res;
}

template<class VALUE>
int TriBoxIsect_Raw(VALUE boxcenter[3],    //
		VALUE boxhalfsize[3],  //
		VALUE triverts[3][3])  //
		{
	return TriTriIsect97::triBoxOverlap(boxcenter, boxhalfsize, triverts);
}

}

#endif
