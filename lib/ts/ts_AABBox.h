/************************
 //  \file   ts_AABBox.h
 //  \brief
 // 
 //  \author czhou
 //  \date   16 juin 2015 
 ***********************/
#ifndef TS_AABBOX_H_
#define TS_AABBOX_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_face.h"
#include "geometry/_tri_tri_intersect.hpp"
//#include "ts_tri_moller.h"

namespace TS {

//VTK_VOXEL (=11)  -- 3D
//VTK_PIXEL (=8)   -- 2D
const Location _ORDER_VTK[8][3] = { { _M, _M, _M }, //
		{ _P, _M, _M }, //
		{ _M, _P, _M }, //
		{ _P, _P, _M }, //
		{ _M, _M, _P }, //
		{ _P, _M, _P }, //
		{ _M, _P, _P }, //
		{ _P, _P, _P } };

template<class TYPE, st DIM>
class AABBox: public Array<Point<TYPE, DIM>, 2> {
public:
	typedef AABBox<TYPE, DIM> Self;
	typedef const AABBox<TYPE, DIM> const_Self;
	typedef st size_type;
	typedef TYPE value_type;
	typedef Point<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Tri* pTri;
	typedef Face<TYPE, DIM> Fac;
	typedef Fac* pFac;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Ver* pVer;
	typedef Segment<TYPE, DIM> Seg;
	typedef Seg* pSeg;

	static const st Dim = DIM;
	static const size_t NumFaces = DIM + DIM;
	static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
protected:
	static const st Max = 1;
	static const st Min = 0;
	static const st X = 0;
	static const st Y = 1;
	static const st Z = 2;

public:
	utPointer bounded;
	OBJ_TYPE type;
public:
	AABBox() {
		reset(EMPTY, nullptr, 0, 0, 0, 0, 0, 0);
	}
	AABBox(OBJ_TYPE t, utPointer bounded,  //
			Float x1, Float y1, Float z1, //
			Float x2, Float y2, Float z2) //
			{
		reset(t, bounded, x1, y1, z1, x2, y2, z2);
	}

	AABBox(pTri bounded) {
		reset(bounded);
	}
	AABBox(pFac bounded) {
		reset(bounded);
	}

	AABBox(pPoi bounded) {
		reset(POINT, bounded);
	}

	AABBox(pSeg bounded) {
		reset(SEGMENT, bounded);
	}

	AABBox(const Self& a, const Self& b) { // union to one box without data
		Self& self = (*this);
		for (st i = 0; i < DIM; ++i) {
			self[Min][i] = MIN(a[Min][i], b[Min][i]);
			self[Max][i] = MAX(a[Max][i], b[Max][i]);
		}
		self.bounded = nullptr;
		self.type = EMPTY;
	}

	AABBox(const Self& a) {
		Self& self = (*this);
		self[Min] = a[Min];
		self[Max] = a[Max];
		self.bounded = a.bounded;
		self.type = a.type;
	}

	Self& operator=(const Self& a) {
		if (this == &a) {
			return *this;
		}
		Self& self = (*this);
		self[Min] = a[Min];
		self[Max] = a[Max];
		self.bounded = a.bounded;
		self.type = a.type;
		return self;
	}

	bool operator<(const Self& rhs) const {
		//compare center point x -> y -> z -> point address
		const Self& self = (*this);
		if ((self[Min][X] + self[Max][X]) / 2.0
				< (rhs[Min][X] + rhs[Max][X]) / 2.0) {
			return true;
		} else if ((self[Min][X] + self[Max][X]) / 2.0
				== (rhs[Min][X] + rhs[Max][X]) / 2.0) {
			if ((self[Min][Y] + self[Max][Y]) / 2.0
					< (rhs[Min][Y] + rhs[Max][Y]) / 2.0) {
				return true;
			} else if ((self[Min][Y] + self[Max][Y]) / 2.0
					== (rhs[Min][Y] + rhs[Max][Y]) / 2.0) {
				if (DIM == 3) {
					if ((self[Min][Z] + self[Max][Z]) / 2.0
							< (rhs[Min][Z] + rhs[Max][Z]) / 2.0) {
						return true;
					}
				}
			}
		}
		return (this < &rhs);
	}


	value_type get(Aix aix, Location loc) const {
		if (aix == _Z && DIM == 2) {
			return 0;
		}
		value_type xv = 0;
		const Self& self = (*this);
		switch (loc) {
		case _M: {
			xv = self[Min][aix];
			break;
		}
		case _C: {
			xv = (self[Min][aix] + self[Max][aix]) * 0.5;
			break;
		}
		case _P: {
			xv = self[Max][aix];
			break;
		}
		default: {
			break;
		}
		}
		return xv;
	}

	value_type get(Location loc, Aix aix) const {
		return get(aix, loc);
	}

	value_type get_d(Aix aix) const {
		const Self& self = (*this);
		return self[Max][aix] - self[Min][aix];
	}

	value_type get_dh(Aix aix) const {
		const Self& self = (*this);
		return (self[Max][aix] - self[Min][aix]) / 2.0;
	}

	Poi get_point(Location lx, Location ly, Location lz) const {
		return Poi(get(_X, lx), get(_Y, ly), get(_Z, lz));
	}

	bool empty() const {
		Self& self = (*this);
		_return_val_if_fail(this->bounded == nullptr, false);
		_return_val_if_fail(self[Min][X] == 0, false);
		_return_val_if_fail(self[Min][Y] == 0, false);
		_return_val_if_fail(self[Max][X] == 0, false);
		_return_val_if_fail(self[Max][Y] == 0, false);
		if (DIM == 3) {
			_return_val_if_fail(self[Min][Z] == 0, false);
			_return_val_if_fail(self[Max][Z] == 0, false);
		}
		return true;
	}
	void reset(OBJ_TYPE t, utPointer bounded,  //
				Float x1, Float y1, Float z1,  // min
				Float x2, Float y2, Float z2); // max
	void reset(pTri bounded);
	void reset(pFac bounded);
	void reset(OBJ_TYPE t, pPoi bounded);
	void reset(OBJ_TYPE t, pSeg bounded);

	bool are_overlapping(const AABBox<TYPE, DIM>& other) const;
	bool are_overlapping(pTri t) const;

	utPointer utp_obj(){
		return this->bounded;
	}
	const utPointer utp_obj() const{
		return this->bounded;
	}
	// show ===========================

	// output =========================
	void output_vtk(const String& fn) const;

}
;

template<class TYPE, st DIM>
void AABBox<TYPE, DIM>::reset(OBJ_TYPE t, utPointer bounded, Float x1, Float y1,
		Float z1, Float x2, Float y2, Float z2) {
	ASSERT(x2 >= x1 && y2 >= y1 && z2 >= z1);
	type = t;
	Self& self = (*this);
	self[Min][X] = x1;
	self[Min][Y] = y1;
	self[Max][X] = x2;
	self[Max][Y] = y2;
	if (DIM == 3) {
		self[Min][Z] = z1;
		self[Max][Z] = z2;
	}
	this->bounded = bounded;
}

template<class TYPE, st DIM>
void AABBox<TYPE, DIM>::reset(OBJ_TYPE t, pPoi pp) {
	if (DIM == 3) {
		reset(t, pp, pp->x(), pp->y(), pp->z(), pp->x(), pp->y(), pp->z());
	} else {
		reset(t, pp, pp->x(), pp->y(), 0, pp->x(), pp->y(), 0);
	}
}
template<class TYPE, st DIM>
void AABBox<TYPE, DIM>::reset(pTri pt) {
	Self& self = (*this);
	OBJ_TYPE t = TRIANGLE;
	Ver& v = (*(pt->get_vertex1()));
	if (DIM == 3) {
		reset(t, pt, v[X], v[Y], v[Z], v[X], v[Y], v[Z]);
	} else {
		reset(t, pt, v[X], v[Y], 0, v[X], v[Y], 0);
	}
	Ver& v2 = (*pt->get_vertex2());
	for (st i = 0; i < DIM; ++i) {
		if (v2[i] > self[Max][i]) {
			self[Max][i] = v2[i];
		}
		if (v2[i] < self[Min][i]) {
			self[Min][i] = v2[i];
		}
	}
	Ver& v3 = (*pt->get_vertex3());
	for (st i = 0; i < DIM; ++i) {
		if (v3[i] > self[Max][i]) {
			self[Max][i] = v3[i];
		}
		if (v3[i] < self[Min][i]) {
			self[Min][i] = v3[i];
		}
	}
}
template<class TYPE, st DIM>
void AABBox<TYPE, DIM>::reset(pFac pt) {
	Self& self = (*this);
	OBJ_TYPE t = FACE;
	Ver& v = (*(pt->get_vertex1()));
	if (DIM == 3) {
		reset(t, pt, v[X], v[Y], v[Z], v[X], v[Y], v[Z]);
	} else {
		reset(t, pt, v[X], v[Y], 0, v[X], v[Y], 0);
	}
	Ver& v2 = (*pt->get_vertex2());
	for (st i = 0; i < DIM; ++i) {
		if (v2[i] > self[Max][i]) {
			self[Max][i] = v2[i];
		}
		if (v2[i] < self[Min][i]) {
			self[Min][i] = v2[i];
		}
	}
	Ver& v3 = (*pt->get_vertex3());
	for (st i = 0; i < DIM; ++i) {
		if (v3[i] > self[Max][i]) {
			self[Max][i] = v3[i];
		}
		if (v3[i] < self[Min][i]) {
			self[Min][i] = v3[i];
		}
	}
}
template<class TYPE, st DIM>
void AABBox<TYPE, DIM>::reset(OBJ_TYPE t, pSeg s) {

	_return_if_fail(s != nullptr);

	Self& self = (*this);
	Ver& ver1 = (*(s->v1));
	Ver& ver2 = (*(s->v2));

	for (st i = 0; i < DIM; ++i) {
		if (ver1[i] > ver2[i]) {
			self[Max][i] = ver1[i];
			self[Min][i] = ver2[i];
		} else {
			self[Max][i] = ver2[i];
			self[Min][i] = ver1[i];
		}
	}
}
template<class TYPE, st DIM>
bool AABBox<TYPE, DIM>::are_overlapping(const AABBox<TYPE, DIM>& other) const {
	//assert(other != nullptr);

	if (this == &other) {
		return true;
	}
	const AABBox<TYPE, DIM> &self = (*this);
	const AABBox<TYPE, DIM> &othe = other;
	for (size_type i = 0; i < DIM; ++i) {
		if (othe[Min][i] > self[Max][i] || othe[Max][i] < self[Min][i]) {
			return false;
		}
	}
	return true;
}

template<class TYPE, st DIM>
void AABBox<TYPE, DIM>::output_vtk(const String& fn) const {
	//uInt NUM_VERTEXES = (DIM == 3) ? 8 : 4;
	FILE *data;
	data = fopen(fn.c_str(), "w");
	if (data == NULL) {
		std::cerr << "!> Open file error! " << fn << " \n";
		exit(-1);
	}
	fprintf(data, "# vtk DataFile Version 3.0\n");
	fprintf(data, "Gird output\n");
	fprintf(data, "ASCII\n");
	fprintf(data, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(data, "POINTS %d float\n", NumVertexes);
	for (uInt i = 0; i < NumVertexes; i++) {
		fprintf(data, "%f %f %f \n",
				get_point(_ORDER_VTK[i][0], _ORDER_VTK[i][1], _ORDER_VTK[i][2]).x(),
				get_point(_ORDER_VTK[i][0], _ORDER_VTK[i][1], _ORDER_VTK[i][2]).y(),
				DIM == 3 ?
						get_point(_ORDER_VTK[i][0], _ORDER_VTK[i][1],
								_ORDER_VTK[i][2]).z() :
						0);
	}
	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", 1, NumVertexes + 1);
	fprintf(data, "%d ", NumVertexes);
	for (int i = 0; i < NumVertexes; i++) {
		fprintf(data, "%d ", i);
	}
	fprintf(data, "\n\n");
	fprintf(data, "CELL_TYPES %d\n", 1);
	fprintf(data, "%d \n", DIM == 3 ? 11 : 8);
	//VTK_VOXEL (=11)  -- 3D
	//VTK_PIXEL (=8)   -- 2D
	fclose(data);
}

template<class TYPE, st DIM>
bool AABBox<TYPE, DIM>::are_overlapping(pTri t) const {
	assert(DIM == 3);
	TYPE bc[3], bh[3], tv[3][3];
	//pVer p1;
	//pVer p2;
	//pVer p3;

	_return_val_if_fail(t != NULL, false);

	bc[0] = this->get(_X, _C);
	bc[1] = this->get(_Y, _C);
	bc[2] = this->get(_Z, _C);
	bh[0] = this->get_dh(_X);
	bh[1] = this->get_dh(_Y);
	bh[2] = this->get_dh(_Z);
	//std::cout<<"bc "<< bc[0] <<" "<< bc[1]<<" "<<bc[2]<<std::endl;
	//std::cout<<"bh "<<bh[0] <<" "<<bh[1]<<" "<<bh[2]<<std::endl;
	auto p1 = t->get_vertex1();
	auto p2 = t->get_vertex2();
	auto p3 = t->get_vertex3();
	//p1->show();
	//p2->show();
	//p3->show();
	for (uInt i = 0; i < DIM; i++) {
		tv[0][i] = (*p1)[i];
		tv[1][i] = (*p2)[i];
		tv[2][i] = (*p3)[i];
	}
	return carpio::TriBoxIsect_Raw(bc, bh, tv);
	//return triBoxOverlap(bc, bh, tv);
}

}

#endif /* TS_AABBOX_H_ */
