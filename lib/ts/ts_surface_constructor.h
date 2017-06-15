/************************
 //  \file   ts_surface_constructor.h
 //  \brief
 // 
 //  \author czhou
 //  \date   26 juin 2015 
 ***********************/
#ifndef TS_SURFACE_CONSTRUCTOR_H_
#define TS_SURFACE_CONSTRUCTOR_H_

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
template<class TYPE>
std::shared_ptr<Surface<TYPE, 2> > ConstructCircle2(          // the surface
		uInt n,                           //the number of triangle
		const TYPE& r) {
	typedef Point<TYPE, 2> Poi;
	typedef std::shared_ptr<Poi> pPoi;
	typedef Vertex<TYPE, 2> Ver;
	typedef std::shared_ptr<Ver> pVer;
	typedef Segment<TYPE, 2> Seg;
	typedef std::shared_ptr<Seg> pSeg;
	typedef Edge<TYPE, 2> Edg;
	typedef std::shared_ptr<Edg> pEdg;
	typedef Triangle<TYPE, 2> Tri;
	typedef std::shared_ptr<Tri> pTri;
	typedef Face<TYPE, 2> Fac;
	typedef std::shared_ptr<Fac> pFac;
	typedef Surface<TYPE, 2> Sur;
	typedef std::shared_ptr<Sur> pSur;

	Vector<pVer> v_vertex;
	Vector<pEdg> v_edge;
	Vector<pFac> v_face;
	pVer pverc(new Ver(0, 0, 0));
	double da = 2 * PI / n;
	v_vertex.push_back(pverc);
	for (uInt i = 0; i < n; i++) {
		TYPE x = r * cos(i * da);
		TYPE y = r * sin(i * da);
		TYPE z = 0;
		pVer pv(new Ver(x, y, z));
		v_vertex.push_back(pv);
	}
	//edge
	for (uInt i = 0; i < n; i++) {
		pEdg pe(new Edg(v_vertex[0], v_vertex[i + 1]));
		pe->attach();
		v_edge.push_back(pe);
	}
	for (uInt i = 1; i < n; i++) {
		pEdg pe(new Edg(v_vertex[i], v_vertex[i + 1]));
		pe->attach();
		v_edge.push_back(pe);
	}
	pEdg pe(new Edg(v_vertex[n], v_vertex[1]));
	pe->attach();
	v_edge.push_back(pe);
	//surface
	pSur sur(new Sur());
	for (uInt i = 0; i < n - 1; i++) {
		pFac pfac(new Fac(v_edge[i], v_edge[i + 1], v_edge[n + i], sur.get()));
		pfac->attach();
		sur->faces.insert(pfac);
	}
	pFac pfac(new Fac(v_edge[n - 1], v_edge[0], v_edge[n + n - 1], sur.get()));
	pfac->attach();
	sur->faces.insert(pfac);
	//for (pVer pv : v_vertex) {
	//	sur->c_vertex.insert(pv);
	//}
	//for (pEdg pe : v_edge) {
	//	sur->c_edge.insert(pe);
	//}
	return sur;
}

template<class TYPE>
std::shared_ptr<Surface<TYPE, 3> > ConstructCircle3(          // the surface
		uInt n,                           //the number of triangle
		const TYPE& r) {
	typedef Point<TYPE, 3> Poi;
	typedef std::shared_ptr<Poi> pPoi;
	typedef Vertex<TYPE, 3> Ver;
	typedef std::shared_ptr<Ver> pVer;
	typedef Segment<TYPE, 3> Seg;
	typedef std::shared_ptr<Seg> pSeg;
	typedef Edge<TYPE, 3> Edg;
	typedef std::shared_ptr<Edg> pEdg;
	typedef Triangle<TYPE, 3> Tri;
	typedef std::shared_ptr<Tri> pTri;
	typedef Face<TYPE, 3> Fac;
	typedef std::shared_ptr<Fac> pFac;
	typedef Surface<TYPE, 3> Sur;
	typedef std::shared_ptr<Sur> pSur;

	Vector<pVer> v_vertex;
	Vector<pEdg> v_edge;
	Vector<pFac> v_face;
	pVer pverc(new Ver(0, 0, 0));
	double da = 2 * PI / n;
	v_vertex.push_back(pverc);
	for (uInt i = 0; i < n; i++) {
		TYPE x = r * cos(i * da);
		TYPE y = r * sin(i * da);
		TYPE z = 0;
		pVer pv(new Ver(x, y, z));
		v_vertex.push_back(pv);
	}
	//edge
	for (uInt i = 0; i < n; i++) {
		pEdg pe(new Edg(v_vertex[0], v_vertex[i + 1]));
		pe->attach();
		v_edge.push_back(pe);
	}
	for (uInt i = 1; i < n; i++) {
		pEdg pe(new Edg(v_vertex[i], v_vertex[i + 1]));
		pe->attach();
		v_edge.push_back(pe);
	}
	pEdg pe(new Edg(v_vertex[n], v_vertex[1]));
	pe->attach();
	v_edge.push_back(pe);
	//surface
	pSur sur(new Sur());
	for (uInt i = 0; i < n - 1; i++) {
		pFac pfac(new Fac(v_edge[i], v_edge[i + 1], v_edge[n + i], sur.get()));
		pfac->attach();
		sur->faces.insert(pfac);
	}
	pFac pfac(new Fac(v_edge[n - 1], v_edge[0], v_edge[n + n - 1], sur.get()));
	pfac->attach();
	sur->faces.insert(pfac);
	//for (pVer pv : v_vertex) {
	//	sur->c_vertex.insert(pv);
	//}
	//for (pEdg pe : v_edge) {
	//	sur->c_edge.insert(pe);
	//}
	return sur;
}

template<class TYPE, st DIM>
void ConstructTriangle(
		Surface<TYPE, DIM>& sur,          // the surface
		const Point<TYPE, DIM> p1, const Point<TYPE, DIM> p2,
		const Point<TYPE, DIM> p3) {
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
	typedef Face<TYPE, DIM> Fac;
	typedef std::shared_ptr<Fac> pFac;
	typedef Surface<TYPE, DIM> Sur;
	typedef std::shared_ptr<Sur> pSur;

	Vertex<TYPE, DIM>* pv1 = new Vertex<TYPE, DIM>(p1.x(), p1.y(), p1.z());
	Vertex<TYPE, DIM>* pv2 = new Vertex<TYPE, DIM>(p2.x(), p2.y(), p2.z());
	Vertex<TYPE, DIM>* pv3 = new Vertex<TYPE, DIM>(p3.x(), p3.y(), p3.z());

	pEdg pe1 = new Edg(pv1, pv2);
	pEdg pe2 = new Edg(pv2, pv3);
	pEdg pe3 = new Edg(pv3, pv1);

	sur.clear();
	pFac pfac = new Fac(pe1, pe2, pe3, &sur);
	sur.faces.insert(pfac);
	sur.c_vertex.insert(pv1);
	sur.c_vertex.insert(pv2);
	sur.c_vertex.insert(pv3);
	sur.c_edge.insert(pe1);
	sur.c_edge.insert(pe2);
	sur.c_edge.insert(pe3);
}
template<class TYPE>
std::shared_ptr<Surface<TYPE, 3> > ConstructCone( //
		uInt n,                           //the number of triangle
		const TYPE& r, const TYPE& zbottom, const TYPE& zpex) {
	typedef Point<TYPE, 3> Poi;
	typedef std::shared_ptr<Poi> pPoi;
	typedef Vertex<TYPE, 3> Ver;
	typedef std::shared_ptr<Ver> pVer;
	typedef Segment<TYPE, 3> Seg;
	typedef std::shared_ptr<Seg> pSeg;
	typedef Edge<TYPE, 3> Edg;
	typedef std::shared_ptr<Edg> pEdg;
	typedef Triangle<TYPE, 3> Tri;
	typedef std::shared_ptr<Tri> pTri;
	typedef Face<TYPE, 3> Fac;
	typedef std::shared_ptr<Fac> pFac;
	typedef Surface<TYPE, 3> Sur;
	typedef std::shared_ptr<Sur> pSur;

	Vector<pVer> v_vertex;
	Vector<pEdg> v_edge;
	Vector<pFac> v_face;
	pVer pverc(new Ver(0, 0, zpex));
	double da = 2 * PI / n;
	v_vertex.push_back(pverc);
	for (uInt i = 0; i < n; i++) {
		TYPE x = r * cos(i * da);
		TYPE y = r * sin(i * da);
		TYPE z = zbottom;
		pVer pv(new Ver(x, y, z));
		v_vertex.push_back(pv);
	}
	//edge
	for (uInt i = 0; i < n; i++) {
		pEdg pe(new Edg(v_vertex[0], v_vertex[i + 1]));
		pe->attach();
		v_edge.push_back(pe);
	}
	for (uInt i = 1; i < n; i++) {
		pEdg pe(new Edg(v_vertex[i], v_vertex[i + 1]));
		pe->attach();
		v_edge.push_back(pe);
	}
	pEdg pe(new Edg(v_vertex[n], v_vertex[1]));
	pe->attach();
	v_edge.push_back(pe);
	//surface
	pSur sur(new Sur);
	for (uInt i = 0; i < n - 1; i++) {
		pFac pfac(new Fac(v_edge[i], v_edge[i + 1], v_edge[n + i], sur.get()));
		pfac->attach();
		sur->faces.insert(pfac);
	}
	pFac pfac(new Fac(v_edge[n - 1], v_edge[0], v_edge[n + n - 1], sur.get()));
	pfac->attach();
	sur->faces.insert(pfac);
	//for (pVer pv : v_vertex) {
	//	sur.c_vertex.insert(pv);
	//}
	//for (pEdg pe : v_edge) {
	//	sur.c_edge.insert(pe);
	//}
	return sur;
}

template<class TYPE>
std::shared_ptr<Surface<TYPE, 3> > ConstructFromGtsFile( //
		const std::string& filename, TYPE dummy = 0.0) {
	typedef TS::Surface<TYPE, 3> Sur;
	typedef std::shared_ptr<TS::Surface<TYPE, 3> > pSur;
	pSur psur(new Sur());
	psur->load_gts_file(filename);
	return psur;
}

}

#endif /* TS_SURFACE_CONSTRUCTOR_H_ */
