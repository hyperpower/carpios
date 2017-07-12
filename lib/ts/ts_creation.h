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

template<class TYPE, st DIM>
class Creation {
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
	static spSur Circle(uInt n, vt r) {
		if (Dim == 2) {
			return Circle2(n, r);
		} else if (Dim == 3) {
			return Circle3(n, r);
		}
		SHOULD_NOT_REACH;
		return nullptr;
	}

	static spSur Circle2(uInt n, vt r) {
		Vector<spVer> v_vertex;
		Vector<spEdg> v_edge;
		Vector<spFac> v_face;
		spVer pverc(new Ver(0, 0, 0));
		double da = 2 * PI / n;
		v_vertex.push_back(pverc);
		for (uInt i = 0; i < n; i++) {
			vt x = r * cos(i * da);
			vt y = r * sin(i * da);
			vt z = 0;
			spVer pv(new Ver(x, y, z));
			v_vertex.push_back(pv);
		}
		//edge
		for (uInt i = 0; i < n; i++) {
			spEdg pe(new Edg(v_vertex[0], v_vertex[i + 1]));
			pe->attach();
			v_edge.push_back(pe);
		}
		for (uInt i = 1; i < n; i++) {
			spEdg pe(new Edg(v_vertex[i], v_vertex[i + 1]));
			pe->attach();
			v_edge.push_back(pe);
		}
		spEdg pe(new Edg(v_vertex[n], v_vertex[1]));
		pe->attach();
		v_edge.push_back(pe);
		//surface
		spSur sur(new Sur());
		for (uInt i = 0; i < n - 1; i++) {
			spFac pfac(
					new Fac(v_edge[i], v_edge[i + 1], v_edge[n + i],
							sur.get()));
			pfac->attach();
			sur->faces.insert(pfac);
		}
		spFac pfac(
				new Fac(v_edge[n - 1], v_edge[0], v_edge[n + n - 1],
						sur.get()));
		pfac->attach();
		sur->faces.insert(pfac);
		return sur;
	}

	static spSur Circle3(uInt n, vt r) {
		Vector<spVer> v_vertex;
		Vector<spEdg> v_edge;
		Vector<spFac> v_face;
		spVer pverc(new Ver(0, 0, 0));
		double da = 2 * PI / n;
		v_vertex.push_back(pverc);
		for (uInt i = 0; i < n; i++) {
			vt x = r * cos(i * da);
			vt y = r * sin(i * da);
			vt z = 0;
			spVer pv(new Ver(x, y, z));
			v_vertex.push_back(pv);
		}
		//edge
		for (uInt i = 0; i < n; i++) {
			spEdg pe(new Edg(v_vertex[0], v_vertex[i + 1]));
			pe->attach();
			v_edge.push_back(pe);
		}
		for (uInt i = 1; i < n; i++) {
			spEdg pe(new Edg(v_vertex[i], v_vertex[i + 1]));
			pe->attach();
			v_edge.push_back(pe);
		}
		spEdg pe(new Edg(v_vertex[n], v_vertex[1]));
		pe->attach();
		v_edge.push_back(pe);
		//surface
		spSur sur(new Sur());
		for (uInt i = 0; i < n - 1; i++) {
			spFac pfac(
					new Fac(v_edge[i], v_edge[i + 1], v_edge[n + i],
							sur.get()));
			pfac->attach();
			sur->faces.insert(pfac);
		}
		spFac pfac(
				new Fac(v_edge[n - 1], v_edge[0], v_edge[n + n - 1],
						sur.get()));
		pfac->attach();
		sur->faces.insert(pfac);
		return sur;
	}

	static spSur TriangleOne(const Poi& p1, const Poi& p2, const Poi& p3) {
		spVer pv1(new Ver(p1.x(), p1.y(), p1.z()));
		spVer pv2(new Ver(p2.x(), p2.y(), p2.z()));
		spVer pv3(new Ver(p3.x(), p3.y(), p3.z()));

		spEdg pe1(new Edg(pv1, pv2));
		pe1->attach();
		spEdg pe2(new Edg(pv2, pv3));
		pe2->attach();
		spEdg pe3(new Edg(pv3, pv1));
		pe3->attach();
		//surface
		spSur sur(new Sur());
		spFac pfac(new Fac(pe1, pe2, pe3, sur.get()));
		pfac->attach();
		sur.faces.insert(pfac);
		return sur;
	}

	static spSur FromGtsFile(const std::string& filename) {
		ASSERT(Dim ==3);
		spSur psur(new Sur());
		psur->load_gts_file(filename);
		return psur;
	}

	static spSur Cone(uInt n,                           //the number of triangle
			const vt& r, const vt& zbottom, const vt& zpex) {
		ASSERT(Dim ==3);
		Vector<spVer> v_vertex;
		Vector<spEdg> v_edge;
		Vector<spFac> v_face;
		spVer pverc(new Ver(0, 0, zpex));
		double da = 2 * PI / n;
		v_vertex.push_back(pverc);
		for (uInt i = 0; i < n; i++) {
			vt x = r * std::cos(i * da);
			vt y = r * sin(i * da);
			vt z = zbottom;
			spVer pv(new Ver(x, y, z));
			v_vertex.push_back(pv);
		}
		//edge
		for (uInt i = 0; i < n; i++) {
			spEdg pe(new Edg(v_vertex[0], v_vertex[i + 1]));
			pe->attach();
			v_edge.push_back(pe);
		}
		for (uInt i = 1; i < n; i++) {
			spEdg pe(new Edg(v_vertex[i], v_vertex[i + 1]));
			pe->attach();
			v_edge.push_back(pe);
		}
		spEdg pe(new Edg(v_vertex[n], v_vertex[1]));
		pe->attach();
		v_edge.push_back(pe);
		//surface
		spSur sur(new Sur);
		for (uInt i = 0; i < n - 1; i++) {
			spFac pfac(
					new Fac(v_edge[i], v_edge[i + 1], v_edge[n + i],
							sur.get()));
			pfac->attach();
			sur->faces.insert(pfac);
		}
		spFac pfac(
				new Fac(v_edge[n - 1], v_edge[0], v_edge[n + n - 1],
						sur.get()));
		pfac->attach();
		sur->faces.insert(pfac);
		return sur;
	}
};

}

#endif /* TS_SURFACE_CONSTRUCTOR_H_ */
