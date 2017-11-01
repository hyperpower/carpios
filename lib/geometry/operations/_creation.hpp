#ifndef _CREATION_HPP_
#define _CREATION_HPP_

#include "../geometry_define.hpp"
#include <array>
#include "../../algebra/arithmetic.hpp"
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
	typedef Segment& ref_Segment;
	typedef const Segment& const_ref_Segment;
	typedef PointChain_<TYPE, DIM> PointChain;

	typedef Polygon_<TYPE> Polygon;
	typedef Contour_<TYPE> Contour;

	typedef TriSurface_<TYPE, DIM> TriSurface;
	typedef TriFace_<TYPE, DIM, TriSurface> TriFace;
	typedef Edge_<TYPE, DIM, TriFace> Edge;
	typedef Vertex_<TYPE, DIM, Edge> Vertex;
	typedef Vertex* pVertex;
	typedef Edge* pEdge;
	typedef TriFace* pTriFace;

	typedef TriSurface_<TYPE, 2> TriSurface2;
	typedef TriSurface_<TYPE, 3> TriSurface3;

	typedef Box_<TYPE, DIM> Box;
	typedef BBox_<TYPE, DIM> BBox;
	typedef BBTree_<BBox> BBTree;

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
			if (pc.size() > 3) {
				// check simple
				if (!(pc.is_simple())) {
					pc.pop_back();
					--i;
				}
			}
			//std::cout<< "i = " << i << "p = "<< p <<std::endl;
		}
	}

	static void Circle(Contour& res, Vt x0, Vt y0, Vt r, int n) {
		ASSERT(Dim == 2);
		ASSERT(n >= 3);
		Float pi = 3.141592653589793238;
		res.clear();
		for (int i = 0; i < n; i++) {
			Vt x = x0 + r * cos(2. * pi / float(n) * i);
			Vt y = y0 + r * sin(2. * pi / float(n) * i);
			res.add(Point(x, y));
		}
	}
	static void Circle(Polygon& res, Vt x0, Vt y0, Vt r, int n) {
		res.clear();
		Contour con;
		Circle(con, x0, y0, r, n);
		res.push_back(con);
	}
	static void Circle(TriSurface& sur, uInt n, Vt r) {
		std::vector<pVertex> v_vertex;
		std::vector<pEdge> v_edge;
		std::vector<pTriFace> v_face;
		pVertex pverc = new Vertex(0, 0, 0);
		double da = 2 * PI / n;
		v_vertex.push_back(pverc);
		for (uInt i = 0; i < n; i++) {
			Vt x = r * cos(i * da);
			Vt y = r * sin(i * da);
			Vt z = 0;
			pVertex pv = new Vertex(x, y, z);
			v_vertex.push_back(pv);
		}
		//edge
		for (uInt i = 0; i < n; i++) {
			pEdge pe = new Edge(v_vertex[0], v_vertex[i + 1]);
			v_edge.push_back(pe);
		}
		for (uInt i = 1; i < n; i++) {
			pEdge pe = new Edge(v_vertex[i], v_vertex[i + 1]);
			v_edge.push_back(pe);
		}
		pEdge pe = new Edge(v_vertex[n], v_vertex[1]);
		v_edge.push_back(pe);
		//surface
		sur.clear();
		for (uInt i = 0; i < n - 1; i++) {
			pTriFace pfac = new TriFace(v_edge[i], v_edge[i + 1], v_edge[n + i],
					&sur);
			sur.insert(pfac);
		}
		pTriFace pfac = new TriFace(v_edge[n - 1], v_edge[0], v_edge[n + n - 1],
				&sur);
		sur.insert(pfac);
	}
	static void TriFaceOne(TriSurface& sur, const Point& a, const Point& b,
			const Point& c) {
		pVertex pv1 = new Vertex(a);
		pVertex pv2 = new Vertex(b);
		pVertex pv3 = new Vertex(c);

		pEdge pe1 = new Edge(pv1, pv2);
		pEdge pe2 = new Edge(pv2, pv3);
		pEdge pe3 = new Edge(pv3, pv1);

		sur.clear();
		pTriFace pfac = new TriFace(pe1, pe2, pe3, &sur);
		sur.insert(pfac);
	}

	static void Cone(TriSurface& sur, uInt n,           //the number of triangle
			const Vt& r, const Vt& zbottom, const Vt& zpex) {
		ASSERT(Dim == 3);
		std::vector<pVertex> v_vertex;
		std::vector<pEdge> v_edge;
		std::vector<pTriFace> v_face;
		pVertex pverc = new Vertex(0, 0, zpex);
		double da = 2 * PI / n;
		v_vertex.push_back(pverc);
		for (uInt i = 0; i < n; i++) {
			Vt x = r * std::cos(i * da);
			Vt y = r * std::sin(i * da);
			Vt z = zbottom;
			pVertex pv = new Vertex(x, y, z);
			v_vertex.push_back(pv);
		}
		//edge
		for (uInt i = 0; i < n; i++) {
			pEdge pe = new Edge(v_vertex[0], v_vertex[i + 1]);
			v_edge.push_back(pe);
		}
		for (uInt i = 1; i < n; i++) {
			pEdge pe = new Edge(v_vertex[i], v_vertex[i + 1]);
			v_edge.push_back(pe);
		}
		pEdge pe = new Edge(v_vertex[n], v_vertex[1]);
		v_edge.push_back(pe);
		//surface
		sur.clear();
		for (uInt i = 0; i < n - 1; i++) {
			pTriFace pfac = new TriFace(v_edge[i], v_edge[i + 1], v_edge[n + i],
					&sur);
			sur.insert(pfac);
		}
		pTriFace pfac = new TriFace(v_edge[n - 1], v_edge[0], v_edge[n + n - 1],
				&sur);
		sur.insert(pfac);
	}

	static void BoundingBoxTree(BBTree& tree, const TriSurface& sur) {
		std::set<BBox> set_box;
		int i = 0;
		for (auto& pf : sur) {
//			auto f = (*pf);    //pFac
			BBox tbox(pf);     //pTri
			set_box.insert(tbox);
			i++;
		}
		tree.rebuild(set_box);
	}

	static std::shared_ptr<BBTree> BoundingBoxTree(const TriSurface& sur) {
		std::shared_ptr<BBTree> sptree(new BBTree());
		BoundingBoxTree(*sptree, sur);
		return sptree;
	}

};
}

#endif

