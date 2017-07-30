//Copyright (C) 2011 by Ivan Fratric
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in
//all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//THE SOFTWARE.

#ifndef _POLYGON_PARTITION_HPP_
#define _POLYGON_PARTITION_HPP_

#include "../geometry_define.hpp"

#include "../objects/_objects.hpp"

#include <list>
#include <set>

namespace carpio {

typedef double tppl_float;

#define TPPL_CCW 1
#define TPPL_CW -1

template<class TYPE, St DIM>
class PolygonPartition_ {
public:
	static const St Dim = DIM;
	typedef TYPE Vt;
	typedef Polygon_<TYPE> Polygon;
	typedef Contour_<TYPE> Contour;
	typedef Point_<TYPE, Dim> Point;
	typedef PointChain_<TYPE, Dim> PointChain;
	typedef Segment_<TYPE, Dim> Segment;
	typedef Operation_<TYPE, Dim> Op;
	typedef Connector_<TYPE, Dim> Connector;

public:
	struct PartitionVertex {
		bool isActive;
		bool isConvex;
		bool isEar;

		Point p;
		tppl_float angle;
		PartitionVertex *previous;
		PartitionVertex *next;

		PartitionVertex() :
				previous(nullptr), next(nullptr) {
		}
	};

	struct MonotoneVertex {
		Point p;
		long previous;
		long next;
	};

	class VertexSorter {
		MonotoneVertex *vertices;
		public:
		VertexSorter(MonotoneVertex *v) :
				vertices(v) {
		}
		bool operator()(long index1, long index2);
	};

	struct Diagonal {
		long index1;
		long index2;
	};

	//dynamic programming state for minimum-weight triangulation
	struct DPState {
		bool visible;
		tppl_float weight;
		long bestvertex;
	};

	//dynamic programming state for convex partitioning
	struct DPState2 {
		bool visible;
		long weight;
		std::list<Diagonal> pairs;
	};

	//edge that intersects the scanline
	struct ScanLineEdge {
		mutable long index;
		Point p1;
		Point p2;

		//determines if the edge is to the left of another edge
		bool operator<(const ScanLineEdge & other) const;

		bool IsConvex(const Point& p1, const Point& p2, const Point& p3) const;
	};

	//standard helper functions
	// ------------------|
	static bool IsConvex(const Point& p1, const Point& p2, const Point& p3) { // this function should be change to CCW
		// p0 = p3
		// p1 = p2
		return Op::IsCCW(p1, p2, p3);
//		tppl_float tmp;
//		tmp = (p3.y() - p1.y()) * (p2.x() - p1.x())
//				- (p3.x() - p1.x()) * (p2.y() - p1.y());
//		if (tmp > 0) {
//			return true;
//		} else {
//			return false;
//
//		}
	}
	bool IsReflex(Point& p1, Point& p2, Point& p3) {
		tppl_float tmp;
		tmp = (p3.y() - p1.y()) * (p2.x() - p1.x())
				- (p3.x() - p1.x()) * (p2.y() - p1.y());
		if (tmp < 0)
			return 1;
		else
			return 0;
	}
//	static bool IsInOn(Point& p1, Point& p2, Point& p3, Point &p) {
//		if (IsConvex(p1, p, p2))
//			return false;
//		if (IsConvex(p2, p, p3))
//			return false;
//		if (IsConvex(p3, p, p1))
//			return false;
//		return true;
//	}

	bool InCone(Point &p1, Point &p2, Point &p3, Point &p) {
		bool convex;

		convex = IsConvex(p1, p2, p3);

		if (convex) {
			if (!IsConvex(p1, p2, p))
				return false;
			if (!IsConvex(p2, p3, p))
				return false;
			return true;
		} else {
			if (IsConvex(p1, p2, p))
				return true;
			if (IsConvex(p2, p3, p))
				return true;
			return false;
		}
	}
	bool InCone(PartitionVertex *v, Point &p) {
		Point p1, p2, p3;

		p1 = v->previous->p;
		p2 = v->p;
		p3 = v->next->p;

		return InCone(p1, p2, p3, p);
	}

	static bool IsIntersect(Point &p11, Point &p12, Point &p21, Point &p22) {
		// same    0
		// touch   1
		// overlap 1
		return Op::IsSegmentIntersect(p11, p12, p21, p22,
				INTERSECT_NORMAL
						| INTERSECT_POINT_POINT
						| INTERSECT_POINT_SEGMENT
						| INTERSECT_POINT_SEGMENT_2);

//		if ((p11.x() == p21.x()) && (p11.y() == p21.y()))
//			return 0;
//		if ((p11.x() == p22.x()) && (p11.y() == p22.y()))
//			return 0;
//		if ((p12.x() == p21.x()) && (p12.y() == p21.y()))
//			return 0;
//		if ((p12.x() == p22.x()) && (p12.y() == p22.y()))
//			return 0;
//
//		Point v1ort, v2ort, v;
//		tppl_float dot11, dot12, dot21, dot22;
//
//		v1ort.x() = p12.y() - p11.y();
//		v1ort.y() = p11.x() - p12.x();
//
//		v2ort.x() = p22.y() - p21.y();
//		v2ort.y() = p21.x() - p22.x();
//
//		v = p21 - p11;
//		dot21 = v.x() * v1ort.x() + v.y() * v1ort.y();
//		v = p22 - p11;
//		dot22 = v.x() * v1ort.x() + v.y() * v1ort.y();
//
//		v = p11 - p21;
//		dot11 = v.x() * v2ort.x() + v.y() * v2ort.y();
//		v = p12 - p21;
//		dot12 = v.x() * v2ort.x() + v.y() * v2ort.y();
//
//		if (dot11 * dot12 > 0)
//			return 0;
//		if (dot21 * dot22 > 0)
//			return 0;
//
//		return 1;
	}



	//helper functions for Triangulate_EC
	void UpdateVertexReflexity(PartitionVertex *v) {
		PartitionVertex *v1 = nullptr, *v3 = nullptr;
		v1 = v->previous;
		v3 = v->next;
		v->isConvex = !IsReflex(v1->p, v->p, v3->p);
	}

	static void UpdateVertex(
			PartitionVertex *v,
			PartitionVertex *vertices,
			long numvertices) {

		PartitionVertex *v1 = nullptr, *v3 = nullptr;
		Point vec1, vec3;

		v1 = v->previous;
		v3 = v->next;

		v->isConvex = IsConvex(v1->p, v->p, v3->p);

		vec1 = Op::Normalize(v1->p - v->p);
		vec3 = Op::Normalize(v3->p - v->p);
		v->angle = vec1.x() * vec3.x() + vec1.y() * vec3.y();

		if (v->isConvex) {
			v->isEar = true;
			for (int i = 0; i < numvertices; i++) {
				if (vertices[i].p == v->p)
					continue;
				if (vertices[i].p == v1->p)
					continue;
				if (vertices[i].p == v3->p)
					continue;

				if (Op::IsInOn(v1->p, v->p, v3->p, vertices[i].p)) {
					v->isEar = false;
					break;
				}
			}
		} else {
			v->isEar = false;
		}
	}


	//helper functions for ConvexPartition_OPT
	void UpdateState(long a, long b, long w, long i, long j,
			DPState2 **dpstates);
	void TypeA(long i, long j, long k, PartitionVertex *vertices,
			DPState2 **dpstates);
	void TypeB(long i, long j, long k, PartitionVertex *vertices,
			DPState2 **dpstates);

	//helper functions for MonotonePartition
	bool Below(Point &p1, Point &p2);
	//	void AddDiagonal(MonotoneVertex *vertices, long *numvertices, long index1,
//			long index2, char *vertextypes,
//			std::set<ScanLineEdge>::iterator *edgeTreeIterators,
//			std::set<ScanLineEdge> *edgeTree, long *helpers);

	//triangulates a monotone polygon, used in Triangulate_MONO
	int TriangulateMonotone(Contour *inPoly, std::list<Contour> *triangles);

	// helper function
	// copy vertices into PartitionVertex
	static PartitionVertex* NewPartitionVertex(const PointChain& pc) {
		ASSERT(pc.size() >= 3);

		long numvertices = pc.size();
//		std::cout << "d size = " << numvertices <<std::endl;
		PartitionVertex* vertices = new PartitionVertex[numvertices];
		long i = 0;
		for (auto& point : pc) {

			vertices[i].isActive = true;
			vertices[i].p = point;
			if (i == (numvertices - 1)) { // last one
				vertices[i].next = &(vertices[0]);
			} else {
				vertices[i].next = &(vertices[i + 1]);
			}
			if (i == 0) {                 //first one
				vertices[i].previous = &(vertices[numvertices - 1]);
			} else {
				vertices[i].previous = &(vertices[i - 1]);
			}
//			std::cout << "i = " << i << "  " << vertices[i].p << std::endl;
			i++;
		}
		return vertices;
	}

public:

	//simple heuristic procedure for removing holes from a list of polygons
	//works by creating a diagonal from the rightmost hole vertex to some visible vertex
	//time complexity: O(h*(n^2)), h is the number of holes, n is the number of vertices
	//space complexity: O(n)
	//params:
	//   inpolys : a list of polygons that can contain holes
	//             vertices of all non-hole polys have to be in counter-clockwise order
	//             vertices of all hole polys have to be in clockwise order
	//   outpolys : a list of polygons without holes
	//returns 1 on success, 0 on failure
	int RemoveHoles(std::list<Contour> *inpolys, std::list<Contour> *outpolys);

	//triangulates a polygon by ear clipping
	//time complexity O(n^2), n is the number of vertices
	//space complexity: O(n)
	//params:
	//   poly : an input polygon to be triangulated
	//          vertices have to be in counter-clockwise order
	//   triangles : a list of triangles (result)
	//returns 1 on success, 0 on failure
	static int EerClipping(const PointChain& poly,
			std::list<PointChain>& triangles) {
		long numvertices;
		PartitionVertex *vertices = nullptr;
		PartitionVertex *ear = nullptr;
		Contour triangle;
		long i, j;
		bool earfound;

		/// if the vertices is less than 3
		/// return
		if (poly.size() < 3)
			return 0;
		if (poly.size() == 3) {
			triangles.push_back(poly);
			return 1;
		}

		/// copy vertices into PartitionVertex
		numvertices = poly.size();
		vertices = NewPartitionVertex(poly);
		// update vertices
		for (i = 0; i < numvertices; i++) {
			UpdateVertex(&(vertices[i]), vertices, numvertices);
		}
//		std::cout << " --------- \n";

		for (i = 0; i < numvertices - 3; i++) {
			earfound = false;
			/// find the most extruded ear
			for (j = 0; j < numvertices; j++) {
				if (!vertices[j].isActive)
					continue;
				if (!vertices[j].isEar)
					continue;
				if (!earfound) {
					earfound = true;
					ear = &(vertices[j]);
				} else {
					if (vertices[j].angle > ear->angle) {
						ear = &(vertices[j]);
					}
				}
			}
//			std::cout << "ear : " << ear->p << std::endl;
			if (!earfound) {
				delete[] vertices;
				return 0;
			}

			triangles.push_back(PointChain(
					ear->previous->p,
					ear->p,
					ear->next->p)
					);
//			std::cout << "T : "<< ear->previous->p << std::endl;
//			std::cout << "  : "<< ear->p << std::endl;
//			std::cout << "  : "<< ear->next->p << std::endl;

			ear->isActive = false;
			ear->previous->next = ear->next;
			ear->next->previous = ear->previous;

			if (i == numvertices - 4)
				break;

			UpdateVertex(ear->previous, vertices, numvertices);
			UpdateVertex(ear->next, vertices, numvertices);
		}
		for (i = 0; i < numvertices; i++) {
			if (vertices[i].isActive) {
				triangles.push_back(PointChain(
						vertices[i].previous->p,
						vertices[i].p,
						vertices[i].next->p)
						);
				break;
			}
		}

		delete[] vertices;

		return 1;
	}

	//triangulates a list of polygons that may contain holes by ear clipping algorithm
	//first calls RemoveHoles to get rid of the holes, and then Triangulate_EC for each resulting polygon
	//time complexity: O(h*(n^2)), h is the number of holes, n is the number of vertices
	//space complexity: O(n)
	//params:
	//   inpolys : a list of polygons to be triangulated (can contain holes)
	//             vertices of all non-hole polys have to be in counter-clockwise order
	//             vertices of all hole polys have to be in clockwise order
	//   triangles : a list of triangles (result)
	//returns 1 on success, 0 on failure
	int Triangulate_EC(std::list<PointChain> *inpolys,
			std::list<PointChain> *triangles) {
		std::list<PointChain> outpolys;
		typename std::list<PointChain>::iterator iter;

		if (!RemoveHoles(inpolys, &outpolys))
			return 0;
		for (iter = outpolys.begin(); iter != outpolys.end(); iter++) {
			if (!Triangulate_EC(&(*iter), triangles))
				return 0;
		}
		return 1;
	}

	//creates an optimal polygon triangulation in terms of minimal edge length
	//time complexity: O(n^3), n is the number of vertices
	//space complexity: O(n^2)
	//params:
	//   poly : an input polygon to be triangulated
	//          vertices have to be in counter-clockwise order
	//   triangles : a list of triangles (result)
	//returns 1 on success, 0 on failure
	int Triangulate_OPT(Contour *poly, std::list<Contour> *triangles);

	//triangulates a polygons by firstly partitioning it into monotone polygons
	//time complexity: O(n*log(n)), n is the number of vertices
	//space complexity: O(n)
	//params:
	//   poly : an input polygon to be triangulated
	//          vertices have to be in counter-clockwise order
	//   triangles : a list of triangles (result)
	//returns 1 on success, 0 on failure
	int Triangulate_MONO(Contour *poly, std::list<Contour> *triangles);

	//triangulates a list of polygons by firstly partitioning them into monotone polygons
	//time complexity: O(n*log(n)), n is the number of vertices
	//space complexity: O(n)
	//params:
	//   inpolys : a list of polygons to be triangulated (can contain holes)
	//             vertices of all non-hole polys have to be in counter-clockwise order
	//             vertices of all hole polys have to be in clockwise order
	//   triangles : a list of triangles (result)
	//returns 1 on success, 0 on failure
	int Triangulate_MONO(std::list<Contour> *inpolys,
			std::list<Contour> *triangles);

	//creates a monotone partition of a list of polygons that can contain holes
	//time complexity: O(n*log(n)), n is the number of vertices
	//space complexity: O(n)
	//params:
	//   inpolys : a list of polygons to be triangulated (can contain holes)
	//             vertices of all non-hole polys have to be in counter-clockwise order
	//             vertices of all hole polys have to be in clockwise order
	//   monotonePolys : a list of monotone polygons (result)
	//returns 1 on success, 0 on failure
	int MonotonePartition(std::list<Contour> *inpolys,
			std::list<Contour> *monotonePolys);

	//partitions a polygon into convex polygons by using Hertel-Mehlhorn algorithm
	//the algorithm gives at most four times the number of parts as the optimal algorithm
	//however, in practice it works much better than that and often gives optimal partition
	//uses triangulation obtained by ear clipping as intermediate result
	//time complexity O(n^2), n is the number of vertices
	//space complexity: O(n)
	//params:
	//   poly : an input polygon to be partitioned
	//          vertices have to be in counter-clockwise order
	//   parts : resulting list of convex polygons
	//returns 1 on success, 0 on failure
	int ConvexPartition_HM(Contour *poly, std::list<Contour> *parts);

	//partitions a list of polygons into convex parts by using Hertel-Mehlhorn algorithm
	//the algorithm gives at most four times the number of parts as the optimal algorithm
	//however, in practice it works much better than that and often gives optimal partition
	//uses triangulation obtained by ear clipping as intermediate result
	//time complexity O(n^2), n is the number of vertices
	//space complexity: O(n)
	//params:
	//   inpolys : an input list of polygons to be partitioned
	//             vertices of all non-hole polys have to be in counter-clockwise order
	//             vertices of all hole polys have to be in clockwise order
	//   parts : resulting list of convex polygons
	//returns 1 on success, 0 on failure
	int ConvexPartition_HM(std::list<Contour> *inpolys,
			std::list<Contour> *parts);

	//optimal convex partitioning (in terms of number of resulting convex polygons)
	//using the Keil-Snoeyink algorithm
	//M. Keil, J. Snoeyink, "On the time bound for convex decomposition of simple polygons", 1998
	//time complexity O(n^3), n is the number of vertices
	//space complexity: O(n^3)
	//   poly : an input polygon to be partitioned
	//          vertices have to be in counter-clockwise order
	//   parts : resulting list of convex polygons
	//returns 1 on success, 0 on failure
	int ConvexPartition_OPT(Contour *poly, std::list<Contour> *parts);
}
;

}

#endif

