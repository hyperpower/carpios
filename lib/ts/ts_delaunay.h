#ifndef TS_DELAUNAY_H_
#define TS_DELAUNAY_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_edge.h"
#include "ts_face.h"
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>

namespace TS {

template<class TYPE, st DIM>
class Delaunay {
public:
	typedef st size_type;
	typedef TYPE vt;
	typedef Point<TYPE, DIM> Poi;
	typedef std::shared_ptr<Poi> spPoi;
	typedef std::shared_ptr<const Poi> spcPoi;
	typedef Vertex<TYPE, DIM> Ver;
	typedef std::shared_ptr<Ver> spVer;
	typedef std::shared_ptr<const Ver> spcVer;
	typedef Segment<TYPE, DIM> Seg;
	typedef std::shared_ptr<Seg> spSeg;
	typedef Edge<TYPE, DIM> Edg;
	typedef std::shared_ptr<Edg> spEdg;
	typedef Triangle<TYPE, DIM> Tri;
	typedef std::shared_ptr<Tri> spTri;
	typedef Face<TYPE, DIM> Fac;
	typedef std::shared_ptr<Fac> spFac;
	typedef Surface<TYPE, DIM> Sur;
	typedef Sur* pSur;
	typedef std::shared_ptr<Sur> spSur;
	typedef List<spSeg> list_spSeg;
	typedef List<spVer> list_spVer;
	typedef List<spPoi> list_spPoi;
	typedef List<spcPoi> list_spcPoi;
	typedef List<spcVer> list_spcVer;
	typedef List<spTri> list_spTri;
	typedef List<spFac> list_spFac;
	typedef List<spSur> list_spSur;
	typedef std::function<void(Fac&)> Fun_Fac;
	static const st Dim;

	pSur _psur;
public:
	Delaunay() {
		_psur = nullptr;
	}

	Delaunay(pSur psur) {
		this->_psur = psur;
	}
	/**
	 * gts_triangle_enclosing:
	 * @klass: the class of the new triangle.
	 * @points: a list of #GtsPoint.
	 * @scale: a scaling factor (must be larger than one).
	 *
	 * Builds a new triangle (including new vertices and edges) enclosing
	 * the plane projection of all the points in @points. This triangle is
	 * equilateral and encloses a rectangle defined by the maximum and
	 * minimum x and y coordinates of the points. @scale is an homothetic
	 * scaling factor. If equal to one, the triangle encloses exactly the
	 * enclosing rectangle.
	 *
	 * Returns: a new #GtsTriangle.
	 */
	spTri triangle_enclosing(list_spcPoi& listv, vt scale) {
		vt xmax, xmin, ymax, ymin;
		//gdouble xo, yo, r;
		//GtsVertex * v1, *v2, *v3;
		//GtsEdge * e1, *e2, *e3;

		if (listv.empty()) {
			return spTri(nullptr);
		}
		auto iter = listv.begin();
		auto& p = (*iter);
		xmax = xmin = p->x();
		ymax = ymin = p->y();
		for (; iter != listv.end(); ++iter) {
			spcPoi p = (*iter);
			if (p->x() > xmax) {
				xmax = p->x();
			} else if (p->x() < xmin) {
				xmin = p->x();
			}
			if (p->y() > ymax) {
				ymax = p->y();
			} else if (p->y() < ymin) {
				ymin = p->y();
			}
		}
		vt xo = (xmax + xmin) / 2.;
		vt yo = (ymax + ymin) / 2.;
		vt r = scale
				* sqrt((xmax - xo) * (xmax - xo) + (ymax - yo) * (ymax - yo));
		if (r == 0.0) {
			r = scale;
		}
		vt SQRT3 = 1.73205080757;
		spVer v1(new Ver(xo + r * SQRT3, yo - r, 0.0));
		spVer v2(new Ver(xo, yo + 2. * r, 0.0));
		spVer v3(new Ver(xo - r * SQRT3, yo - r, 0.0));
		spEdg e1(new Edg(v1, v2));
		spEdg e2(new Edg(v2, v3));
		spEdg e3(new Edg(v3, v1));
		return spTri(new Tri(e1, e2, e3));
	}

	/**
	 * gts_delaunay_add_vertex:
	 * @surface: a #GtsSurface.
	 * @v: a #GtsVertex.
	 * @guess: %NULL or a #GtsFace belonging to @surface to be used as an initial
	 * guess for point location.
	 *
	 * Adds vertex @v to the Delaunay triangulation defined by
	 * @surface. If @v is not contained in the convex hull bounding
	 * @surface, @v is not added to the triangulation.
	 *
	 * Returns: %NULL is @v has been successfully added to @surface or was
	 * already contained in @surface, @v if @v is not contained in the
	 * convex hull bounding surface or a #GtsVertex having the same x and
	 * y coordinates as @v.
	 */
	spVer add_vertex(Sur& sur, spVer v, spFac guess) {
		spFac f(nullptr);

		if (!(f = point_locate(v, sur, guess))){
			return v;
		}
		return nullptr;
		//return gts_delaunay_add_vertex_to_face(surface, v, f);
	}

	/**
	 * gts_point_locate:
	 * @p: a #GtsPoint.
	 * @surface: a #GtsSurface.
	 * @guess: %NULL or a face of @surface close to @p.
	 *
	 * Locates the face of the planar projection of @surface containing
	 * @p. The planar projection of @surface must define a connected set
	 * of triangles without holes and bounded by a convex boundary. The
	 * algorithm is randomized and performs in O(n^1/3) expected time
	 * where n is the number of triangles of @surface.
	 *
	 * If a good @guess is given the point location can be significantly faster.
	 *
	 * Returns: a #GtsFace of @surface containing @p or %NULL if @p is not
	 * contained within the boundary of @surface.
	 */
	spFac point_locate (spPoi p,
				    spFac guess)
	{
	  spFac fr;
	  spPoi o;

	  //g_return_val_if_fail (p != NULL, NULL);
	  //g_return_val_if_fail (surface != NULL, NULL);
	  //g_return_val_if_fail (guess == NULL ||
	  //gts_face_has_parent_surface (guess, surface), NULL);
	  guess->has_parent_surface(_psur);

	  if (guess == nullptr){
	    guess = closest_face (_psur, p);
	  }else{
		  if (guess->orientation() > 0.0){
			  return nullptr;
		  }
	  }
	  if (guess == nullptr){
	    return nullptr;
	  }

	  //spPoi o = guess->barycenter();
	  // fr = _point_locate (o, p, guess, _psur);

	  return fr;
	}

};

}

#endif
