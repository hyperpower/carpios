/************************
 //  \file   ts_face.h
 //  \brief
 // 
 //  \author czhou
 //  \date   15 juin 2015 
 ***********************/
#ifndef TS_FACE_H_
#define TS_FACE_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_triangle.h"
#include "ts_edge.h"

namespace TS {

template<class TYPE, st DIM> class Surface;

template<class TYPE, st DIM>
class Face: public Triangle<TYPE, DIM>, public std::enable_shared_from_this<
		Face<TYPE, DIM> > {
public:
	typedef Triangle<TYPE, DIM> base_class;
	typedef Face<TYPE, DIM> self_class;
	typedef st size_type;
	typedef Point<TYPE, DIM> Poi;
	typedef std::shared_ptr<Poi> spPoi;
	typedef Vertex<TYPE, DIM> Ver;
	typedef std::shared_ptr<Ver> spVer;
	typedef Segment<TYPE, DIM> Seg;
	typedef std::shared_ptr<Seg> spSeg;
	typedef Edge<TYPE, DIM> Edg;
	typedef std::shared_ptr<Edg> spEdg;
	typedef Triangle<TYPE, DIM> Tri;
	typedef std::shared_ptr<Tri> spTri;
	typedef Face<TYPE, DIM> Fac;
	typedef std::shared_ptr<Fac> spFac;
	typedef std::shared_ptr<const Fac> const_spFac;
	typedef Surface<TYPE, DIM> Sur;
	typedef std::shared_ptr<Sur> spSur;
	typedef Sur* pSur;
	typedef List<spSeg> list_spSeg;
	typedef List<spVer> list_spVer;
	typedef List<spTri> list_spTri;
	typedef List<spFac> list_spFac;
	typedef List<Sur*> list_pSur;
public:
	list_pSur surfaces;
public:
	Face(spEdg a, spEdg b, spEdg c, pSur sur) :
			base_class(a, b, c) {
		surfaces.push_back(sur);
	}

	void attach() {
		spFac f = get_this();
		this->e1->faces.push_back(f);
		this->e2->faces.push_back(f);
		this->e3->faces.push_back(f);
	}

	void attach_surface(pSur ps){
		surfaces.push_back(ps);
	}

	spFac get_this() {
		return this->shared_from_this();
	}
	const_spFac get_this() const {
		return this->shared_from_this();
	}

	/**
	 * gts_face_has_parent_surface:
	 * @f: a #GtsFace
	 * @s: a #GtsSurface.
	 *
	 * Returns: %TRUE if @f belongs to @s, %FALSE otherwise.
	 */
	bool has_parent_surface(const Sur* ps) const {
		_return_val_if_fail(ps != nullptr, false);
		for (pSur s : this->surfaces) {
			if (ps == s) {
				return true;
			}
		}
		return false;
	}

	list_spFac get_neighbor_faces(const Sur* ps) const {
		ASSERT(ps != nullptr);
		list_spFac res;
		for (st i = 0; i < 3; i++) {
			spEdg e = this->edge(i);
			for (auto iter = e->begin_face(); iter != e->end_face(); iter++) {
				spFac spf = *iter;
				if (spf != this->get_this() && spf->has_parent_surface(ps)) {
					res.push_back(spf);
				}
			}
		}
		return std::move(res);
	}



	void show() const {
		std::ios::fmtflags f(std::cout.flags());
		std::cout.setf(std::ios::right);
		this->e1->show();
		this->e2->show();
		this->e3->show();
		std::cout << "    -> sur = " << surfaces.size() << "\n";
		std::cout.setf(f);
	}

};

}

#endif /* TS_FACE_H_ */
