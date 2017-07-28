/************************
 //  \file   ts_face.h
 //  \brief
 //
 //  \author czhou
 //  \date   15 juin 2015
 ***********************/
#ifndef TS_FACE_H_
#define TS_FACE_H_

#include "_vertex.hpp"
#include "_edge.hpp"

namespace carpio {

template<class TYPE, St DIM> class TriSurface_;

template<class TYPE, St DIM>
class TriFace_ {
	static const St Dim = DIM;
	typedef TriFace_<TYPE, DIM> self_class;
	typedef TriFace_<TYPE, DIM> TriFace;
	typedef TriFace* pTriFace;
	typedef St size_type;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Point;
	typedef Point* pPoi;
	typedef Vertex_<TYPE, DIM> Vertex;
	typedef Vertex* pVertex;
	typedef Edge_<TYPE, DIM> Edge;
	typedef Edge* pEdge;
	typedef const Edge* const_pEdge;
	typedef TriSurface_<TYPE, DIM> TriSurface;
	typedef TriSurface* pTriSurface;
	typedef const TriSurface* const_pTriSurface;

	typedef std::list<pTriSurface> list_pTriSurface;
	typedef std::list<pTriFace> list_pTriFace;
public:
	pEdge e1;
	pEdge e2;
	pEdge e3;

	list_pTriSurface surfaces;
public:
	TriFace_(pEdge a, pEdge b, pEdge c, pTriSurface psur) :
			e1(a), e2(b), e3(c) {
		surfaces.push_back(psur);

		e1->faces.push_back(this);
		e2->faces.push_back(this);
		e3->faces.push_back(this);

	}

	void attach_surface(pTriSurface ps) {
		surfaces.push_back(ps);
	}

	/**
	 * gts_face_has_parent_surface:
	 * @f: a #GtsFace
	 * @s: a #GtsSurface.
	 *
	 * Returns: %TRUE if @f belongs to @s, %FALSE otherwise.
	 */
	bool has_parent_surface(const_pTriSurface ps) const {
		_return_val_if_fail(ps != nullptr, false);
		for (pTriSurface s : this->surfaces) {
			if (ps == s) {
				return true;
			}
		}
		return false;
	}

	list_pTriFace get_neighbor_faces(const_pTriSurface ps) const {
		ASSERT(ps != nullptr);
		list_pTriFace res;
		for (St i = 0; i < 3; i++) {
			pEdge e = this->edge(i);
			for (auto iter = e->begin_face(); iter != e->end_face(); iter++) {
				pTriFace spf = *iter;
				if (spf != this && spf->has_parent_surface(ps)) {
					res.push_back(spf);
				}
			}
		}
		return std::move(res);
	}

	pEdge operator[](const St& idx) {
		ASSERT(idx >= 0 && idx < 3);
		pEdge l[] = { e1, e2, e3 };
		return l[idx];
	}

	const_pEdge operator[](const St& idx) const {
		ASSERT(idx >= 0 && idx < 3);
		pEdge l[] = { e1, e2, e3 };
		return l[idx];
	}

	pEdge edge(const St& idx) {
		return this->operator [](idx);
	}

	const_pEdge edge(const St& idx) const{
		return this->operator [](idx);
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

