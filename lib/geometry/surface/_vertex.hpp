/*
 *  _vertex.h
 *
 *  Created on: May 30, 2015
 *  Author: zhou
 */

#ifndef _VERTEX_HPP_
#define _VERTEX_HPP_

#include "../geometry.hpp"
#include "_edge.hpp"

namespace carpio {

//template<class TYPE, St DIM> class TriSurface_;
template<class TYPE, St DIM> class Edge_;

template<class TYPE, St DIM>
class Vertex_: public Point_<TYPE, DIM> {
public:
	static const St Dim = DIM;
	typedef Point_<TYPE, DIM> base_class;
	typedef Vertex_<TYPE, DIM> self_class;
	typedef St size_type;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef Vertex_<TYPE, DIM> Ver;
	typedef Ver* pVer;
	typedef Edge_<TYPE, DIM> Edg;
	typedef Edge* pEdg;

	typedef std::list<pEdg> list_pEdg;
	typedef std::list<pVer> list_pVer;
public:
	list_pEdg _ledges;
public:
	Vertex_(const base_class& poi) :
			base_class(poi), _ledges() {
	}

	Vertex_(const Vt&x, const Vt& y, const Vt& z = 0) :
			base_class(x, y, z), _ledges() {
	}

	void attach(pEdg pe){
		_ledges.push_back(pe);
	}

	inline bool vertex_is_unattached() {
		if (_ledges.empty())
			return true;
		return false;
	}

	/**
	 * vertex_replace:
	 * @with: another #GtsVertex.
	 *
	 * Replaces vertex this with vertex @with. this and @with must be
	 * different.  All the #GtsSegment which have @v has one of their
	 * vertices are updated.  The segments list of vertex @v is freed and
	 * @v->segments is set to %NULL.
	 */
	void vertex_replace(const Ver& with) {
		_IF_FALSE_RETRUN(this != &with);

		for (auto iter = _ledges.begin(); iter != _ledges.end(); ++iter) {
			pEdge& s = (*iter);
			if (s->v1 != with && s->v2 != with)
				with->_ledges.push_front(s);
			if (s->v1 == this)
				s->v1 = with;
			if (s->v2 == this)
				s->v2 = with;
		}
		this->_ledges.clear();
	}
	/**
	 * vertices_are_connected:
	 * @v2: another #GtsVertex.
	 *
	 * Returns: if @v1 and @v2 are the vertices of the same #GtsSegment
	 * this segment else %NULL.
	 */
	pEdg vertices_are_connected(pVer v2) {
		for (auto iter = this->_ledges->begin(); iter != this->_ledges->end();
				++iter) {
			pEdge& s = (*iter);
			if (s->v1 == v2 || s->v2 == v2)
				return s;
		}
		return nullptr;
	}
	/**
	 * vertices_from_segments:
	 * @segments: a list of Segment.
	 *
	 * Returns: a list of Vertex, vertices of a Segment in list segments.
	 * Each element in the list is unique (no duplicates).
	 */
	void vertices_from_segments(list_pVer& lpver) const {

	}

	void show() const {
		std::ios::fmtflags f(std::cout.flags());
		std::cout.setf(std::ios::right);
		std::cout << "ver ( " << this->at(0) << " , " << this->at(1);
		if (DIM == 3) {
			std::cout << " , " << this->at(2) << " )";
		} else {
			std::cout << " )";
		}
		std::cout << " -> seg = " << _ledges.size() << "\n";
		std::cout.setf(f);
	}
}
;

}

#endif /* _VERTEX_H_ */
