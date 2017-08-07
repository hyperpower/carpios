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

template<class TYPE, St DIM, class EDGE>
class Vertex_: public Point_<TYPE, DIM> {
public:
	static const St Dim = DIM;
	typedef Point_<TYPE, DIM> base_class;
	typedef Vertex_<TYPE, DIM, EDGE> self_class;
	typedef St size_type;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef Vertex_<TYPE, DIM, EDGE> Ver;
	typedef Ver* pVer;
	typedef EDGE Edg;
	typedef Edg* pEdg;

	typedef std::list<pEdg> list_pEdg;
	typedef std::list<pVer> list_pVer;

	typedef typename std::list<pEdg>::iterator iterator_pEdg;
	typedef typename std::list<pEdg>::const_iterator const_iterator_pEdg;

public:
	list_pEdg _ledges;
	public:
	Vertex_(const base_class& poi) :
			base_class(poi), _ledges() {
	}

	Vertex_(const Vt&x, const Vt& y, const Vt& z = 0) :
			base_class(x, y, z), _ledges() {
	}

	void attach(pEdg pe) {
		if (!has_edge(pe)) {
			_ledges.push_back(pe);
		}
	}

	void detach(pEdg pe) {
		if (has_edge(pe)) {
			_ledges.remove(pe);
		}
	}

	list_pEdg& edges() {
		return _ledges;
	}

	const list_pEdg& edges() const {
		return _ledges;
	}

	inline bool is_unattached() {
		if (_ledges.empty())
			return true;
		return false;
	}

	bool has_edge(pEdg pe) {
		_RETURN_VAL_IF_FAIL(pe != nullptr, false);
		for (auto& e : _ledges) {
			if (e == pe) {
				return true;
			}
		}
		return false;
	}

	bool has_one_edge(pEdg pe) {
		_RETURN_VAL_IF_FAIL(pe != nullptr, false);
		if (this->_ledges.size() != 1) {
			return false;
		}
		/// compare pointer address
		return (*(this->_ledges.begin())) == pe;
	}

	iterator_pEdg begin_edge() {
		return this->_ledges.begin();
	}
	const_iterator_pEdg begin_edge() const {
		return this->_ledges.begin();
	}
	iterator_pEdg end_edge() {
		return this->_ledges.end();
	}
	const_iterator_pEdg end_edge() const {
		return this->_ledges.end();
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
			pEdg& s = (*iter);
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
	 * is_connected:
	 * @v2: another #Vertex.
	 *
	 * Returns: if @v2 is connect to one of edges in ledges
	 * this segment else %NULL.
	 */
	pEdg is_connected(pVer v2) {
		for (auto& e : this->_ledges) {
			if (e->v1 == v2 || e->v2 == v2)
				return e;
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

	pVer other_vertex_on_edge(pEdg pe) {
		ASSERT(this->has_edge(pe));
		if (pe->vertex(0) == this) {
			return pe->vertex(1);
		} else {
			return pe->vertex(0);
		}
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
	/**
	 * IsConnected:
	 * @v1: one #Vertex
	 * @v2: another #Vertex.
	 *
	 * Returns: if @v1 and @v2 have the pedge in the list, return true;
	 * this segment else %NULL.
	 */
	static pEdg IsConnected(pVer v1, pVer v2) {
		for (auto& e : v1->_ledges) {
			if (v2->has_edge(e)) {
				return e;
			}
		}
		return nullptr;
	}

	static bool IsClosed(const std::list<pVer>& lver) {
		for (auto iter = lver.begin(); iter != lver.end(); ++iter) {
			auto itern = std::next(iter, 1);
			if(itern == lver.end()){
				itern = lver.begin();
			}
			pVer pv = *iter;
			pVer pvn = *itern;
			if(IsConnected(pv, pvn) == nullptr){
				return false;
			}
		}
		return true;
	}

}
;

}

#endif /* _VERTEX_H_ */
