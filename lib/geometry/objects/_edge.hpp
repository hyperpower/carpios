/*
 * ts_edge.h
 *
 *  Created on: May 31, 2015
 *      Author: zhou
 */

#ifndef _EDGE_HPP_
#define _EDGE_HPP_

#include "_vertex.hpp"
#include <iomanip>
#include <list>
#include <set>
#include <functional>

namespace carpio {

//template<class TYPE, St DIM> class TriFace_;
//template<class TYPE, St DIM> class TriSurface_;

template<class TYPE, St DIM, class EDGE>
class Vertex_;

template<class TYPE, St DIM, class FACE>
class Edge_ {
public:
	static const St Dim = DIM;
	typedef Edge_<TYPE, DIM, FACE> Self;
	typedef St size_type;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef Vertex_<TYPE, DIM, Self> Ver;
	typedef Ver* pVer;
	typedef const Ver* const_pVer;
	typedef Edge_<TYPE, DIM, FACE> Edg;
	typedef Edg* pEdg;
	typedef FACE Fac;
	typedef Fac* pFac;
	typedef const Fac* const_pFac;

	typedef std::list<pEdg> list_pEdge;
	typedef std::list<pVer> list_pVertex;
	typedef std::list<pFac> list_pTriFace;

	typedef typename std::list<pFac>::iterator iterator_face;
	typedef typename std::list<pFac>::const_iterator const_iterator_face;

public:
	pVer v1;
	pVer v2;

	list_pTriFace faces;

	Any _any_data;

public:
	Edge_(pVer a, pVer b) :
			v1(a), v2(b) {
		v1->attach(this);
		v2->attach(this);
	}

	inline bool is_unattached() {
		if (faces.empty())
			return true;
		return false;
	}

	iterator_face begin_face() {
		return faces.begin();
	}
	const_iterator_face begin_face() const {
		return faces.begin();
	}
	iterator_face end_face() {
		return faces.end();
	}
	const_iterator_face end_face() const {
		return faces.end();
	}

	pVer operator[](const St& idx) {
		ASSERT(idx >= 0 && idx < 2);
		return idx == 0 ? v1 : v2;
	}

	const_pVer operator[](const St& idx) const {
		ASSERT(idx >= 0 && idx < 2);
		return idx == 0 ? v1 : v2;
	}

	pVer vertex(const St& idx) {
		ASSERT(idx >= 0 && idx < 2);
		return idx == 0 ? v1 : v2;
	}

	const_pVer vertex(const St& idx) const {
		ASSERT(idx >= 0 && idx < 2);
		return idx == 0 ? v1 : v2;
	}

	bool has_vertex(const_pVer pv) const {
		_RETURN_VAL_IF_FAIL(pv != nullptr, false);
		return v1 == pv || v2 == pv;
	}

	bool has_face(const_pFac pf) const {
		_RETURN_VAL_IF_FAIL(pf != nullptr, false);
		for (const_pFac s : this->faces) {
			if (s == pf) {
				return true;
			}
		}
		return false;
	}

	bool has_one_face(const_pFac pf) const {
		_RETURN_VAL_IF_FAIL(pf != nullptr, false);
		if (this->faces.size() != 1) {
			return false;
		}
		/// compare pointer address
		return (*(this->faces.begin())) == pf;
	}

	void detach(pFac pf) {
		if (has_face(pf)) {
			faces.remove(pf);
		}
	}

	void attach(pFac pf) {
		if (!has_face(pf)) {
			faces.push_back(pf);
		}
	}

	void reverse() {
		pVer tmp;
		tmp = v1;
		v1 = v2;
		v2 = tmp;
	}

	Poi tangent() const {
		return (*v2) - (*v1);
	}

	St size_faces() const {
		return faces.size();
	}

	Any& data() {
		return _any_data;
	}

	const Any& data() const{
		return _any_data;
	}

	// show =====================================
	void show() const {
		std::ios::fmtflags f(std::cout.flags());
		std::cout.setf(std::ios::right);
		if (DIM == 3) {
			std::cout << std::setprecision(5);
			std::cout << std::fixed;
			std::cout << "edg (" << this->v1->at(0) << " " << this->v1->at(1)
					<< " " << this->v1->at(2) << ")";
			std::cout << "->(" << this->v2->at(0) << " " << this->v2->at(1)
					<< " " << this->v2->at(2) << ") ";
			std::cout << " -> tri " << faces.size() << "\n";
		} else if (DIM == 2) {
			std::cout << "edg (" << this->v1->at(0) << " " << this->v1->at(1)
					<< " " << ")";
			std::cout << "->(" << this->v2->at(0) << " " << this->v2->at(1)
					<< " " << ") ";
			std::cout << " -> tri " << faces.size() << "\n";
		}
		std::cout.setf(f);
	}

	bool is_connected(pVer pv) const {
		_RETURN_VAL_IF_FAIL(pv != nullptr, false);
		return (pv == this->v1) || (pv == this->v2);
	}

	bool is_connected_2(const Edg& e) const {
		return v2 == (e.v1) || v2 == (e.v2);
	}
	bool is_connected_1(const Edg& e) const {
		return v1 == (e.v1) || v1 == (e.v2);
	}

	static bool IsSameV(  //
			const Edg& e, //
			const Ver& ver1, //
			const Ver& ver2) { //
		return ((*(e.v1)) == ver1 && (*(e.v2)) == ver2)
				|| ((*(e.v1)) == ver2 && (*(e.v2)) == ver1);
	}

	static bool IsSameDirectionV(  //
				const Edg& e, //
				const Ver& ver1, //
				const Ver& ver2) { //
		if (IsSameV(e,ver1, ver2)){
			return ((*(e.v1)) == ver1 && (*(e.v2)) == ver2);
		}
		return false;
	}
	static bool IsSameA(  //
			const Edg& e, //
			const Ver& ver1, //
			const Ver& ver2) { //
		return ((e.v1) == &ver1 && (e.v2) == &ver2)
				|| ((e.v1) == &ver2 && (e.v2) == &ver1);
	}
	static bool IsSameV( // comapare value
			const Edg& e1,   //
			const Edg& e2) { //
		return ((*(e1.v1)) == (*(e2.v1)) && (*(e1.v2)) == (*(e2.v2)))
				|| ((*(e1.v1)) == (*(e2.v2)) && (*(e1.v2)) == (*(e2.v1)));
	}
	static bool IsSameA( // comapare adress
			const Edg& e1,   //
			const Edg& e2) { //
		return (((e1.v1)) == ((e2.v1)) && ((e1.v2)) == ((e2.v2)))
				|| (((e1.v1)) == ((e2.v2)) && ((e1.v2)) == ((e2.v1)));
	}
	/// compare pointer address
	static bool IsConnected(
			const Edg& e1,
			const Edg& e2) {
		return ((e1.v1) == (e2.v1) || (e1.v1) == (e2.v2)
				|| (e1.v2) == (e2.v1) || (e1.v2) == (e2.v2));
	}

	/// compare value
	static bool IsTouched(  //
			const Edg& s1, //
			const Edg& s2) {
		return (*(s1.v1) == *(s2.v1) || *(s1.v1) == *(s2.v2)
				|| *(s1.v2) == *(s2.v1) || *(s1.v2) == *(s2.v2));
	}

	/// for each vertex in list of edges
	template<class Container>
	static void ForEachVertex(
			Container& edges,
			std::function<void(pVer&)> fun) {
		typename Container::value_type dummy;
		_ForEachVertex(edges, fun, dummy);
	}

	template<class Container>
	static void _ForEachVertex(
			Container& edges,
			std::function<void(pVer&)> fun,
			pEdg dummy) {
		std::set<pVer> tmp;
		for (auto& edg : edges) {
			for (int i = 0; i < 2; ++i) {
				pVer v = edg->vertex(i);
				if (tmp.find(v) == tmp.end()) {
					fun(v);
					tmp.insert(v);
				}
			}
		}
	}

	template<class Container>
	static void _ForEachVertex(
			Container& edges,
			std::function<void(pVer&)> fun,
			Edg dummy) {
		std::set<pEdg> tmp;
		for (auto& edg : edges) {
			for (int i = 0; i < 2; ++i) {
				pVer v = edg.vertex(i);
				if (tmp.find(v) == tmp.end()) {
					fun(v);
					tmp.insert(v);
				}
			}
		}
	}

};

template<class TYPE, St DIM, class FACE, class DATA>
class EdgeD_: public Edge_<TYPE, DIM, FACE> {
public:
	static const St Dim = DIM;
	typedef Edge_<TYPE, DIM, FACE> Base;
	typedef EdgeD_<TYPE, DIM, FACE, DATA> Self;
	typedef St size_type;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef Vertex_<TYPE, DIM, Self> Ver;
	typedef Ver* pVer;
	typedef const Ver* const_pVer;
	typedef EdgeD_<TYPE, DIM, FACE, DATA> Edg;
	typedef Edg* pEdg;
	typedef FACE Fac;
	typedef Fac* pFac;
	typedef const Fac* const_pFac;

protected:
	DATA _data;
	public:
	EdgeD_(pVer a, pVer b) :
			Base(a, b) {
	}
	const DATA& data() const {
		return _data;
	}
	DATA& data() {
		return _data;
	}
};

}

#endif /* _EDGE_HPP_ */

