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

	bool is_connected(pVer pv) const{
		_RETURN_VAL_IF_FAIL(pv != nullptr, false);
		return (pv == this->v1) || (pv == this->v2);
	}

	bool is_connected_2(const Edg& e) const {
		return v2 == (e.v1) || v2 == (e.v2);
	}
	bool is_connected_1(const Edg& e) const {
		return v1 == (e.v1) || v1 == (e.v2);
	}

	static bool IsConnected(  //
			const Edg& s, //
			const Ver& e1, //
			const Ver& e2) { //
		return ((*(s.v1)) == e1 && (*(s.v2)) == e2)
				|| ((*(s.v1)) == e2 && (*(s.v2)) == e1);
	}
	static bool IsIdentical(     //
			const Edg& s1,   //
			const Edg& s2) { //
		return ((*(s1.v1)) == (*(s2.v1)) && (*(s1.v2)) == (*(s2.v2)))
				|| ((*(s1.v1)) == (*(s2.v2)) && (*(s1.v2)) == (*(s2.v1)));
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

};

}

#endif /* _EDGE_HPP_ */

