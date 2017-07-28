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

namespace carpio {

template<class TYPE, St DIM> class TriFace_;
template<class TYPE, St DIM> class TriSurface_;

template<class TYPE, St DIM>
class Edge_ {
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
	typedef Edg* pEdg;
	typedef TriFace_<TYPE, DIM> Fac;
	typedef Fac* pFac;

	typedef std::list<pEdg> list_pEdge;
	typedef std::list<pVer> list_pVertex;
	typedef std::list<pFac> list_pTriFace;

	typedef typename std::list<pFac>::iterator iterator_face;
	typedef typename std::list<pFac>::const_iterator const_iterator_face;

public:
	pVer v1;
	pVer v2;

	list_pTriFace faces;

	Edge_(pVer a, pVer b) :
			v1(a), v2(b) {
		v1->attach(this);
		v2->attach(this);
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

};


}

#endif /* _EDGE_HPP_ */

