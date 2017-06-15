/*
 * ts_edge.h
 *
 *  Created on: May 31, 2015
 *      Author: zhou
 */

#ifndef TS_EDGE_H_
#define TS_EDGE_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_face.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include <iomanip>

namespace TS {

template<class TYPE, st DIM> class Triangle;
template<class TYPE, st DIM> class Face;
template<class TYPE, st DIM> class Surface;

template<class TYPE, st DIM>
class Edge: public Segment<TYPE, DIM>, public std::enable_shared_from_this<
		Edge<TYPE, DIM> > {
public:
	typedef Segment<TYPE, DIM> base_class;
	typedef Edge<TYPE, DIM> self_class;
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
	typedef Surface<TYPE, DIM> Sur;
	typedef std::shared_ptr<Sur> spSur;
	typedef List<spSeg> list_spSeg;
	typedef List<spVer> list_spVer;
	typedef List<spFac> list_spFac;
	typedef List<spTri> list_spTri;

	typedef typename list_spFac::iterator iterator_face;
	typedef typename list_spFac::const_iterator const_iterator_face;
public:
	list_spFac faces;

	Edge(spVer a, spVer b) :
			base_class(a, b) {
	}

	void attach(){
		this->v1->attach(get_this());
		this->v2->attach(get_this());
	}

	spEdg get_this() {
		return this->shared_from_this();
	}

	iterator_face begin_face(){
		return faces.begin();
	}
	const_iterator_face begin_face() const{
		return faces.begin();
	}
	iterator_face end_face(){
		return faces.end();
	}
	const_iterator_face end_face() const{
		return faces.end();
	}

	// show =====================================
	void show() const;

};

template<class TYPE, st DIM>
void Edge<TYPE, DIM>::show() const {
	std::ios::fmtflags f(std::cout.flags());
	std::cout.setf(std::ios::right);
	if (DIM == 3) {
		std::cout << std::setprecision(5);
		std::cout << std::fixed;
		std::cout << "edg (" << this->v1->at(0) << " " << this->v1->at(1) << " "
				<< this->v1->at(2) << ")";
		std::cout << "->(" << this->v2->at(0) << " " << this->v2->at(1) << " "
				<< this->v2->at(2) << ") ";
		std::cout << " -> tri " << faces.size() << "\n";
	} else if (DIM == 2) {
		std::cout << "edg (" << this->v1->at(0) << " " << this->v1->at(1) << " "
				<< ")";
		std::cout << "->(" << this->v2->at(0) << " " << this->v2->at(1) << " "
				<< ") ";
		std::cout << " -> tri " << faces.size() << "\n";
	}
	std::cout.setf(f);
}

}

#endif /* TS_TS_EDGE_H_ */
