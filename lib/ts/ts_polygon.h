/*
 * ts_polygon.h
 *
 *  Created on: Jul 19, 2017
 *      Author: zhou
 */

#ifndef _TS_TS_POLYGON_H_
#define _TS_POLYGON_H_

#include "ts_define.h"
#include "ts_vertex.h"

namespace TS {

template<class TYPE, st DIM> class Segment;
template<class TYPE, st DIM> class Surface;
template<class TYPE, st DIM> class Edge;

template<class TYPE, st DIM>
class Poly {
public:
	typedef Poly<TYPE, DIM> Self;
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
	typedef Surface<TYPE, DIM> Sur;
	typedef std::shared_ptr<Sur> spSur;
	typedef List<spEdg> list_spEdg;
	typedef List<spVer> list_spVer;
	typedef List<spTri> list_spTri;
	typedef List<spSur> list_spSur;
	typedef typename list_spVer::iterator iterator;
	typedef typename list_spVer::const_iterator const_iterator;
	static const st Dim = DIM;
protected:
	list_spVer _lv;
	bool _hole;
public:
	//constructors/destructors
	Poly() :
			_lv(), _hole(false) {
	}
	~Poly() {

	}

	Poly(const Self &src) {
		this->_hole = src._hole;
		this->_lv = src._lv;
	}
	Self& operator=(const Self& src) {
		this->_hole = src._hole;
		this->_lv = src._lv;
		return *this;
	}

	st size() {
		return _lv.size();
	}

	void push_back(spVer ver){
		_lv.push_back(ver);
	}

};

#endif /* LIB_TS_TS_POLYGON_H_ */
