/*
 * _contour.hpp
 *
 *  Created on: Jul 5, 2017
 *      Author: zhou
 */

#ifndef _CONTOUR_HPP_
#define _CONTOUR_HPP_


#include "geometry_define.hpp"
#include "_point.hpp"
#include "_relation.hpp"
#include <array>

namespace carpio {

// the new polygon class ======================================================
template<class TYPE>
class Contour_ {
public:
	typedef Contour_<TYPE> Self;
	typedef Point_<TYPE, 2> Poi;
	typedef Point_<TYPE, 2>& ref_Poi;
	typedef const Point_<TYPE, 2>& const_ref_Poi;
	typedef typename std::vector<Poi>::iterator iterator;
	typedef typename std::vector<Poi>::const_iterator const_iterator;
	typedef typename std::vector<Poi>::size_type St;
	typedef Segment_<TYPE, 2> Segment;
	typedef Segment_<TYPE, 2>& ref_Segment;
	typedef TYPE Vt;

	typedef ArrayListT<Poi> ArrP;
	typedef Operation_<TYPE, 2> Op;

protected:

	// Set of points conforming the external contour
	std::vector<Poi> _vertices;
	// Holes of the contour. They are stored as the indexes of the holes in a polygon class
	std::vector<St> _holes;
	// is the contour an external contour? (i.e., is it not a hole?)
	bool _external;
	// this will be false, before calling the function conterclockwise
	bool _precomputedCC;
	// is count clock wise
	bool _CC;
public:
	Contour_() :
			_vertices(), _holes(), _external(true), _precomputedCC(false), _CC(
					false) {
	}

	Self& operator=(const Self&a) {
		if (this == &a) {
			return *this;
		} else {
			this->_vertices = a._vertices;
			this->_holes = a._holes;
			this->_external = a._external;
			this->_precomputedCC = a._precomputedCC;
			this->_CC = a._CC;
		}
		return *this;
	}

	Vt area() const {
		if (empty()) {
			return 0.0;
		}
		Vt s = 0.0;
		for (St i = 1; i < _vertices.size() - 1; i++) {
			s = s + Op::Cro(_vertices[i + 1], _vertices[i], _vertices[0]); // det to cro
		}
		return std::abs(s) / 2.0;
	}
	bool empty() const {
		if (_vertices.size() == 0) {
			return true;
		} else {
			return false;
		}
	}

	/** Get the p-th vertex of the external contour */
	Poi& vertex(unsigned p) {
		return _vertices[p];
	}
	const Poi& vertex(unsigned p) const {
		return _vertices[p];
	}
	Segment segment(unsigned p) const {
		return (p == nvertices() - 1) ?
				Segment(_vertices.back(), _vertices.front()) :
				Segment(_vertices[p], _vertices[p + 1]);
	}
	/** Number of vertices and edges */
	unsigned nvertices() const {
		return _vertices.size();
	}
	unsigned nedges() const {
		return _vertices.size();
	}
	/** Get the bounding box */
	void boundingbox(Poi& min, Poi& max) {
		min.x() = -std::numeric_limits<double>::max();
		min.y() = -std::numeric_limits<double>::max();
		max.x() = max.y() = std::numeric_limits<double>::max();
		iterator i = begin();
		while (i != end()) {
			if (i->x() < min.x())
				min.x() = i->x();
			if (i->x() > max.x())
				max.x() = i->x();
			if (i->_y < min.y())
				min.y() = i->_y;
			if (i->_y > max.y())
				max.y() = i->_y;
			++i;
		}
	}
	/** Return if the contour is counterclockwise oriented */
	bool counterclockwise() {
		if (_precomputedCC)
			return _CC;
		_precomputedCC = true;
		double area = 0.0;
		for (unsigned int c = 0; c < nvertices() - 1; c++)
			area += vertex(c).x() * vertex(c + 1).y()
					- vertex(c + 1).x() * vertex(c).y();
		area += vertex(nvertices() - 1).x() * vertex(0).y()
				- vertex(0).x() * vertex(nvertices() - 1).y();
		return _CC = area >= 0.0;
	}
	/** Return if the contour is clockwise oriented */
	bool clockwise() {
		return !counterclockwise();
	}
	void changeOrientation() {
		reverse(_vertices.begin(), _vertices.end());
		_CC = !_CC;
	}
	void setClockwise() {
		if (counterclockwise())
			changeOrientation();
	}
	void setCounterClockwise() {
		if (clockwise())
			changeOrientation();
	}

	void move(Vt x, Vt y) {
		for (St i = 0; i < _vertices.size(); i++) {
			_vertices[i].x() += x;
			_vertices[i].y() += y;
		}
	}
	void add(const Poi& s) {
		_vertices.push_back(s);
	}
	void erase(iterator i) {
		_vertices.erase(i);
	}
	void clear() {
		_vertices.clear();
		_holes.clear();
	}
	iterator begin() {
		return _vertices.begin();
	}
	iterator end() {
		return _vertices.end();
	}
	Poi& back() {
		return _vertices.back();
	}
	const Poi& back() const {
		return _vertices.back();
	}
	void addHole(unsigned ind) {
		_holes.push_back(ind);
	}
	unsigned nholes() const {
		return _holes.size();
	}
	unsigned hole(unsigned p) const {
		return _holes[p];
	}
	bool external() const {
		return _external;
	}
	void setExternal(bool e) {
		_external = e;
	}
	Vt max(int dim) const { //
		ASSERT(_vertices.size() > 0);
		Float max = _vertices[0].val(dim);
		for (St i = 1; i < _vertices.size(); i++) {
			if (_vertices[i].val(dim) > max) {
				max = _vertices[i].val(dim);
			}
		}
		return max;
	}

	Vt min(int dim) const { //
		ASSERT(_vertices.size() > 0);
		Float min = _vertices[0].val(dim);
		for (St i = 1; i < _vertices.size(); i++) {
			if (_vertices[i].val(dim) < min) {
				min = _vertices[i].val(dim);
			}
		}
		return min;
	}

	St find_closest_vertex(const Poi& p) const {
		ASSERT(_vertices.size() > 0);
		St idx = 0;
		Vt mindis = Op::Distance(_vertices[0], p);
		for (St i = 1; i < _vertices.size(); i++) {
			Vt dis = Op::Distance(_vertices[i], p);
			if (dis < mindis) {
				idx = i;
				mindis = dis;
			}
		}
		return idx;
	}

	St find_closest_vertex(const Vt& x, const Vt& y) const {
		Poi p(x, y);
		return this->find_closest_vertex(p);
	}

	/*
	 * special function
	 * aix = x, y
	 * v   = (a value on coordinate)
	 * l_seg_idx (return)
	 *
	 * Find all the segments across x=v, y=v
	 */
	void find_seg_across(std::list<St>& l_seg_idx, int aix, Vt val) {
		l_seg_idx.clear();
		St i = 0;
		int flag = GEL(val, vertex(i).val(aix));
		int nf;
		for (++i; i < this->size_vertexs(); i++) {
			Vt pv = vertex(i).val(aix);
			nf = GEL(val, pv);
			if (flag != nf) {
				l_seg_idx.push_back(i - 1);
				flag = nf;
			}
		}
		nf = GEL(val, vertex(0).val(aix));
		if (nf != flag) {
			l_seg_idx.push_back(size_vertexs() - 1);
		}
	}

	inline St size_vertexs() const {
		return _vertices.size();
	}
	inline St size_segments() const {  //
		return _vertices.size();
	}

};
template<class TYPE>
std::ostream& operator<<(std::ostream& o, Contour_<TYPE>& c) {
	o << c.nvertices() << std::endl;
	typename Contour_<TYPE>::iterator i = c.begin();
	while (i != c.end()) {
		o << "    (" << i->x() << ", " << i->y() << ")" << std::endl;
		++i;
	}
	return o;
}





}


#endif
