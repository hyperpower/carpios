/*
 * _point_chain.hpp
 *
 *  Created on: Jul 6, 2017
 *      Author: zhou
 */

#ifndef _POINT_CHAIN_HPP_
#define _POINT_CHAIN_HPP_

#include "geometry_define.hpp"
#include "_point.hpp"
#include "_contour.hpp"
#include "../algebra/array_list.hpp"
#include "_segment.hpp"
#include "_polygon.hpp"
#include <array>
#include <vector>
#include <limits>
#include <list>
#include <fstream>
#include <algorithm>


namespace carpio {

template<class TYPE, St DIM>
class PointChain_ {
public:
	static const St Dim = DIM;
	typedef Point_<TYPE, DIM> Point;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Segment_<TYPE, DIM>& ref_Segment;
	typedef PointChain_<TYPE, DIM>& PointChain;
	typedef TYPE vt;
	typedef typename std::list<Point>::iterator iterator;
protected:
	/** Linked point chain */
	std::list<Point> l;
	bool _closed; // is the chain closed, that is, is the first point is linked with the last one?
public:
	PointChain_() :
			l(), _closed(false) {
	}
	void init(const Segment& s) {
		l.push_back(s.ps());
		l.push_back(s.pe());
	}
	bool LinkSegment(const Segment& s) {
		if (s.ps() == l.front()) {
			if (s.pe() == l.back())
				_closed = true;
			else
				l.push_front(s.pe());
			return true;
		}
		if (s.pe() == l.back()) {
			if (s.ps() == l.front())
				_closed = true;
			else
				l.push_back(s.ps());
			return true;
		}
		if (s.pe() == l.front()) {
			if (s.ps() == l.back())
				_closed = true;
			else
				l.push_front(s.ps());
			return true;
		}
		if (s.ps() == l.back()) {
			if (s.pe() == l.front())
				_closed = true;
			else
				l.push_back(s.pe());
			return true;
		}
		return false;
	}
	bool LinkPointChain(PointChain& chain) {
		if (chain.l.front() == l.back()) {
			chain.l.pop_front();
			l.splice(l.end(), chain.l);
			return true;
		}
		if (chain.l.back() == l.front()) {
			l.pop_front();
			l.splice(l.begin(), chain.l);
			return true;
		}
		if (chain.l.front() == l.front()) {
			l.pop_front();
			std::reverse(chain.l.begin(), chain.l.end());
			l.splice(l.begin(), chain.l);
			return true;
		}
		if (chain.l.back() == l.back()) {
			l.pop_back();
			reverse(chain.l.begin(), chain.l.end());
			l.splice(l.end(), chain.l);
			return true;
		}
		return false;
	}
	bool closed() const {
		return _closed;
	}
	iterator begin() {
		return l.begin();
	}
	iterator end() {
		return l.end();
	}
	void clear() {
		l.clear();
	}
	unsigned int size() const {
		return l.size();
	}
};

template<class TYPE, St DIM>
class Connector_ {
public:
	static const St Dim = DIM;
	typedef Point_<TYPE, DIM> Point;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Segment_<TYPE, DIM>& ref_Segment;
	typedef TYPE vt;
	typedef Polygon_<TYPE> Polygon;
	//typedef ArrayListT<Point> ArrP;
	typedef PointChain_<TYPE, DIM> PointChain;
	typedef Contour_<TYPE> Contour;
	typedef typename std::list<PointChain>::iterator iterator;
protected:
	std::list<PointChain> openPolygons;
	std::list<PointChain> closedPolygons;
public:
	Connector_() :
			openPolygons(), closedPolygons() {
	}
	~Connector_() {
	}
	void add(const Segment& s) {
		iterator j = openPolygons.begin();
		while (j != openPolygons.end()) {
			if (j->LinkSegment(s)) {
				if (j->closed())
					closedPolygons.splice(closedPolygons.end(), openPolygons,
							j);
				else {
					iterator k = j;
					for (++k; k != openPolygons.end(); k++) {
						if (j->LinkPointChain(*k)) {
							openPolygons.erase(k);
							break;
						}
					}
				}
				return;
			}
			j++;
		}
		// The segment cannot be connected with any open polygon
		openPolygons.push_back(PointChain());
		openPolygons.back().init(s);
	}

	iterator begin() {
		return closedPolygons.begin();
	}
	iterator end() {
		return closedPolygons.end();
	}
	void clear() {
		closedPolygons.clear();
		openPolygons.clear();
	}
	unsigned int size() const {
		return closedPolygons.size();
	}

	void toPolygon(Polygon& p) {

		for (iterator it = begin(); it != end(); it++) {
			p.push_back(Contour());
		Contour& contour = p.back();
			for (typename PointChain::iterator it2 = it->begin(); it2 != it->end();
					it2++)
				contour.add(*it2);
		}
	}
};

}



#endif
