/*
 * _point_chain.hpp
 *
 *  Created on: Jul 6, 2017
 *      Author: zhou
 */

#ifndef _POINT_CHAIN_HPP_
#define _POINT_CHAIN_HPP_

#include "../geometry_define.hpp"
#include "_point.hpp"
#include "_contour.hpp"
#include "../operations/_operation.hpp"
#include "../../algebra/array_list.hpp"
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
	typedef Segment& ref_Segment;
	typedef PointChain_<TYPE, DIM>& PointChain;
	typedef TYPE vt;
	typedef typename std::list<Point>::iterator iterator;
	typedef typename std::list<Point>::const_iterator const_iterator;
	typedef Operation_<TYPE, DIM> Op;
	typedef Point value_type;
	protected:
	/** Linked point chain */
	std::list<Point> _lpoints;
	bool _closed; // is the chain closed, that is, is the first point is linked with the last one?
public:
	PointChain_() :
			_lpoints(), _closed(false) {
	}
	PointChain_(const Point& p1) :
			_closed(false) {
		_lpoints.push_back(p1);
	}
	PointChain_(const Point& p1, const Point& p2) :
			_closed(false) {
		_lpoints.push_back(p1);
		_lpoints.push_back(p2);
	}
	PointChain_(const Point& p1, const Point& p2, const Point& p3, bool close =
			true) :
			_closed(close) {
		_lpoints.push_back(p1);
		_lpoints.push_back(p2);
		_lpoints.push_back(p3);
	}
	PointChain_(const Point& p1, const Point& p2, const Point& p3,
			const Point& p4, bool close = true) :
			_closed(close) {
		_lpoints.push_back(p1);
		_lpoints.push_back(p2);
		_lpoints.push_back(p3);
		_lpoints.push_back(p4);
	}

	template<class Container_Point>
	PointChain_(const Container_Point& ver,
			bool closed = true) : _closed(closed){
		std::copy(ver.begin(), ver.end(), std::back_inserter(_lpoints));
	}

	void init(const Segment& s) {
		_lpoints.push_back(s.ps());
		_lpoints.push_back(s.pe());
	}
	bool LinkSegment(const Segment& s) {
		if (s.ps() == _lpoints.front()) {
			if (s.pe() == _lpoints.back())
				_closed = true;
			else
				_lpoints.push_front(s.pe());
			return true;
		}
		if (s.pe() == _lpoints.back()) {
			if (s.ps() == _lpoints.front())
				_closed = true;
			else
				_lpoints.push_back(s.ps());
			return true;
		}
		if (s.pe() == _lpoints.front()) {
			if (s.ps() == _lpoints.back())
				_closed = true;
			else
				_lpoints.push_front(s.ps());
			return true;
		}
		if (s.ps() == _lpoints.back()) {
			if (s.pe() == _lpoints.front())
				_closed = true;
			else
				_lpoints.push_back(s.pe());
			return true;
		}
		return false;
	}
	bool LinkPointChain(PointChain& chain) {
		if (chain._lpoints.front() == _lpoints.back()) {
			chain._lpoints.pop_front();
			_lpoints.splice(_lpoints.end(), chain._lpoints);
			return true;
		}
		if (chain._lpoints.back() == _lpoints.front()) {
			_lpoints.pop_front();
			_lpoints.splice(_lpoints.begin(), chain._lpoints);
			return true;
		}
		if (chain._lpoints.front() == _lpoints.front()) {
			_lpoints.pop_front();
			std::reverse(chain._lpoints.begin(), chain._lpoints.end());
			_lpoints.splice(_lpoints.begin(), chain._lpoints);
			return true;
		}
		if (chain._lpoints.back() == _lpoints.back()) {
			_lpoints.pop_back();
			std::reverse(chain._lpoints.begin(), chain._lpoints.end());
			_lpoints.splice(_lpoints.end(), chain._lpoints);
			return true;
		}
		return false;
	}
	void set_close() {
		_closed = true;
	}
	bool closed() const {
		return _closed;
	}
	bool empty() const {
		return _lpoints.empty();
	}
	iterator begin() {
		return _lpoints.begin();
	}
	iterator end() {
		return _lpoints.end();
	}

	const_iterator begin() const {
		return _lpoints.begin();
	}
	const_iterator end() const {
		return _lpoints.end();
	}

	void push_back(const Point& p) {
		_lpoints.push_back(p);
	}
	void pop_back() {
		_lpoints.pop_back();
	}
	void clear() {
		_lpoints.clear();
	}
	unsigned int size() const {
		return _lpoints.size();
	}

	bool is_simple() const {
		return Op::IsSimple(_lpoints.begin(), _lpoints.end(), closed());
	}

	double perimeter() const {
		if (empty()) {
			return 0;
		}
		double res = 0;
		const_iterator iter_end = std::prev(_lpoints.end(), 1);
		for (const_iterator iter = _lpoints.begin();
				iter != iter_end; ++iter) {
			const_iterator iter1 = std::next(iter);
			res += Op::Distance(*iter, *iter1);
		}
		if (closed()) {
			res += Op::Distance(*iter_end, *(_lpoints.begin()));
		}
		return res;
	}

	void show() const {
		int count = 0;
		for (auto& p : _lpoints) {
			std::cout << "i = " << count << " " << p << std::endl;
			count++;
		}
	}

}
;

template<class TYPE, St DIM>
class Connector_ {
public:
	static const St Dim = DIM;
	typedef Point_<TYPE, DIM> Point;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Segment& ref_Segment;
	typedef TYPE vt;
	typedef Polygon_<TYPE> Polygon;
	//typedef ArrayListT<Point> ArrP;
	typedef PointChain_<TYPE, DIM> PointChain;
	typedef Contour_<TYPE> Contour;
	typedef typename std::list<PointChain>::iterator iterator;
	protected:
	std::list<PointChain> _lopen;
	std::list<PointChain> _lclosed;
	public:
	Connector_() :
			_lopen(), _lclosed() {
	}
	~Connector_() {
	}
	void add(const Segment& s) {
		iterator j = _lopen.begin();
		while (j != _lopen.end()) {
			if (j->LinkSegment(s)) {
				if (j->closed())
					_lclosed.splice(_lclosed.end(), _lopen, j);
				else {
					iterator k = j;
					for (++k; k != _lopen.end(); k++) {
						if (j->LinkPointChain(*k)) {
							_lopen.erase(k);
							break;
						}
					}
				}
				return;
			}
			j++;
		}
		// The segment cannot be connected with any open polygon
		_lopen.push_back(PointChain());
		_lopen.back().init(s);
	}

	iterator begin() {
		return _lclosed.begin();
	}
	iterator end() {
		return _lclosed.end();
	}
	void clear() {
		_lclosed.clear();
		_lopen.clear();
	}
	unsigned int size() const {
		return _lclosed.size();
	}

	void toPolygon(Polygon& p) {
		for (iterator it = begin(); it != end(); it++) {
			p.push_back(Contour());
			Contour& contour = p.back();
			for (typename PointChain::iterator it2 = it->begin();
					it2 != it->end(); it2++)
				contour.add(*it2);
		}
	}
};

}

#endif
