#ifndef _POLYGON_HPP_
#define _POLYGON_HPP_

#include "geometry_define.hpp"
#include "_point.hpp"
#include "../algebra/array_list.hpp"
#include "_segment.hpp"
#include <array>
#include <vector>
#include <limits>
#include <list>
#include <fstream>

namespace carpio {

template<typename TYPE>
class Polygon_ {
public:
	typedef Point_<TYPE, 2> Point;
	typedef Point_<TYPE, 2>& ref_Point;
	typedef const Point_<TYPE, 2>& const_ref_Point;
	typedef Segment_<TYPE, 2> Segment;
	typedef Segment_<TYPE, 2>& ref_Segment;
	typedef TYPE vt;
	typedef ArrayListT<Point> ArrP;

public:
	/*
	 *  Contructor
	 */
	Polygon_() :
			_arrp() {
	}
	Polygon_(const Polygon_ &a) {
		this->_arrp = a._arrp;
	}
	Polygon_(const ArrP &a) {
		ASSERT(a.size() >= 3);
		this->_arrp = a;
		_trim_same_points();
	}
	void reconstruct(const ArrP & a) {
		ASSERT(a.size() >= 3);
		this->_arrp = a;
		_trim_same_points();
	}

	Polygon_& operator=(const Polygon_ &a) {
		if (this == &a) {
			return *this;
		} else {
			this->_arrp = a._arrp;
		}
		return *this;
	}
	vt area() const {
		if (empty()) {
			return 0.0;
		}
		vt s = 0.0;
		for (St i = 1; i < _arrp.size() - 1; i++) {
			s = s + Cro(_arrp[i + 1], _arrp[i], _arrp[0]); // det to cro
		}
		return Abs(s) / 2.0;
	}
	void clear() {
		_arrp.resize(0);
	}
	bool empty() const {
		if (_arrp.size() == 0) {
			return true;
		} else {
			return false;
		}
	}
	void show(const std::string& name = "") {
		std::cout << "Polygon  ---- " << name << "\n";
		for (St i = 0; i < _arrp.size(); i++) {
			std::cout << "> AP[ " << i << " ]=( " << _arrp[i].x() << " , "
					<< _arrp[i].y() << " )" << "\n";
		}
	}
	void reverse() {
		if (empty()) {
			return;
		} else {
			_arrp.reverse();
		}
	}

	inline St size_vertexs() const {
		return _arrp.size();
	}
	inline St size_segments() const {  //
		return _arrp.size();
	}
	const_ref_Point vertex(St i) const {
		ASSERT(i < _arrp.size());
		return _arrp[i];
	}
	ref_Point vertex(St i) {
		ASSERT(i < _arrp.size());
		return _arrp[i];
	}
	Segment segment(St i) const {
		return Segment(_arrp[i], _arrp[(i + 1) % _arrp.size()]);
	}

	vt max_x() const { //
		ASSERT(_arrp.size() > 0);
		Float max = _arrp[0].x();
		for (St i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].x() > max) {
				max = _arrp[i].x();
			}
		}
		return max;
	}
	vt min_x() const { //
		ASSERT(_arrp.size() > 0);
		Float min = _arrp[0].x();
		for (St i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].x() < min) {
				min = _arrp[i].x();
			}
		}
		return min;
	}
	vt max_y() const {
		ASSERT(_arrp.size() > 0);
		Float max = _arrp[0].y();
		for (St i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].y() > max) {
				max = _arrp[i].y();
			}
		}
		return max;
	}
	vt min_y() const { //
		ASSERT(_arrp.size() > 0);
		Float min = _arrp[0].y();
		for (St i = 1; i < _arrp.size(); i++) {
			if (_arrp[i].y() < min) {
				min = _arrp[i].y();
			}
		}
		return min;
	}

	St find_closest_vertex(const Point& p) const {
		ASSERT(_arrp.size() > 0);
		St idx = 0;
		vt mindis = Distance(_arrp[0], p);
		for (St i = 1; i < _arrp.size(); i++) {
			vt dis = Distance(_arrp[i], p);
			if (dis < mindis) {
				idx = i;
				mindis = dis;
			}
		}
		return idx;
	}

	St find_closest_vertex(const vt& x, const vt& y) const {
		Point p(x, y);
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
	void find_seg_across(std::list<St>& l_seg_idx, Axes aix, vt val) {
		l_seg_idx.clear();
		St i = 0;
		int flag = GEL(val, vertex(i).val(aix));
		int nf;
		for (++i; i < size_vertexs(); i++) {
			vt pv = vertex(i).val(aix);
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

	void find_seg_connect_to_vertex(std::list<St>& l_seg_idx,
			St ver_idx) const {
		St n = this->size_vertexs();
		ASSERT(ver_idx < n);
		l_seg_idx.clear();
		St prev = ver_idx - 1;
		if (ver_idx == 0) {
			prev = n - 1;
		}
		l_seg_idx.push_back(prev);
		l_seg_idx.push_back(ver_idx);
	}
protected:
	void _trim_same_points() {
		for (int i = 0; i < _arrp.size() - 1; i++) {
			if (_arrp[i] == _arrp[i + 1]) {
				_arrp.erase(i);
				i--;
			}
		}
		if (_arrp[0] == _arrp[_arrp.size() - 1]) {
			_arrp.pop_back();
		}
	}

	/*
	 * data
	 */
	ArrP _arrp;
};
/*
 * Function out of class
 */
template<typename VALUE>
void CreatCircle(Polygon_<VALUE>& s, VALUE x0, VALUE y0, VALUE r, int n) {
	ASSERT(n >= 3);
	Float pi = 3.141592653589793238;
	typedef typename Polygon_<VALUE>::ArrP ArrPoint;
	typedef typename Polygon_<VALUE>::Point Poi;
	ArrPoint arrp;
	for (int i = 0; i < n; i++) {
		Float x = x0 + r * cos(2. * pi / float(n) * i);
		Float y = y0 + r * sin(2. * pi / float(n) * i);
		arrp.push_back(Poi(x, y));
	}
	s.reconstruct(arrp);
}

template<typename VALUE>
void CreatCube(Polygon_<VALUE>& s, VALUE x0, VALUE y0, VALUE x1, VALUE y1) {
	ASSERT(x1 > x0);
	ASSERT(y1 > y0);
	typedef typename Polygon_<VALUE>::ArrP ArrPoint;
	typedef typename Polygon_<VALUE>::Point Poi;
	ArrPoint arrp;
	VALUE dx = x1 - x0;
	VALUE dy = y1 - y0;
	arrp.push_back(Poi(x0, y0));
	arrp.push_back(Poi(x0 + dx, y0));
	arrp.push_back(Poi(x1, y1));
	arrp.push_back(Poi(x0, y0 + dy));
	s.reconstruct(arrp);
}

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

	typedef Operation_<TYPE, 2> Operation;

protected:

	/** Set of points conforming the external contour */
	std::vector<Poi> _vertices;
	/** Holes of the contour. They are stored as the indexes of the holes in a polygon class */
	std::vector<St> _holes;
	bool _external; // is the contour an external contour? (i.e., is it not a hole?)
	bool _precomputedCC; // this will be false, before calling the function conterclockwise
	bool _CC;             // is count clock wise
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
			s = s
					+ Operation::Cro(_vertices[i + 1], _vertices[i],
							_vertices[0]); // det to cro
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
		Vt mindis = Operation::Distance(_vertices[0], p);
		for (St i = 1; i < _vertices.size(); i++) {
			Vt dis = Operation::Distance(_vertices[i], p);
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

template<class TYPE>
class Polygons_ {
public:
	typedef Contour_<TYPE> Contour;
	typedef typename std::vector<Contour>::iterator iterator;
	typedef Point_<TYPE, 2> Poi;
	typedef Point_<TYPE, 2>& ref_Poi;

protected:
	/** Set of contours conforming the polygon */
	std::vector<Contour> contours;
public:
	Polygons_() :
			contours() {
	}
	Polygons_(const std::string& filename) {
		std::ifstream f(filename.c_str());
		if (!(f.is_open())) {
			std::cerr << "Error opening " << filename << '\n';
			exit(1);
		}
		f >> *this;
		if (!(f.eof()))
			std::cerr << "An error reading file " << filename << " happened\n";
	}
	/** Get the p-th contour */
	Contour& contour(unsigned p) {
		return contours[p];
	}
	Contour& operator[](unsigned int p) {
		return contours[p];
	}
	/** Number of contours */
	St ncontours() const {
		return contours.size();
	}
	/** Number of vertices */
	St nvertices() const {
		St nv = 0;
		for (St i = 0; i < ncontours(); i++)
			nv += contours[i].nvertices();
		return nv;
	}
	/** Get the bounding box */
	void boundingbox(Poi& min, Poi& max) {
		min.x() = min.y() = std::numeric_limits<double>::max();
		max.x() = max.y() = -std::numeric_limits<double>::max();
		Poi mintmp;
		Poi maxtmp;
		for (unsigned int i = 0; i < ncontours(); i++) {
			contours[i].boundingbox(mintmp, maxtmp);
			if (mintmp.x() < min.x())
				min.x() = mintmp.x();
			if (maxtmp.x() > max.x())
				max.x() = maxtmp.x();
			if (mintmp.y() < min.y())
				min.y() = mintmp.y();
			if (maxtmp.y() > max.y())
				max.y() = maxtmp.y();
		}
	}

	void move(double x, double y) {
		for (St i = 0; i < contours.size(); i++)
			contours[i].move(x, y);
	}

	void push_back(const Contour& c) {
		contours.push_back(c);
	}
	Contour& back() {
		return contours.back();
	}
	const Contour& back() const {
		return contours.back();
	}
	void pop_back() {
		contours.pop_back();
	}
	void erase(iterator i) {
		contours.erase(i);
	}
	void clear() {
		contours.clear();
	}

	iterator begin() {
		return contours.begin();
	}
	iterator end() {
		return contours.end();
	}
	void computeHoles() {
		SHOULD_NOT_REACH;
	}
};
template<class TYPE>
std::ostream& operator<<(std::ostream& o, Polygons_<TYPE>& p) {
	o << p.ncontours() << std::endl;
	for (unsigned int i = 0; i < p.ncontours(); i++)   // write the contours
		o << p.contour(i);
	for (unsigned int i = 0; i < p.ncontours(); i++) { // write the holes of every contour
		if (p.contour(i).nholes() > 0) {
			o << i << ": ";
			for (unsigned int j = 0; j < p.contour(i).nholes(); j++)
				o << p.contour(i).hole(j)
						<< (j == p.contour(i).nholes() - 1 ? '\n' : ' ');
		}
	}
	return o;
}
template<class TYPE>
std::istream& operator>>(std::istream& is, Polygons_<TYPE>& p) {
	typedef Contour_<TYPE> Contour;
	typedef typename Contour_<TYPE>::Poi Poi;
	// read the contours
	int ncontours;
	double px,py;
	is >> ncontours;
	for (int i = 0; i < ncontours; i++) {
		int npoints;
		is >> npoints;
		p.push_back(Contour());
		Contour& contour = p.back();
		for (int j = 0; j < npoints; j++) {
			is >> px >> py;
			if (j > 0 && px == contour.back().x() && py == contour.back().y())
				continue;
			if (j == npoints - 1 && px == contour.vertex(0).x()
					&& py == contour.vertex(0).y())
				continue;
			contour.add(Poi(px, py));
		}
		if (contour.nvertices() < 3) {
			p.pop_back();
			continue;
		}
	}
	// read holes information
	int contourId;
	char aux;
	std::string restOfLine;
	while (is >> contourId) {
		is >> aux; // read the character :
		if (aux != ':')
			break;
		std::getline(is, restOfLine);
		std::istringstream iss(restOfLine);
		int hole;
		while (iss >> hole) {
			p[contourId].addHole(hole);
			p[hole].setExternal(false);
		}
		if (!iss.eof())
			break;
	}
	return is;
}

}

#endif
