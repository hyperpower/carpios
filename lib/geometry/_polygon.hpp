#ifndef _POLYGON_HPP_
#define _POLYGON_HPP_

#include "geometry_define.hpp"
#include "_operation.hpp"
#include "_point.hpp"
#include "_contour.hpp"
#include "_bbox.hpp"
#include "../algebra/array_list.hpp"
#include "_segment.hpp"
#include <array>
#include <vector>
#include <limits>
#include <list>
#include <fstream>

namespace carpio {

template<typename TYPE, St DIM>
class Operation_;

template<class TYPE>
class Contour_;

template<class TYPE>
class Polygon_ {
public:
	typedef Contour_<TYPE> Contour;
	typedef typename std::vector<Contour>::iterator iterator;
	typedef typename std::vector<Contour>::const_iterator const_iterator;
	typedef Point_<TYPE, 2> Point;
	typedef Point_<TYPE, 2>& ref_Point;

	typedef Segment_<TYPE, 2> Segment;
	typedef Segment_<TYPE, 2>& ref_Segment;

protected:
	/** Set of contours conforming the polygon */
	std::vector<Contour> contours;
public:
	Polygon_() :
			contours() {
	}
	Polygon_(const std::string& filename) {
		std::ifstream f(filename.c_str());
		if (!(f.is_open())) {
			std::cerr << "Error opening " << filename << '\n';
			exit(1);
		}
		f >> *this;
		if (!(f.eof()))
			std::cerr << "An error reading file " << filename << " happened\n";
	}
	void read_file(const std::string& filename) {
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
		ASSERT(p < contours.size());
		return contours[p];
	}
	Contour& operator[](unsigned int p) {
		ASSERT(p < contours.size());
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
	void boundingbox(Point& min, Point& max) {
		min.x() = min.y() = std::numeric_limits<double>::max();
		max.x() = max.y() = -std::numeric_limits<double>::max();
		Point mintmp;
		Point maxtmp;
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

	const_iterator begin() const {
		return contours.begin();
	}
	const_iterator end() const {
		return contours.end();
	}

	void computeHoles() {
		SHOULD_NOT_REACH;
	}

};
template<class TYPE>
std::ostream& operator<<(std::ostream& o, Polygon_<TYPE>& p) {
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
std::istream& operator>>(std::istream& is, Polygon_<TYPE>& p) {
	typedef Contour_<TYPE> Contour;
	typedef typename Contour_<TYPE>::Point Point;
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
			contour.add(Point(px, py));
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

///*
// * Function out of class
// */
//template<typename VALUE>
//void CreatCircle(Polygon_<VALUE>& s, VALUE x0, VALUE y0, VALUE r, int n) {
//	ASSERT(n >= 3);
//	Float pi = 3.141592653589793238;
//	typedef typename Polygon_<VALUE>::ArrP ArrPoint;
//	typedef typename Polygon_<VALUE>::Point Poi;
//	ArrPoint arrp;
//	for (int i = 0; i < n; i++) {
//		Float x = x0 + r * cos(2. * pi / float(n) * i);
//		Float y = y0 + r * sin(2. * pi / float(n) * i);
//		arrp.push_back(Poi(x, y));
//	}
//	s.reconstruct(arrp);
//}
//
//template<typename VALUE>
//void CreatCube(Polygon_<VALUE>& s, VALUE x0, VALUE y0, VALUE x1, VALUE y1) {
//	ASSERT(x1 > x0);
//	ASSERT(y1 > y0);
//	typedef typename Polygon_<VALUE>::ArrP ArrPoint;
//	typedef typename Polygon_<VALUE>::Point Poi;
//	ArrPoint arrp;
//	VALUE dx = x1 - x0;
//	VALUE dy = y1 - y0;
//	arrp.push_back(Poi(x0, y0));
//	arrp.push_back(Poi(x0 + dx, y0));
//	arrp.push_back(Poi(x1, y1));
//	arrp.push_back(Poi(x0, y0 + dy));
//	s.reconstruct(arrp);
//}

}

#endif
