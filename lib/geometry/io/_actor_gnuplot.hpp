/*
 * _actor_gnuplot.hpp
 *
 *  Created on: Jul 5, 2017
 *      Author: zhou
 */

#ifndef _ACTOR_GNUPLOT_GEOMETRY_HPP_
#define _ACTOR_GNUPLOT_GEOMETRY_HPP_

#include "../geometry_define.hpp"
#include <array>
#include "../objects/_objects.hpp"
#include "../../io/gnuplot.h"
#include <memory>
#include <cmath>

namespace carpio {

template<typename TYPE, St DIM>
class GPA_Geometry_ {
public:
	static const St Dim = DIM;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Point;
	typedef Point_<TYPE, DIM>& ref_Point;
	typedef const Point_<TYPE, DIM>& const_ref_Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Segment_<TYPE, DIM>& ref_Segment;
	typedef const Segment_<TYPE, DIM>& const_ref_Segment;

	typedef Contour_<TYPE> Contour;
	typedef PointChain_<TYPE, DIM> PointChain;
	typedef Polygon_<TYPE> Polygon;

	typedef std::shared_ptr<Gnuplot_actor> spActor;
	typedef std::list<spActor> list_spActor;

public:
	static spActor Lines(const Contour& contour, int base_idx = 0,
			int color_idx = -1) {
		spActor actor = spActor(new Gnuplot_actor());
		actor->command() = "using 1:2:3 title \"\" ";
		actor->style() = "with lines lc variable";
		if (contour.empty()) {
			actor->data().push_back("");
			return actor;
		}
		for (St i = 0; i < contour.size_vertexs(); ++i) {
			const Point& p = contour.v(i);
			if (color_idx >= 0) {
				actor->data().push_back(ToString(p.x(), p.y(), color_idx, " "));
			} else {
				actor->data().push_back(
						ToString(p.x(), p.y(), i + base_idx, " "));
			}
		}
		const Point& pstart = contour.v(0);
		if (color_idx >= 0) {
			actor->data().push_back(
					ToString(pstart.x(), pstart.y(), color_idx, " "));
		} else {
			actor->data().push_back(
					ToString(pstart.x(), pstart.y(), base_idx, " "));
		}
		actor->data().push_back("");
		return actor;
	}

	static spActor Lines(const PointChain& points, int base_idx = 0,
			int color_idx = -1) {
		spActor actor = spActor(new Gnuplot_actor());
		actor->command() = "using 1:2:3 title \"\" ";
		actor->style() = "with lines lc variable";
		if (points.empty()) {
			actor->data().push_back("");
			return actor;
		}
		int i = 0;
		for (auto& p : points) {
			if (color_idx >= 0) {
				actor->data().push_back(ToString(p.x(), p.y(), color_idx, " "));
			} else {
				actor->data().push_back(
						ToString(p.x(), p.y(), i + base_idx, " "));
			}
			i++;
		}
		// If point chain is closed, the first point should be add.
		if (points.closed()) {
			const Point& pstart = *(points.begin());
			if (color_idx >= 0) {
				actor->data().push_back(
						ToString(pstart.x(), pstart.y(), color_idx, " "));
			} else {
				actor->data().push_back(
						ToString(pstart.x(), pstart.y(), base_idx, " "));
			}
		}
		actor->data().push_back("");
		return actor;
	}

	static spActor Points(const Point& p, int color_idx = -1) {
		spActor actor = spActor(new Gnuplot_actor());
		actor->command() = "using 1:2:3 title \"\" ";
		actor->style() = "with points lc variable";

		if (color_idx >= 0) {
			actor->data().push_back(ToString(p.x(), p.y(), color_idx, " "));
		} else {
			actor->data().push_back(ToString(p.x(), p.y(), 0, " "));
		}

		actor->data().push_back("");
		return actor;
	}

};

}

#endif /* _ACTOR_GNUPLOT_HPP_ */
