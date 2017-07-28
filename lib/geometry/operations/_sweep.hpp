/*
 * _sweep.hpp
 *
 *  Created on: Jul 6, 2017
 *      Author: zhou
 */

#ifndef _SWEEP_HPP_
#define _SWEEP_HPP_

#include "../geometry_define.hpp"
#include "_operation.hpp"
#include "_polygon_boolean.hpp"
#include <array>
#include <sstream>
#include <set>

namespace carpio {

template<typename TYPE, St DIM>
class SweepEvent_ {
public:
	typedef SweepEvent_<TYPE, DIM> SweepEvent;
	typedef Point_<TYPE, DIM> Point;
	typedef Segment_<TYPE, DIM> Segment;
	typedef Operation_<TYPE, DIM> Op;

	// point associated with the event
	Point p;
	// is the point the left endpoint of the segment (p, other->p)?
	bool left;
	// Polygon to which the associated segment belongs to
	PolygonType pl;
	// Event associated to the other endpoint of the segment
	SweepEvent *other;
	/**  Does the segment (p, other->p) represent an inside-outside transition in the polygon for a vertical ray from (p.x, -infinite) that crosses the segment? */
	bool inOut;
	EdgeType type;
	// Only used in "left" events. Is the segment (p, other->p) inside the other polygon?
	bool inside;
	// Only used in "left" events. Position of the event (line segment) in S
	typename std::set<SweepEvent*>::iterator* poss;

	/** Class constructor */
	SweepEvent_(const Point& pp, bool isleft, PolygonType apl, SweepEvent* o,
			EdgeType t = NORMAL) :
			p(pp), left(isleft), pl(apl), other(o), type(t), poss(0) {
	}
	/** Class destructor */
	~SweepEvent_() {
		delete poss;
	}
	/** Return the line segment associated to the SweepEvent */
	Segment segment() {
		return Segment(p, other->p);
	}
	/** Is the line segment (p, other->p) below point x */
	bool below(const Point& x) const {
		return (left) ?
				Op::SignedArea(p, other->p, x) > 0 :
				Op::SignedArea(other->p, p, x) > 0;
	}
	/** Is the line segment (p, other->p) above point x */
	bool above(const Point& x) const {
		return !below(x);
	}

	void show() const {
		const char* namesEventTypes[] = { " (NORMAL) ", " (NON_CONTRIBUTING) ",
				" (SAME_TRANSITION) ", " (DIFFERENT_TRANSITION) " };
		std::cout << " Point: " << p;
		if (other == nullptr) {
			std::cout << " Other point: " << "NULL";
		} else {
			std::cout << " Other point: " << other->p;
		}
		std::cout << (left ? " (Left) " : " (Right) ")
				<< (inside ? " (Inside) " : " (Outside) ")
				<< (inOut ? " (In-Out) " : " (Out-In) ") << "Type: "
				<< namesEventTypes[type] << " Polygon: "
				<< (pl == SUBJECT ? " (SUBJECT)" : " (CLIPPING)") << std::endl;
	}
};

template<typename TYPE, St DIM>
struct SweepEventComp_: public std::binary_function<SweepEvent_<TYPE, DIM>*,
		SweepEvent_<TYPE, DIM>*, bool> {
	typedef SweepEvent_<TYPE, DIM> SweepEvent;
	typedef SweepEventComp_<TYPE, DIM> SweepEventComp;

	typedef Operation_<TYPE, DIM> Op;
	bool operator()(SweepEvent* e1, SweepEvent* e2) {
		// Different x-coordinate
		// The event with lower x-coordinate is processed first
		if (e1->p.x() > e2->p.x()) {
			return true;
		}
		if (e2->p.x() > e1->p.x()) {
			return false;
		}
		// Different points, but same x-coordinate.
		// The event with lower y-coordinate is processed first
		if (e1->p != e2->p) {
			return e1->p.y() > e2->p.y();
		}
		// Same point, but one is a left endpoint and the other a right endpoint.
		// The right endpoint is processed first
		if (e1->left != e2->left) {
			return e1->left;
		}
		// Same point, both events are left endpoints or both are right endpoints.
		// The event associate to the bottom segment is processed first
		return e1->above(e2->other->p);
	}
};

template<typename TYPE, St DIM>
struct SegmentComp_: public std::binary_function<SweepEvent_<TYPE, DIM>*,
		SweepEvent_<TYPE, DIM>*, bool> {
	typedef SweepEvent_<TYPE, DIM> SweepEvent;
	typedef Operation_<TYPE, DIM> Op;
	typedef SweepEventComp_<TYPE, DIM> SweepEventComp;
	bool operator()(SweepEvent* e1, SweepEvent* e2) {
		if (e1 == e2) {
			return false;
		}
		if (Op::SignedArea(e1->p, e1->other->p, e2->p) != 0
				|| Op::SignedArea(e1->p, e1->other->p, e2->other->p) != 0) {
			// Segments are not collinear
			// If they share their left endpoint use the right endpoint to sort
			if (e1->p == e2->p)
				return e1->below(e2->other->p);

			// Different points
			SweepEventComp comp;
			if (comp(e1, e2)) // has the line segment associated to e1 been inserted into S after the line segment associated to e2 ?
				return e2->above(e1->p);
			// The line segment associated to e2 has been inserted into S after the line segment associated to e1
			return e1->below(e2->p);
		}
		// Segments are collinear
		if (e1->pl != e2->pl) {
			return e1->pl < e2->pl;
		}
		// Just a consistent criterion is used
		if (e1->p == e2->p) {
			return e1 < e2;
		}
		SweepEventComp comp;
		return comp(e1, e2);
	}
};

}

#endif /* LIB_GEOMETRY__SWEEP_HPP_ */
