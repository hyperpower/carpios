/*
 * _polygon_boolean.hpp
 *
 *  Created on: Jun 19, 2017
 *      Author: zhou
 */

#ifndef _POLYGON_BOOLEAN_HPP_
#define _POLYGON_BOOLEAN_HPP_

#include "geometry_define.hpp"
#include "_point.hpp"
#include "_point_chain.hpp"
#include "../algebra/array_list.hpp"
#include "_segment.hpp"
#include "_operation.hpp"
#include "_polygon.hpp"
#include "_sweep.hpp"

#include <array>
#include <vector>
#include <limits>
#include <list>
#include <set>
#include <fstream>
#include <queue>

#define _DEBUG_

namespace carpio {

enum BooleanOpType {
	INTERSECTION, UNION, DIFFERENCE, XOR
};
enum EdgeType {
	NORMAL, NON_CONTRIBUTING, SAME_TRANSITION, DIFFERENT_TRANSITION
};
enum PolygonType {
	SUBJECT, CLIPPING
};

template<typename TYPE, St DIM>
class SweepEvent_;

template<typename TYPE, St DIM>
class SweepEventComp_;

template<typename TYPE, St DIM>
class SegmentComp_;

template<typename TYPE, St DIM>
class Connector_;

template<typename TYPE>
class Clip_ {
public:
	static const St Dim = 2;
	typedef Clip_<TYPE> Clip;
	typedef Polygon_<TYPE> Polygon;
	typedef SweepEvent_<TYPE, Dim> SweepEvent;
	typedef SweepEventComp_<TYPE, Dim> SweepEventComp;
	typedef Point_<TYPE, Dim> Point;
	typedef Segment_<TYPE, Dim> Segment;
	typedef Operation_<TYPE, Dim> Op;
	typedef Connector_<TYPE, Dim> Connector;
	typedef SegmentComp_<TYPE, Dim> SegmentComp;
public:
	/** @brief Event Queue */
	std::priority_queue<SweepEvent*, std::vector<SweepEvent*>, SweepEventComp> eq;
	/** @brief It holds the events generated during the computation of the boolean operation */
	std::deque<SweepEvent> eventHolder;
	/** @brief Polygon 1 */
	Polygon& subject;
	/** @brief Polygon 2 */
	Polygon& clipping;
	/** To compare events */
	SweepEventComp sec;
	/** @brief Number of intersections (for statistics) */
	int nint;

	/** Class constructor */
	Clip_(Polygon& sp, Polygon& cp) :
			eq(), eventHolder(), subject(sp), clipping(cp), sec(), nint(0) {
	}

	/** Compute the boolean operation */
	void compute(BooleanOpType op, Polygon& result) {
		// Test 1 for trivial result case -----------------
		// At least one of the polygons is empty
		if (1 == _trivial_1(op, result)) {
			return;
		}
		// Test 2 for trivial result case -----------------
		Point minsubj,maxsubj,minclip,maxclip;
		subject.boundingbox(minsubj, maxsubj);
		clipping.boundingbox(minclip, maxclip);
		if (minsubj.x() > maxclip.x() || minclip.x() > maxsubj.x()
				|| minsubj(_Y_) > maxclip(_Y_) || minclip(_Y_) > maxsubj(_Y_)) {
			// the bounding boxes do not overlap
			if (op == DIFFERENCE)
				result = subject;
			if (op == UNION) {
				result = subject;
				for (St i = 0; i < clipping.ncontours(); i++)
					result.push_back(clipping.contour(i));
			}
			return;
		}
		// Boolean operation is not trivial ---------------

		// Insert all the endpoints associated to the line segments into the event queue
		for (St i = 0; i < subject.ncontours(); i++)
			for (St j = 0; j < subject.contour(i).nvertices(); j++)
				this->_process_segment(subject.contour(i).segment(j), SUBJECT);
		for (St i = 0; i < clipping.ncontours(); i++)
			for (St j = 0; j < clipping.contour(i).nvertices(); j++)
				this->_process_segment(clipping.contour(i).segment(j),
						CLIPPING);

		Connector connector; // to connect the edge solutions
		std::set<SweepEvent*, SegmentComp> S; // Status line
		typename std::set<SweepEvent*, SegmentComp>::iterator it,sli,prev,next;
		SweepEvent* e;
		const double MINMAXX = std::min(maxsubj.x(), maxclip.x()); // for optimization 1

		while (!eq.empty()) {
			e = eq.top();
			eq.pop();
#ifdef _DEBUG_
			std::cout << "Process event: ";
			e->show();
#endif
			// optimization 1
			if ((op == INTERSECTION && (e->p.x() > MINMAXX))
					|| (op == DIFFERENCE && e->p.x() > maxsubj.x())) {
				connector.toPolygon(result);
				return;
			}
			if ((op == UNION && (e->p.x() > MINMAXX))) {
				// add all the non-processed line segments to the result
				if (!e->left)
					connector.add(e->segment());
				while (!eq.empty()) {
					e = eq.top();
					eq.pop();
					if (!e->left)
						connector.add(e->segment());
				}
				connector.toPolygon(result);
				return;
			}
			// end of optimization 1

			if (e->left) { // the line segment must be inserted into S
				it = S.insert(e).first;
				e->poss = new typename std::set<SweepEvent*>::iterator(it);
				next = prev = it;
				(prev != S.begin()) ? --prev : prev = S.end();

				// Compute the inside and inOut flags
				if (prev == S.end()) { // there is not a previous line segment in S?
					e->inside = e->inOut = false;
				} else if ((*prev)->type != NORMAL) {
					if (prev == S.begin()) { // e overlaps with prev
						e->inside = true; // it is not relevant to set true or false
						e->inOut = false;
					} else { // the previous two line segments in S are overlapping line segments
						sli = prev;
						sli--;
						if ((*prev)->pl == e->pl) {
							e->inOut = !(*prev)->inOut;
							e->inside = !(*sli)->inOut;
						} else {
							e->inOut = !(*sli)->inOut;
							e->inside = !(*prev)->inOut;
						}
					}
				} else if (e->pl == (*prev)->pl) { // previous line segment in S belongs to the same polygon that "e" belongs to
					e->inside = (*prev)->inside;
					e->inOut = !(*prev)->inOut;
				} else { // previous line segment in S belongs to a different polygon that "e" belongs to
					e->inside = !(*prev)->inOut;
					e->inOut = (*prev)->inside;
				}

#ifdef _DEBUG_
				std::cout << "Status line after insertion: " << std::endl;
				for (typename std::set<SweepEvent*, SegmentComp>::const_iterator it2 =
						S.begin(); it2 != S.end(); it2++)
					(*it2)->show();
#endif

				// Process a possible intersection between "e" and its next neighbor in S
				if ((++next) != S.end())
					_possiblei_intersection(e, *next);

				// Process a possible intersection between "e" and its previous neighbor in S
				if (prev != S.end())
					_possiblei_intersection(*prev, e);
			} else { // the line segment must be removed from S
				next = prev = sli = *(e->other->poss); // S.find (e->other);

				// Get the next and previous line segments to "e" in S
				++next;
				(prev != S.begin()) ? --prev : prev = S.end();

				// Check if the line segment belongs to the Boolean operation
				switch (e->type) {
				case (NORMAL):
					switch (op) {
					case (INTERSECTION):
						if (e->other->inside)
							connector.add(e->segment());
						break;
					case (UNION):
						if (!e->other->inside)
							connector.add(e->segment());
						break;
					case (DIFFERENCE):
						if (((e->pl == SUBJECT) && (!e->other->inside))
								|| (e->pl == CLIPPING && e->other->inside))
							connector.add(e->segment());
						break;
					case (XOR):
						connector.add(e->segment());
						break;
					}
					break;
				case (SAME_TRANSITION):
					if (op == INTERSECTION || op == UNION)
						connector.add(e->segment());
					break;
				case (DIFFERENT_TRANSITION):
					if (op == DIFFERENCE)
						connector.add(e->segment());
					break;
				}
				// delete line segment associated to e from S and check for intersection between the neighbors of "e" in S
				S.erase(sli);
				if (next != S.end() && prev != S.end())
					_possiblei_intersection(*prev, *next);
			}

		} // end while

	}
	/** Number of intersections found (for statistics) */
	int nInt() const {
		return nint;
	}

protected:
	int _trivial_1(BooleanOpType op, Polygon& result) {
		// Test 1 for trivial result case
		// At least one of the polygons is empty
		if (subject.ncontours() * clipping.ncontours() == 0) {
			if (op == DIFFERENCE)
				result = subject;
			if (op == UNION)
				result = (subject.ncontours() == 0) ? clipping : subject;
			return 1;
		} else {
			return 0;
		}
	}

	/** @brief Compute the events associated to segment s, and insert them into pq and eq */
	void _process_segment(const Segment& s, PolygonType pl) {
		// if the two edge endpoints are equal, the segment is dicarded
		if (s.ps() == s.pe()) {
			// in the future this can be done as preprocessing to avoid "polygons" with less than 3 edges
			return;
		}
		SweepEvent* e1 = this->_store(SweepEvent(s.ps(), true, pl, 0));
		SweepEvent* e2 = this->_store(SweepEvent(s.pe(), true, pl, e1));
		e1->other = e2;

		if (e1->p.x() < e2->p.x()) {
			e2->left = false;
		} else if (e1->p.x() > e2->p.x()) {
			e1->left = false;
		} else if (e1->p.y() < e2->p.y()) { // the line segment is vertical. The bottom endpoint is the left endpoint
			e2->left = false;
		} else {
			e1->left = false;
		}
		eq.push(e1);
		eq.push(e2);
	}
	/** @brief Process a posible intersection between the segment associated to the left events e1 and e2 */
	void _possiblei_intersection(SweepEvent *e1, SweepEvent *e2) {
		// you can uncomment these two lines if self-intersecting polygons are not allowed
		//	if ((e1->pl == e2->pl) )
		//		return false;

		Point ip1,ip2;  // intersection points
		int nintersections;

		if (!(nintersections = Op::FindIntersection(e1->segment(),
				e2->segment(), ip1, ip2)))
			return;

		if ((nintersections == 1)
				&& ((e1->p == e2->p) || (e1->other->p == e2->other->p)))
			return; // the line segments intersect at an endpoint of both line segments

		if (nintersections == 2 && e1->pl == e2->pl)
			return; // the line segments overlap, but they belong to the same polygon

		// The line segments associated to e1 and e2 intersect
		nint += nintersections;

		if (nintersections == 1) {
			// if ip1 is not an endpoint of the line segment associated to e1 then divide "e1"
			if (e1->p != ip1 && e1->other->p != ip1)
				_divide_segment(e1, ip1);
			// if ip1 is not an endpoint of the line segment associated to e2 then divide "e2"
			if (e2->p != ip1 && e2->other->p != ip1)
				_divide_segment(e2, ip1);
			return;
		}

		// The line segments overlap
		std::vector<SweepEvent *> sortedEvents;
		if (e1->p == e2->p) {
			sortedEvents.push_back(0);
		} else if (sec(e1, e2)) {
			sortedEvents.push_back(e2);
			sortedEvents.push_back(e1);
		} else {
			sortedEvents.push_back(e1);
			sortedEvents.push_back(e2);
		}
		if (e1->other->p == e2->other->p) {
			sortedEvents.push_back(0);
		} else if (sec(e1->other, e2->other)) {
			sortedEvents.push_back(e2->other);
			sortedEvents.push_back(e1->other);
		} else {
			sortedEvents.push_back(e1->other);
			sortedEvents.push_back(e2->other);
		}

		if (sortedEvents.size() == 2) { // are both line segments equal?
			e1->type = e1->other->type = NON_CONTRIBUTING;
			e2->type = e2->other->type =
					(e1->inOut == e2->inOut) ?
							SAME_TRANSITION : DIFFERENT_TRANSITION;
			return;
		}
		if (sortedEvents.size() == 3) { // the line segments share an endpoint
			sortedEvents[1]->type = sortedEvents[1]->other->type =
					NON_CONTRIBUTING;
			if (sortedEvents[0])      // is the right endpoint the shared point?
				sortedEvents[0]->other->type =
						(e1->inOut == e2->inOut) ?
								SAME_TRANSITION : DIFFERENT_TRANSITION;
			else
				// the shared point is the left endpoint
				sortedEvents[2]->other->type =
						(e1->inOut == e2->inOut) ?
								SAME_TRANSITION : DIFFERENT_TRANSITION;
			_divide_segment(
					sortedEvents[0] ? sortedEvents[0] : sortedEvents[2]->other,
					sortedEvents[1]->p);
			return;
		}
		if (sortedEvents[0] != sortedEvents[3]->other) { // no line segment includes totally the other one
			sortedEvents[1]->type = NON_CONTRIBUTING;
			sortedEvents[2]->type =
					(e1->inOut == e2->inOut) ?
							SAME_TRANSITION : DIFFERENT_TRANSITION;
			_divide_segment(sortedEvents[0], sortedEvents[1]->p);
			_divide_segment(sortedEvents[1], sortedEvents[2]->p);
			return;
		}
		// one line segment includes the other one
		sortedEvents[1]->type = sortedEvents[1]->other->type = NON_CONTRIBUTING;
		_divide_segment(sortedEvents[0], sortedEvents[1]->p);
		sortedEvents[3]->other->type =
				(e1->inOut == e2->inOut) ?
						SAME_TRANSITION : DIFFERENT_TRANSITION;
		_divide_segment(sortedEvents[3]->other, sortedEvents[2]->p);
	}
	/** @brief Divide the segment associated to left event e, updating pq and (implicitly) the status line */
	void _divide_segment(SweepEvent *e, const Point& p) {
		// "Right event" of the "left line segment" resulting from dividing e
		// (the line segment associated to e)
		SweepEvent *r = _store(SweepEvent(p, false, e->pl, e, e->type));
		// "Left event" of the "right line segment" resulting from dividing e
		// (the line segment associated to e)
		SweepEvent *l = _store(
				SweepEvent(p, true, e->pl, e->other, e->other->type));
		if (sec(l, e->other)) { // avoid a rounding error. The left event would be processed after the right event
			std::cout << "Oops" << std::endl;
			e->other->left = true;
			l->left = false;
		}
		if (sec(e, r)) { // avoid a rounding error. The left event would be processed after the right event
			std::cout << "Oops2" << std::endl;
			//		cout << *e << endl;
		}
		e->other->other = l;
		e->other = r;
		eq.push(l);
		eq.push(r);
	}
	/** @brief Store the SweepEvent e into the event holder, returning the address of e */
	SweepEvent *_store(const SweepEvent& e) {
		eventHolder.push_back(e);
		return &eventHolder.back();
	}

}
;

}

#endif /* _POLYGON_BOOLEAN_HPP_ */
