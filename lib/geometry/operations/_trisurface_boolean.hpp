/*
 * _polygon_boolean.hpp
 *
 *  Created on: Jun 19, 2017
 *      Author: zhou
 */

#ifndef _TRISURFACE_BOOLEAN_HPP_
#define _TRISURFACE_BOOLEAN_HPP_

#include "geometry/geometry_define.hpp"
#include "geometry/objects/_objects.hpp"
#include "algebra/array_list.hpp"
#include "_operation.hpp"
#include "_intersection.hpp"
#include "_sweep.hpp"

//#include "io/plotly.h"
//#include "geometry/io/_actor_plotly.hpp"

#include <array>
#include <vector>
#include <limits>
#include <list>
#include <set>
#include <fstream>
#include <queue>

//#define _DEBUG_

namespace carpio {

template<typename TYPE>
class TriSurface_Boolean_ {
public:
	static const St Dim = 3;
	typedef Polygon_<TYPE> Polygon;
	typedef Point_<TYPE, Dim> Point;
	typedef Segment_<TYPE, Dim> Segment;
	typedef Operation_<TYPE, Dim> Op;
	typedef Creation_<TYPE, Dim> Cr;
	typedef Intersection_<TYPE, Dim> Isc;

	typedef TriSurface_<TYPE, Dim> TriSurface;
	typedef TriSurface* pTriSurface;
	typedef const TriSurface* const_pTriSurface;

	typedef typename TriSurface::Fac TriFace;
	typedef typename TriSurface::pFac pTriFace;

	typedef typename TriFace::Tag TagTriFace;
	typedef typename TriSurface::Edg Edge;
	typedef Edge* pEdge;
	typedef typename TriSurface::Ver Vertex;
	typedef typename TriSurface::pVer pVertex;

	typedef Point_<TYPE, 2> Point2;

	typedef Box_<TYPE, Dim> Box;
	typedef BBox_<TYPE, Dim> BBox;
	typedef BBTree_<BBox> BBTree;

	typedef ArrayListV<short> ArrFlag;
	typedef std::list<pEdge> list_pEdge;
	typedef std::list<pVertex> list_pVertex;

protected:
	pTriSurface _subject;
	pTriSurface _clipping;

public:
	TriSurface_Boolean_(TriSurface& sub, TriSurface& clip) {
		_subject = &sub;
		_clipping = &clip;
	}

	/** Compute the boolean operation */
	int compute(BooleanOpType op, TriSurface& result) {
		/// trivial case 1
		if (1 == trivial_1(op, result)) {
			return _SUCCESS;
		}
		/// Crate tree
		BBTree tsub;
		BBTree tclip;
		Cr::BoundingBoxTree(tsub, *_subject);
		Cr::BoundingBoxTree(tclip, *_clipping);

		if (!Isc::Check(tsub, tclip, 1)) {
			/// trivial case 2
			/// not overlap
			if (op == DIFFERENCE)
				TriSurface::Copy(_subject, &result);
			if (op == UNION) {
				TriSurface::Copy(_subject, &result);
				for (auto& pfc : *_clipping) {
					result.insert(pfc);
				}
			}
			return _SUCCESS;
		}

		/// normal case
		std::list<pTriFace> lsub, lclip;
		std::list<std::array<Point, 2> > lsegres;

		_intersection_list(tsub, tclip, lsub, lclip, lsegres);
		std::cout << "list sub     = " << lsub.size() << "\n";
		std::cout << "list clip    = " << lclip.size() << "\n";
		std::cout << "list segres  = " << lsegres.size() << "\n";
		std::cout << "intersection list =========\n";

		std::list<pEdge> ledge;
		_new_edge_list(lsub, lclip, lsegres, ledge);
		std::cout << "new edge list =============\n";

		/// recontruct face
		std::list<pVertex> lvend;
		for (auto& pfsub : lsub) {
//			_reconstruct_subface(pfsub, lsub, lclip, ledge);
			std::list<pTriFace> lonclip;
			std::list<pEdge> lonedge;
			_choose_clip_face_on_subject(pfsub, lsub, lclip, ledge, lonclip, lonedge);
			_find_ends_vertex(lonedge, lvend);
			std::cout<< "ends p = " << lvend.size() <<std::endl;

			break;
		}
//		typedef PlotlyActor_Geometry_<TYPE, 3> PAG;
//		Plotly p;
//		auto actor = PAG::Surface(*_subject);
//		p.add(actor);
//		actor->set_opacity(0.6);
//		auto actor2 = PAG::SurfaceWireFrame(*_clipping);
//		p.add(actor2);
//		actor2->set_opacity(0.8);
//
//		for (auto& ap : lsegres) {
//			auto actorap = PAG::AsSegment(ap[0], ap[1]);
//			p.add(actorap);
//		}
//
//		auto actorpoints = PAG::ScalarPoints(lvend);
//		p.add(actorpoints);
//		auto actor3 = PAG::Surface(result);
//		//p.add(actor3);
//		//actor3->set_opacity(0.5);
//		p.plot();

	}
protected:
	void _intersection_list(
			BBTree& tsub,
			BBTree& tclip,
			std::list<pTriFace>& lsub,
			std::list<pTriFace>& lclip,
			std::list<std::array<Point, 2> >& lsegres
			) {
		for (typename BBTree::iterator iter = tsub.begin();
				iter != tsub.end();
				++iter) {
			if (iter->is_leaf()) {
				BBox& bsub = (iter->box());
				/// define function
				typename BBTree::Fun_flag fun =
						[&bsub, &lsub, &lclip, &lsegres](
								typename BBTree::pNode pn, bool& flag) {
							if (!(Isc::Check(pn->box(), bsub))) {
								flag = false;
							} else {
								flag = true;
							}
							if(pn->is_leaf() && Isc::Check(pn->box(), bsub)) {
								std::tuple<bool, bool, Point, Point> res;
								BBox& bclip = pn->box();
								bool isfind = Isc::Find(bsub, bclip, res, 0.0);
								if(isfind) {
									if(!(std::get<1>(res))) { //not coplane
										pTriFace pfclip = any_cast<pTriFace>(bclip.get_obj());
										pTriFace pfsub = any_cast<pTriFace>(bsub.get_obj());
										Point& p1 = std::get<2>(res);
										Point& p2 = std::get<3>(res);
										if(p1 != p2) {
											std::array<Point, 2> apoint { {p1 , p2}};
											lsegres.push_back(apoint);
											lsub.push_back(pfsub);
											lclip.push_back(pfclip);
										}
									}
								}
							}
						};

				tclip.PreOrder_flag(fun);
			}
		}
	}

	void _new_edge_list(
			std::list<pTriFace>& lsub,
			std::list<pTriFace>& lclip,
			std::list<std::array<Point, 2> >& lsegres,
			std::list<pEdge>& ledg) {
		ledg.clear();
		std::set<pVertex> used;
		auto iter_sub = lsub.begin();
		auto iter_clip = lclip.begin();
		auto iter_arr = lsegres.begin();
		for (; iter_sub != lsub.end();) {
			auto fsub = *iter_sub;
			auto fclip = *iter_clip;
			auto arr = *iter_arr;

			// add vertex on face
			for (int i = 0; i < 3; ++i) {
				used.insert(fsub->vertex(i));
				used.insert(fclip->vertex(i));
			}

			Point& p1 = arr[0];
			Point& p2 = arr[1];
			pVertex v1 = nullptr, v2 = nullptr;
			for (auto& pvold : used) {
//				std::cout<< "p1 = " << p1 <<" " << *pvold << std::endl;
//				std::cout<< "p2 = " << p2 <<" " << *pvold << std::endl;
				if (Op::Distance(p1, *(pvold)) < 1e-10) { // found
					v1 = pvold;
				}
				if (Op::Distance(p2, *(pvold)) < 1e-10) { // found
					v2 = pvold;
				}
			}
			if (v1 == nullptr) {
				v1 = new Vertex(p1);
				used.insert(v1);
			}
			if (v2 == nullptr) {
				v2 = new Vertex(p2);
				used.insert(v2);
			}
			pEdge pe = new Edge(v1, v2);
			/// direction
			auto nsub = fsub->normal();
			auto nclip = fclip->normal();
			auto te = pe->tangent();
			double svolume = Op::TripleScalar(nsub, nclip, te);
			if (svolume < 0) {
				pe->reverse();
			}
			ledg.push_back(pe);
			++iter_sub;
			++iter_clip;
			++iter_arr;
		}
	}

	struct IntersectUnit {
		pTriFace face_subject;
		pTriFace face_clipping;
		pEdge edge;
		int group;

		IntersectUnit(pTriFace ps, pTriFace pc, pEdge pe) {
			face_subject = ps;
			face_clipping = pc;
			edge = pe;
			group = 0;
		}

		struct CompareFaceClipping {
			bool operator()(
					const IntersectUnit& a,
					const IntersectUnit& b) const {
				return a.face_clipping < b.face_clipping;
			}
		};

		struct CompareEdge {
			bool operator()(
					const IntersectUnit& a,
					const IntersectUnit& b) const {
				return a.edge < b.edge;
			}
		};
	};

	int _choose_clip_face_on_subject(
			pTriFace pfsub,
			const std::list<pTriFace>& lsub,
			const std::list<pTriFace>& lclip,
			const std::list<pEdge>& ledge,
			std::list<pTriFace>& lonclip,
			std::list<pEdge>& lonedge) {
		auto iter_sub = lsub.begin();
		auto iter_clip = lclip.begin();
		auto iter_edg = ledge.begin();
		/// get edges on subject face
		for (; iter_sub != lsub.end();) {
			pTriFace pfs = *iter_sub;
			pTriFace pfc = *iter_clip;
			pEdge edg = *iter_edg;
			if (pfs == pfsub) {
				// on subject face
				lonclip.push_back(pfc);
				lonedge.push_back(edg);
			}
			++iter_sub;
			++iter_clip;
			++iter_edg;
		}
		ASSERT(lonclip.size() == lonedge.size());
		return _SUCCESS;
	}

	int _find_ends_vertex(
			std::list<pEdge>& ledge,
			std::list<pVertex>& lvend) {
		std::set<pEdge> smain;
		for (auto& e : ledge) {
			smain.insert(e);
		}
		// find two end vertex
		std::function<void(pVertex&)> fun_find_end_vertex =
				[&smain, &lvend](pVertex pv) {
					int count = 0;
					for(auto it = pv->begin_edge(); it!= pv->end_edge();++it) {
						pEdge e = *it;
						if(smain.find(e) != smain.end()) {
							++count;
						}
					}
					if(count == 1) {
						lvend.push_back(pv);
					}
					ASSERT(count != 0 && count <=2);
				};

		Edge::ForEachVertex(smain, fun_find_end_vertex);

	}

	int _reconstruct_subface(
			pTriFace pfsub,
			std::list<pTriFace>& lsub,
			std::list<pTriFace>& lclip,
			std::list<pEdge>& ledge) {
		auto iter_sub = lsub.begin();
		auto iter_clip = lclip.begin();
		auto iter_edg = ledge.begin();
		/// get edges on subject face
		std::list<IntersectUnit> l_on;
		for (; iter_sub != lsub.end();) {
			pTriFace pfs = *iter_sub;
			pTriFace pfc = *iter_clip;
			pEdge edg = *iter_edg;
			if (pfs == pfsub) {
				// on subject face
				l_on.push_back(IntersectUnit(pfs, pfc, edg));
			}
			++iter_sub;
			++iter_clip;
			++iter_edg;
		}
		/// edge on the subject face should be connected as a point chain
		/// this step will connect Edges as vertex list
		/// 1 regroup edges
		int groups = _regroup(l_on);
		std::cout << "Groups = " << groups << "\n";

		//std::cout << "group size = " << group.size() << std::endl;
		for (int i = 0; i < groups; i++) {
			std::list<pVertex> lver;
			//_to_vertex_list(l_on, lver, i);
			// gver[i] = lver;
			std::list<pEdge> ledge;
			_sort_edges(l_on, i, lver);
		}

		///

//		_reconstruct_subface_one(pfsub, lclip_on, ledge_on);

	}

	void _projection_on_face_plane(
			pTriFace pfsub,
			std::list<pVertex> lver,
			std::list<IntersectUnit>& l_on,
			pVertex theone, int idxgroup = 0) {
		/// find other two vertices

		/// the first three points are the vertex of the triangle

	}

	int _regroup(
			std::list<IntersectUnit>& l_on) {
		/// regroup by edge connect relation
		typedef std::map<pEdge, IntersectUnit*> Map;
		typedef std::pair<pEdge, IntersectUnit*> Pair;
		Map smain;
		Map ssec;
		for (auto& unit : l_on) {
			smain.insert(Pair(unit.edge, &unit));
			ssec.insert(Pair(unit.edge, &unit));
		}
		int flag = 0;
		while (smain.size() > 0) {
			std::cout << "flag = " << flag << std::endl;
			int count = 0;
			for (auto& unit : smain) {
				if (count == 0) {
					unit.second->group = flag;
					ssec.erase(unit.first);
				}
				for (auto& unit_other : ssec) {
					if (Edge::IsConnected(*(unit.first), *(unit_other.first))) {
						unit.second->group = flag;
						ssec.erase(unit_other.first);
						break;
					}
				}
				count++;
			}
			smain = ssec;
			flag++;
		}
		return flag;  // return the size of groups

	}

	void _to_vertex_list(
			std::list<IntersectUnit>& l_on,
			std::list<pVertex>& lver,
			int idxgroup = 0) {
		typedef std::map<pEdge, IntersectUnit*> Map;
		typedef std::pair<pEdge, IntersectUnit*> Pair;
		Map smain;
		for (auto& unit : l_on) {
			if (unit.group == idxgroup) {
				smain.insert(Pair(unit.edge, &unit));
			}
		}
		lver.clear();
		pEdge pbegin = smain.begin()->first;
		lver.push_back(pbegin->vertex(0));

		while (smain.size() > 0) {
			pVertex v = lver.back();
			std::cout << "ver " << *v << "\n";
			short found = 0;
			for (auto iter = v->begin_edge(); iter != v->end_edge(); ++iter) {
				pEdge edge = *(iter);
				if (smain.find(edge) != smain.end()) {
					//found
					found = 1;
					pVertex vo = v->other_vertex_on_edge(edge);
					smain.erase(edge);
					lver.push_back(vo);
					break;
				}
			}
			ASSERT(found == 1);
		}
	}

	int _choose_the_one_vertex(
			pTriFace pfsub,
			std::list<IntersectUnit>& l_on,
			pVertex theone, int idxgroup = 0) {
		/// find orientation
		ArrFlag orif(3, 0);
		for (int i = 0; i < 3; i++) {
			pVertex v = pfsub->vertex(i);
			for (auto& unit : l_on) {
				if (unit.group == idxgroup) {
					int side = Op::OnWhichSide3(
							*v,
							*(unit.face_clipping->vertex(0)),
							*(unit.face_clipping->vertex(1)),
							*(unit.face_clipping->vertex(2)));
					if (side > 0) {
						orif[i] = 1;
						break;
					}
				}
			}
		}
		/// count 1
		int count1 = orif.count_equal(1);
		if (count1 == 1) {
			for (int i = 0; i < 3; i++) {
				if (orif[i] == 1) {
					theone = pfsub->vertex(i);
					return count1;
				}
			}
		}
		if (count1 == 2) {
			for (int i = 0; i < 3; i++) {
				if (orif[i] == 0) {
					theone = pfsub->vertex(i);
					return count1;
				}
			}
		}
		if (count1 == 3 || count1 == 0) {
			theone = pfsub->vertex(0);
			return 3;
		}
	}

	int _reconstract_contour_on_triface(
			pTriFace pf,
			std::list<pEdge>& ledge,
			std::array<short, 3>& orif, short flag) {
		int count = 0;
		for (auto& f : orif) {
			if (f == flag) {
				count++;
			}
		}
		if (count == 1) {
			/// there are two vertices that connect to edges
			/// which are not belong to ledge

		}
	}

	int trivial_1(BooleanOpType op, TriSurface& result) {
		if (_subject->size_face() * _clipping->size_face() == 0) {
			if (op == DIFFERENCE)
				TriSurface::Copy(_subject, &result);
			if (op == UNION)
				if (_subject->empty()) {
					TriSurface::Copy(_clipping, &result);
				} else {
					TriSurface::Copy(_subject, &result);
				}
			return 1;
		} else {
			return 0;
		}
	}

}
;

}
#endif
