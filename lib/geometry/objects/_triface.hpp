/************************
 //  \file   ts_face.h
 //  \brief
 //
 //  \author czhou
 //  \date   15 juin 2015
 ***********************/
#ifndef _TRIFACE_HPP_
#define _TRIFACE_HPP_

#include "_vertex.hpp"
#include "_edge.hpp"

namespace carpio {

struct TagTriFace: public TagGeometry {
	TagTriFace() {
	}
	;
};

template<class TYPE, St DIM, class SURFACE>
class TriFace_ {
public:
	static const St Dim = DIM;
	typedef TagTriFace Tag;
	typedef TriFace_<TYPE, DIM, SURFACE> Self;
	typedef Self Fac;
	typedef Fac* pFac;
	typedef St size_type;
	typedef TYPE Vt;
	typedef Point_<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef Edge_<TYPE, DIM, Self> Edg;
	typedef Edg* pEdg;
	typedef const Edg* const_pEdg;
	typedef Vertex_<TYPE, DIM, Edg> Ver;
	typedef Ver* pVer;
	typedef const Ver* const_pVer;
	typedef SURFACE Sur;
	typedef Sur* pSur;
	typedef const Sur* const_pSur;

	typedef std::list<pSur> list_pSur;
	typedef std::list<pFac> list_pFac;

	typedef Operation_<TYPE, DIM> Op;

public:
	pEdg e1;
	pEdg e2;
	pEdg e3;

	list_pSur surfaces;

	Any _any_data;

public:
	TriFace_(pEdg a, pEdg b, pEdg c, pSur psur) :
			e1(a), e2(b), e3(c) {
		if (!has_parent_surface(psur)) {
			surfaces.push_back(psur);
		}

		e1->faces.push_back(this);
		e2->faces.push_back(this);
		e3->faces.push_back(this);
	}

	TriFace_(pVer v1, pVer v2, pVer v3, pSur psur) {
		e1 = new Edg(v1, v2);
		e2 = new Edg(v2, v3);
		e3 = new Edg(v3, v1);

		if (!has_parent_surface(psur)) {
			surfaces.push_back(psur);
		}

		e1->faces.push_back(this);
		e2->faces.push_back(this);
		e3->faces.push_back(this);
	}

	void attach(pSur ps) {
		if (!has_parent_surface(ps)) {
			surfaces.push_back(ps);
		}
	}

	void detach(pSur ps) {
		if (has_parent_surface(ps)) {
			surfaces.remove(ps);
		}
	}

	inline bool is_unattached() {
		if (surfaces.empty())
			return true;
		return false;
	}

	/**
	 * gts_face_has_parent_surface:
	 * @f: a #GtsFace
	 * @s: a #GtsSurface.
	 *
	 * Returns: %TRUE if @f belongs to @s, %FALSE otherwise.
	 */
	bool has_parent_surface(const_pSur ps) const {
		_RETURN_VAL_IF_FAIL(ps != nullptr, false);
		for (auto& s : this->surfaces) {
			if (ps == s) {
				return true;
			}
		}
		return false;
	}

	bool has_one_parent_surface(const_pSur ps) const {
		_RETURN_VAL_IF_FAIL(ps != nullptr, false);
		if (this->surfaces.size() != 1) {
			return false;
		}
		/// compare pointer address
		return (*(this->surfaces.begin())) == ps;
	}

	list_pFac get_neighbor_faces(const_pSur ps) const {
		ASSERT(ps != nullptr);
		list_pFac res;
		for (St i = 0; i < 3; i++) {
			const_pEdg e = this->edge(i);
			for (auto iter = e->begin_face(); iter != e->end_face(); iter++) {
				pFac spf = *iter;
				if (spf != this && spf->has_parent_surface(ps)) {
					res.push_back(spf);
				}
			}
		}
		return std::move(res);
	}

	pEdg operator[](const St& idx) {
		ASSERT(idx >= 0 && idx < 3);
		pEdg l[] = { e1, e2, e3 };
		return l[idx];
	}

	const_pEdg operator[](const St& idx) const {
		ASSERT(idx >= 0 && idx < 3);
		pEdg l[] = { e1, e2, e3 };
		return l[idx];
	}

	pEdg edge(const St& idx) {
		return this->operator [](idx);
	}

	const_pEdg edge(const St& idx) const {
		return this->operator [](idx);
	}

	pVer vertex1() const {
		return e1->v1;
	}
	pVer vertex2() const {
		if (e1->v1 == e2->v1) {
			return e2->v2;
			// case 1 ----------
			//     v2 *#
			//       /
			//   e2 /
			//     /
			// v1 /
			//    ----------*
			//   v1   e1   v2

		} else if (e1->v2 == e2->v2) {
			// case 2 ----------
			//         \ v1
			//          \
				//           \ e2
			//            \
				//             \  v2
			//    ----------*#
			//   v1   e1   v2
			return e1->v2;
		} else if (e1->v1 == e2->v2) {
			// case 3 ----------
			//     v1 #
			//       /
			//   e2 /
			//     /
			// v2 /
			//   *----------*
			//   v1   e1   v2
			return e2->v1;
		} else if (e1->v2 == e2->v1) {
			// case 2 ----------
			//         * v1
			//          \
				//           \ e2
			//            \
				//             \  v2
			//    ----------*#
			//   v1   e1   v2
			return e1->v2;
		}
		SHOULD_NOT_REACH;
		return nullptr;
		// return e1->v2;
	}
	pVer vertex3() const {
		if (e1->v1 == e2->v1) {
			return e1->v2;
		} else if (e1->v2 == e2->v2) {
			return e2->v1;
		} else if (e1->v1 == e2->v2) {
			return e1->v2;
		} else if (e1->v2 == e2->v1) {
			return e2->v2;
		}
		SHOULD_NOT_REACH;
		return nullptr;
		//return (e1->v1 == e2->v1) || (e1->v2 == e2->v1) ? e2->v2 : e2->v1;
	}
	pVer vertex(St i) const {
		assert(i < 3);
		if (i == 0) {
			return vertex1();
		}
		if (i == 1) {
			return vertex2();
		}
		if (i == 2) {
			return vertex3();
		}
		return nullptr;
	}
	Vt vertex(St i, int a) const {
		assert(i < 3);
		if (i == 0) {
			return vertex1()->val(a);
		}
		if (i == 1) {
			return vertex2()->val(a);
		}
		if (i == 2) {
			return vertex3()->val(a);
		}
		return 0;
	}

	bool has_vertex(pVer pv) const {
		return pv == e1->v1 || pv == e1->v2
				|| pv == e2->v1 || pv == e2->v2
				|| pv == e3->v1 || pv == e3->v2;
	}

	bool has_edge(pEdg pe) const {
		return pe == e1 || pe == e2 || pe == e3;
	}

	Poi centroid() const {
		return Op::Centroid(
				*(this->vertex(0)),
				*(this->vertex(1)),
				*(this->vertex(2)));
	}

	Poi normal() const {
		return Op::Normal(
				*(this->vertex(0)),
				*(this->vertex(1)),
				*(this->vertex(2)));
	}

	Vt signed_area() const {
		ASSERT(Dim == 2);
		return Op::SignedArea(
				*(this->vertex(0)),
				*(this->vertex(1)),
				*(this->vertex(2)));
	}

	bool is_inon(const Poi& poi) const {
		ASSERT(Dim == 2);
		return Op::IsInOn(
				*(this->vertex(0)),
				*(this->vertex(1)),
				*(this->vertex(2)), poi);
	}

	typename Op::Box box() const {
		return Op::BoundingBox(
				*(this->vertex(0)),
				*(this->vertex(1)),
				*(this->vertex(2)));
	}

	Any& data() {
		return _any_data;
	}

	const Any& data() const {
		return _any_data;
	}

	/**
	 * revert:
	 * @t: a #GtsTriangle.
	 *
	 * Changes the orientation of triangle @t, turning it inside out.
	 */
	void revert() {
		pEdg e;

		e = this->e1;
		this->e1 = this->e2;
		this->e2 = e;
	}

	pVer opposite_vertex(pEdg pe) {
		ASSERT(has_edge(pe));
		for (St i = 0; i < 3; ++i) {
			pVer pv = vertex(i);
			if (!pe->has_vertex(pv)) {
				return pv;
			}
		}
		return nullptr;
	}

	pEdg opposite_edge(pVer pv) {
		ASSERT(has_vertex(pv));
		for (St i = 0; i < 3; ++i) {
			pEdg pe = edge(i);
			if (!(pe->has_vertex(pv))) {
				return pe;
			}
		}
		return nullptr;
	}

	void show() const {
		std::ios::fmtflags f(std::cout.flags());
		std::cout.setf(std::ios::right);
		this->e1->show();
		this->e2->show();
		this->e3->show();
		std::cout << "    -> sur = " << surfaces.size() << "\n";
		std::cout.setf(f);
	}

	/**
	 * triangles_common_edge:
	 * @t1: a #GtsTriangle.
	 * @t2: a #GtsTriangle.
	 *
	 * Returns: a #GtsEdge common to both @t1 and @t2 or %NULL if @t1 and @t2
	 * do not share any edge.
	 */
	static pEdg GetCommon_pEdg( //
			pFac t1,  //
			pFac t2) {
		_RETURN_VAL_IF_FAIL(t1 != nullptr, nullptr);
		_RETURN_VAL_IF_FAIL(t2 != nullptr, nullptr);

		if (t1->e1 == t2->e1 || t1->e1 == t2->e2 || t1->e1 == t2->e3)
			return t1->e1;
		if (t1->e2 == t2->e1 || t1->e2 == t2->e2 || t1->e2 == t2->e3)
			return t1->e2;
		if (t1->e3 == t2->e1 || t1->e3 == t2->e2 || t1->e3 == t2->e3)
			return t1->e3;
		return pEdg(nullptr);
	}

	/**
	 * gts_triangles_are_compatible:
	 * @t1: a #GtsTriangle.
	 * @t2: a #GtsTriangle.
	 * @e: a #GtsEdge used by both @t1 and @t2.
	 *
	 * Checks if @t1 and @t2 have compatible orientations i.e. if @t1 and
	 * @t2 can be part of the same surface without conflict in the surface
	 * normal orientation.
	 *
	 * Returns: %TRUE if @t1 and @t2 are compatible, %FALSE otherwise.
	 */
	static bool AreCompatible( //
			pFac t1,  //
			pFac t2,  //
			pEdg e) {
		pEdg e1 = nullptr, e2 = nullptr;

		_RETURN_VAL_IF_FAIL(t1 != nullptr, false);
		_RETURN_VAL_IF_FAIL(t2 != nullptr, false);
		_RETURN_VAL_IF_FAIL(e != nullptr, false);

		if (t1->e1 == e)
			e1 = t1->e2;
		else if (t1->e2 == e)
			e1 = t1->e3;
		else if (t1->e3 == e)
			e1 = t1->e1;
		else
			assert(false);
		if (t2->e1 == e)
			e2 = t2->e2;
		else if (t2->e2 == e)
			e2 = t2->e3;
		else if (t2->e3 == e)
			e2 = t2->e1;
		else
			assert(false);
		if (e1->v1 == e2->v1 || e1->v1 == e2->v2 || e1->v2 == e2->v1
				|| e1->v2 == e2->v2)
			return false;
		return true;
	}

};

}

#endif /* TS_FACE_H_ */

