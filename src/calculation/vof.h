#ifndef _VOF_H_
#define _VOF_H_

#include "../carpio_define.hpp"
#include "domain/domain.hpp"
#include "geometry/geometry.hpp"
#include "algebra/matrix.hpp"
#include "algebra/arithmetic.hpp"

#include <functional>
#include <math.h>

namespace carpio {
/*
 * The vof in the node
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class Vof_face_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Vof_face_<COO_VALUE, VALUE, DIM> Self;
	typedef Vof_face_<COO_VALUE, VALUE, DIM>& ref_Self;

	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pGrid;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;

	typedef Line_<vt> Line;
	typedef Line_<vt>* pLine;
	typedef const Line_<vt>* const_pLine;

	typedef Plane_<vt> Plane;
	typedef Plane_<vt>* pPlane;
	typedef const Plane_<vt>* const_pPlane;

	typedef Segment_<cvt, 2> Segment;
	typedef Segment_<cvt, 2>* pSegment;
	typedef const Segment_<cvt, 2>* const_pSegment;

	typedef Point_<cvt, 2> Point2D;
	typedef Point_<cvt, 2>* pPoint2D;

protected:
	/*
	 * Data
	 */
	pLine _pe2;  //equation  2D
	pPlane _pe3; //          3D
	//
	pSegment _ps2; //the segment 2D
	//       _ps3;               3D
public:
	Vof_face_() {
		_set_all_null();
	}
	Vof_face_(vt alpha, vt A, vt B, vt C = 0) {
		if (Dim == 2) {
			_pe2 = new Line(A, B, alpha);
			_pe3 = nullptr;
			_ps2 = nullptr;
		} else {
			_pe2 = nullptr;
			_pe3 = new Plane(A, B, C, alpha);
			_ps2 = nullptr;
		}
	}
	~Vof_face_() {
		if (_pe2 != nullptr) {
			delete _pe2;
			_pe2 = nullptr;
		}
		if (_pe3 != nullptr) {
			delete _pe3;
			_pe3 = nullptr;
		}
		if (_ps2 != nullptr) {
			delete _ps2;
			_ps2 = nullptr;
		}
	}
	bool has_segment() const {
		if (Dim == 2) {
			if (_ps2 == nullptr) {
				return false;
			} else {
				if (_ps2->empty()) {
					return false;
				} else {
					return true;
				}
			}
		} else {
			return false; //unfinish
		}
	}
	void set(const Line& line) {
		if (_pe2 != nullptr) {
			this->_pe2->reconstruct(line.a(), line.b(), line.alpha());
		} else {
			this->_pe2 = new Line(line);
		}
	}
	void set(const Segment& seg) {
		if (_ps2 != nullptr) {
			this->_ps2->reconstruct(seg);
		} else {
			this->_ps2 = new Segment(seg);
		}
	}
	bool has_equ() const {
		if (Dim == 2) {
			if (_pe2 != nullptr) {
				return true;
			} else {
				return false;
			}
		} else {
			if (_pe3 != nullptr) {
				return true;
			} else {
				return false;
			}
		}
	}
	bool has_seg() const {
		if (Dim == 2) {
			if (_ps2 != nullptr) {
				return true;
			} else {
				return false;
			}
		} else {
			if (_ps2 != nullptr) { //!!!!!!!
				return true;
			} else {
				return false;
			}
		}
	}
	pLine get_pLine() {
		if (Dim == 2) {
			return this->_pe2;
		} else {
			return nullptr;
		}
	}
	const_pLine get_pLine() const {
		if (Dim == 2) {
			return this->_pe2;
		} else {
			return nullptr;
		}
	}
	pSegment get_pSeg() {
		if (Dim == 2) {
			return this->_ps2;
		} else {
			return nullptr;
		}
	}
	const_pSegment get_pSeg() const {
		if (Dim == 2) {
			return this->_ps2;
		} else {
			return nullptr;
		}
	}

	bool empty() {
		if (_pe2 == nullptr && _pe3 == nullptr) {
			return true;
		} else {
			if (Dim == 2) {
				return (_pe2 == nullptr) ? true : false;
			}
		}
		SHOULD_NOT_REACH;
		return false;
	}
	/*
	 * calculate face
	 */

protected:
	void _set_all_null() {
		_pe2 = nullptr;
		_pe3 = nullptr;
		_ps2 = nullptr;
	}

};
/*
 * the vof class
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class Vof_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const St NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Vof_<COO_VALUE, VALUE, DIM> Self;
	typedef Vof_<COO_VALUE, VALUE, DIM>& ref_Self;

	typedef Grid_<COO_VALUE, VALUE, DIM> Grid;
	typedef Grid_<COO_VALUE, VALUE, DIM> *pGrid;
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pGrid;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	typedef const Node_<COO_VALUE, VALUE, DIM> *const_pNode;
	typedef Vof_face_<COO_VALUE, VALUE, DIM> Voff;
	typedef Vof_face_<COO_VALUE, VALUE, DIM>* pVoff;
	typedef const Vof_face_<COO_VALUE, VALUE, DIM>* const_pVoff;
	typedef Vof_face_<COO_VALUE, VALUE, DIM>& ref_Voff;

	typedef Stencil_<COO_VALUE, VALUE, DIM, DIM> Stencil;
	typedef MatrixS<vt, 3, 3> Matrix3_3;

	typedef Line_<vt> Line;
	typedef Line_<vt>* pLine;
	typedef const Line_<vt>* const_pLine;

	typedef Segment_<cvt, 2> Segment;
	typedef Segment_<cvt, 2>* pSegment;

	typedef Point_<cvt, 2> Point2D;
	typedef Point_<cvt, 2>* pPoint2D;
protected:
	/*
	 * Data
	 */
	pGrid _pg;

	St _c_idx;  //idx c on data
	St _f_idx;  //idx vof_face on utPoint

public:
	/*
	 * constructor
	 */
	Vof_(pGrid pg, St ci, St fi) :
			_pg(pg), _c_idx(ci), _f_idx(fi) {
		_new_faces();
	}
	~Vof_() {
		_delete_faces();
	}
	/*
	 * get
	 */
	St& c_idx() {
		return _c_idx;
	}
	const St& c_idx() const {
		return _c_idx;
	}
	St& f_idx() {
		return _f_idx;
	}
	const St& f_idx() const {
		return _f_idx;
	}
	pGrid get_pgrid() {
		return _pg;
	}
	const_pGrid get_pgrid() const {
		return _pg;
	}
	pVoff get_pface(pNode pn) {
		ASSERT(pn != nullptr);
		ASSERT(pn->data != nullptr);
		return CAST(pVoff, pn->data->utp(_f_idx));
	}
	const_pVoff get_pface(const_pNode pn) const {
		ASSERT(pn != nullptr);
		ASSERT(pn->data != nullptr);
		return CAST(const_pVoff, pn->data->utp(_f_idx));
	}
	vt get_c(pNode pn) {
		ASSERT(pn != nullptr);
		return pn->cdva(_c_idx);
	}

	/*
	 * set initial color function
	 */
	void set_color(const Shape_<cvt, Dim>& shape) {
		for (typename Grid::iterator_leaf iter = _pg->begin_leaf();
				iter != _pg->end_leaf(); ++iter) {
			pNode pn = iter.get_pointer();
			if (pn != nullptr) {
				Shape2D sn, res;
				CreatCube(sn, pn->p(_M_, _X_), pn->p(_M_, _Y_), pn->p(_P_, _X_),
						pn->p(_P_, _Y_));
				Intersect(sn, shape, res);
				vt rv = res.volume();
				vt sv = sn.volume();
				if (res.empty()) { //node is all out
					pn->cd(this->_c_idx) = 0.0;
				} else if (Abs(rv - sv) < 1e-8) { // all in
					pn->cd(this->_c_idx) = 1.0;
				} else {  //intersect
					pn->cd(this->_c_idx) = rv / sv;
				}
			}
		}
	}
	/*
	 *  calculate a and b
	 *
	 *  stencil to matrix
	 */
	void get_matrix(Matrix3_3& mat, Stencil& sten) {
		ASSERT(Dim == 2);
		// the stencil is the same as matrix, which should be 3x3
		typename Stencil::const_pNode pnc = sten.center_pnode();
		ASSERT(pnc != nullptr);
		// set broundry condition (symmetric)
		mat(1, 1) = pnc->cdva(_c_idx);
		// x m
		typename Stencil::pNode pnt = nullptr;
		pnt = sten.get_pnode(-1, 0);
		if (pnt == nullptr) {
			mat(0, 1) = mat(1, 1);
		} else {
			mat(0, 1) = pnt->cdva(_c_idx);
		}
		// x p
		pnt = nullptr;
		pnt = sten.get_pnode(1, 0);
		if (pnt == nullptr) {
			mat(2, 1) = mat(1, 1);
		} else {
			mat(2, 1) = pnt->cdva(_c_idx);
		}
		// y m
		pnt = nullptr;
		pnt = sten.get_pnode(0, -1);
		if (pnt == nullptr) {
			mat(1, 0) = mat(1, 1);
		} else {
			mat(1, 0) = pnt->cdva(_c_idx);
		}
		// y p
		pnt = nullptr;
		pnt = sten.get_pnode(0, 1);
		if (pnt == nullptr) {
			mat(1, 2) = mat(1, 1);
		} else {
			mat(1, 2) = pnt->cdva(_c_idx);
		}
		//corner mm xy
		pnt = nullptr;
		pnt = sten.get_pnode(-1, -1);
		if (pnt == nullptr) {
			mat(0, 0) = (mat(1, 0) + mat(0, 1)) * 0.5;
		} else {
			mat(0, 0) = pnt->cdva(_c_idx);
		}
		//corner pm xy
		pnt = nullptr;
		pnt = sten.get_pnode(1, -1);
		if (pnt == nullptr) {
			mat(2, 0) = (mat(1, 0) + mat(2, 1)) * 0.5;
		} else {
			mat(2, 0) = pnt->cdva(_c_idx);
		}
		//corner pp xy
		pnt = nullptr;
		pnt = sten.get_pnode(1, 1);
		if (pnt == nullptr) {
			mat(2, 2) = (mat(1, 2) + mat(2, 1)) * 0.5;
		} else {
			mat(2, 2) = pnt->cdva(_c_idx);
		}
		//corner mp xy
		pnt = nullptr;
		pnt = sten.get_pnode(-1, 1);
		if (pnt == nullptr) {
			mat(0, 2) = (mat(1, 2) + mat(0, 1)) * 0.5;
		} else {
			mat(0, 2) = pnt->cdva(_c_idx);
		}
	}

	void cal_norm_Young(const Matrix3_3 &C, vt& mx, vt& my) {
		vt cpm, cps, cpp;
		vt csm, csp;
		vt cmm, cms, cmp;
		cpm = C(2, 0);
		cps = C(2, 1);
		cpp = C(2, 2);
		csm = C(1, 0);
		csp = C(1, 2);
		cmm = C(0, 0);
		cms = C(0, 1);
		cmp = C(0, 2);
		mx = cpp + 2.0 * cps + cpm - cmp - 2.0 * cms - cmm;
		my = cpp + 2.0 * csp + cmp - cpm - 2.0 * csm - cmm;
	}

	/*
	 *  construct face
	 *  C get line equation
	 */

	bool _has_face(pNode pn) {
		if (pn->data->utp(_f_idx) == nullptr) {
			return false;
		} else {
			return true;
		}
	}
	void _new_face(pNode pn) {
		ASSERT(pn != nullptr);
		if (!_has_face(pn)) {
			pn->data->utp(_f_idx) = new Voff();
		}
	}
	void _new_face(pNode pn, Line line) {
		ASSERT(pn != nullptr);
		if (!_has_face(pn)) {
			pn->data->utp(_f_idx) = new Voff(line.alpha(), line.a(), line.b());
		}
	}
	void _delete_face(pNode pn) {
		ASSERT(pn != nullptr);
		if (_has_face(pn)) {
			pVoff pvf = CAST(pVoff, pn->data->utp(_f_idx));
			delete pvf;
			pn->data->utp(_f_idx) = nullptr;
		}
	}
	void _new_faces() {
		for (typename Grid::iterator_leaf iter = _pg->begin_leaf();
				iter != _pg->end_leaf(); ++iter) {
			_new_face(iter.get_pointer());
		}
	}
	void _delete_faces() {
		for (typename Grid::iterator_leaf iter = _pg->begin_leaf();
				iter != _pg->end_leaf(); ++iter) {
			_delete_face(iter.get_pointer());
		}
	}
	void _clear_face(pNode pn) {
		ASSERT(pn != nullptr);
		if (_has_face(pn)) {
			pVoff pvf = CAST(pVoff, pn->data->utp(_f_idx));
			delete pvf;
			pn->data->utp(_f_idx) = new Voff();
		}
	}
	void _clear_faces() {
		for (typename Grid::iterator_leaf iter = _pg->begin_leaf();
				iter != _pg->end_leaf(); ++iter) {
			_clear_face(iter.get_pointer());
		}
	}
	void _construct_face_2d(pNode pn) {
		ASSERT(pn != nullptr);
		ASSERT(pn->data != nullptr);
		// 1 find stencil
		Stencil sten(pn, _X_, 1, 1, _Y_, 1, 1);
		// 2 calculate normal vector (a, b)
		Matrix3_3 mat;
		this->get_matrix(mat, sten);
		vt x, y;
		this->cal_norm_Young(mat, x, y);
		// 3 get eqution, save to vof_face
		Line line = _cal_line(-x, -y, pn->cdva(_c_idx));
		line.show();
		pVoff pvff = this->get_pface(pn);
		if (pvff != nullptr) {
			pvff->set(line);
		} else {
			_new_face(pn, line);
		}
	}
	void construct_faces() {
		for (typename Grid::iterator_leaf iter = _pg->begin_leaf();
				iter != _pg->end_leaf(); ++iter) {
			pNode pn = iter.get_pointer();
			if (IsInRange(0.0, get_c(pn), 1.0, _oo_)) {
				_construct_face_2d(pn);
			}
		}
	}
	void _construct_seg_2d(pNode pn) {
		// make sure pn has line
		ASSERT(pn != nullptr);
		ASSERT(pn->data != nullptr);
		pVoff pvff = this->get_pface(pn);
		ASSERT(pvff != nullptr);
		if (!pvff->has_equ()) {
			_construct_face_2d(pn);
		}
		pLine pl = pvff->get_pLine();
		ASSERT(pn != nullptr);
		Segment s = this->_cal_interct_points((*pl), 1.0, 1.0);
		s.scale(pn->d(_X_), pn->d(_Y_));
		s.transfer(pn->p(_M_, _X_), pn->p(_M_, _Y_));
		pvff->set(s);
	}

	void construct_segments() {
		for (typename Grid::iterator_leaf iter = _pg->begin_leaf();
				iter != _pg->end_leaf(); ++iter) {
			pNode pn = iter.get_pointer();
			if (IsInRange(0.0, get_c(pn), 1.0, _oo_)) {
				_construct_seg_2d(pn);
			}
		}
	}

protected:
	/**
	 * \brief   known a,b in ax+by=alpha and C, calculate alpha \n
	 *          no matter what a and b are, they will be change to abs(a) and abs(b)
	 *          return alpha, ax+by=alpha, a>b>0;
	 * \param   Float a a in ax+by=alpha
	 * \param   Float b b in ax+by=alpha
	 * \param   Float C the color function
	 * \return  alpha
	 */
	vt _cal_alpha(vt a, vt b, vt C) {
		vt c1, c2, alpha;
		vt absa = (a < 0) ? (-a) : a;
		vt absb = (b < 0) ? (-b) : b;
		vt m, n;
		n = (absa >= absb) ? (absa) : absb;
		m = (absa <= absb) ? (absa) : absb;
		c1 = m / 2 / n;
		c2 = 1 - c1;
		if (C >= 0 && C <= c1) {
			alpha = sqrt(2 * C * m * n);
		} else if (C > c1 && C < c2) {
			alpha = (2 * C * n + m) / 2;
		} else { //(C>=c2 && C<=1)
			alpha = m + n - sqrt(2 * (1 - C) * m * n);
		}
		return alpha;
	}

	/**
	 * \brief   Calculate the direction of the line normal vetor
	 *          <pre>
	 *                          ^y
	 *                          |
	 *                          |
	 *                    2     |     1
	 *                          |
	 *                          |
	 *          --+------+------+------+------+--> x
	 *                          |
	 *                          |
	 *                    3     |     4
	 *                          |
	 *                          |
	 *          </pre>
	 * \param   Float a a in (a,b)
	 * \param   Float b b in (a,b)
	 * \return  the case
	 */
	inline int which_case_4(Float a, Float b) {
		const int cVOF_whichcase4[2][2] = { { 1, 4 }, { 2, 3 } };
		return cVOF_whichcase4[a >= 0 ? 0 : 1][b >= 0 ? 0 : 1];
	}
	/**
	 * \brief   Calculate the direction of the line normal vetor
	 *          <pre>
	 *                          ^y
	 *                +         |         +
	 *                  +    3  |  2    +
	 *                    +     |     +
	 *                      +   |   +
	 *                  4     + | +      1
	 *          --+------+------+------+------+--> x
	 *                        + | +
	 *                  5   +   |   +    8
	 *                    +     |     +
	 *                  +     6 | 7     +
	 *                +         |         +
	 *          </pre>
	 * \param   Float a a in ax+by=C
	 * \param   Float b b in ax+by=C
	 * \return  the case
	 */
	inline int _which_case_8(Float a, Float b) {
		static const int cVOF_whichcase8[2][2][2] = { { { 2, 1 }, { 7, 8 } }, {
				{ 3, 4 }, { 6, 5 } } };
		return cVOF_whichcase8[a >= 0 ? 0 : 1][b >= 0 ? 0 : 1][
				abs(b) >= abs(a) ? 0 : 1];
	}

	/**
	 * \brief   known a,b in ax+by=alpha and C, calculate alpha \n
	 *
	 *          return alpha, ax+by=alpha;
	 * \param   Float a a in ax+by=alpha
	 * \param   Float b b in ax+by=alpha
	 * \param   Float C the color function
	 * \return  Line
	 */
	Line _cal_line(vt a, vt b, vt C) {
		if (a == 0.0 && b == 0.0) {
			a = 2 * SMALL;
		} else if (IsEqual(a, 0.0) && IsEqual(b, 0.0)) {
			a = 2 * SMALL * a / Abs(a);
			b = 2 * SMALL * b / Abs(b);
		}
		Line res;
		Float alpha;
		Float absa = Abs(a);
		Float absb = Abs(b);
		Float m, n;

		n = (absa >= absb) ? (absa) : absb;
		m = (absa <= absb) ? (absa) : absb;
		alpha = _cal_alpha(a, b, C);
		switch (_which_case_8(a, b)) {
		case 1: {
			res.reconstruct(n, m, alpha);
			break;
		}
		case 2: {
			res.reconstruct(m, n, alpha);
			break;
		}
		case 3: {
			res.reconstruct(-m, n, alpha - m);
			break;
		}
		case 4: {
			res.reconstruct(-n, m, alpha - n);
			break;
		}
		case 5: {
			res.reconstruct(-n, -m, alpha - m - n);
			break;
		}
		case 6: {
			res.reconstruct(-m, -n, alpha - n - m);
			break;
		}
		case 7: {
			res.reconstruct(m, -n, alpha - n);
			break;
		}
		case 8: {
			res.reconstruct(n, -m, alpha - m); //
			break;
		}
		}
		return res;
	}

	Segment _cal_interct_points(const Line &l, vt c1, vt c2) {
		vt m1 = l.a();
		vt m2 = l.b();
		vt alpha = l.alpha();

		vt am1, am2, aalpha;
		int lcase = which_case_4(m1, m2);
		switch (lcase) {
		case 1: {
			am1 = m1;
			am2 = m2;
			aalpha = alpha;
			break;
		}
		case 2: {
			am1 = -m1;
			am2 = m2;
			aalpha = alpha - m1 * c1;
			break;
		}
		case 3: {
			am1 = -m1;
			am2 = -m2;
			aalpha = alpha - m1 * c1 - m2 * c2;
			break;
		}
		case 4: {
			am1 = m1;
			am2 = -m2;
			aalpha = alpha - m2 * c2;
			break;
		}
		}

		ASSERT(!(IsZero(am1) && IsZero(am2)));
		//This step produce error;
		if (IsZero(am1)) {
			am1 = SMALL;
		}
		if (IsZero(am2)) {
			am2 = SMALL;
		}

		//The first special case
		if (aalpha < 0) {
			return Segment();
		}
		//The second special case
		if ((aalpha - am1 * c1) / am2 > c2 && (aalpha - am2 * c2) / am1 > c1) {
			return Segment();
		}

		//The normal case
		Point2D p1, p2;
		if (StepFun(aalpha - c1 * am1) == 1) {
			//Point F
			p1.x() = c1;
			p1.y() = (aalpha - am1 * c1) / am2;
		} else {
			//Point E
			p1.x() = aalpha / am1;
			p1.y() = 0.0;
		}
		if (StepFun(aalpha - c2 * am2) == 1) {
			//Point G
			p2.x() = (aalpha - am2 * c2) / am1;
			p2.y() = c2;
		} else {
			//Point E
			p2.x() = 0.0;
			p2.y() = aalpha / am2;
		}
		if (p1 == p2) {
			std::cerr << " >! Intersect Point Equal. \n"
					<< "(VOF.h -->> Segment2D calInterctPoints(const Line2D &l, Float c1, Float c2)"
					<< std::endl;
			return Segment();
		} else {
			//The region on the leftside of segment equal ONE
			switch (lcase) {
			case 1: {
				return Segment(p2, p1);
				break;
			}
			case 2: {
				p1.x() = -p1.x() + c1;
				p2.x() = -p2.x() + c1;
				return Segment(p1, p2);
				break;
			}
			case 3: {
				p1.x() = -p1.x() + c1;
				p2.x() = -p2.x() + c1;
				p1.y() = -p1.y() + c2;
				p2.y() = -p2.y() + c2;
				return Segment(p2, p1);
				break;
			}
			case 4: {
				p1.y() = -p1.y() + c2;
				p2.y() = -p2.y() + c2;
				return Segment(p1, p2);
				break;
			}
			default: {
				std::cerr << "Return case error. "
						<< "(VOF.h -->> Segment2D calInterctPoints(const Line2D &l, Float c1, Float c2)"
						<< std::endl;

			}
			}
			return Segment();
		}
	}

};

}
#endif
