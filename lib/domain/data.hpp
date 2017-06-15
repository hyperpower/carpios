#ifndef DATA_H_
#define DATA_H_

#include "domain_define.hpp"
#include "algebra/array_list.hpp"

#include <array>

namespace carpio {

template<typename VALUE, St DIM>
class Data_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM); //

	typedef VALUE vt;
	typedef Data_<VALUE, DIM> Self;

	typedef void (*pfunction)(Self*, utPointer);

protected:
	int _idx;
	ArrayListV<vt> _center;
	ArrayListV<vt> _face[NumFaces];
	//ArrayListT<utPointer> _face_untype[NumFaces];
	ArrayListV<vt> _vertex[NumVertexes];
	ArrayListT<utPointer> _untype;
public:
	Data_() {
		_idx = 0;
	}
	Data_(const St& nc, const St& nf, const St& nv, utPointer utp) :
			_center(nc), _untype(1) {
		_idx = 0;
		for (St i = 0; i < NumFaces; ++i) {
			_face[i].reconstruct(nf);
			//_face_untype[i].reconstruct(0);
		}
		for (St i = 0; i < NumVertexes; ++i) {
			_vertex[i].reconstruct(nv);
		}
		_untype[0] = utp;
	}
	Data_(const St& nc, const St& nf, const St& nv, const St& nutp) :
			_center(nc), _untype(nutp) {
		_idx = 0;
		for (St i = 0; i < NumFaces; ++i) {
			_face[i].reconstruct(nf);
			//_face_untype[i].reconstruct(0);
		}
		for (St i = 0; i < NumVertexes; ++i) {
			_vertex[i].reconstruct(nv);
		}
		_untype.assign(nullptr);
	}
	Data_(const St& nc, const St& nf, const St& nv, const St& nutp,
			const St& nfutp) :
			_center(nc), _untype(nutp) {
		_idx = 0;
		for (St i = 0; i < NumFaces; ++i) {
			_face[i].reconstruct(nf);
			//_face_untype[i].reconstruct(nfutp);
		}
		for (St i = 0; i < NumVertexes; ++i) {
			_vertex[i].reconstruct(nv);
		}
		_untype.assign(nullptr);
	}
	Data_(const Self& s) :
			_center(s._center.size()), _untype(s._untype.size()) {
		_idx = s._idx;
		for (St i = 0; i < NumFaces; ++i) {
			_face[i] = s._face[i];
			//_face_untype[i] = s._face_untype[i];
		}
		for (St i = 0; i < NumVertexes; ++i) {
			_vertex[i] = s._vertex[i];
		}
		_center = s._center;
		_untype = s._untype;
	}
	void reconstruct(const St& nc, const St& nf, const St& nv, const St& nutp) {
		_center.reconstruct(nc);
		for (St i = 0; i < NumFaces; ++i) {
			_face[i].reconstruct(nf);
		}
		for (St i = 0; i < NumVertexes; ++i) {
			_vertex[i].reconstruct(nv);
			//_face_untype[i].reconstruct(0);
		}
		_untype.reconstruct(nutp);
		_untype.assign(nullptr);
	}

	inline int& idx() {
		return _idx;
	}
	inline const int& idx() const {
		return _idx;
	}
	inline vt& center(const St& i) {
		ASSERT(i < _center.size());
		return _center[i];
	}

	inline const vt& center(const St& i) const {
		ASSERT(i < _center.size());
		return _center[i];
	}
	inline bool has_center(const St& i) const {
		return i < _center.size();
	}

	inline ArrayListV<vt>& face(const Direction& d) {
		St idx_face = FaceDirectionInOrder(d);
		return _face[idx_face];
	}

	inline const ArrayListV<vt>& face(const Direction& d) const {
		St idx_face = FaceDirectionInOrder(d);
		return _face[idx_face];
	}

	inline ArrayListV<vt>& vertex(const Direction d, const St& i) {
		ASSERT(i < this->NumVertexes);
		return _vertex[i];
	}

	inline const ArrayListV<vt>& vertex(const Direction d, const St& i) const {
		ASSERT(i < this->NumVertexes);
		return _vertex[i];
	}

	inline utPointer& utp(const St& i) {
		ASSERT(i < _untype.size());
		return _untype[i];
	}

	inline const_utPointer& utp(const St& i) const {
		ASSERT(i < _untype.size());
		return _untype[i];
	}

	/*
	 *  resize
	 */
	void resize(const St& nc, const St& nf, const St& nv, const St& nutp) {
		this->resize_center(nc);
		for (St i = 0; i < NumFaces; ++i) {
			this->resize_face(i, nf);
		}
		for (St i = 0; i < NumVertexes; ++i) {
			this->resize_vertex(i, nv);
		}
		this->resize_utp(nutp);
	}
	void resize_center(St len) {
		ASSERT(len >= 0);
		this->_center.resize(len);
	}
	void resize_face(St idx, St len) {
		ASSERT(idx >= 0 && idx < NumFaces);
		this->_face[idx].resize(len);
	}
	void resize_vertex(St idx, St len) {
		ASSERT(idx >= 0 && idx < NumVertexes);
		this->_vertex[idx].resize(len);
	}
	void resize_utp(St len) {
		ASSERT(len >= 0);
		St len_o = this->_untype.size();
		this->_untype.resize(len);
		for(St i = len_o; i< len;i++){
			this->_untype[i] = nullptr;
		}
	}

	bool empty() const {
		bool res = true;
		res = res && (_center.size() == 0);
		for (int i = 0; i < NumFaces; ++i) {
			res = res && (_face[i].size() == 0);
		}
		for (int i = 0; i < NumVertexes; ++i) {
			res = res && (_vertex[i].size() == 0);
		}
		return res;
	}

	void show_info() const {
		std::cout << "center data:" << this->_center.size() << "\n";
		std::cout << "face data   :" << "\n";
		for (int i = 0; i < NumFaces; ++i) {
			std::cout << "    " << i << "       :" << _face[i].size() << "\n";
		}
		for (int i = 0; i < NumVertexes; ++i) {
			std::cout << "    " << i << "       :" << _vertex[i].size() << "\n";
		}
	}

};

template<typename COO_VALUE, typename VALUE, int DIM>
class PData_ {  //Point Data
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM); //

	static const int Flag_Invalid = 0;
	static const int Flag_Center = 1;
	static const int Flag_Face = 2;
	static const int Flag_Vertex = 3;

	typedef VALUE vt;
	typedef COO_VALUE cvt;
	typedef PData_<COO_VALUE, VALUE, DIM> Self;
	typedef PData_<COO_VALUE, VALUE, DIM>& ref_Self;
	typedef const PData_<COO_VALUE, VALUE, DIM>& const_ref_Self;
protected:
	cvt _p[Dim];
	ArrayListV<int> _flag;
	ArrayListV<vt> _val;
	ArrayListV<St> _idx;
	/*
	 * assign p
	 */
	void _assign_p(const cvt& x, const cvt& y = 0, const cvt& z = 0) {
		_p[0] = x;
		if (Dim >= 2) {
			_p[1] = y;
		}
		if (Dim >= 3) {
			_p[2] = z;
		}
	}
	void _assign_p_origin() {
		for (St i = 0; i < Dim; ++i) {
			_p[i] = 0.0;
		}
	}
	void _copy_p(const cvt _lp[]) {
		for (St i = 0; i < Dim; ++i) {
			_p[i] = _lp[i];
		}
	}
	bool _is_same_size() const {
		return (_val.size() == _flag.size() && _idx.size() == _flag.size());
	}
public:
	/*
	 *  constructor
	 */
	PData_() :
			_flag(), _val(), _idx() {
		_assign_p_origin();
	}
	PData_(St n) :
			_flag(n), _val(n), _idx(n) {
		_flag.assign(0);
		_assign_p_origin();
	}
	PData_(St n, cvt x, cvt y = 0, cvt z = 0) :
			_flag(n), _val(n), _idx(n) {
		_flag.assign(0);
		_assign_p(x, y, z);
	}
	PData_(ArrayListV<St> arridx, cvt x, cvt y = 0, cvt z = 0) :
			_flag(arridx.size()), _val(arridx.size()), _idx(arridx.size()) {
		_flag.assign(0);
		_idx = arridx;
		_assign_p(x, y, z);
	}
	PData_(const Self& _pd) :
			_flag(_pd._flag.size()), _val(_pd._val.size()), _idx(
					_pd._idx.size()) {
		_flag = _pd._flag;
		_val = _pd._val;
		_idx = _pd._idx;
		_copy_p(_pd._p);
	}
	ref_Self& operator=(const const_ref_Self& _pd) {
		if (this == &_pd) {
			return *this;
		}
		_flag = _pd._flag;
		_val = _pd._val;
		_idx = _pd._idx;
		_copy_p(_pd._p);
		return *this;
	}
	/*
	 *  I/O
	 */
	St size() const {
		return _flag.size();
	}

	inline cvt& x() {
		return _p[0];
	}

	inline const cvt& x() const {
		return _p[0];
	}

	inline cvt& y() {
		ASSERT(Dim >= 2);
		return _p[1];
	}

	inline const cvt& y() const {
		ASSERT(Dim >= 2);
		return _p[1];
	}

	inline cvt& z() {
		ASSERT(Dim >= 3);
		return _p[2];
	}

	inline const cvt& z() const {
		ASSERT(Dim >= 3);
		return _p[2];
	}

	inline cvt& p(Axes axes) {  //point
		switch (axes) {
		case _X_: {
			return x();
			break;
		}
		case _Y_: {
			return y();
			break;
		}
		case _Z_: {
			return z();
			break;
		}
		default:
			SHOULD_NOT_REACH;
			return x();
		}
	}
	inline const cvt& p(Axes axes) const {  //point const
		switch (axes) {
		case _X_: {
			return x();
			break;
		}
		case _Y_: {
			return y();
			break;
		}
		case _Z_: {
			return z();
			break;
		}
		default:
			SHOULD_NOT_REACH;
			return x();
		}
	}

	inline void set_point(cvt x, cvt y = 0, cvt z = 0) {
		_assign_p(x, y, z);
	}

	inline vt& val(St n) {
		ASSERT(is_valid(n));
		return _val[n];
	}

	inline const vt& val(St n) const {
		ASSERT(is_valid(n));
		return _val[n];
	}

	inline St& idx(St n) {
		return _idx[n];
	}

	inline const St& idx(St n) const {
		return _idx[n];
	}

	inline int& flag(St n) {
		return _flag[n];
	}

	inline const int& flag(St n) const {
		return _flag[n];
	}

	inline ArrayListV<vt>& arr_val() {
		return _val;
	}

	inline const ArrayListV<vt>& arr_val() const {
		return _val;
	}

	inline ArrayListV<St>& arr_idx() {
		return _idx;
	}

	inline const ArrayListV<St>& arr_idx() const {
		return _idx;
	}

	inline void set_val(St i, const vt& val, St idx, int flag) {
		_val[i] = val;
		_idx[i] = idx;
		_flag[i] = flag;
	}

	inline bool is_valid(St i) const {
		if (_flag[i] == Flag_Invalid) {
			return false;
		} else {
			return true;
		}
	}

	inline bool has_valid_val() const {
		for (St i = 0; i < _flag.size(); ++i) {
			if (_flag[i] != Flag_Invalid) {
				return true;
			}
		}
		return false;
	}

	inline St count_valid_val() const {
		St res = 0;
		for (St i = 0; i < _flag.size(); ++i) {
			if (_flag[i] != Flag_Invalid) {
				res++;
			}
		}
		return res;
	}

	inline St count_invalid_val() const {
		return size() - count_valid_val();
	}

	void set_all_center() {
		for (St i = 0; i < _flag.size(); ++i) {
			_flag[i] = Flag_Center;
		}
	}

	/*
	 *  show
	 */

	void show_point(std::string name = "") const {
		if (name != "") {
			std::cout << name << "\n";
		}
		for (St i = 0; i < Dim; ++i) {
			std::cout << "p" << i << " = " << _p[i] << "\n";
		}
	}
};

}

#endif /* DATA_H_ */
