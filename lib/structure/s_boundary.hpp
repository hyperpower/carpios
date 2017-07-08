#ifndef _S_BOUNDARY_HPP
#define _S_BOUNDARY_HPP

#include "s_define.hpp"
#include <map>
#include <utility>

namespace structure {

class BoundaryCondition {
public:
	typedef Vt vt;
	typedef Vt cvt;
	typedef std::function<vt(vt, cvt, cvt, cvt)> Fun;
	typedef BoundaryCondition Self;
	static const int _BC1_ = 1;
	static const int _BC2_ = 2;
	static const int _BC3_ = 3; // periodic
protected:
	// data
	// 1 Bondary conditon type
	// 2 function
	int _type;
	Fun _function;
	Fun _function2;
	Fun _function3;
public:
	// Constructor
	BoundaryCondition() {
		// default boundary condition is symmetric boundary condition
		_type = _BC2_;
		Fun df = [](vt t, cvt x, cvt y, cvt z) {return 0;};
		_function = df;
		_function2 = df;
		_function3 = df;
	}
	/*
	 * this constructor should not used to BC2
	 */
	BoundaryCondition(int type, Fun fun) :
			_type(type), _function(fun), _function2(fun), _function3(fun) {
	}
	BoundaryCondition(int type, Fun fun_x, Fun fun_y, Fun fun_z) :
			_type(type), _function(fun_x), _function2(fun_y), _function3(fun_z) {
	}
	BoundaryCondition(const Self& self) :
			_type(self._type), _function(self._function), _function2(
					self._function2), _function3(self._function3) {
	}
	// get
	int get_type() const {
		return _type;
	}
	vt get_val(vt t, cvt x, cvt y, cvt z) const {
		return _function(t, x, y, z);
	}
	/*
	 * for the vector type boundary condition
	 * (x,y,z) is the location
	 * Axes    is axes
	 */
	vt get_val(vt t, cvt x, cvt y, cvt z, St a) const {
		switch (a) {
		case _X_: {
			return _function(t, x, y, z);
			break;
		}
		case _Y_: {
			return _function2(t, x, y, z);
			break;
		}
		case _Z_: {
			return _function3(t, x, y, z);
			break;
		}
		}
		SHOULD_NOT_REACH;
		return 0;
	}
	// set
	void set_function(Fun fun, short a = _X_) {
		switch (a) {
		case _X_: {
			_function = fun;
			break;
		}
		case _Y_: {
			_function2 = fun;
			break;
		}
		case _Z_: {
			_function3 = fun;
			break;
		}
		}
	}
	void set_function(Fun fun_x, Fun fun_y, Fun fun_z) {
		_function = fun_x;
		_function2 = fun_y;
		_function3 = fun_z;
	}
	void set_default_1_bc(const vt& val) {
		_type = _BC1_;
		Fun f = [val](vt t, cvt x, cvt y, cvt z) {return val;};
		_function = f;
		_function2 = f;
		_function3 = f;
	}
	void set_default_1_bc(Fun fun) {
		_type = _BC1_;
		Fun f = fun;
		_function = f;
		_function2 = f;
		_function3 = f;
	}
	void set_default_2_bc(const vt& val) {
		_type = _BC2_;
		Fun f = [val](vt t,cvt x, cvt y, cvt z) {return val;};
		_function = f;
		_function2 = f;
		_function3 = f;
	}
	void set_default_2_bc(Fun fun) {
		_type = _BC2_;
		Fun f = fun;
		_function = f;
		_function2 = f;
		_function3 = f;
	}

};

class BoundaryIndex {
protected:
	struct BCID {
		St shape_idx;
		St seg_idx;
		std::string val_name;
	};

	struct BCID_compare {
		typedef BCID BCid;
		bool operator()(const BCid& lhs, const BCid& rhs) const {
			if (lhs.shape_idx < rhs.shape_idx) {
				return true;
			} else if (lhs.shape_idx == rhs.shape_idx) {
				if (lhs.seg_idx < rhs.seg_idx) {
					return true;
				} else if (lhs.seg_idx == rhs.seg_idx) {
					return lhs.val_name < rhs.val_name;
				}
			}
			return false;
		}

	};
public:
	typedef Vt cvt;
	typedef Vt vt;
	typedef BoundaryCondition BC;
	typedef std::shared_ptr<BoundaryCondition> spBC;
	typedef std::shared_ptr<const BoundaryCondition> spcBC;
	//typedef BoundaryCondition* pBoundaryCondition;
	//typedef const BoundaryCondition<cvt, vt>* const_pBoundaryCondition;
	//typedef BCID_compare_<cvt,vt> BCID_compare;
	typedef std::pair<const BCID, spBC> BCNode;
	typedef std::list<BCNode> list_BCNode;
	typedef std::map<BCID, spBC, BCID_compare> BCMap;
protected:
	BCMap _BCmap;
	spcBC _pdefault_BC;

	typedef typename BCMap::iterator iterator;
	typedef typename BCMap::const_iterator const_iterator;

public:
	//constructor
	BoundaryIndex() :
			_BCmap() {
		_pdefault_BC = spcBC(new BC());
	}
	~BoundaryIndex() {
	}
	//
	void insert(St si, St segi, const std::string& vi, spBC pbc) {
		//
		BCID key = { si, segi, vi };
		//key.seg_idx = segi;
		//key.shape_idx = si;
		//key.val_idx = vi;
		BCNode bcn(key, pbc);
		_BCmap.insert(bcn);
	}

	spcBC find(St si, St segi, const std::string& vali) const {
		BCID key;
		key.seg_idx = segi;
		key.shape_idx = si;
		key.val_name = vali;
		const_iterator it = _BCmap.find(key);
		if (it != _BCmap.end()) {
			// found
			return (it->second);
		} else {
			// not found
			return _pdefault_BC;
		}
	}

	void copy(const std::string& vn1, const std::string& vn2) {
		list_BCNode lbc;
		for (const_iterator it = _BCmap.begin(); it != _BCmap.end(); ++it) {
			const BCID& key = it->first;
			if (key.val_name == vn1) {
				BCID nkey = { key.shape_idx, key.seg_idx, vn2 };
				spBC spbc = it->second;
				BCNode bcn(nkey, spbc);
				lbc.push_back(bcn);
			}
		}
		for (auto term : lbc) {
			this->_BCmap.insert(term);
		}
	}

	void show() const {
		for (auto term : this->_BCmap) {
			std::cout << term.first.shape_idx << " " << term.first.seg_idx
					<< " " << term.first.val_name;
			std::cout << " " << term.second->get_type() << "\n";
		}
	}

};





}
#endif
