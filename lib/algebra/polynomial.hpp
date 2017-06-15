/*
 * Expression.h
 *
 *  Created on: Feb 2, 2015
 *      Author: zhou
 */

#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include "algebra_define.hpp"
#include <utility>
#include <map>
#include <stdint.h>
#include <iomanip>
//#include "utility/hash_table.hpp"

namespace carpio {

template<class T>
struct IsZero_: std::unary_function<const T&, bool> {
	bool operator()(const T& number) const {
		return (number == 0);
	}
};

template<class COE, class TERM, class EXP,                                    //
		class ISZERO_COE = IsZero_<COE>,                                     //
		class ISZERO_EXP = IsZero_<EXP>,                                     //
		class COMPARE_TERM = std::less<TERM>,                                 //
		class COMPARE_EXP = std::less<EXP> >                                 //
class Polynomial_ {
protected:
	typedef TERM T;
	typedef EXP E;
	typedef COE C;
	typedef std::pair<TERM, EXP> Key;

	typedef ISZERO_COE is_zero_coe;
	typedef ISZERO_EXP is_zero_exp;
	ISZERO_COE _is_zero_coe;
	ISZERO_EXP _is_zero_exp;

	struct _compare {
		COMPARE_TERM _compare_term;
		COMPARE_EXP _compare_exp;

		ISZERO_EXP __is_zero_exp;
		bool operator()(const Key& lhs, const Key& rhs) const {
			if (__is_zero_exp(lhs.second) && __is_zero_exp(rhs.second)) {
				return false;
			}
			if (_compare_term(lhs.first, rhs.first)) {
				return true;
			} else if (!_compare_term(rhs.first, lhs.first)) {
				return _compare_exp(lhs.second, rhs.second);
			}
			return false;
		}
	};

	typedef std::map<Key, C, _compare> Map;
	Map _map;
public:

	typedef Polynomial_<COE, TERM, EXP, ISZERO_COE, ISZERO_EXP, COMPARE_TERM,
			COMPARE_EXP> Self;
	typedef Polynomial_<COE, TERM, EXP, ISZERO_COE, ISZERO_EXP, COMPARE_TERM,
			COMPARE_EXP>& ref_Self;
	typedef const Polynomial_<COE, TERM, EXP, ISZERO_COE, ISZERO_EXP,
			COMPARE_TERM, COMPARE_EXP>& const_ref_Self;
	typedef typename Map::const_reference const_reference;
	typedef typename Map::pointer pointer;
	typedef typename Map::const_pointer const_pointer;
	typedef typename Map::iterator iterator;
	typedef typename Map::const_iterator const_iterator;

	typedef typename Map::difference_type difference_type;
	typedef typename Map::size_type size_type;

	class Term: public std::pair<const Key, C> {
	public:
		Term() {
		}
		Term(const Key& key, const C& v) :
				std::pair<const Key, C>(key, v) {
		}
		Term(const std::pair<const Key, C>& p) :
				std::pair<const Key, C>(p.first, p.second) {
		}
		Term(const COE& coe, const TERM& term, const EXP& exp) :
				std::pair<const Key, C>(std::pair<TERM, EXP>(term, exp), coe) {
		}
		const COE& coe() const {
			return this->second;
		}
		COE& coe() {
			return this->second;
		}
		const TERM& term() const {
			return this->first.first;
		}
		TERM& term() {
			return this->first.first;
		}
		const TERM& exp() const {
			return this->first.second;
		}
		TERM& exp() {
			return this->first.second;
		}
    };

public:
	Polynomial_() :
			_map() {
	}
	Polynomial_(const_ref_Self rhs) :
			_map(rhs._map) {
		//std::cout<<"copy c"<<std::endl;
	}
	ref_Self operator=(const_ref_Self rhs) {
		//std::cout<<"operator="<<std::endl;
		this->_map = rhs._map;
		return *this;
	}
	/*
	 *  iterator
	 */
	iterator begin() {
		return _map.begin();
	}
	const_iterator begin() const {
		return _map.begin();
	}
	iterator end() {
		return _map.end();
	}
	const_iterator end() const {
		return _map.end();
	}
	/*
	 * capacity
	 */
	bool empty() const {
		return _map.empty();
	}
	size_type size() const {
		return _map.size();
	}
	C& operator[](const Key& k) {
		return _map[k];
	}
	const C& operator[](const Key& k) const {
		return _map[k];
	}

	static bool IsConstant(const iterator& iter) {
		ISZERO_EXP iz;
		return iz(iter->first.second);
	}
	static bool IsConstant(const const_iterator& iter) {
		ISZERO_EXP iz;
		return iz(iter->first.second);
	}

	static bool IsZeroCoe(const iterator& iter) {
		ISZERO_COE ic;
		return ic(iter->second);
	}

	static C GetCoe(const iterator& iter) {
		return iter->second;
	}
	static C GetCoe(const_iterator& iter) {
		return iter->second;
	}

	static E GetExp(const iterator& iter) {
		return iter->first.second;
	}

	static E GetExp(const_iterator& iter) {
		return iter->first.second;
	}

	static T GetTerm(const iterator& iter) {
		return iter->first.first;
	}

	static T GetTerm(const_iterator& iter) {
		return iter->first.first;
	}

protected:
	typedef std::pair<iterator, bool> _ret;
public:

	_ret insert(const std::pair<const Key, C> &x) {
		iterator it = _map.find(x.first);
		if (it != _map.end()) {
			it->second = it->second + x.second;
			return _ret(it, true);
		}
		return _map.insert(x);
	}
	_ret insert(const Term &x) {
		iterator it = _map.find(x.first);
		if (it != _map.end()) {
			it->second = it->second + x.second;
			return _ret(it, true);
		}
		return _map.insert(x);
	}

	_ret insert(const C& c, const T& t, const E& e) {
		const Key k(t, e);
		std::pair<const Key, C> p(k, c);
		return this->insert(p);
	}

	iterator find(const std::pair<const Key, C> &x) {
		return _map.find(x.first);
	}
	const_iterator find(const std::pair<const Key, C> &x) const {
		return _map.find(x.first);
	}

	iterator find(const Term &x) {
		return _map.find(x.first);
	}
	const_iterator find(const Term &x) const {
		return _map.find(x.first);
	}

	void erase(iterator& iter) {
		_map.erase(iter);
	}

	void clear() {
		_map.clear();
	}

	void plus(const_ref_Self poly) {
		for (Polynomial_::const_iterator iter = poly.begin();
				iter != poly.end(); ++iter) {
			this->insert((*iter));
		}
	}

	void minus(const_ref_Self poly) {
		for (Polynomial_::const_iterator iter = poly.begin();
				iter != poly.end(); ++iter) {
			const Key& k = iter->first;
			const C& v = iter->second;
			Term term(k, -v);       //minus
			this->insert(term);
		}
	}
	void times(const COE& rhs) {   //overload operator*
		for (Polynomial_::iterator iter = this->begin(); iter != this->end();
				++iter) {
			iter->second *= rhs;
		}
	}
	void divide(const COE& rhs) {  //overload operator/
		for (Polynomial_::iterator iter = this->begin(); iter != this->end();
				++iter) {
			iter->second = iter->second / rhs;
		}
	}

	void trim_zero() {
		iterator it_b = _map.begin();
		iterator it_e = _map.end();
		iterator it_t;

		while (it_b != it_e) {
			if (_is_zero_coe(it_b->second)) {  // Criteria checking here
				it_t = it_b;            // Keep a reference to the iter
				++it_b;                 // Advance in the map
				_map.erase(it_t);       // Erase it !!!
			} else {
				++it_b;                 // Just move on ...
			}
		}
	}

	void merge_const() {
		iterator it_b = _map.begin();
		iterator it_e = _map.end();
		iterator it_t;
		size_type flag = 0;
		std::pair<Key, C> ct;
		while (it_b != it_e) {
			if (_is_zero_exp(it_b->first.second)) {  // Criteria checking here
				if (flag == 0) {
					ct.first.first = it_b->first.first;
					ct.first.second = it_b->first.second;
					ct.second = it_b->second;
				} else {
					ct.second = ct.second + it_b->second;
				}
				flag++;
				it_t = it_b;            // Keep a reference to the iter
				++it_b;                 // Advance in the map
				_map.erase(it_t);       // Erase it !!!
			} else {
				++it_b;                 // Just move on ...
			}
		}
		if (flag > 0) {
			_map.insert(ct);
		}
	}

	void concise() {
		trim_zero();
		merge_const();
	}

	void show() const {
		// all the valuable should overload <<
		std::cout << "size = " << this->size() << "\n";
		for (Polynomial_::const_iterator iter = this->begin();
				iter != this->end(); ++iter) {
			std::cout << "  " << std::setw(5) << iter->second << "  ";
			std::cout << "  " << std::setw(5) << iter->first.first << "  ";
			std::cout << "  " << std::setw(5) << iter->first.second << "";
			std::cout << "\n";
		}
	}
};

//Polynomial operator-(const Polynomial &, const Polynomial &);
//Polynomial operator+(const Polynomial &, const Polynomial &);
//Polynomial operator*(const Float&, const Polynomial&);

}

#endif /* ALGEBRA_EXPRESSION_H_ */
