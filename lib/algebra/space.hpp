/************************
 //  \file   Space.h
 //  \brief
 //
 //  \author czhou
 //  \date   20 janv. 2015
 ***********************/
#ifndef SPACE_H_
#define SPACE_H_

#include <iostream>
#include <assert.h>

#include "../type_define.hpp"

#include "array_list.hpp"

namespace carpio {

template<typename T, St DIM>
class SpaceT {
  public:
    static const St Dim = DIM;
    // type definitions===================
    typedef T value_t;
    typedef SpaceT<value_t, Dim> self;
    typedef SpaceT<value_t, Dim>* pself;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef St size_type;
    //typedef St difference_type;

    typedef typename ArrayListT<value_t>::iterator iterator;
    typedef typename ArrayListT<value_t>::const_iterator const_iterator;
    //static const int DIM = Dim;
  private:
    std::array<size_type, Dim> m_len;
    ArrayListT<value_t> m_mp;
  public:
    //constructor==========================
    SpaceT();
    SpaceT(const self& a);
    SpaceT(size_type iLen, size_type = 0, size_type = 0);

    void reconstruct(size_type iLen, size_type = 0, size_type = 0);
    //=============================================
    SpaceT<T, DIM>& operator=(const SpaceT<T, DIM>& a);
    //=============================================
    ~SpaceT() {
    }

    //Capacity=====================================
    size_type size() const {
        return m_mp.size();
    }
    size_type size_i() const {
        return m_len[0];
    }
    size_type size_j() const {
        return Dim >= 2 ? m_len[1] : 0;
    }
    size_type size_k() const {
        return Dim >= 3 ? m_len[2] : 0;
    }
    bool empty() const {
        return m_mp.empty();
    }
    /*
     *  iterator
     */
    iterator begin() {
        return m_mp.begin();
    }
    const_iterator begin() const {
        return m_mp.begin();
    }
    iterator end() {
        return m_mp.end();
    }
    const_iterator end() const {
        return m_mp.end();
    }
    //Element access===============================
    size_type to_1d_idx(size_type i, size_type = 0, size_type = 0) const;

    reference operator()(size_type i, size_type = 0, size_type = 0);
    const_reference operator()(size_type i,
                               size_type = 0,
                               size_type = 0) const;
    reference at(size_type i, size_type = 0, size_type = 0);
    const_reference at(size_type i,
                       size_type = 0,
                       size_type = 0) const;

    reference at_1d(size_type i);
    const_reference at_1d(size_type i) const;

    T get(size_type i,
          size_type = 0,
          size_type = 0);
    void set(const T& value,
             size_type i,
             size_type = 0,
             size_type = 0);
    void assign(const T& value);
    //element access===============================

    T* getpValue(size_type i,
                 size_type = 0,
                 size_type = 0);  //not good

    inline bool check_idx(size_type dim, size_type idx) const {
        ASSERT(dim < Dim);
        if (idx >= 0 && idx < m_len[dim]) {
            return true;
        } else {
            return false;
        }
    }
    inline bool check_idx_ijk(size_type i, size_type j, size_type k) const {
        return check_idx(0, i) && ((Dim >= 2) ? check_idx(1, j) : true)
               && ((Dim >= 3) ? check_idx(2, k) : true);
    }

    inline size_type count_equal(const T& nd) const { //overload ==
        return m_mp.count_equal(nd);
    }
};

template<typename T, St Dim>
SpaceT<T, Dim>::SpaceT() {
    m_len.fill(0);
}
template<typename T, St Dim>
SpaceT<T, Dim>::SpaceT(const SpaceT<T, Dim>& a) {
    this->m_len = a.m_len;
    this->m_mp = a.m_mp;
}
template<typename T, St Dim>
SpaceT<T, Dim>::SpaceT(size_type iLen, size_type jLen, size_type kLen) {
    this->m_len[0] = iLen;
    if (Dim >= 2) {
        ASSERT(iLen > 0 && jLen > 0);
        this->m_len[1] = jLen;
        this->m_mp.reconstruct(iLen * jLen);
    }
    if (Dim >= 3) {
        ASSERT(iLen > 0 && jLen > 0 && kLen > 0);
        this->m_len[2] = kLen;
        this->m_mp.reconstruct(iLen * jLen * kLen);
    }
}
template<typename T, St Dim>
void SpaceT<T, Dim>::reconstruct(size_type iLen, size_type jLen,
                                 size_type kLen) {
    size_type Len = 0;
    this->m_len[0] = iLen;
    if (Dim == 1) {
        Len = iLen;
    }
    if (Dim >= 2) {
        ASSERT(iLen > 0 && jLen > 0);
        this->m_len[1] = jLen;
        Len = iLen * jLen;
    }
    if (Dim >= 3) {
        ASSERT(iLen > 0 && jLen > 0 && kLen > 0);
        this->m_len[2] = kLen;
        Len = iLen * jLen * kLen;
    }
    this->m_mp.reconstruct(Len);
}
template<typename T, St DIM>
SpaceT<T, DIM>& SpaceT<T, DIM>::operator=(const SpaceT<T, DIM>& a) {
    if (this == &a) {
        return *this;
    }
    this->m_len = a.m_len;
    this->m_mp = a.m_mp;
    return *this;
}
template<typename T, St DIM>
typename SpaceT<T, DIM>::size_type SpaceT<T, DIM>::to_1d_idx(
    SpaceT<T, DIM>::size_type i, SpaceT<T, DIM>::size_type j,
    SpaceT<T, DIM>::size_type k) const {
    ASSERT(i < this->m_len[0]);
    if (Dim >= 2)
        ASSERT(j < this->m_len[1]);
    if (Dim >= 3)
        ASSERT(k < this->m_len[2]);
    std::array<size_type, Dim> inp;
    inp[0] = i;
    if (Dim >= 2) {
        inp[1] = j;
    }
    if (Dim >= 3) {
        inp[2] = k;
    }
    size_type idx = 0;
    for (size_type ii = 0; ii < Dim; ii++) {
        size_type b = 1;
        for (size_type jj = ii + 1; jj < Dim; jj++) {
            b *= m_len[jj];
        }
        idx += b * inp[ii];
    }
    return idx;
}

template<typename T, St DIM>
T& SpaceT<T, DIM>::at(SpaceT<T, DIM>::size_type i, SpaceT<T, DIM>::size_type j,
                      SpaceT<T, DIM>::size_type k) {
    size_type idx = to_1d_idx(i, j, k);
    return m_mp[idx];
}
template<typename T, St DIM>
const T& SpaceT<T, DIM>::at(SpaceT<T, DIM>::size_type i,
                            SpaceT<T, DIM>::size_type j, SpaceT<T, DIM>::size_type k) const {
    size_type idx = to_1d_idx(i, j, k);
    return m_mp[idx];
}
template<typename T, St DIM>
typename SpaceT<T, DIM>::reference SpaceT<T, DIM>::at_1d(
    SpaceT<T, DIM>::size_type i) {
    return m_mp[i];
}
template<typename T, St DIM>
typename SpaceT<T, DIM>::const_reference SpaceT<T, DIM>::at_1d(
    SpaceT<T, DIM>::size_type i) const {
    return m_mp[i];
}
template<typename T, St DIM>
T& SpaceT<T, DIM>::operator()(SpaceT<T, DIM>::size_type i,
                              SpaceT<T, DIM>::size_type j, SpaceT<T, DIM>::size_type k) {
    return at(i, j, k);
}
template<typename T, St DIM>
const T& SpaceT<T, DIM>::operator()(SpaceT<T, DIM>::size_type i,
                                    SpaceT<T, DIM>::size_type j, SpaceT<T, DIM>::size_type k) const {
    return at(i, j, k);
}
template<typename T, St DIM>
T SpaceT<T, DIM>::get(SpaceT<T, DIM>::size_type i, SpaceT<T, DIM>::size_type j,
                      SpaceT<T, DIM>::size_type k) {
    return at(i, j, k);
}
template<typename T, St DIM>
void SpaceT<T, DIM>::set(const T& value, SpaceT<T, DIM>::size_type i,
                         SpaceT<T, DIM>::size_type j, SpaceT<T, DIM>::size_type k) {
    this->at(i, j, k) = value;
}
template<typename T, St DIM>
void SpaceT<T, DIM>::assign(const T& value) {
    m_mp.assign(value);
}
//template<typename T, St DIM>
//T* SpaceT<T, DIM>::getpValue(SpaceT<T, DIM>::size_type i,
//		SpaceT<T, DIM>::size_type j, SpaceT<T, DIM>::size_type k) {
//	return &at(i, j, k);
//}

}

#endif /* SPACE_H_ */
