/***********************************************************************
 * Header:
 *    vector_list
 * Summary:
 *    This class contains the notion of a vector: a bucket to hold
 *    data for the user. This is just a simple version of std::vector
 *
 *    This will contain the class definition of:
 *        Vector         : A class that holds stuff
 *
 * Author
 *    Chengsi Zhou revise the code from Br. Helfrich
 ************************************************************************/

#ifndef _VECTOR_LIST_HPP_
#define _VECTOR_LIST_HPP_

#include "algebra_define.hpp"
#include <iostream>

namespace carpio {

/************************************************
 * CONTAINER
 * A class that holds stuff
 ***********************************************/
template<class T>
class VectorList_ {
  public:
    // type definitions
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef T& reference;
    typedef const T& const_reference;
    typedef St size_type;
    typedef St difference_type;

  private:
    T* _data;       // dynamically allocated array of T
    St _num;        // how many items are currently in the Vector?
    St _cap;        // how many items can I put on the Vector before full?

  public:
    // default constructor : empty and kinda useless
    VectorList_() :
        _num(0), _cap(0), _data(nullptr) {
    }

    // copy constructor : copy it
    VectorList_(const VectorList_ & rhs) throw (const char*);

    // non-default constructor : pre-allocate
    VectorList_(St cap) throw (const char*);

    ~VectorList_() {
        if (_data != nullptr) {
            delete[] _data;
            _data = nullptr;
        }
    }

    // is the container currently empty
    inline bool empty() const {
        return _num == 0;
    }

    // remove all the items from the vector
    void clear();

    // how many items are currently in the vector?
    St size() const {
        return _num;
    }

    void resize(St n) {
        ASSERT(n >= 0);
        if (n < _cap) {
            _num = n;
        } else {
            _cap = n;
            T* old = _data;
            _allocate();
            for (int i = 0; i < _num; i++)
                _data[i] = old[i];
            delete[] old;
        }
    }

    // add an item to the vector
    void push_back(const T& t) throw (const char*);

    // remove an item from a vector
    void pop_back() throw (const char*) {
        SHOULD_NOT_REACH;
    }

    St capacity() const {
        return _cap;
    }

    VectorList_& operator =(const VectorList_<T>& lhs);

    // Index ===================================
    reference at(St index) {
        ASSERT(index >= 0 && index < _num);
        return _data[index];
    }
    const_reference at(St index) const {
        ASSERT(index >= 0 && index < _num);
        return _data[index];
    }
    reference operator [](St index) {
        ASSERT(index >= 0 && index < _num);
        return _data[index];
    }
    const_reference operator [](St index) const {
        ASSERT(index >= 0 && index < _num);
        return _data[index];
    }
    reference operator()(St index) {
        ASSERT(index >= 0 && index < _num);
        return _data[index];
    }
    const_reference operator()(St index) const {
        ASSERT(index >= 0 && index < _num);
        return _data[index];
    }

    // Iterator =================================

    // return an iterator to the beginning of the list
    iterator begin() {
        return _data;
    }
    const_iterator begin() const {
        return _data;
    }
    // return an iterator to the end of the list
    iterator end() {
        return _data + _num;
    }
    const_iterator end() const {
        return _data + _num;
    }

    iterator rbegin() {
        return _data + (_num - 1);
    }
    const_iterator rbegin() const {
        return _data + (_num - 1);
    }
    iterator rend() {
        return _data - 1;
    }
    const_iterator rend() const {
        return _data - 1;
    }

  protected:
    void _allocate();
};

template<class T>
void VectorList_<T>::_allocate() {
    // attempt to allocate
    try {
        _data = new T[_cap];
    } catch (std::bad_alloc&) {
        throw "ERROR: Unable to allocate buffer";
    }
}

/*******************************************
 * VECTOR :: COPY CONSTRUCTOR
 *******************************************/
template<class T>
VectorList_<T>::VectorList_(const VectorList_<T>& rhs) throw (const char*) {
    //*this = rhs;
    ASSERT(rhs._cap >= 0);
    // do nothing if there is nothing to do
    if (rhs._cap == 0) {
        _cap = _num = 0;
        _data = nullptr;
        return;
    }
    _allocate();
    // copy over the _cap and size
    ASSERT(rhs._num >= 0 && rhs._num <= rhs._cap);
    _cap = rhs._cap;
    _num = rhs._num;
    // copy the items over one at a time using the assignment operator
    for (int i = 0; i < _num; i++)
        _data[i] = rhs._data[i];
    // the rest needs to be filled with the default value for T
    for (int i = _num; i < _cap; i++)
        _data[i] = T();
}

/**********************************************
 * VECTOR : NON-DEFAULT CONSTRUCTOR
 * Preallocate the vector to "n"
 **********************************************/
template<class T>
VectorList_<T>::VectorList_(St cap) throw (const char*) {
    ASSERT(cap >= 0);
    // do nothing if there is nothing to do
    if (cap == 0) {
        this->_cap = this->_num = 0;
        this->_data = nullptr;
        return;
    }
    // copy over the stuff
    this->_cap = cap;
    this->_num = 0;
    _allocate();
    // initialize the container by calling the default constructor
    //for (int i = 0; i < _num; i++)
    //	_data[i] = T();
}

template<class T>
VectorList_<T>& VectorList_<T>::operator =(const VectorList_<T>& copy) {
    this->_num = copy._num;
    this->_cap = copy._cap;
    _allocate();
    for (int i = 0; i < _num; i++)
        _data[i] = copy._data[i];
    for (int i = _num; i < _cap; i++)
        _data[i] = T();
    return *this;
}

/***************************************************
 * CONTAINER :: PUSH_BACK
 * Insert an item on the end of the container
 **************************************************/
template<class T>
void VectorList_<T>::push_back(const T& t) throw (const char*) {
    if (_num == 0) {
        _cap = 1;
        _allocate();
    } else if (_cap == _num) {
        T* old = _data;
        _cap *= 2;
        _allocate();
        for (St i = 0; i < _num; i++)
            _data[i] = old[i];
        delete[] old;
    }
    _data[_num++] = t;
}

template<class T>
void VectorList_<T>::clear() {
    _num = 0;
}

}

#endif
