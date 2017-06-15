
#ifndef _MATRIX_SMALL_HPP_
#define _MATRIX_SMALL_HPP_

#include <iostream>
#include <iterator>

#include "algebra_define.hpp"

namespace carpio {

template<typename T, St DIM1, St DIM2>
class MatrixS {
public:
    // type definitions===================
    typedef MatrixS<T, DIM1, DIM2> self;
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef St size_type;
    typedef St difference_type;

    T elems[DIM1 > 0 ? DIM1 : 1][DIM2 > 0 ? DIM2 : 1];

public:
    //constructor==========================
    MatrixS() {
    }
    MatrixS(const T& a) {
        this->assign(a);
    }
    MatrixS(const MatrixS<T, DIM1, DIM2>& a) {
        typedef size_type st;
        for (st i = 0; i < DIM1; i++) {
            for (st j = 0; j < DIM2; j++) {
                elems[i][j] = a[i][j];
            }
        }
    }
    //=============================================
    self& operator=(const self &a) {
        if (this == &a) {
            return *this;
        } else {
            typedef size_type st;
            for (st i = 0; i < DIM1; i++) {
                for (st j = 0; j < DIM2; j++) {
                    this->elems[i][j] = a[i][j];
                }
            }
        }
        return *this;
    }
    //=============================================
    ~MatrixS() {
    }
    //Capacity=====================================
    static inline size_type size() {
        return DIM1 * DIM2;
    }
    static inline size_type size_i() {
        return DIM1;
    }
    static inline size_type size_j() {
        return DIM2;
    }
    //Element access===============================
    reference operator()(size_type i, size_type j) {
        ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
        return elems[i][j];
    }
    const_reference operator()(size_type i, size_type j) const {
        ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
        return elems[i][j];
    }
    reference at(size_type i, size_type j) {
        ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
        return elems[i][j];
    }
    const_reference at(size_type i, size_type j) const {
        ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
        return elems[i][j];
    }
    T get(size_type i, size_type j) {
        ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
        return elems[i][j];
    }
    //T* getpValue(size_type i, size_type j) {
    //	ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
    //	return elems[i][j];
    //}
    void set(size_type i, size_type j, const T& value) {
        ASSERT_MSG(i < DIM1 && j < DIM2, "out of range");
        elems[i][j] = value;
    }
    void set_row(size_type i, const T& value) {
        ASSERT_MSG(i < DIM1, "out of range");
        typedef size_type st;
        for (st j = 0; j < DIM2; j++) {
            this->elems[i][j] = value;
        }
    }
    void set_col(size_type j, const T& value) {
        ASSERT_MSG(j < DIM2, "out of range");
        typedef size_type st;
        for (st i = 0; i < DIM2; i++) {
            this->elems[i][j] = value;
        }
    }
    void assign(const T& value) {
        typedef size_type st;
        for (st i = 0; i < DIM1; i++) {
            for (st j = 0; j < DIM2; j++) {
                this->elems[i][j] = value;
            }
        }
    }
    inline void ones() {
        assign(1);
    }
    inline void zeros() {
        assign(0);
    }
    //
    void show() const {
        std::cout << "> MatrixS " << DIM1 << " x " << DIM2 << "\n";
        std::cout << "> ";
        for (int i = 0; i < DIM1; i++) {
            for (int j = 0; j < DIM2; j++) {
                std::cout << std::scientific << this->elems[i][j] << "  ";
            }
            std::cout << std::endl;
            std::cout << "> ";
        }
        std::cout << "< ----------\n";
    }
};

typedef MatrixS<Float, 2, 2> Matrix2x2;
typedef MatrixS<Float, 3, 3> Matrix3x3;
typedef MatrixS<Float, 4, 4> Matrix4x4;

template<typename T, St DIM1, St DIM2>
static inline void zeros(MatrixS<Float, DIM1, DIM2>& m) {
    m.assign(0);
}

template<typename T, St DIM1, St DIM2>
static inline void ones(MatrixS<Float, DIM1, DIM2>& m) {
    m.assign(1);
}

/*
 * \brief calculate the determinant of a 2x2 matrix.
 *
 *        Adapted from:
 *        Matrix Inversion
 *        by Richard Carling
 *        from "Graphics Gems", Academic Press, 1990
 */
template<typename T>
static inline T det2x2(const T& a, const T& b, const T& c, const T& d) {
    return a * d - b * c;
}

template<typename T>
inline T determinant(const MatrixS<T, 2, 2>& m) {
    T a, b, c, d;
    a = m(0, 0);
    b = m(0, 1);
    c = m(1, 0);
    d = m(1, 1);
    return det2x2(a, b, c, d);
}

/*
 * \brief calculate the determinant of a 3x3 matrix
 *        in the form
 *
 *        | a1,  b1,  c1 |
 *        | a2,  b2,  c2 |
 *        | a3,  b3,  c3 |
 *
 *        Adapted from:
 *        Matrix Inversion
 *        by Richard Carling
 *        from "Graphics Gems", Academic Press, 1990
 */
template<typename T>
static inline T det3x3(const T& a1, const T& a2, const T& a3, const T& b1,
                       const T& b2, const T& b3, const T& c1, const T& c2, const T& c3) {
    T ans3;
    ans3 = a1 * det2x2(b2, b3, c2, c3) - b1 * det2x2(a2, a3, c2, c3)
           + c1 * det2x2(a2, a3, b2, b3);
    return ans3;
}

template<typename T>
inline T determinant(const MatrixS<T, 3, 3>& m) {
    const T& a1 = m(0, 0), b1 = m(0, 1), c1 = m(0, 2);
    const T& a2 = m(1, 0), b2 = m(1, 1), c2 = m(1, 2);
    const T& a3 = m(2, 0), b3 = m(2, 1), c3 = m(2, 2);
    return det3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
}

/**
 * \brief matrix_determinant:
 * \param m Matrix 4 x 4 .
 *
 * \returns: the value of det(\p m).
 */
template<typename T>
inline T determinant(const MatrixS<T, 4, 4>& m) {
    T ans4;
    const T& a1 = m(0, 0), b1 = m(0, 1), c1 = m(0, 2), d1 = m(0, 3);
    const T& a2 = m(1, 0), b2 = m(1, 1), c2 = m(1, 2), d2 = m(1, 3);
    const T& a3 = m(2, 0), b3 = m(2, 1), c3 = m(2, 2), d3 = m(2, 3);
    const T& a4 = m(3, 0), b4 = m(3, 1), c4 = m(3, 2), d4 = m(3, 3);

    ans4 = a1 * det3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4)
           - b1 * det3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4)
           + c1 * det3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4)
           - d1 * det3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);

    return ans4;
}

template<typename T>
inline void adjoint(const MatrixS<T, 4, 4>& m, MatrixS<T, 4, 4>& ma) {
    const T& a1 = m(0, 0), b1 = m(0, 1), c1 = m(0, 2), d1 = m(0, 3);
    const T& a2 = m(1, 0), b2 = m(1, 1), c2 = m(1, 2), d2 = m(1, 3);
    const T& a3 = m(2, 0), b3 = m(2, 1), c3 = m(2, 2), d3 = m(2, 3);
    const T& a4 = m(3, 0), b4 = m(3, 1), c4 = m(3, 2), d4 = m(3, 3);
    /* row column labeling reversed since we transpose rows & columns */
    // the matrix without transpose is called cofactor matrix
    // the matrix with    transpose is called adjoint  matrix
    ma(0, 0) = det3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
    ma(1, 0) = -det3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
    ma(2, 0) = det3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
    ma(3, 0) = -det3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);

    ma(0, 1) = -det3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
    ma(1, 1) = det3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
    ma(2, 1) = -det3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
    ma(3, 1) = det3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);

    ma(0, 2) = det3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
    ma(1, 2) = -det3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
    ma(2, 2) = det3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
    ma(3, 2) = -det3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);

    ma(0, 3) = -det3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
    ma(1, 3) = det3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
    ma(2, 3) = -det3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
    ma(3, 3) = det3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
}

template<typename T>
inline void adjoint(const MatrixS<T, 3, 3>& m, MatrixS<T, 3, 3>& ma) {
    /* row column labeling reversed since we transpose rows & columns */
    // the matrix without transpose is called cofactor matrix
    // the matrix with    transpose is called adjoint  matrix
    ma(0, 0) = (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1));
    ma(0, 1) = (m(2, 1) * m(0, 2) - m(0, 1) * m(2, 2));
    ma(0, 2) = (m(0, 1) * m(1, 2) - m(1, 1) * m(0, 2));
    ma(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2));
    ma(1, 1) = (m(0, 0) * m(2, 2) - m(2, 0) * m(0, 2));
    ma(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2));
    ma(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1));
    ma(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1));
    ma(2, 2) = (m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0));
}

template<typename T>
inline int inverse(const MatrixS<T, 4, 4>& m, MatrixS<T, 4, 4>& m_inv) {
    typedef MatrixS<T, 4, 4> matrix;
    T det = determinant(m);
    if (det == 0.) {
        return -1;
    }
    adjoint(m, m_inv);
    for (typename matrix::size_type i = 0; i < 4; ++i) {
        for (typename matrix::size_type j = 0; j < 4; ++j) {
            m_inv(i, j) /= det;
        }
    }
    return 1;
}

template<typename T>
inline int inverse(const MatrixS<T, 3, 3>& m, MatrixS<T, 3, 3>& m_inv) {
    typedef MatrixS<T, 3, 3> matrix;
    T det = determinant(m);
    if (det == 0.) {
        return -1;
    }
    adjoint(m, m_inv);
    for (typename matrix::size_type i = 0; i < m.iLen(); ++i) {
        for (typename matrix::size_type j = 0; j < m.jLen(); ++j) {
            m_inv(i, j) /= det;
        }
    }
    return 1;
}

template<typename T>
void transpose(const MatrixS<T, 4, 4> m, MatrixS<T, 4, 4> mt) {
    mt(0, 0) = m(0, 0);
    mt(1, 0) = m(0, 1);
    mt(2, 0) = m(0, 2);
    mt(3, 0) = m(0, 3);
    mt(0, 1) = m(1, 0);
    mt(1, 1) = m(1, 1);
    mt(2, 1) = m(1, 2);
    mt(3, 1) = m(1, 3);
    mt(0, 2) = m(2, 0);
    mt(1, 2) = m(2, 1);
    mt(2, 2) = m(2, 2);
    mt(3, 2) = m(2, 3);
    mt(0, 3) = m(3, 0);
    mt(1, 3) = m(3, 1);
    mt(2, 3) = m(3, 2);
    mt(3, 3) = m(3, 3);
}

template<typename T>
void transpose(const MatrixS<T, 3, 3> m, MatrixS<T, 3, 3> mt) {
    mt(0, 0) = m(0, 0);
    mt(1, 0) = m(0, 1);
    mt(2, 0) = m(0, 2);
    mt(0, 1) = m(1, 0);
    mt(1, 1) = m(1, 1);
    mt(2, 1) = m(1, 2);
    mt(0, 2) = m(2, 0);
    mt(1, 2) = m(2, 1);
    mt(2, 2) = m(2, 2);
}

template<typename T>
void transpose(const MatrixS<T, 2, 2> m, MatrixS<T, 2, 2> mt) {
    mt(0, 0) = m(0, 0);
    mt(1, 0) = m(0, 1);
    mt(0, 1) = m(1, 0);
    mt(1, 1) = m(1, 1);
}

}
#endif
