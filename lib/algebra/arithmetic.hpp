/************************
 //  \file   Arithmetic.h
 //  \brief
 //
 //  \author czhou
 //  \date   12 févr. 2014
 ***********************/
#ifndef ARITHMETIC_H_
#define ARITHMETIC_H_

#include "../type_define.hpp"

#include "array_list.hpp"
#include <math.h>
#include <iostream>

namespace carpio {

const Float PI = 3.14159265358979323846;
const Float EXP = 2.71828182845904523536;

template<typename TYPE>
inline int StepFun(TYPE x) {
    return (x <= 0) ? 0 : 1;
}

inline int Sign(Float x) {
    if (x < 0.0) {
        return -1;
    } else if (x > 0.0) {
        return 1;
    } else {
        return 0;
    }
}

enum Range {
    _oo_, _oc_, _co_, _cc_,
};

template<typename TYPE>
inline bool IsInRange(const TYPE &down, const TYPE &value, const TYPE &up,
                      const Range &range) {
    switch (range) {
    case _oo_:
        return (down < value && value < up);
    case _oc_:
        return (down < value && value <= up);
    case _co_:
        return (down <= value && value < up);
    case _cc_:
        return (down <= value && value <= up);
    }
    return false;
}
template<class ST>
class IdxRange {
protected:
    ArrayListV<ST> _arr;
public:
    typedef typename ArrayListV<ST>::iterator iterator;
    typedef typename ArrayListV<ST>::const_iterator const_iterator;
    IdxRange(const ST& down, const ST& up) {
        ASSERT(down < up);
        _arr.reconstruct(up - down);
        _arr.assign_forward(down, 1);
    }
    //iterater
    iterator begin() {
        return _arr.begin();
    }
    const_iterator begin() const {
        return _arr.begin();
    }
    iterator end() {
        return _arr.end();
    }
    const_iterator end() const {
        return _arr.end();
    }

};

// round to n digit of a
// roundto(3.145,1)=3.1    roundto(3.145,2)=3.15
inline Float RoundTo(Float a, int n) {
    return round(a * pow(10, n)) / pow(10, n);
}
// this function return a^2+b^2
inline Float SqareSum(Float &a, Float &b) {
    return a * a + b * b;
}

// this function return (a+b)*(a+b)
inline Float SumSqare(Float &a, Float &b) {
    return (a + b) * (a + b);
}

inline int CountSignificanceDigit(Float a) {
    for (int i = 0; i < 30; i++) {
        if (RoundTo(a, i) == a) {
            return i;
        }
    }
    return -1;
}

inline Float VolumeOfCircularTruncatedCone(Float R, Float r, Float h) {
    return PI * h * (R * R + r * r + R * r) / 3;
}

template<class TYPE>
inline TYPE Max(TYPE a, TYPE b, bool (*Comp_ge)(const TYPE &, const TYPE &)) {
    return Comp_ge(a, b) ? a : b;
}

template<class TYPE>
inline TYPE Max(TYPE a, TYPE b, TYPE c,
                bool (*Comp_ge)(const TYPE &, const TYPE &)) {
    TYPE tmp = Comp_ge(a, b) ? a : b;
    return Comp_ge(tmp, c) ? tmp : c;
}

template<class TYPE>
inline TYPE Min(TYPE a, TYPE b, bool (*Comp_le)(const TYPE &, const TYPE &)) {
    return Comp_le(a, b) ? a : b;
}

template<class TYPE>
inline TYPE Min(TYPE a, TYPE b, TYPE c,
                bool (*Comp_le)(const TYPE &, const TYPE &)) {
    TYPE tmp = Comp_le(a, b) ? a : b;
    return Comp_le(tmp, c) ? tmp : c;
}

template<class TYPE>
inline TYPE Mid(TYPE a, TYPE b, TYPE c,
                bool (*Comp_ge)(const TYPE &, const TYPE &)) {
    int idx = Comp_ge(a, b) ? 1 : 2;
    if (idx == 1)
        idx = Comp_ge(a, c) ? 1 : 3;
    else
        idx = Comp_ge(b, c) ? 2 : 3;
    if (idx == 1)
        return Comp_ge(b, c) ? b : c;
    else if (idx == 2)
        return Comp_ge(a, c) ? a : c;
    else
        return Comp_ge(a, b) ? a : b;
}

template<class TYPE>
inline void Sort(const TYPE &a, const TYPE &b, const TYPE &c,  //
                 bool (*Comp)(const TYPE &, const TYPE &), //
                 TYPE &big, TYPE &mid, TYPE &small) {
    int idx = Comp(a, b) ? 1 : 2;
    if (idx == 1)
        idx = Comp(a, c) ? 1 : 3;
    else
        idx = Comp(b, c) ? 2 : 3;
    if (idx == 1) {
        big = a;
        if (Comp(b, c)) {
            mid = b;
            small = c;
        } else {
            mid = c;
            small = b;
        }
        return;
    }
    if (idx == 2) {
        big = b;
        if (Comp(a, c)) {
            mid = a;
            small = c;
        } else {
            mid = c;
            small = a;
        }
        return;
    }
    big = c;
    if (Comp(a, b)) {
        mid = a;
        small = b;
    } else {
        mid = b;
        small = a;
    }
}

template<class TYPE>
bool CompGreat(const TYPE &a, const TYPE &b) {
    return a > b;
}

template<class TYPE>
bool CompLess(const TYPE &a, const TYPE &b) {
    return a < b;
}

template<class TYPE>
inline void Sort(const TYPE &a, const TYPE &b, const TYPE &c,  //
                 bool (*Comp)(const TYPE &, const TYPE &), //
                 int &big, int &mid, int &small) {
    int idx = Comp(a, b) ? 0 : 1;
    if (idx == 0)
        idx = Comp(a, c) ? 0 : 2;
    else
        idx = Comp(b, c) ? 1 : 2;
    if (idx == 0) {
        big = 0;
        if (Comp(b, c)) {
            mid = 1;
            small = 2;
        } else {
            mid = 2;
            small = 1;
        }
        return;
    }
    if (idx == 1) {
        big = 1;
        if (Comp(a, c)) {
            mid = 0;
            small = 2;
        } else {
            mid = 2;
            small = 0;
        }
        return;
    }
    big = 2;
    if (Comp(a, b)) {
        mid = 0;
        small = 1;
    } else {
        mid = 1;
        small = 0;
    }
}

template<class TYPE>
inline void Swap(TYPE &a, TYPE &b) //
{
    TYPE tmp = a;
    a = b;
    b = tmp;
}

template<class TYPE>
inline void SortIncrease(TYPE &a, TYPE &b, TYPE &c) //
{
    if (b < a) {
        swap(a, b);
    }
    if (c < a) {
        swap(a, c);
    }
    if (c < b) {
        swap(b, c);
    }
}

template<class TYPE>
int SolveQuadraticEquation(const TYPE &a, const TYPE &b, const TYPE &c,
                           Float &x1, Float &x2) {
    Float discri = 0;
    int numroot;
    discri = b * b - 4.0 * a * c;
    if (discri == 0) {
        numroot = 1;
    } else if (discri > 0) {
        numroot = 2;
    } else {
        numroot = 0;
    }

    if (numroot == 2) {
        x1 = (-b - sqrt(discri)) / 2 / a;
        x2 = (-b + sqrt(discri)) / 2 / a;
        return 2;
    } else if (numroot == 1) {
        x1 = -b / 2 / a;
        x2 = x1;
        return 1;
    } else {
        return 0;
    }
}

template<class TYPE>
int SolveCubicEquation(const TYPE &a, const TYPE &b, const TYPE &c,
                       const TYPE &d, Float &x1, Float &x2, Float &x3) {
    ASSERT(a != 0);
    Float A = b * b - 3.0 * a * c;
    Float B = b * c - 9.0 * a * d;
    Float C = c * c - 3.0 * b * d;
    Float discri = B * B - 4.0 * A * C;
    //case 1 has three equal real roots
    if (A == 0 && B == 0) {
        x1 = -b / 3.0 / a;
        x2 = x1;
        x3 = x1;
        return 1;
    }
    if (discri > 0) {
        Float Y1 = A * b + 1.5 * a * (-B + sqrt(discri));
        Float Y2 = A * b + 1.5 * a * (-B - sqrt(discri));
        Float cuberY1 = Y1 < 0 ? -pow(-Y1, 1.0 / 3.0) : pow(Y1, 1.0 / 3.0);
        Float cuberY2 = Y2 < 0 ? -pow(-Y2, 1.0 / 3.0) : pow(Y2, 1.0 / 3.0);
        x1 = (-b - cuberY1 - cuberY2) / 3.0 / a;
        //ignore complex roots
        x2 = x1;
        x3 = x1;
        return 2;
    }
    if (discri == 0) {
        Float K = B / A;
        x1 = -b / a + K;
        x2 = K / 2.0;
        x3 = x2;
        SortIncrease(x1, x2, x3);
        return 3;
    }
    if (discri < 0) {
        Float T = (2.0 * A * b - 3.0 * a * B) / (2.0 * pow(A, 1.5));
        Float sita3 = acos(T) / 3.0;
        x1 = (-b - 2.0 * sqrt(A) * cos(sita3)) / (3.0 * a);
        x2 = (-b + sqrt(A) * (cos(sita3) + sqrt(3.0) * sin(sita3))) / (3.0 * a);
        x3 = (-b + sqrt(A) * (cos(sita3) - sqrt(3.0) * sin(sita3))) / (3.0 * a);
        SortIncrease(x1, x2, x3);
        return 4;
    }
    return -1;
}

inline Float Rand(Float r1, Float r2) {
    ASSERT(r1 != r2);
    Float rnum1 = rand() % 100;
    Float rmax = std::max(r1, r2);
    Float rmin = std::min(r1, r2);
    return rmin + (rmax - rmin) / 100 * rnum1;
}

template<class TYPE>
inline TYPE Max(const TYPE& a, const TYPE& b) {
    return a >= b ? a : b;
}
template<class TYPE>
inline TYPE Max(const TYPE& a, const TYPE& b, const TYPE& c) {
    return Max(Max(a, b), c);
}
template<class TYPE>
inline TYPE Max(const TYPE a[], St len) {
    TYPE max = a[0];
    for (St i = 1; i < len; ++i) {
        if (max < a[i]) {
            max = a[i];
        }
    }
    return max;
}
template<class TYPE>
inline TYPE Min(const TYPE& a, const TYPE& b) {
    return a <= b ? a : b;
}
template<class TYPE>
inline TYPE Abs(const TYPE& s) {
    return s < 0 ? -s : s;
}
template<class TYPE>
inline bool IsEqual(const TYPE& a, const TYPE& b) {
    return Abs(a - b) < SMALL;
}
template<class TYPE>
inline bool IsZero(const TYPE& a) {
    return Abs(a - 0.0) < SMALL;
}
// Greater, equal or less
template<class TYPE>
inline int GEL(const TYPE& a, const TYPE& v) {
    //greater equal or less
    if (v < a) {
        return -1;
    } else if (v == a) {
        return 0;
    } else {
        return 1;
    }
}

template<class CVT, class VT>
struct Interp_ {

    ///
    template<class DIM, class POI, class ARR>
    static VT Linear_all(const DIM& dim, const POI& p, const POI& p1,
                         const POI& p2,
                         const ARR& av) {
        if (dim == 1) {
            return Linear(p[0], p1[0], p2[0], av[0], av[1]);
        }
        if (dim == 2) {
            return Bilinear(
                       p[0], p[1],
                       p1[0], p1[1], //
                       p2[0], p2[1], //
                       av[0], av[1], av[2], av[3]);
        }
        if (dim == 3) {
            return Trilinear(     //
                       p[0], p[1], p[2],  //
                       p1[0], p1[1], p1[2],  //
                       p2[0], p2[1], p2[2],  //
                       av[0], av[1], av[2], av[3], av[4], av[5], av[6], av[7]);
        }
    }

    /// Linear interpolate.
    /// Linear reference
    /// u is a normalized u
    static VT Linear_r(CVT u, const VT& a, const VT& b) {
        return ((a * (1.0 - u)) + (b * u));
    }

    static VT Linear(CVT x, CVT x1, CVT x2, const VT& a, const VT& b) {
        CVT u = (x - x1) / (x2 - x1);
        return Linear_r(u, a, b);
    }

    /// Bilinear interpolate.
    /* Corners are:
     *   1b -----2c
     * ^  |       |
     * |  |       |
     * v  |       |
     *   0a -----3d
     *     u->
     */
    static VT Bilinear_r(CVT u, CVT v, const VT& a, const VT& b, const VT& c,
                         const VT& d) {
        return ((a * (1 - u) * (1 - v)) + (d * (1 - u) * v) + (b * u * (1 - v))
                + (c * u * v));
    }

    static VT Bilinear(     //
        CVT x, CVT y,   //
        CVT x1, CVT y1, //
        CVT x2, CVT y2, //
        const VT& a, const VT& b, const VT& c, const VT& d) {
        CVT u = (x - x1) / (x2 - x1);
        CVT v = (y - y1) / (y2 - y1);
        return Bilinear_r(u, v, a, b, c, d);
    }

    /// Trilinear interpolate.
    static VT Trilinear_r(CVT u, CVT v, CVT w, const VT& a, const VT& b,
                          const VT& c, const VT& d, const VT& e, const VT& f, const VT& g,
                          const VT& h) {
        VT t0(Bilinear_r(u, v, a, b, c, d));
        VT t1(Bilinear_r(u, v, e, f, g, h));

        return (((1 - w) * t0) + (w * t1));
    }

    static VT Trilinear(
        //
        CVT x, CVT y,
        CVT z, //
        CVT x1, CVT y1,
        CVT z1, //
        CVT x2, CVT y2,
        CVT z2, //
        const VT& a, const VT& b, const VT& c, const VT& d, const VT& e,
        const VT& f, const VT& g, const VT& h) {
        CVT u = (x - x1) / (x2 - x1);
        CVT v = (y - y1) / (y2 - y1);
        CVT w = (z - z1) / (z2 - z1);
        return Trilinear_r(u, v, w, a, b, c, d, e, f, g, h);
    }

};

} //namespace Larus

#endif /* ARITHMETIC_H_ */
