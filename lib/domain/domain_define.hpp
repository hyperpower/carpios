#ifndef DOMAIN_DEFINE_H_
#define DOMAIN_DEFINE_H_

#include "../type_define.hpp"
#include "algebra/arithmetic.hpp"

#include <string>
#include <sstream>

namespace carpio {

/*enum Direction {
 //D_X==================
 D_MMX_M = 8 + 0, //001 000
 D_MPX_M = 8 + 2, //001 010
 D_PMX_M = 8 + 4, //001 100
 D_PPX_M = 8 + 6, //001 110
 //
 D_MMX_P = 8 + 1, //001 001
 D_MPX_P = 8 + 3, //001 011
 D_PMX_P = 8 + 5, //001 101
 D_PPX_P = 8 + 7, //001 111

 //D_Y==================
 D_MYM_M = 16 + 0, //010 000
 D_MYP_M = 16 + 1, //010 001
 D_PYM_M = 16 + 4, //010 100
 D_PYP_M = 16 + 5, //010 101
 //
 D_MYM_P = 16 + 2, //010 010
 D_MYP_P = 16 + 3, //010 011
 D_PYM_P = 16 + 6, //010 110
 D_PYP_P = 16 + 7, //010 111

 //D_Z==================
 D_ZMM_M = 32 + 0, //100 000
 D_ZMP_M = 32 + 1, //100 001
 D_ZPM_M = 32 + 2, //100 010
 D_ZPP_M = 32 + 3, //100 011
 //
 D_ZMM_P = 32 + 4, //100 100
 D_ZMP_P = 32 + 5, //100 101
 D_ZPM_P = 32 + 6, //100 110
 D_ZPP_P = 32 + 7, //100 111

 //D_ZX==================
 D_ZMX_MM = 40 + 0, //101 000
 D_ZMX_MP = 40 + 1, //101 001
 D_ZMX_PM = 40 + 4, //101 100
 D_ZMX_PP = 40 + 5, //101 101
 //
 D_ZPX_MM = 40 + 2, //101 010
 D_ZPX_MP = 40 + 3, //101 011
 D_ZPX_PM = 40 + 6, //101 110
 D_ZPX_PP = 40 + 7, //101 111

 //D_ZY==================
 D_ZYM_MM = 48 + 0, //110 000
 D_ZYM_MP = 48 + 2, //110 010
 D_ZYM_PM = 48 + 4, //110 100
 D_ZYM_PP = 48 + 6, //110 110
 //
 D_ZYP_MM = 48 + 1, //110 001
 D_ZYP_MP = 48 + 3, //110 011
 D_ZYP_PM = 48 + 5, //110 101
 D_ZYP_PP = 48 + 7, //110 111

 //D_YX==================
 D_MYX_MM = 24 + 0, //011 000
 D_MYX_MP = 24 + 1, //011 001
 D_MYX_PM = 24 + 2, //011 010
 D_MYX_PP = 24 + 3, //011 011
 //
 D_PYX_MM = 24 + 4, //011 100
 D_PYX_MP = 24 + 5, //011 101
 D_PYX_PM = 24 + 6, //011 110
 D_PYX_PP = 24 + 7, //011 111

 //D_ZYX==================
 D_ZYX_MMM = 56 + 0, //111 000
 D_ZYX_MMP = 56 + 1, //111 001
 D_ZYX_MPM = 56 + 2, //111 010
 D_ZYX_MPP = 56 + 3, //111 011
 D_ZYX_PMM = 56 + 4, //111 100
 D_ZYX_PMP = 56 + 5, //111 101
 D_ZYX_PPM = 56 + 6, //111 110
 D_ZYX_PPP = 56 + 7, //111 111
 };*/

typedef unsigned short Direction;

static const Direction _XM_ = 8; //001 000
static const Direction _XP_ = 9; //001 001
static const Direction _YM_ = 16; //010 000
static const Direction _YP_ = 18; //010 010
static const Direction _ZM_ = 32; //100 000
static const Direction _ZP_ = 36; //100 100

static const int _E_ = 0;
static const int _F_C_ = 1;
static const int _C_F_ = 2;






enum Plane {
	_XY_ = 24, _YZ_ = 48, _ZX_ = 40,
};


inline bool GetBit(const Direction &d, St i) {
	const Direction ARR_BIT[6] = { 1, 2, 4, 8, 16, 32 };
	return (d & ARR_BIT[i]) == ARR_BIT[i];
}

inline bool IsP(const Orientation& ori) {
	return (ori == _P_);
}
inline bool IsM(const Orientation& ori) {
	return (ori == _M_);
}
inline bool IsC(const Orientation& ori) {
	return (ori == _C_);
}
inline Orientation Opposite(const Orientation& ori) {
	if (IsC(ori)) {
		return ori;
	} else if (IsP(ori)) {
		return _M_;
	} else {
		return _P_;
	}
}

inline Direction ToDirection(const Plane &p, const Orientation &o1,
		const Orientation &o2) {
	ASSERT(o1 != _C_);
	ASSERT(o2 != _C_);
	switch (p) {
	case _XY_:
		return p + o1 + (o2 << 1);
	case _YZ_:
		return p + (o1 << 1) + (o2 << 2);
	case _ZX_:
		return p + (o1 << 2) + o2;
	}
	SHOULD_NOT_REACH;
	return 0;
}

inline Direction ToDirection(const Orientation &x, const Orientation &y,
		const Orientation &z) {
	ASSERT(x != _C_);
	ASSERT(y != _C_);
	ASSERT(z != _C_);
	return 56 + x + (y << 1) + (z << 2);
}

inline unsigned short HI(const Direction &d) {
	return d >> 3;
}

inline unsigned short LO(const Direction &d) {
	return d & 7;
}
inline Direction Opposite(const Direction& ori) {
	Direction res = ori;
	unsigned short hi = HI(ori);
	//unsigned short hi = LO(ori);
	return (res ^ hi);
}
static const short COUNT_1[8] = { 0, 1, 1, 2, 1, 2, 2, 3 };

inline bool IsFaceDirection(const Direction &d) {
	return COUNT_1[(d >> 3)] == 1;
}
inline bool IsCornerDirection(const Direction &d) {
	return COUNT_1[(d >> 3)] == 2;
}
inline bool IsVertexDirection(const Direction &d) {
	return COUNT_1[(d >> 3)] == 3;
}

inline Direction ToFaceDirection(const Orientation &o, const Axes& a) {
	ASSERT(o != _C_);
	St lo = (o == _P_) ? 7 : 0;
	St hi = 0;
	switch (a) {
	case _X_: {
		hi = 1;
		break;
	}
	case _Y_: {
		hi = 2;
		break;
	}
	case _Z_: {
		hi = 4;
		break;
	}
	}
	St nlo = lo & hi;
	return (hi << 3) | nlo;
}

inline Direction ToCornerDirection( //
		const Orientation &o1, const Axes& a1, //
		const Orientation &o2, const Axes& a2 //
		) {
	ASSERT(o1 != _C_);
	ASSERT(o2 != _C_);
	Direction d1 = ToFaceDirection(o1, a1);
	Direction d2 = ToFaceDirection(o2, a2);
	return d1 | d2;
}

inline void FaceDirectionToOrientationAndAxes(const Direction &d,
		Orientation &o, Axes &a) {
	ASSERT(IsFaceDirection(d));
	unsigned short hi = d >> 3;
	if ((hi & 1) == 1) {
		a = _X_;
        o = (GetBit(d, 0)) ? _P_ : _M_;
		return;
	}
	if ((hi & 2) == 2) {
		a = _Y_;
		o = (GetBit(d, 1)) ? _P_ : _M_;
		return;
	}
	if ((hi & 4) == 4) {
		a = _Z_;
		o = (GetBit(d, 2)) ? _P_ : _M_;
		return;
	}
	SHOULD_NOT_REACH;
}
inline void CornerDirectionToOrientationAndAxes(const Direction &d,
		Orientation &o1, Axes &a1, Orientation &o2, Axes& a2) {
	ASSERT(IsCornerDirection(d));
	unsigned short hi = d >> 3;
	if ((hi & 1) == 1 && (hi & 2) == 2) {
		a1 = _X_;
		a2 = _Y_;
		o1 = (GetBit(d, 0)) ? _P_ : _M_;
		o2 = (GetBit(d, 1)) ? _P_ : _M_;
		return;
	}
	if ((hi & 2) == 2 && (hi& 4)== 4) {
		a1 = _Y_;
		a2 = _Z_;
		o1 = (GetBit(d, 1)) ? _P_ : _M_;
		o2 = (GetBit(d, 2)) ? _P_ : _M_;
		return;
	}
	if ((hi & 4) == 4 && (hi & 1) == 1) {
		a1 = _Z_;
		a2 = _X_;
		o1 = (GetBit(d, 2)) ? _P_ : _M_;
		o2 = (GetBit(d, 0)) ? _P_ : _M_;
		return;
	}
	SHOULD_NOT_REACH;
}
inline Axes FaceDirectionToAxes(const Direction &d) {
	ASSERT(IsFaceDirection(d));
	Axes a;
	unsigned short hi = d >> 3;
	if ((hi & 1) == 1) {
		a = _X_;
		return a;
	}
	if ((hi & 2) == 2) {
		a = _Y_;
		return a;
	}

	if ((hi & 4) == 4) {
		a = _Z_;
		return a;
	}
	SHOULD_NOT_REACH;
	return a;
}

/*
 * Does Direction on axes active
 */
inline bool IsDirectionOn(const Direction &d, const Axes& a) {
	unsigned short hi = d >> 3;
	switch (a) {
	case _X_:
		return ((hi & 1) == 1);
	case _Y_:
		return ((hi & 2) == 2);
	case _Z_:
		return ((hi & 4) == 4);
	}
	return false;
}

inline bool IsFacePDirection(const Direction &d) {
	unsigned short hi = d >> 3;
	unsigned short low = d & 7;
	return (COUNT_1[hi] == 1) && (hi & low) != 0;
}

inline bool IsXYDirection(const Direction &d) {
	return (d >> 3) == 3;
}

inline bool IsYZDirection(const Direction &d) {
	return (d >> 3) == 6;
}

inline bool IsZXDirection(const Direction &d) {
	return (d >> 3) == 5;
}


inline bool IsXYZDirection(const Direction &d) {
	return (d >> 3) == 7;
}

inline Orientation ToOrientation(const Direction &d, const Axes &a) {
	if (IsDirectionOn(d, a)) {
		switch (a) {
		case _X_:
			return (GetBit(d, 0) == true) ? _P_ : _M_;
		case _Y_:
			return (GetBit(d, 1) == true) ? _P_ : _M_;
		case _Z_:
			return (GetBit(d, 2) == true) ? _P_ : _M_;
		}
	} else {
		return _C_;
	}
}

inline Direction FaceDirectionInOrder(const St& i) {
	ASSERT(i < 6);
	static const Direction ARR_FD[] = { 8, 9, 16, 18, 32, 36 };
	return ARR_FD[i];
}
inline St FaceDirectionInOrder(const Direction& dir) {
	ASSERT(IsFaceDirection(dir));
	static const Direction ARR_FD[] = { 8, 9, 16, 18, 32, 36 };
	for (St i = 0; i < 6; i++) {
		if (ARR_FD[i] == dir) {
			return i;
		}
	}
	SHOULD_NOT_REACH;
	return 1;
}
inline St FaceDirectionInOrder(const Axes& a, const Orientation& o){
	// axes and orientation construct a direction
	Direction d = ToFaceDirection(o, a);
	return FaceDirectionInOrder(d);
}

inline Direction XYDirectionInOrder(const St& i) {
	ASSERT(i < 4);
	return 24 + i;
}

inline Direction YZDirectionInOrder(const St& i) {
	ASSERT(i < 4);
	return 48 + i * 2;
}

inline Direction ZXDirectionInOrder(const St& i) {
	ASSERT(i < 4);
	static const Direction ARR_ZXD[] = { 40 + 0, //101 000
	40 + 1, //101 001
	40 + 4, //101 100
	40 + 5, //101 101
			};
	return ARR_ZXD[i];
}

inline Direction XYZDirectionInOrder(const St& i) {
	ASSERT(i < 8);
	return 56 + i;
}

inline Direction DirectionInOrder(const St& i) {
	ASSERT(i < 26);
	if (i < 6)
		return FaceDirectionInOrder(i);
	if (i < 10)
		return XYDirectionInOrder(i - 6);
	if (i < 14)
		return YZDirectionInOrder(i - 10);
	if (i < 18)
		return ZXDirectionInOrder(i - 14);
	return XYZDirectionInOrder(i - 18);
}

inline std::string ToString(const Axes& a) {
	if (a == _X_) {
		return "x";
	}
	if (a == _Y_) {
		return "y";
	}
	return "z";
}
inline std::string ToString(const Orientation& a) {
	if (a == _M_) {
		return "m";
	}
	if (a == _P_) {
		return "p";
	}
	return "c";
}
inline std::string ToString(const Direction& d) {
	std::stringstream sst;
	std::stringstream sst2;
	unsigned short hi = HI(d);
	unsigned short lo = LO(d);
	if ((hi & 1) == 1) {
		sst << "X";
		sst2 << (((lo & 1) == 1) ? "P" : "M");
	}
	if ((hi & 2) == 2) {
		sst << "Y";
		sst2 << (((lo & 2) == 2) ? "P" : "M");
	}
	if ((hi & 4) == 4) {
		sst << "Z";
		sst2 << (((lo & 4) == 4) ? "P" : "M");
	}
	sst << "_" << sst2.str();
	return sst.str();
}



inline Axes VerticalAxes2D(const Axes& a) {
	ASSERT(a != _Z_);
	return a == _X_ ? _Y_ : _X_;
}

inline Axes VerticalAxes1(const Axes& a) {
	switch (a) {
	case _X_: {
		return _Y_;
	}
	case _Y_: {
		return _Z_;
	}
	case _Z_: {
		return _X_;
	}
	default: {
		SHOULD_NOT_REACH;
		return _X_;
	}
	}
}
inline Axes VerticalAxes2(const Axes& a) {
	switch (a) {
	case _X_: {
		return _Z_;
	}
	case _Y_: {
		return _X_;
	}
	case _Z_: {
		return _Y_;
	}
	default: {
		SHOULD_NOT_REACH;
		return _X_;
	}
	}
}


//default type
//typedef double CooValueType;
//typedef double ValueType;

}

#endif /* GRIDDEF_H_ */
