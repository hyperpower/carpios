#ifndef TYPEDEFINE_HPP
#define TYPEDEFINE_HPP

#include "stddef.h"
#include <assert.h>


#define ASSERT(expr) assert(expr)
#define ASSERT_MSG(expr, msg) assert((expr)&&(msg))

#define SHOULD_NOT_REACH assert((false)&&(" >! Should not reach"))
#define CAST(type, p)           ((type)p)
#define CAST_REF(type, p)       (*((type)p))
#define _IF_TRUE_RETRUN(expr)   if(expr){return;};
#define _IF_FALSE_RETRUN(expr)  if(false==(expr)){return;};
#define SMALL (1e-10)
//return code
#define _SUCCESS   0
#define _ERROR     1
#define _WARNING   2


#define _RETURN_VAL_IF_FAIL(expr,val)  {       \
		if (!(expr))                           \
			return (val);                      \
	    }

#define _RETURN_IF_FAIL(expr)  {               \
		if (!(expr))                           \
			return ;                           \
	    }

namespace carpio {
// value type
typedef unsigned int St; //size type
typedef int Int;
typedef unsigned int uInt;
typedef double Float;
typedef void* utPointer;
typedef const void* const_utPointer;


enum Axes {
	_X_ = 0, //
	_Y_ = 1, //
	_Z_ = 2, //
};
inline Axes ToAxes(const St& i) {
	ASSERT(i >= 0 && i < 3);
	switch (i) {
	case 0:
		return _X_;
	case 1:
		return _Y_;
	case 2:
		return _Z_;
	default:
		ASSERT_MSG(false, "Error input i");
	}
	SHOULD_NOT_REACH;
	return _X_;
}

template<typename C>
struct IsIterable
{
    typedef char true_type;
    typedef long false_type;

    template<class T>
    static true_type  is_beg_iterable(
    		int i,
            typename T::const_iterator = C().begin());
    template<class T>
    static false_type is_beg_iterable(...);

    enum { value = sizeof(is_beg_iterable<C>(0)) == sizeof(true_type) };
};

/*
 * geometry
 */
enum Orientation {
	_M_ = 0, //
	_P_ = 1, //
	_C_ = 2, //
};

}




#endif
