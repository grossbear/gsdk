////////////////////////////////////////////////////////////////////////////////
// mathdefs.h
//
// Math library defines
////////////////////////////////////////////////////////////////////////////////

#ifndef _MATH_LIB_DEFINES_H_
#define _MATH_LIB_DEFINES_H_

// Primary functions that has SSE implementation
#define MATH_PRIM_SSE

// Library uses it's own conversion from floating types to integer types
//#define MATH_NO_FPU

// Library uses ILM Function converting float type to half type
#define MATH_ILM_HALF_CONVERSION

// Build classes and functions instantiation using following data types
#define MATH_FLOAT_INST
//#define MATH_DOUBLE_INST
//#define MATH_LONG_DOUBLE_INST
//#define MATH_FIXED_INST
//#define MATH_FIXED64_INST
#define MATH_HALF_INST

#define MATH_STD_SQRT 0
#define MATH_FAST_SQRT
#define MATH_TABLE_SQRT
#define MATH_SERIES_SQRT

#define MATH_STD_TRIGONOMETRY
#define MATH_TABLE_TRIGONOMETRY
#define MATH_SERIES_TRIGONOMETRY


#endif //_MATH_LIB_DEFINES_H_
