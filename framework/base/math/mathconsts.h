////////////////////////////////////////////////////////////////////////////////
//  mathconsts.h
//
//  Mathematical Constants
//
////////////////////////////////////////////////////////////////////////////////

#ifndef _MATHLIB_CONSTS_
#define _MATHLIB_CONSTS_

////////////////////////////////////////////////////////////////////////////////

// Constants defines
#define CONST_PI            3.1415926535897932384626433832795
#define CONST_2PI           6.2831853071795864769252867665590
#define CONST_DEG_IN_RAD    0.0174532925199432957692369076848
#define CONST_RAD_IN_DEG    57.295779513082320876798154814114
#define CONST_PI_BY_2       1.5707963267948966192313216916398
#define CONST_PI_BY_4       0.7853981633974483096156608458119

#define CONST_1_BY_PI       0.31830988618379067153776752674502
#define CONST_2_BY_PI       0.63661977236758134307553505349006
#define CONST_1_BY_2PI      0.15915494309189533576888376337251

#define CONST_SQRT_2        1.4142135623730950488016887242097
#define CONST_1_BY_SQRT_2   0.7071067811865475244008443621048
#define CONST_2_BY_SQRT_PI  1.1283791670955125738961589031215
#define CONST_SQRT_3        1.7320508075688772935274463415059

#define CONST_E             2.7182818284590452353602874713527
#define CONST_LOG2_E        1.4426950408889634073599246810019

#define CONST_LOG10_E       0.4342944819032518276511289189166
#define CONST_LN_2          0.6931471805599453094172321214581
#define CONST_LN_10         2.3025850929940456840179914546844
#define CONST_PHI           1.6180339887498948482045888343656


////////////////////////////////////////////////////////////////////////////////
#ifdef _WIN32
#define MATH_INT64_MAX  0x7fffffffffffffff
#define MATH_INT64_MIN  0x8000000000000000
#else //__GNUC__
#define MATH_INT64_MAX  0x7fffffffffffffffLL
#define MATH_INT64_MIN  0x8000000000000000LL
#endif

#define MATH_INT_MAX    0x7fffffff
#define MATH_INT_MIN    0x80000000

#define MATH_DOUBLE_MAX 1.7976931348623158e+308
#define MATH_FLOAT_MAX  3.402823466e+38F

#define FLOAT_EPS       1e-6
#define DOUBLE_EPS      1e-12

#ifdef __cplusplus
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
struct CMConst
{
    constexpr static const TREAL MATH_PI          = TREAL(CONST_PI);
    constexpr static const TREAL MATH_2PI         = TREAL(CONST_2PI);
    constexpr static const TREAL MATH_PI_BY_2     = TREAL(CONST_PI_BY_2);
    constexpr static const TREAL MATH_PI_BY_4     = TREAL(CONST_PI_BY_4);
    constexpr static const TREAL MATH_1_BY_PI     = TREAL(CONST_1_BY_PI);
    constexpr static const TREAL MATH_2_BY_PI     = TREAL(CONST_2_BY_PI);
    constexpr static const TREAL MATH_1_BY_2PI    = TREAL(CONST_1_BY_2PI);

    constexpr static const TREAL MATH_DEG_IN_RAD  = TREAL(CONST_DEG_IN_RAD);
    constexpr static const TREAL MATH_RAD_IN_DEG  = TREAL(CONST_RAD_IN_DEG);

    constexpr static const TREAL MATH_SQRT_2      = TREAL(CONST_SQRT_2);
    constexpr static const TREAL MATH_1_BY_SQRT_2 = TREAL(CONST_1_BY_SQRT_2);
    constexpr static const TREAL MATH_1_BY_SQRT_PI= TREAL(CONST_2_BY_SQRT_PI);
    constexpr static const TREAL MATH_SQRT_3      = TREAL(CONST_SQRT_3);

    constexpr static const TREAL MATH_E           = TREAL(CONST_E);
    constexpr static const TREAL MATH_LOG2_E      = TREAL(CONST_LOG2_E);
    constexpr static const TREAL MATH_LOG10_E     = TREAL(CONST_LOG10_E);
    constexpr static const TREAL MATH_LN_2        = TREAL(CONST_LN_2);
    constexpr static const TREAL MATH_LN_10       = TREAL(CONST_LN_10);
    constexpr static const TREAL MATH_PHI         = TREAL(CONST_PHI);
};
////////////////////////////////////////////////////////////////////////////////

/*
// For C++ ver. 98
// Don't know if ever be used.
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_PI             = TREAL(CONST_PI);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_2PI            = TREAL(CONST_2PI);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_PI_BY_2        = TREAL(CONST_PI_BY_2);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_PI_BY_4        = TREAL(CONST_PI_BY_4);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_1_BY_PI        = TREAL(CONST_1_BY_PI);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_2_BY_PI        = TREAL(CONST_2_BY_PI);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_1_BY_2PI       = TREAL(CONST_1_BY_2PI);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_DEG_IN_RAD     = TREAL(CONST_DEG_IN_RAD);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_RAD_IN_DEG     = TREAL(CONST_RAD_IN_DEG);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_SQRT_2         = TREAL(CONST_SQRT_2);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_1_BY_SQRT_2    = TREAL(CONST_1_BY_SQRT_2);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_1_BY_SQRT_PI   = TREAL(CONST_2_BY_SQRT_PI);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_SQRT_3         = TREAL(CONST_SQRT_3);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_E              = TREAL(CONST_E);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_LOG2_E         = TREAL(CONST_LOG2_E);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_LOG10_E        = TREAL(CONST_LOG10_E);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_LN_2           = TREAL(CONST_LN_2);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_LN_10          = TREAL(CONST_LN_10);
////////////////////////////////////////////////////////////////////////////////
template <typename TREAL>
const TREAL CMConst<TREAL>::MATH_PHI            = TREAL(CONST_PHI);
////////////////////////////////////////////////////////////////////////////////
*/
#endif //__cplusplus
#endif //_MATHLIB_CONSTS_
