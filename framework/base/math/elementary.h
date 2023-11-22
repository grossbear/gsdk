////////////////////////////////////////////////////////////////////////////////
// elementary.h
//
// Math elementary functions declaration
////////////////////////////////////////////////////////////////////////////////

#ifndef __ELEMANTARYFUNC_H__
#define __ELEMANTARYFUNC_H__

#include "mtypes.h"

///////////////////////////////////////////////////////////////////////////////
// Library functions declarations
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Square root compute functions
///////////////////////////////////////////////////////////////////////////////

// Table square root function
float m_tsqrt(float x);

// Fast reverse square root function
float m_rfsqrt(float f);

// Approximation square root function
float m_asqrt(float x);

// 32-bit fixed point square root function
uint32t m_isqrt(uint32t x);

// 64-bit fixed point square root function
uint64t m_isqrt(uint64t x);


///////////////////////////////////////////////////////////////////////////////
// Trigonometry functions
///////////////////////////////////////////////////////////////////////////////

// Table sinus function
float m_tsinf(float angle);

// Float table cosinus function
float m_tcosf(float angle);

// Sinus and cosinus precalculated table function
void m_tsincosf(float angle, float & sin_val, float & cos_val);

// 32-bit fixed point table sinus function
tfixed32 m_tsinx(const tfixed32 & x);

// 32-bit fixed point table cosinus function
tfixed32 m_tcosx(const tfixed32 & x);

// Sinus and cosinus precalculated table function
void m_tsincosx(const tfixed32 & angle,
                 tfixed32 & sin_val, tfixed32 & cos_val);



///////////////////////////////////////////////////////////////////////////////
// Template inline functions declarations
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Power functions

// Square root of x
template <typename TReal>
TReal msqrt(TReal x);

// Base raised to the power of exponent
template <typename TReal>
TReal mpow(TReal base, TReal exp);


///////////////////////////////////////////////////////////////////////////////
// Trigonometric functions

// Compute cosine function
template <typename TReal>
TReal mcos(TReal x);

// Compute sine function
template <typename TReal>
TReal msin(TReal x);

// Compute tangent function
template <typename TReal>
TReal mtan(TReal x);

// Compute arc cosine function
template <typename TReal>
TReal macos(TReal x);

// Compute arc sine function
template <typename TReal>
TReal masin(TReal x);

// Compute arc tangent function
template <typename TReal>
TReal matan(TReal x);

// Compute arc tangent with two parameters function
template <typename TReal>
TReal matan2(TReal x, TReal y);




///////////////////////////////////////////////////////////////////////////////
// Inline template functions definitions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Power functions

///////////////////////////////////////////////////////////////////////////////
// Square root of x
template <typename TReal>
inline TReal msqrt(TReal x)
{
    return sqrt(x);
}

///////////////////////////////////////////////////////////////////////////////
// Base raised to the power of exponent
template <typename TReal>
inline TReal mpow(TReal base, TReal exp)
{
    return pow(base, exp);
}


///////////////////////////////////////////////////////////////////////////////
// Trigonometric functions

// Compute cosine function
template <typename TReal>
inline TReal mcos(TReal x)
{
    return cos(x);
}

// Compute sine function
template <typename TReal>
inline TReal msin(TReal x)
{
    return sin(x);
}

// Compute tangent function
template <typename TReal>
inline TReal mtan(TReal x)
{
    return tan(x);
}

// Compute arc cosine function
template <typename TReal>
inline TReal macos(TReal x)
{
    return acos(x);
}

// Compute arc sine function
template <typename TReal>
inline TReal masin(TReal x)
{
    return asin(x);
}

// Compute arc tangent function
template <typename TReal>
inline TReal matan(TReal x)
{
    return atan(x);
}

// Compute arc tangent with two parameters function
template <typename TReal>
inline TReal matan2(TReal x, TReal y)
{
    return atan2(x,y);
}




///////////////////////////////////////////////////////////////////////////////
// Inline template functions specialization
///////////////////////////////////////////////////////////////////////////////

#ifdef MATH_FIXED_INST
///////////////////////////////////////////////////////////////////////////////
// Power functions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Square root of x
inline tfixed32 msqrt(tfixed32 x)
{
    uint32t ui_root = m_isqrt((uint32t) x);
    tfixed32 fx_root(ui_root);
    return fx_root;
}

///////////////////////////////////////////////////////////////////////////////
// Base raised to the power of exponent
inline tfixed32 mpow(tfixed32 base, tfixed32 exp)
{
    tfixed32 fx_pow = (tfixed32) pow((float) base, (float) exp);
    return fx_pow;
}
#endif //MATH_FIXED_INST


#ifdef MATH_FIXED64_INST
///////////////////////////////////////////////////////////////////////////////
// Power functions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Square root of x
inline tfixed64 msqrt(tfixed64 x)
{
    uint64t ui_root = m_isqrt((uint64t) x);
    tfixed64 fx_root(ui_root);
    return fx_root;
}

///////////////////////////////////////////////////////////////////////////////
// Base raised to the power of exponent
inline tfixed64 mpow(tfixed64 base, tfixed64 exp)
{
    tfixed64 fx_pow = (tfixed64) pow((float) base, (float) exp);
    return fx_pow;
}
#endif //MATH_FIXED64_INST


#ifdef MATH_HALF_INST
///////////////////////////////////////////////////////////////////////////////
// Power functions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Square root of x
inline thalf msqrt(thalf x)
{
    thalf h_root = (thalf) sqrt((float) x);
    return h_root;
}

///////////////////////////////////////////////////////////////////////////////
// Base raised to the power of exponent
inline thalf mpow(thalf base, thalf exp)
{
    thalf h_pow = (thalf) pow((float) base, (float) exp);
    return h_pow;
}
#endif //MATH_HALF_INST




//--#define MATH_FLOAT_INST
//--#define MATH_DOUBLE_INST
//--#define MATH_LONG_DOUBLE_INST
//--#define MATH_FIXED_INST
//--#define MATH_FIXED64_INST
//--#define MATH_HALF_INST


//TReal mmod(TReal x, TReal y);
//TReal mpow(TReal x, TReal y);
//TReal mfmod(TReal x, TReal y);
//TReal mexp(TReal x);

//log10
//loge
//log2

//lerp
//cos_interp
//square_interp

//sin
//cos
//tan
//cot

//asin
//acos
//atan
//acot

//Trigonometric functions
/*
cos Compute cosine (function)
sin Compute sine (function)
tan Compute tangent (function)
acos    Compute arc cosine (function)
asin    Compute arc sine (function)
atan    Compute arc tangent (function)
atan2   Compute arc tangent with two parameters (function)
*/

//Hyperbolic functions
/*
cosh    Compute hyperbolic cosine (function)
sinh    Compute hyperbolic sine (function)
tanh    Compute hyperbolic tangent (function)
acosh   Compute area hyperbolic cosine (function)
asinh   Compute area hyperbolic sine (function)
atanh   Compute area hyperbolic tangent (function)
*/


//Exponential and logarithmic functions
/*
exp Compute exponential function (function)
frexp   Get significand and exponent (function)
ldexp   Generate value from significand and exponent (function)
log Compute natural logarithm (function)
log10   Compute common logarithm (function)
modf    Break into fractional and integral parts (function)
exp2    Compute binary exponential function (function)
expm1   Compute exponential minus one (function)
ilogb   Integer binary logarithm (function)
log1p   Compute logarithm plus one (function)
log2    Compute binary logarithm (function)
logb    Compute floating-point base logarithm (function)
scalbn  Scale significand using floating-point base exponent (function)
scalbln Scale significand using floating-point base exponent (long) (function)
*/

//Power functions
/*
pow Raise to power (function)
sqrt    Compute square root (function)
cbrt    Compute cubic root (function)
hypot   Compute hypotenuse (function)
*/

#endif //__ELEMANTARYFUNC_H__
