////////////////////////////////////////////////////////////////////////////////
// primary.h
//
// Math primary functions declarations
////////////////////////////////////////////////////////////////////////////////


#ifndef __PRIMARYFUNC_H__
#define __PRIMARYFUNC_H__

#include "platform.h"
#include "mathdefs.h"
#include "mathconsts.h"
#include "intreal.h"
#include "mtypes.h"


///////////////////////////////////////////////////////////////////////////////
// Functions declarations
///////////////////////////////////////////////////////////////////////////////

// Converting 32-bit integer to float
float   mitof(int32t i);
// Converting float to 32-bit integer
int32t  mftoi(float f);

// Converting 64-bit integer to double
double  mitod(int64t i);
// Converitng double to 64-bit integer
int64t  mdtoi(double d);

// Checking if value is not a number
bool misnan(float f);
// Checking if value is not a number
bool misnan(double d);

// Check if value is signaling NaN
bool msnan(float f);
// Check if value is signaling NaN
bool msnan(double d);

// Check if value is quiet NaN
bool mqnan(float f);
// Check if value is quiet NaN
bool mqnan(double d);

// Check if value is infinity
bool misinf(float f);
// Check if value is infinity
bool misinf(double d);

// Check if value is positive infinity
bool mpinf(float f);
// Check if value is positive infinity
bool mpinf(double d);

// Check if value is negative infinity
bool mninf(float f);
// Check if value is negative infinity
bool mninf(double d);

// Lomont compare function
// Fast function to compare two floating point numbers
bool mlcmp(float af, float bf, int32t max_diff);


#ifdef MATH_HALF_INST
// Converting 32-bit integer to half type
thalf   mitoh(int32t x);
// Converting half to 32-bit integer
int32t  mhtoi(thalf x);

// Converting float to half type
thalf   mftoh(float x);
// Converting half to float
float   mhtof(thalf x);

// Checking if value is not a number
bool misnan(thalf & x);
// Check if value is signaling NaN
bool msnan(thalf & x);
// Check if value is quiet NaN
bool mqnan(thalf & x);

// Check if value is infinity
bool misinf(thalf & x);
// Check if value is positive infinity
bool mpinf(thalf & x);
// Check if value is negative infinity
bool mninf(thalf & x);
///////////////////////////////////////////////////////////////////////////////

#endif // MATH_HALF_INST


///////////////////////////////////////////////////////////////////////////////
// Template functions declarations
///////////////////////////////////////////////////////////////////////////////

// Calculates the floor of a value x
template <typename TReal>
TReal mfloor(TReal x);

// Calculates the ceiling of a value x
template <typename TReal>
TReal mceil(TReal x);

// Finds fraction part of number x
template <typename TReal>
TReal mfrc(TReal x);

// Normalizing Angle To [-PI,PI]
template <typename TReal>
TReal mnorma(TReal rad);

// Inverting value sign
template <typename TReal>
TReal minvert(TReal x);

// Reverse value
template <typename TReal>
TReal mreverse(TReal x);




///////////////////////////////////////////////////////////////////////////////
// Template inline functions declaration
///////////////////////////////////////////////////////////////////////////////

// Pack values from [-1,1] to [0,1]
template <typename TReal>
TReal mpack01(TReal x);

// Unpack values from [0,1] to [-1,1]
template <typename TReal>
TReal munpack01(TReal x);

// Getting maximum value
template <typename T>
T mmax(T x, T y);

// Getting minimum value
template <typename T>
T mmin(T x, T y);

// Check if number is a primary number
template <typename T>
bool misprim(T n);

// Absolute value functions
template <typename T>
T mabs(T x);

// Negate value
template <typename T>
T mnegate(T x);

// Is value greater than 0
template <typename T>
bool mgre0(T x);

// Is value greater or equals zero
template <typename TReal>
bool mgreq0(TReal x);

// Is value equals zero
template <typename T>
bool mis0(T x);

// Is value is less or equals zero
template <typename T>
bool mlesseq0(T x);

// Is value less than 0
template <typename T>
bool mless0(T x);

// Clamp Float Value To 0
template <typename T>
T mclamp0(T x);

// Clamp value to 1
template <typename T>
T mclamp1(T x);

// Clamp value to [0,1]
template <typename T>
T mclamp01(T x);

// Clamping value between min and max values
template <typename T>
T mclamp(const T & min, const T & max, const T & val);

// Is number is power of 2
template <typename T>
bool mispow2(T n);

// Sgn
template <typename T>
T msgn(T x);




///////////////////////////////////////////////////////////////////////////////
// Inline template functions definition
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Pack values from [-1,1] to [0,1]
template <typename TReal>
inline TReal mpack01(TReal x)
{
    static const TReal half_val(0.5);
    return x * half_val + half_val;
}

///////////////////////////////////////////////////////////////////////////////
// Unpack values from [0,1] to [-1,1]
template <typename TReal>
inline TReal munpack01(TReal x)
{
    static const TReal two_val(2.0);
    static const TReal one_val(1.0);
    return x * two_val - one_val;
}

///////////////////////////////////////////////////////////////////////////////
// Absolute value functions
template <typename T>
inline T mabs(T i)
{
    return (i >> (sizeof(T) - 1)) ? ((~i) + 1) : (i);
}

///////////////////////////////////////////////////////////////////////////////
// Negate value
template <typename T>
inline T mnegate(T i)
{
    return -i;
}

///////////////////////////////////////////////////////////////////////////////
// Is value greater than 0
template <typename T>
inline bool mgre0(T i)
{
    return i > 0;
}

///////////////////////////////////////////////////////////////////////////////
// Is value greater or equals zero
template <typename T>
inline bool mgreq0(T i)
{
    return i >= 0;
}

///////////////////////////////////////////////////////////////////////////////
// Is value equals zero
template <typename T>
inline bool mis0(T i)
{
    return i == 0;
}

///////////////////////////////////////////////////////////////////////////////
// Is value is less or equals zero
template <typename T>
inline bool mlesseq0(T i)
{
    return i <= 0;
}

///////////////////////////////////////////////////////////////////////////////
// Is value less than 0
template <typename T>
inline bool mless0(T i)
{
    return i < 0;
}

///////////////////////////////////////////////////////////////////////////////
// Clamp Value To 0
template <typename T>
inline T mclamp0(T i)
{
    return (i > 0) ? (0) : (i);
}

///////////////////////////////////////////////////////////////////////////////
// Clamp value to 1
template <typename T>
inline T mclamp1(T i)
{
    return (i > 1) ? (1) : (i);
}

///////////////////////////////////////////////////////////////////////////////
// Clamp value to [0,1]
template <typename T>
inline T mclamp01(T i)
{
    return (mless0(i)) ? (0) : (mclamp1(i));
}

///////////////////////////////////////////////////////////////////////////////
// Template function clamping value between min and max values
template <typename T>
inline T  mclamp(const T & min, const T & max, const T & val)
{
    return (val < min) ? (min) : ((val > max) ? max : val);
}

///////////////////////////////////////////////////////////////////////////////
// Getting Maximum Value
template <typename T>
T mmax(T x, T y)
{
    return (x > y) ? (x) : (y);
}
///////////////////////////////////////////////////////////////////////////////
// Getting Minimum Value
template <typename T>
T mmin(T x, T y)
{
    return (x < y) ? (x) : (y);
}


///////////////////////////////////////////////////////////////////////////////
// Is number is power of 2
template <typename T>
inline bool mispow2(T n)
{
    return !(n & (n - 1));
}

///////////////////////////////////////////////////////////////////////////////
// Sgn
template <typename T>
inline T msgn(T x)
{
    return (x > 0) ? (1) : ((x < 0) ? (-1) : (0));
}




///////////////////////////////////////////////////////////////////////////////
// Template inline functions specialization
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Absolute value functions
inline float  mabs(float f)
{
    return *(float*) & (*(int32t*) &f &= MATH_INT_MAX);
}
///////////////////////////////////////////////////////////////////////////////
inline double mabs(double d)
{
    return *(double*) & (*(int64t*) &d &= MATH_INT64_MAX);
}

///////////////////////////////////////////////////////////////////////////////
// Negate value
inline float  mnegate(float f)
{
    return *(float*) & (*(int32t*) &f ^= MATH_INT_MIN);
}
///////////////////////////////////////////////////////////////////////////////
inline double mnegate(double d)
{
    return *(double*) & (*(int64t*) &d ^= MATH_INT64_MIN);
}


///////////////////////////////////////////////////////////////////////////////
// Is value near zero (no template)
inline bool  mnear0(float f)
{
    return mabs(f) < FLOAT_EPS;
}
///////////////////////////////////////////////////////////////////////////////
inline bool  mnear0(double d)
{
    return mabs(d) < DOUBLE_EPS;
}

///////////////////////////////////////////////////////////////////////////////
// Is value greater than 0
inline bool mgre0(float f)
{
    return *(int32t*) &f > 0;
}
///////////////////////////////////////////////////////////////////////////////
inline bool mgre0(double d)
{
    return *(int64t*) &d > 0;
}

///////////////////////////////////////////////////////////////////////////////
// Is value greater or equals zero
inline bool mgreq0(float f)
{
    return *(uint32t*) &f <= 0x80000000UL;
}
///////////////////////////////////////////////////////////////////////////////
inline bool mgreq0(double d)
{
    return *(uint64t*) &d <= 0x8000000000000000UL;
}

///////////////////////////////////////////////////////////////////////////////
// Is value equals zero
inline bool mis0(float f)
{
    return (*(int32t*) &f & MATH_INT_MAX) == 0;
}
///////////////////////////////////////////////////////////////////////////////
inline bool mis0(double d)
{
    return (*(int64t*) &d & MATH_INT64_MAX) == 0;
}

///////////////////////////////////////////////////////////////////////////////
// Is value is less or equals zero
inline bool mlesseq0(float f)
{
    return *(int32t*) &f <= 0;
}
///////////////////////////////////////////////////////////////////////////////
inline bool mlesseq0(double d)
{
    return *(int64t*) &d <= 0;
}

///////////////////////////////////////////////////////////////////////////////
// Is value less than 0
inline bool mless0(float f)
{
    return *(int32t*) &f < 0;
}
///////////////////////////////////////////////////////////////////////////////
inline bool mless0(double d)
{
    return *(int64t*) &d < 0;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Clamp Float Value To 0
inline float mclamp0(float f)
{
    int32t s = (*(int32t*) &f) >> 31;
    s = ~s;
    *(int32t*) &f &= s;

    return f;
}
///////////////////////////////////////////////////////////////////////////////
inline double mclamp0(double d)
{
    int64t s = (*(int64t*) &d) >> 63;
    s = ~s;
    *(int64t*) &d &= s;

    return d;
}

///////////////////////////////////////////////////////////////////////////////
// Clamp value to 1
inline float mclamp1(float f)
{
    return *(int32t*) &f > 0x3f800000 ? 1.0f : f;
}
///////////////////////////////////////////////////////////////////////////////
inline double mclamp1(double d)
{
    return *(int64t*) &d > 0x3ff0000000000000 ? 1.0 : d;
}

///////////////////////////////////////////////////////////////////////////////
// Clamp value to [0,1]
inline float mclamp01(float f)
{
    if(mless0(f))
        return 0.0f;

    return mclamp1(f);
}
///////////////////////////////////////////////////////////////////////////////
inline double mclamp01(double d)
{
    if(mless0(d))
        return 0.0;

    return mclamp1(d);
}


///////////////////////////////////////////////////////////////////////////////
// Sgn
inline float msgn(float x)
{
    return (mless0(x)) ? (-1.0f) : ((mgre0(x)) ? (1.0f) : (0.0f));
}
///////////////////////////////////////////////////////////////////////////////
inline double msgn(double x)
{
    return (mless0(x)) ? (-1.0) : ((mgre0(x)) ? (1.0) : (0.0));
}


///////////////////////////////////////////////////////////////////////////////
#ifdef  MATH_LONG_DOUBLE_INST
inline long double mabs(long double d)
{
    return (d < 0.0) ? (-d) : (d);
}
#endif

///////////////////////////////////////////////////////////////////////////////
#ifdef  MATH_FIXED_INST
inline tfixed32 mabs(tfixed32 x)
{
    return (*(int32t*)x) < 0 ? (-x) : (x);
}
#endif

///////////////////////////////////////////////////////////////////////////////
#ifdef  MATH_FIXED64_INST
inline tfixed64 mabs(tfixed64 x)
{
    return (*(int64t*)x) < 0 ? (-x) : (x);
}
#endif

///////////////////////////////////////////////////////////////////////////////
#ifdef  MATH_HALF_INST
inline thalf mabs(thalf h)
{
    return *(thalf*) & (*(uint16t*) &h &= MATH_INT16_MAX);
}

inline bool mless0(thalf h)
{
    return *(uint16t*) &h & MATH_INT16_MIN;
}
#endif


#endif //__PRIMARYFUNC_H__
