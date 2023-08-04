///////////////////////////////////////////////////////////////////////////////
// primary.cpp
//
// Math primary functions definition
///////////////////////////////////////////////////////////////////////////////

//#include <assert.h>
#include <math.h>
#include "mathdefs.h"
#include "mathconsts.h"
#include "primary.h"
#include "intreal.h"

///////////////////////////////////////////////////////////////////////////////
// Bias constant used for fast conversions between int and float. First element
// in INTORFLOAT union is integer -- We are storing biased exponent 23,
// mantissa .1, which is equivalent to 1.5 x 2^23.
static const INTFLOAT  biasflt = {((23 + 127) << 23) + (1 << 22)};
// 2^23 = 8388608

///////////////////////////////////////////////////////////////////////////////
// Converting functions from integer to float type
int32t mftoi(float f)
{
    ASSERT(f >= -4194304.0f && f <= 4194304.0f);

    INTFLOAT fi;
    fi.f = f + biasflt.f;
    return fi.i - biasflt.i;
}

///////////////////////////////////////////////////////////////////////////////
// Converting functions from float to integer
float mitof(int32t i)
{
    ASSERT(i >= -4194304 && i <= 4194304);

    INTFLOAT fi;
    fi.i = i + biasflt.i;
    return fi.f - biasflt.f;
}

///////////////////////////////////////////////////////////////////////////////
// Bias constant used for fast conversions between int and double. First
// element in INTORDOUBLE union is integer -- We are storing biased
// exponent 52, mantissa .1, which is equivalent to 1.5 x 2^52.
static const INTDOUBLE biasdbl = {(uint64t(52 + 1023) << 52) + (uint64t(1) << 51)};
// 2^52 = 4503599627370496

///////////////////////////////////////////////////////////////////////////////
// Converting functions from integer to double type
int64t mdtoi(double d)
{
    ASSERT(d >= -2251799813685248.0 && d <= 2251799813685248.0);

    INTDOUBLE di;
    di.d = d + biasdbl.d;
    return di.i - biasdbl.i;
}
///////////////////////////////////////////////////////////////////////////////
// Converting functions from double to integer type
double mitod(int64t i)
{
    ASSERT(i >= -2251799813685248 && i <= 2251799813685248);

    INTDOUBLE di;
    di.i = i + biasdbl.i;
    return di.d - biasdbl.d;
}

///////////////////////////////////////////////////////////////////////////////
// Making float number from three components
float mmakef(int32t sign, int32t exp, int32t mant)
{
    int32t i = (sign << 31) | ((exp & 255) << 23) | (mant & ((1 << 23) - 1));
    float f = *(float*)&i;

    return f;
}
///////////////////////////////////////////////////////////////////////////////
// Splitting float number to three components
void msplitf(float f, int32t &sign, int32t &exp, int32t &mant)
{
    int32t i = *(int32t*)&f;

    sign = (i >> 31) != 0 ? 1 : 0;
    exp =  (i >> 23) & 255;
    mant = i & ((1 << 23) - 1);
}

///////////////////////////////////////////////////////////////////////////////
// Making double number
double mmaked(int64t sign, int64t exp, int64t mant)
{
    int64t i = (sign << 63) | ((exp & 2047) << 52) |
                                (mant & ((int64t(1) << 52) - 1));
    double d = *(double*)&i;

    return d;
}
///////////////////////////////////////////////////////////////////////////////
// Splitting double number
void msplitd(double d, int64t &sign, int64t &exp, int64t &mant)
{
    int64t i = *(int64t*)&d;

    sign = (i >> 63) != 0 ? 1 : 0;
    exp  = (i >> 52) & 2047;
    mant = i & ((int64t(1) << 52) - 1);
}


///////////////////////////////////////////////////////////////////////////////
// Checking if value is not a number
int32t misnan(float f)
{
    int32t s,e,m;
    msplitf(f,s,e,m);

    return (e == 255) && (m != 0);
}
///////////////////////////////////////////////////////////////////////////////
// Checking if value is not a number
int64t misnan(double d)
{
    int64t s,e,m;
    msplitd(d,s,e,m);

    return (e == 2047) && (m != 0);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is signaling NaN
int32t msnan(float f)
{
    int32t s,e,m;
    msplitf(f,s,e,m);

    return (e == 255) && (m != 0) && (m >> 22);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is signaling NaN
int64t msnan(double d)
{
    int64t s,e,m;
    msplitd(d,s,e,m);

    return (e == 2047) && (m != 0) && (m >> 51);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is quiet NaN
int32t mqnan(float f)
{
    int32t s,e,m;
    msplitf(f,s,e,m);

    return (e == 255) && (m != 0) && !(m >> 22);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is quiet NaN
int64t mqnan(double d)
{
    int64t s,e,m;
    msplitd(d,s,e,m);

    return (e == 2047) && (m != 0) && !(m >> 51);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is infinity
int32t misinf(float f)
{
    int32t s,e,m;
    msplitf(f,s,e,m);

    return (e == 255) && (m == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is infinity
int64t misinf(double d)
{
    int64t s,e,m;
    msplitd(d,s,e,m);

    return (e == 2047) && (m == 0);
}


///////////////////////////////////////////////////////////////////////////////
// Check if value is positive infinity
int32t mpinf(float f)
{
    int32t s,e,m;
    msplitf(f,s,e,m);

    return (s == 0) && (e == 255) && (m == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is positive infinity
int64t mpinf(double d)
{
    int64t s,e,m;
    msplitd(d,s,e,m);

    return (s == 0) && (e == 2047) && (m == 0);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is negative infinity
int32t mninf(float f)
{
    int32t s,e,m;
    msplitf(f,s,e,m);

    return (s > 0) && (e == 255) && (m == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is negative infinity
int64t mninf(double d)
{
    int64t s,e,m;
    msplitd(d,s,e,m);

    return (s > 0) && (e == 2047) && (m == 0);
}

///////////////////////////////////////////////////////////////////////////////
// Lomont compare function
// Fast function to compare two floating point numbers
bool mlcmp(float af, float bf, int32t max_diff)
{
    int32t ai;
    int32t bi;

    ai = *(int32t *) & af;
    bi = *(int32t *) & bf;

    int32t test = (ai^bi) >> 31;

    ASSERT((test==0) || (0xffffffff == (uint32t)test));

    int32t diff = (((0x80000000 - ai) & (test)) | (ai & (~test))) - bi;
    int32t v1 = max_diff + diff;
    int32t v2 = max_diff - diff;
    return (v1 | v2) >= 0;
}

///////////////////////////////////////////////////////////////////////////////
// Calculates the floor of a value x
float   mfloor(float x)
{
    return floorf(x);
}
///////////////////////////////////////////////////////////////////////////////
double  mfloor(double x)
{
    return floor(x);
}

///////////////////////////////////////////////////////////////////////////////
// Calculates the ceiling of a value x
float   mceil(float x)
{
    return ceilf(x);
}
///////////////////////////////////////////////////////////////////////////////
double  mceil(double x)
{
    return ceil(x);
}


#ifdef MATH_LONG_DOUBLE_INST
///////////////////////////////////////////////////////////////////////////////
// Calculates the floor of a value x
long double mfloor(long double x)
{
    return floor(x);
}
///////////////////////////////////////////////////////////////////////////////
// Calculates the ceiling of a value x
long double mceil(long double x)
{
    return ceilf(x);
}
#endif //MATH_LONG_DOUBLE_INST

#ifdef MATH_FIXED_INST
///////////////////////////////////////////////////////////////////////////////
// Calculates the floor of a value x
tfixed32 mfloor(tfixed32 x)
{
    float f_val = (float)x;
    f_val = floorf(f_val);
    tfixed32 res_fixed(f_val);

    return res_fixed;
}
///////////////////////////////////////////////////////////////////////////////
// Calculates the ceiling of a value x
tfixed32   mceil(tfixed32 x)
{
    float f_val = (float)x;
    f_val = mceilf(f_val);
    tfixed32 res_fixed(f_val);

    return res_fixed;
}
#endif //MATH_FIXED_INST

#ifdef MATH_FIXED64_INST
///////////////////////////////////////////////////////////////////////////////
// Calculates the floor of a value x
tfixed64   mfloor(tfixed64 x)
{
    double f_val = (double)x;
    f_val = floor(f_val);
    tfixed64 res_fixed(f_val);

    return res_fixed;
}
///////////////////////////////////////////////////////////////////////////////
// Calculates the ceiling of a value x
tfixed64   mceil(tfixed64 x)
{
    double f_val = (double)x;
    f_val = floor(f_val);
    tfixed64 res_fixed(f_val);

    return res_fixed;
}
#endif //MATH_FIXED64_INST

#ifdef MATH_HALF_INST
///////////////////////////////////////////////////////////////////////////////
// Calculates the floor of a value x
thalf   mfloor(thalf h)
{
    float f_val = (float)h;
    f_val = floorf(f_val);
    thalf res_half(f_val);

    return res_half;
}
///////////////////////////////////////////////////////////////////////////////
// Calculates the ceiling of a value x
thalf   mceil(thalf h)
{
    float f_val = (float)h;
    f_val = ceilf(f_val);
    thalf res_half(f_val);

    return res_half;
}
#endif //MATH_HALF_INST




////////////////////////////////////////////////////////////////////////////////
// Tempalate functions definitions
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Pack values from [-1,1] to [0,1]
template <typename T>
void mpack01(marray2<T> & out, const marray2<T> & in)
{
    out.elem1 = mpack01(in.elem1);
    out.elem2 = mpack01(in.elem2);
}
////////////////////////////////////////////////////////////////////////////////
template <typename T>
void mpack01(marray3<T> & out, const marray3<T> & in)
{
    out.elem1 = mpack01(in.elem1);
    out.elem2 = mpack01(in.elem2);
    out.elem3 = mpack01(in.elem3);
}
////////////////////////////////////////////////////////////////////////////////
template <typename T>
void mpack01(marray4<T> & out, const marray4<T> & in)
{
    out.elem1 = mpack01(in.elem1);
    out.elem2 = mpack01(in.elem2);
    out.elem3 = mpack01(in.elem3);
    out.elem4 = mpack01(in.elem4);
}

////////////////////////////////////////////////////////////////////////////////
// Unpack values from [0,1] to [-1,1]
template <typename T>
void munpack01(marray2<T> & out, const marray2<T> & in)
{
    out.elem1 = munpack01( in.elem1 );
    out.elem2 = munpack01( in.elem2 );
}
////////////////////////////////////////////////////////////////////////////////
template <typename T>
void munpack01(marray3<T> & out, const marray3<T> & in)
{
    out.elem1 = munpack01( in.elem1 );
    out.elem2 = munpack01( in.elem2 );
    out.elem3 = munpack01( in.elem3 );
}
////////////////////////////////////////////////////////////////////////////////
template <typename T>
void munpack01(marray4<T> & out, const marray4<T> & in)
{
    out.elem1 = munpack01( in.elem1 );
    out.elem2 = munpack01( in.elem2 );
    out.elem3 = munpack01( in.elem3 );
    out.elem4 = munpack01( in.elem4 );
}

///////////////////////////////////////////////////////////////////////////////
// Finds fraction part of number x
template <typename TReal>
TReal mfrc(TReal x)
{
    TReal intval = mfloor(mabs(x));

    if(mless0(x))
        return intval + x;
    else
        return mabs(intval - x);
}

////////////////////////////////////////////////////////////////////////////////
// Normalizing angle to [-PI,PI]
template <typename TReal>
TReal mnorma(TReal rad)
{
    TReal alpha = rad + CMConst<TReal>::MATH_PI;
    alpha = alpha * CMConst<TReal>::MATH_1_BY_2PI;
    alpha = mfrc(alpha);

    return alpha * CMConst<TReal>::MATH_2PI - CMConst<TReal>::MATH_PI;
}

////////////////////////////////////////////////////////////////////////////////
// Check if number is a primary number
template <typename T>
bool misprim(T n)
{
    if(n <= 1) return false;
    //if( !((n>>3) & 0x0111b) )
    if(!(n>>3))
    {
        //0100 - 4
        //0110 - 6
        return (n ^ 0b100) && (n ^ 0b110);
    }

    return n%2 && n%3 && n%5 && n%7;
}




////////////////////////////////////////////////////////////////////////////////
// Functions specilization
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
// Template functions instantiation
////////////////////////////////////////////////////////////////////////////////

#ifdef MATH_FLOAT_INST
////////////////////////////////////////////////////////////////////////////////
// Float type functions instantiation

////////////////////////////////////////////////////////////////////////////////
// Pack values from [-1,1] to [0,1]
template void mpack01(marray2<float> & out, const marray2<float> & in);
template void mpack01(marray3<float> & out, const marray3<float> & in);
template void mpack01(marray4<float> & out, const marray4<float> & in);

////////////////////////////////////////////////////////////////////////////////
// Unpack values from [0,1] to [-1,1]
template void munpack01(marray2<float> & out, const marray2<float> & in);
template void munpack01(marray3<float> & out, const marray3<float> & in);
template void munpack01(marray4<float> & out, const marray4<float> & in);

///////////////////////////////////////////////////////////////////////////////
// Finds fraction part of number x
template float mfrc(float x);

////////////////////////////////////////////////////////////////////////////////
// Normalizing angle to [-PI,PI]
template float mnorma(float rad);

////////////////////////////////////////////////////////////////////////////////
#endif //MATH_FLOAT_INST


#ifdef MATH_DOUBLE_INST
////////////////////////////////////////////////////////////////////////////////
// Double type functions instantiation

////////////////////////////////////////////////////////////////////////////////
// Pack values from [-1,1] to [0,1]
template void mpack01(marray2<double> & out, const marray2<double> & in);
template void mpack01(marray3<double> & out, const marray3<double> & in);
template void mpack01(marray4<double> & out, const marray4<double> & in);

////////////////////////////////////////////////////////////////////////////////
// Unpack values from [0,1] to [-1,1]
template void munpack01(marray2<double> & out, const marray2<double> & in);
template void munpack01(marray3<double> & out, const marray3<double> & in);
template void munpack01(marray4<double> & out, const marray4<double> & in);

///////////////////////////////////////////////////////////////////////////////
// Finds fraction part of number x
template double mfrc(double x);

////////////////////////////////////////////////////////////////////////////////
// Normalizing angle to [-PI,PI]
template double mnorma(double rad);

////////////////////////////////////////////////////////////////////////////////
#endif //MATH_DOUBLE_INST


#ifdef MATH_LONG_DOUBLE_INST
////////////////////////////////////////////////////////////////////////////////
// Long double type functions instantiation

////////////////////////////////////////////////////////////////////////////////
// Pack values from [-1,1] to [0,1]
template void mpack01(marray2<long double> & out, const marray2<long double> & in);
template void mpack01(marray3<long double> & out, const marray3<long double> & in);
template void mpack01(marray4<long double> & out, const marray4<long double> & in);

////////////////////////////////////////////////////////////////////////////////
// Unpack values from [0,1] to [-1,1]
template void munpack01(marray2<long double> & out, const marray2<long double> & in);
template void munpack01(marray3<long double> & out, const marray3<long double> & in);
template void munpack01(marray4<long double> & out, const marray4<long double> & in);

///////////////////////////////////////////////////////////////////////////////
// Finds fraction part of number x
template long double mfrc(long double x);

////////////////////////////////////////////////////////////////////////////////
// Normalizing angle to [-PI,PI]
template long double mnorma(long double rad);

////////////////////////////////////////////////////////////////////////////////
#endif //MATH_LONG_DOUBLE_INST


#ifdef MATH_FIXED_INST
////////////////////////////////////////////////////////////////////////////////
// Fixed type functions instantiation

////////////////////////////////////////////////////////////////////////////////
// Pack values from [-1,1] to [0,1]
template void mpack01(marray2<tfixed32> & out, const marray2<tfixed32> & in);
template void mpack01(marray3<tfixed32> & out, const marray3<tfixed32> & in);
template void mpack01(marray4<tfixed32> & out, const marray4<tfixed32> & in);

////////////////////////////////////////////////////////////////////////////////
// Unpack values from [0,1] to [-1,1]
template void munpack01(marray2<tfixed32> & out, const marray2<tfixed32> & in);
template void munpack01(marray3<tfixed32> & out, const marray3<tfixed32> & in);
template void munpack01(marray4<tfixed32> & out, const marray4<tfixed32> & in);

///////////////////////////////////////////////////////////////////////////////
// Finds fraction part of number x
template tfixed32 mfrc(tfixed32 x);

////////////////////////////////////////////////////////////////////////////////
// Normalizing angle to [-PI,PI]
template tfixed32 mnorma(tfixed32 rad);

////////////////////////////////////////////////////////////////////////////////
#endif //MATH_FIXED_INST


#ifdef MATH_FIXED64_INST
////////////////////////////////////////////////////////////////////////////////
// Fixed type functions instantiation

////////////////////////////////////////////////////////////////////////////////
// Pack values from [-1,1] to [0,1]
template void mpack01(marray2<tfixed64> & out, const marray2<tfixed64> & in);
template void mpack01(marray3<tfixed64> & out, const marray3<tfixed64> & in);
template void mpack01(marray4<tfixed64> & out, const marray4<tfixed64> & in);

////////////////////////////////////////////////////////////////////////////////
// Unpack values from [0,1] to [-1,1]
template void munpack01(marray2<tfixed64> & out, const marray2<tfixed64> & in);
template void munpack01(marray3<tfixed64> & out, const marray3<tfixed64> & in);
template void munpack01(marray4<tfixed64> & out, const marray4<tfixed64> & in);

///////////////////////////////////////////////////////////////////////////////
// Finds fraction part of number x
template tfixed64 mfrc(tfixed64 x);

////////////////////////////////////////////////////////////////////////////////
// Normalizing angle to [-PI,PI]
template tfixed64 mnorma(tfixed64 rad);

////////////////////////////////////////////////////////////////////////////////
#endif //MATH_FIXED64_INST


#ifdef MATH_HALF_INST
////////////////////////////////////////////////////////////////////////////////
// Fixed type functions instantiation

////////////////////////////////////////////////////////////////////////////////
// Pack values from [-1,1] to [0,1]
template void mpack01(marray2<thalf> & out, const marray2<thalf> & in);
template void mpack01(marray3<thalf> & out, const marray3<thalf> & in);
template void mpack01(marray4<thalf> & out, const marray4<thalf> & in);

////////////////////////////////////////////////////////////////////////////////
// Unpack values from [0,1] to [-1,1]
template void munpack01(marray2<thalf> & out, const marray2<thalf> & in);
template void munpack01(marray3<thalf> & out, const marray3<thalf> & in);
template void munpack01(marray4<thalf> & out, const marray4<thalf> & in);

///////////////////////////////////////////////////////////////////////////////
// Finds fraction part of number x
template thalf mfrc(thalf x);

////////////////////////////////////////////////////////////////////////////////
// Normalizing angle to [-PI,PI]
template thalf mnorma(thalf rad);

////////////////////////////////////////////////////////////////////////////////
#endif //MATH_HALF_INST


////////////////////////////////////////////////////////////////////////////////
// Check if number is a primary number
template bool misprim<uint8t>(uint8t n);
template bool misprim<uint16t>(uint16t n);
template bool misprim<uint32t>(uint32t n);
template bool misprim<uint64t>(uint64t n);
////////////////////////////////////////////////////////////////////////////////





/*
///////////////////////////////////////////////////////////////////////////////
// Clamping value
float mclamp(const float & min, const float & max, const float & f)
{
    float fval = f;
    fval -= min;
    fval /= (max - min);
    fval = mclamp01(fval);
    fval *= (max - min);
    fval += min;

    ASSERT((min <= fval) && (fval <= max));

    return fval;
}

///////////////////////////////////////////////////////////////////////////////
double mclamp(const double & min, const double & max, const double & d)
{
    double dval = d;
    dval -= min;
    dval /= (max - min);
    dval = mclamp01(dval);
    dval *= (max - min);
    dval += min;

    ASSERT((min <= dval) && (dval <= max));

    return dval;
}
///////////////////////////////////////////////////////////////////////////////
*/

/*
///////////////////////////////////////////////////////////////////////////////////////
// Higher Power Of 2
///////////////////////////////////////////////////////////////////////////////////////
M_API int8t mhpow2(int8t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n + 1;

    return (n>>7) ? 0 : n;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API uint8t mhpow2(uint8t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;

    return n + 1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API int16t mhpow2(int16t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;
    n = n + 1;

    return (n>>15) ? 0 : n;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API uint16t mhpow2(uint16t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;

    return n + 1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API int32t mhpow2(int32t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;
    n = n | n>>16;
    n = n + 1;

    return (n>>31) ? 0 : n;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API uint32t mhpow2(uint32t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;
    n = n | n>>16;

    return n + 1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API int64t mhpow2(int64t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;
    n = n | n>>16;
    n = n | n>>32;
    n = n + 1;

    return (n>>63) ? 0 : n;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API uint64t mhpow2(uint64t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;
    n = n | n>>16;
    n = n | n>>32;

    return n + 1;
}

///////////////////////////////////////////////////////////////////////////////////////
// Lower Power Of 2
///////////////////////////////////////////////////////////////////////////////////////
M_API int8t mlpow2(int8t n)
{
    if(n & 0x80)
        return 0;

    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;

    return (n + 1)>>1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API uint8t mlpow2(uint8t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;

    return (n + 1)>>1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API int16t mlpow2(int16t n)
{
    if(n & 0x8000)
        return 0;

    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;

    return (n + 1)>>1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API uint16t mlpow2(uint16t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;

    return (n + 1)>>1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API int32t mlpow2(int32t n)
{
    if(n & 0x80000000)
        return 0;

    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;
    n = n | n>>16;

    return (n + 1)>>1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API uint32t mlpow2(uint32t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;
    n = n | n>>16;

    return (n + 1)>>1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API int64t mlpow2(int64t n)
{
    if(n & 0x8000000000000000)
        return 0;

    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;
    n = n | n>>16;
    n = n | n>>32;

    return (n + 1)>>1;
}
///////////////////////////////////////////////////////////////////////////////////////
M_API uint64t mlpow2(uint64t n)
{
    n = n - 1;
    n = n | n>>1;
    n = n | n>>2;
    n = n | n>>4;
    n = n | n>>8;
    n = n | n>>16;
    n = n | n>>32;

    return (n + 1)>>1;
}




///////////////////////////////////////////////////////////////////////////////////////
template bool misprim<uint8t>(uint8t n);
///////////////////////////////////////////////////////////////////////////////////////
template bool misprim<uint16t>(uint16t n);
///////////////////////////////////////////////////////////////////////////////////////
template bool miste <typename Ttype>
bool misprim(Ttype n)
{
    if(n <= 1) return false;
    //if( !((n>>3) & 0x0111b) )+
    if(!(n>>3))
    {
        //0100 - 4
        //0110 - 6
        return (n ^ 0b100) && (n ^ 0b110);
    }

    return n%2 && n%3 && n%5 && n%7;
}
prim<uint32t>(uint32t n);
///////////////////////////////////////////////////////////////////////////////////////
template bool misprim<uint64t>(uint64t n);
///////////////////////////////////////////////////////////////////////////////////////
*/
