///////////////////////////////////////////////////////////////////////////////
// primary.cpp
//
// Math primary functions definitions
///////////////////////////////////////////////////////////////////////////////

#include "primary.h"
#include <math.h>

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
// Check if value is not a number
bool misnan(float f)
{
    INTFLOAT flt;
    flt.f = f;

    uint32t exp = flt.bits.exponent;
    uint32t mant = flt.bits.mantissa;

    return (exp == 255) && (mant != 0);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is not a number
bool misnan(double d)
{
    INTDOUBLE dbl;
    dbl.d = d;

    uint64t exp = dbl.bits.exponent;
    uint64t mant = dbl.bits.mantissa;

    return (exp == 2047) && (mant != 0);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is signaling NaN
bool msnan(float f)
{
    INTFLOAT flt;
    flt.f = f;

    uint32t exp = flt.bits.exponent;
    uint32t mant = flt.bits.mantissa;

    return (exp == 255) && (mant != 0) && (mant >> 22);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is signaling NaN
bool msnan(double d)
{
    INTDOUBLE dbl;
    dbl.d = d;

    uint64t exp = dbl.bits.exponent;
    uint64t mant = dbl.bits.mantissa;

    return (exp == 2047) && (mant != 0) && (mant >> 51);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is quiet NaN
bool mqnan(float f)
{
    INTFLOAT flt;
    flt.f = f;

    uint32t exp = flt.bits.exponent;
    uint32t mant = flt.bits.mantissa;

    return (exp == 255) && (mant != 0) && !(mant >> 22);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is quiet NaN
bool mqnan(double d)
{
    INTDOUBLE dbl;
    dbl.d = d;

    uint64t exp = dbl.bits.exponent;
    uint64t mant = dbl.bits.mantissa;

    return (exp == 2047) && (mant != 0) && !(mant >> 51LL);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is infinity
bool misinf(float f)
{
    INTFLOAT flt;
    flt.f = f;

    uint32t exp = flt.bits.exponent;
    uint32t mant = flt.bits.mantissa;

    return (exp == 255) && (mant == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is infinity
bool misinf(double d)
{
    INTDOUBLE dbl;
    dbl.d = d;

    uint64t exp = dbl.bits.exponent;
    uint64t mant = dbl.bits.mantissa;

    return (exp == 2047) && (mant == 0);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is positive infinity
bool mpinf(float f)
{
    INTFLOAT flt;
    flt.f = f;

    uint32t sign = flt.bits.sign;
    uint32t exp = flt.bits.exponent;
    uint32t mant = flt.bits.mantissa;

    return (sign == 0) && (exp == 255) && (mant == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is positive infinity
bool mpinf(double d)
{
    INTDOUBLE dbl;
    dbl.d = d;

    uint64t sign = dbl.bits.sign;
    uint64t exp = dbl.bits.exponent;
    uint64t mant = dbl.bits.mantissa;

    return (sign == 0) && (exp == 2047) && (mant == 0);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is negative infinity
bool mninf(float f)
{
    INTFLOAT flt;
    flt.f = f;

    uint32t sign = flt.bits.sign;
    uint32t exp = flt.bits.exponent;
    uint32t mant = flt.bits.mantissa;

    return (sign > 0) && (exp == 255) && (mant == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is negative infinity
bool mninf(double d)
{
    INTDOUBLE dbl;
    dbl.d = d;

    uint64t sign = dbl.bits.sign;
    uint64t exp = dbl.bits.exponent;
    uint64t mant = dbl.bits.mantissa;

    return (sign > 0) && (exp == 2047) && (mant == 0);
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


#ifdef MATH_HALF_INST

extern uint16t m_float_to_half_amd(float x);
extern float m_half_to_float_amd(uint16t x);

///////////////////////////////////////////////////////////////////////////////
// Converting 32-bit integer to half type
thalf   mitoh(int32t x)
{
    float f_val = mitof(x);

    thalf out;
    uint16t * half_ptr = (uint16t *) out;
    * half_ptr = m_float_to_half_amd(f_val);

    return out;
}
///////////////////////////////////////////////////////////////////////////////
// Converting half to 32-bit integer
int32t  mhtoi(thalf x)
{
    uint16t * half_ptr = (uint16t *) x;
    float f_val = m_half_to_float_amd( * half_ptr );
    int32t i_val = mftoi(f_val);

    return i_val;
}

///////////////////////////////////////////////////////////////////////////////
// Converting float to half type
thalf   mftoh(float x)
{
    thalf out;
    uint16t * half_ptr = (uint16t *) out;
    * half_ptr = m_float_to_half_amd(x);

    return out;
}
///////////////////////////////////////////////////////////////////////////////
// Converting half to float
float   mhtof(thalf x)
{
    uint16t * half_ptr = (uint16t *) x;
    float out = m_half_to_float_amd( * half_ptr );

    return out;
}

///////////////////////////////////////////////////////////////////////////////
// Checking if value is not a number
bool misnan(thalf & x)
{
    INTHALF hlf;
    hlf.h = * (uint16t *) x;

    uint16t exp = hlf.bits.exponent;
    uint16t mant = hlf.bits.mantissa;

    return (exp == 255) && (mant != 0);

}
///////////////////////////////////////////////////////////////////////////////
// Check if value is signaling NaN
bool msnan(thalf & x)
{
    INTHALF hlf;
    hlf.h = * (uint16t *) x;

    uint16t exp = hlf.bits.exponent;
    uint16t mant = hlf.bits.mantissa;

    return (exp == 255) && (mant != 0) && !(mant >> 9);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is quiet NaN
bool mqnan(thalf & x)
{
    INTHALF hlf;
    hlf.h = * (uint16t *) x;

    uint16t exp = hlf.bits.exponent;
    uint16t mant = hlf.bits.mantissa;

    return (exp == 255) && (mant != 0) && (mant >> 9);
}

///////////////////////////////////////////////////////////////////////////////
// Check if value is infinity
bool misinf(thalf & x)
{
    INTHALF hlf;
    hlf.h = * (uint16t *) x;

    uint16t exp = hlf.bits.exponent;
    uint16t mant = hlf.bits.mantissa;

    return (exp == 31) && (mant == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is positive infinity
bool mpinf(thalf & x)
{
    INTHALF hlf;
    hlf.h = * (uint16t *) x;

    uint16t sign = hlf.bits.sign;
    uint16t exp = hlf.bits.exponent;
    uint16t mant = hlf.bits.mantissa;

    return (sign == 0) && (exp == 31) && (mant == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Check if value is negative infinity
bool mninf(thalf & x)
{
    INTHALF hlf;
    hlf.h = * (uint16t *) x;

    uint16t sign = hlf.bits.sign;
    uint16t exp = hlf.bits.exponent;
    uint16t mant = hlf.bits.mantissa;

    return (sign == 1) && (exp == 31) && (mant == 0);
}

#endif // MATH_HALF_INST




////////////////////////////////////////////////////////////////////////////////
// Tempalate functions definitions
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Calculates the floor of a value x
template <typename TReal>
TReal mfloor(TReal x)
{
    return floor(x);
}

////////////////////////////////////////////////////////////////////////////////
// Calculates the ceiling of a value x
template <typename TReal>
TReal mceil(TReal x)
{
    return mceil(x);
}

////////////////////////////////////////////////////////////////////////////////
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
// Normalizing Angle To [-PI,PI]
template <typename TReal>
TReal mnorma(TReal rad)
{
    TReal alpha = rad + CMConst<TReal>::MATH_PI;
    alpha = alpha * CMConst<TReal>::MATH_1_BY_2PI;
    alpha = mfrc(alpha);
    TReal norm_angle = alpha * CMConst<TReal>::MATH_2PI;
    norm_angle = norm_angle - CMConst<TReal>::MATH_PI;

    return norm_angle;
}

////////////////////////////////////////////////////////////////////////////////
// Inverting value sign
template <typename TReal>
TReal minvert(TReal x)
{
    return -x;
}

////////////////////////////////////////////////////////////////////////////////
// Reverse value y = 1/x
template <typename TReal>
TReal mreverse(TReal x)
{
    TReal val_one(1.f);
    return val_one/x;
}




////////////////////////////////////////////////////////////////////////////////
// Template functions specialization
////////////////////////////////////////////////////////////////////////////////

#ifdef MATH_FIXED_INST
////////////////////////////////////////////////////////////////////////////////
// Calculates the floor of a value x
tfixed32 mfloor(tfixed32 x)
{
    float flt = (float) x;
    flt = floor(flt);
    tfixed32 fx(flt);

    return fx;
}

////////////////////////////////////////////////////////////////////////////////
// Calculates the ceiling of a value x
tfixed32 mceil(tfixed32 x)
{

    float flt = (float) x;
    flt = floor(flt);
    tfixed32 fx(flt);

    return fx;
}
#endif //MATH_FIXED_INST
////////////////////////////////////////////////////////////////////////////////


#ifdef MATH_FIXED64_INST
////////////////////////////////////////////////////////////////////////////////
// Calculates the floor of a value x
tfixed64 mfloor(tfixed64 x)
{
    double dbl = (double) x;
    dbl = floor(dbl);
    tfixed64 fx(dbl);

    return fx;
}

////////////////////////////////////////////////////////////////////////////////
// Calculates the ceiling of a value x
tfixed64 mceil(tfixed64 x)
{
    float flt = (float) x;
    flt = floor(flt);
    tfixed32 fx(flt);

    return fx;
}
#endif //MATH_FIXED64_INST
////////////////////////////////////////////////////////////////////////////////

#ifdef MATH_HALF_INST
////////////////////////////////////////////////////////////////////////////////
// Calculates the floor of a value x
thalf mfloor(thalf x)
{
    float flt = (float) x;
    flt = floor(flt);
    thalf hf(flt);

    return hf;
}

////////////////////////////////////////////////////////////////////////////////
// Calculates the ceiling of a value x
thalf mceil(thalf x)
{
    float flt = (float) x;
    flt = floor(flt);
    thalf hf(flt);

    return hf;
}
#endif //MATH_HALF_INST
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
// Template functions instantiation
////////////////////////////////////////////////////////////////////////////////


#ifdef MATH_FLOAT_INST
////////////////////////////////////////////////////////////////////////////////
// Float type functions instantiation

// Calculates the floor of a value x
template float mfloor(float x);
// Calculates the ceiling of a value x
template float mceil(float x);
// Finds fraction part of number x
template float mfrc(float x);
// Normalizing angle to [-PI,PI]
template float mnorma(float rad);
// Inverting value sign
template float minvert(float x);
// Reverse value
template float mreverse(float x);
#endif //MATH_FLOAT_INST


#ifdef MATH_DOUBLE_INST
////////////////////////////////////////////////////////////////////////////////
// Double type functions instantiation

// Calculates the floor of a value x
template double mfloor(double x);
// Calculates the ceiling of a value x
template double mceil(double x);
// Finds fraction part of number x
template double mfrc(double x);
// Normalizing angle to [-PI,PI]
template double mnorma(double rad);
// Inverting value sign
template double minvert(double x);
// Reverse value
template double mreverse(double x);
#endif //MATH_DOUBLE_INST


#ifdef MATH_LONG_DOUBLE_INST
////////////////////////////////////////////////////////////////////////////////
// Long double type functions instantiation

// Calculates the floor of a value x
template long double mfloor(long double x);
// Calculates the ceiling of a value x
template long double mceil(long double x);
// Finds fraction part of number x
template long double mfrc(long double x);
// Normalizing angle to [-PI,PI]
template long double mnorma(long double rad);
// Inverting value sign
template long double minvert(long double x);
// Reverse value
template long double mreverse(long double x);
#endif //MATH_LONG_DOUBLE_INST


#ifdef MATH_FIXED_INST
////////////////////////////////////////////////////////////////////////////////
// Fixed type functions instantiation

// Calculates the floor of a value x
template tfixed32 mfloor(tfixed32 x);
// Calculates the ceiling of a value x
template tfixed32 mceil(tfixed32 x);
// Finds fraction part of number x
template tfixed32 mfrc(tfixed32 x);
// Normalizing angle to [-PI,PI]
template tfixed32 mnorma(tfixed32 rad);
// Inverting value sign
template tfixed32 minvert(tfixed32 x);
// Reverse value
template tfixed32 mreverse(tfixed32 x);
#endif //MATH_FIXED_INST


#ifdef MATH_FIXED64_INST
////////////////////////////////////////////////////////////////////////////////
// Fixed type functions instantiation

// Calculates the floor of a value x
template tfixed64 mfloor(tfixed64 x);
// Calculates the ceiling of a value x
template tfixed64 mceil(tfixed64 x);
// Finds fraction part of number x
template tfixed64 mfrc(tfixed64 x);
// Normalizing angle to [-PI,PI]
template tfixed64 mnorma(tfixed64 rad);
// Inverting value sign
template tfixed64 minvert(tfixed64 x);
// Reverse value
template tfixed64 mreverse(tfixed64 x);
#endif //MATH_FIXED64_INST


#ifdef MATH_HALF_INST
////////////////////////////////////////////////////////////////////////////////
// Fixed type functions instantiation

// Calculates the floor of a value x
template thalf mfloor(thalf x);
// Calculates the ceiling of a value x
template thalf mceil(thalf x);
// Finds fraction part of number x
template thalf mfrc(thalf x);
// Normalizing angle to [-PI,PI]
template thalf mnorma(thalf rad);
// Inverting value sign
template thalf minvert(thalf x);
// Reverse value
template thalf mreverse(thalf x);
#endif //MATH_HALF_INST


/*
-- for integer types
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
*/
/*
///////////////////////////////////////////////////////////////////////////////
// Set to zero if value is near zero (no template)
inline float   malign0(float x)
{
    if((x > -FLOAT_EPS) && (x < FLOAT_EPS))
        return 0.0f;
    else
        return x;
}
///////////////////////////////////////////////////////////////////////////////
inline double  malign0(double x)
{
    if((x > -DOUBLE_EPS) && (x < DOUBLE_EPS))
        return 0.0;
    else
        return x;
}
*/
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
