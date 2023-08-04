///////////////////////////////////////////////////////////////////////////////
//  thalf.cpp
//
//  Half precision floating point class and conversion functions
//
//
///////////////////////////////////////////////////////////////////////////////

#include "thalf.h"

//-------------------------------------------------------------------------
// Limits
//
// Visual C++ will complain if HALF_MIN, HALF_NRM_MIN etc. are not float
// constants, but at least one other compiler (gcc 2.96) produces incorrect
// results if they are.
//-------------------------------------------------------------------------
/*
#if (defined _WIN32 || defined _WIN64) && defined _MSC_VER

    #define HALF_MIN    5.96046448e-08f // Smallest positive half

    #define HALF_NRM_MIN    6.10351562e-05f // Smallest positive normalized half

    #define HALF_MAX    65504.0f    // Largest positive half

    #define HALF_EPSILON    0.00097656f // Smallest positive e for which
                            // half (1.0 + e) != half (1.0)
#else

    #define HALF_MIN    5.96046448e-08  // Smallest positive half

    #define HALF_NRM_MIN    6.10351562e-05  // Smallest positive normalized half

    #define HALF_MAX    65504.0     // Largest positive half

    #define HALF_EPSILON    0.00097656  // Smallest positive e for which
                            // half (1.0 + e) != half (1.0)
#endif


#define HALF_MANT_DIG   11      // Number of digits in mantissa
// (significand + hidden leading 1)

#define HALF_DIG    2       // Number of base 10 digits that
// can be represented without change

#define HALF_RADIX  2       // Base of the exponent

#define HALF_MIN_EXP    -13     // Minimum negative integer such that
// HALF_RADIX raised to the power of
// one less than that integer is a
// normalized half

#define HALF_MAX_EXP    16      // Maximum positive integer such that
// HALF_RADIX raised to the power of
// one less than that integer is a
// normalized half

#define HALF_MIN_10_EXP -4      // Minimum positive integer such
// that 10 raised to that power is
// a normalized half

#define HALF_MAX_10_EXP 4       // Maximum positive integer such
// that 10 raised to that power is
// a normalized half

*/
//---------------------------------------------------------------------------
//
// Implementation --
//
// Representation of a float:
//
//  We assume that a float, f, is an IEEE 754 single-precision
//  floating point number, whose bits are arranged as follows:
//
//      31 (msb)
//      |
//      | 30     23
//      | |      |
//      | |      | 22                    0 (lsb)
//      | |      | |                     |
//      X XXXXXXXX XXXXXXXXXXXXXXXXXXXXXXX
//
//      s e        m
//
//  S is the sign-bit, e is the exponent and m is the significand.
//
//  If e is between 1 and 254, f is a normalized number:
//
//              s    e-127
//      f = (-1)  * 2      * 1.m
//
//  If e is 0, and m is not zero, f is a denormalized number:
//
//              s    -126
//      f = (-1)  * 2      * 0.m
//
//  If e and m are both zero, f is zero:
//
//      f = 0.0
//
//  If e is 255, f is an "infinity" or "not a number" (NAN),
//  depending on whether m is zero or not.
//
//  Examples:
//
//      0 00000000 00000000000000000000000 = 0.0
//      0 01111110 00000000000000000000000 = 0.5
//      0 01111111 00000000000000000000000 = 1.0
//      0 10000000 00000000000000000000000 = 2.0
//      0 10000000 10000000000000000000000 = 3.0
//      1 10000101 11110000010000000000000 = -124.0625
//      0 11111111 00000000000000000000000 = +infinity
//      1 11111111 00000000000000000000000 = -infinity
//      0 11111111 10000000000000000000000 = NAN
//      1 11111111 11111111111111111111111 = NAN
//
// Representation of a half:
//
//  Here is the bit-layout for a half number, h:
//
//      15 (msb)
//      |
//      | 14  10
//      | |   |
//      | |   | 9        0 (lsb)
//      | |   | |        |
//      X XXXXX XXXXXXXXXX
//
//      s e     m
//
//  S is the sign-bit, e is the exponent and m is the significand.
//
//  If e is between 1 and 30, h is a normalized number:
//
//              s    e-15
//      h = (-1)  * 2     * 1.m
//
//  If e is 0, and m is not zero, h is a denormalized number:
//
//              S    -14
//      h = (-1)  * 2     * 0.m
//
//  If e and m are both zero, h is zero:
//
//      h = 0.0
//
//  If e is 31, h is an "infinity" or "not a number" (NAN),
//  depending on whether m is zero or not.
//
//  Examples:
//
//      0 00000 0000000000 = 0.0
//      0 01110 0000000000 = 0.5
//      0 01111 0000000000 = 1.0
//      0 10000 0000000000 = 2.0
//      0 10000 1000000000 = 3.0
//      1 10101 1111000001 = -124.0625
//      0 11111 0000000000 = +infinity
//      1 11111 0000000000 = -infinity
//      0 11111 1000000000 = NAN
//      1 11111 1111111111 = NAN
//
// Conversion:
//
//  Converting from a float to a half requires some non-trivial bit
//  manipulations.  In some cases, this makes conversion relatively
//  slow, but the most common case is accelerated via table lookups.
//
//  Converting back from a half to a float is easier because we don't
//  have to do any rounding.  In addition, there are only 65536
//  different half numbers; we can convert each of those numbers once
//  and store the results in a table.  Later, all conversions can be
//  done using only simple table lookups.
//
//---------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// Float to half conversion
// Including zeros, denormalized numbers and exponent overflows
uint16t m_float_to_half_ilm(float f)
{

    ///////////////////////////////////////////////////////////////////////////
    // Our floating point number, f, is represented by the bit
    // pattern in integer i.  Disassemble that bit pattern into
    // the sign, s, the exponent, e, and the significand, m.
    // Shift s into the position where it will go in in the
    // resulting half number.
    // Adjust e, accounting for the different exponent bias
    // of float and half (127 versus 15).
    ///////////////////////////////////////////////////////////////////////////

    int32t i = *(int32t*)&f;

    register int32t  s =  (i >> 16)  & 0x00008000;
    register int32t  e =  ((i >> 23) & 0x000000ff) - (127 - 15);
    register int32t  m =   i         & 0x007fffff;


    // Now reassemble s, e and m into a half:
    if (e <= 0)
    {
        if (e < -10)
        {
            ///////////////////////////////////////////////////////////////////
            // E is less than -10.  The absolute value of f is
            // less than HALF_MIN (f may be a small normalized
            // float, a denormalized float or a zero).
            //
            // We convert f to a half zero with the same sign as f.
            ///////////////////////////////////////////////////////////////////

            return s;
        }

        ///////////////////////////////////////////////////////////////////////
        // E is between -10 and 0.  F is a normalized float
        // whose magnitude is less than HALF_NRM_MIN.
        //
        // We convert f to a denormalized half.
        ///////////////////////////////////////////////////////////////////////

        m = (m | 0x00800000) >> (1 - e);

        ///////////////////////////////////////////////////////////////////////
        // Round to nearest, round "0.5" up.
        //
        // Rounding may cause the significand to overflow and make
        // our number normalized.  Because of the way a half's bits
        // are laid out, we don't have to treat this case separately;
        // the code below will handle it correctly.
        ///////////////////////////////////////////////////////////////////////

        if (m &  0x00001000)
            m += 0x00002000;

        ///////////////////////////////////////////////////////////////////////
        // Assemble the half from s, e (zero) and m.
        ///////////////////////////////////////////////////////////////////////

        return s | (m >> 13);
    }
    else if (e == 0xff - (127 - 15))
    {
        if (m == 0)
        {
            ///////////////////////////////////////////////////////////////////
            // F is an infinity; convert f to a half
            // infinity with the same sign as f.
            ///////////////////////////////////////////////////////////////////

            return s | 0x7c00;
        }
        else
        {
            ///////////////////////////////////////////////////////////////////
            // F is a NAN; we produce a half NAN that preserves
            // the sign bit and the 10 leftmost bits of the
            // significand of f, with one exception: If the 10
            // leftmost bits are all zero, the NAN would turn
            // into an infinity, so we have to set at least one
            // bit in the significand.
            ///////////////////////////////////////////////////////////////////

            m >>= 13;
            return s | 0x7c00 | m | (m == 0);
        }
    }
    else
    {
        ///////////////////////////////////////////////////////////////////////
        // E is greater than zero.  F is a normalized float.
        // We try to convert f to a normalized half.
        ///////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////
        // Round to nearest, round "0.5" up
        ///////////////////////////////////////////////////////////////////////

        if (m &  0x00001000)
        {
            m += 0x00002000;

            if (m & 0x00800000)
            {
                m =  0;     // overflow in significand,
                e += 1;     // adjust exponent
            }
        }

        ///////////////////////////////////////////////////////////////////////
        // Handle exponent overflow
        ///////////////////////////////////////////////////////////////////////

        if (e > 30)
        {
            //--overflow ();        // Cause a hardware floating point overflow;

            return s | 0x7c00;  // if this returns, the half becomes an
        }                       // infinity with the same sign as f.

        ///////////////////////////////////////////////////////////////////////
        // Assemble the half from s, e and m.
        ///////////////////////////////////////////////////////////////////////

        return s | (e << 10) | (m >> 13);
    }
}

///////////////////////////////////////////////////////////////////////////////
// Conversion function from ATI SDK
uint16t m_float_to_half_amd(float f)
{
    uint32t x = *(uint32t *)&f;
    uint32t sign = (uint16t)(x >> 31);
    uint32t mantissa;
    int32t exp;
    uint16t half;

    //Get mantissa
    mantissa = x & ((1 << 23) - 1);

    //Get exponent bits
    exp =  (int32t)((x >> 23) & 0xFF) - 127;

    if (exp > 16)
    {
        // Largest possible half float number
        exp = 16;
        mantissa = (1 << 23) - 1;
    }
    else if (exp <= -15)
    {
        // All float denorm and 0 values map to half 0
        mantissa = 0;
        exp = -15;
    }

    half = (((uint16t)sign) << 15 )  | (((uint16t)(exp + 15 )) << 10) |
             (uint16t)(mantissa >> 13);

    return half;
}

///////////////////////////////////////////////////////////////////////////////
// Converting half type to float type functions
float m_half_to_float_ilm(uint16t h_val)
{
    int32t sign = (h_val >> 15) & 0x00000001;
    int32t exp  = (h_val >> 10) & 0x0000001f;
    int32t mant =  h_val        & 0x000003ff;

    if (exp == 0)
    {
        if (mant == 0)
        {
            // Plus or minus zero
            sign = sign << 31;
            return *(float*)&sign; //(s << 31)
        }
        else
        {
            // Denormalized number -- renormalize it
            while (!(mant & 0x00000400))
            {
                mant <<= 1;
                exp  -=  1;
            }

            exp  += 1;
            mant &= ~0x00000400;
        }
    }
    else if (exp == 31)
    {
        if (mant == 0)
        {
            // Positive or negative infinity
            sign = ((sign << 31) | 0x7f800000);
            return *(float*)&sign; //((s << 31) | 0x7f800000)
        }
        else
        {
            // Nan -- preserve sign and significand bits
            sign = ((sign << 31) | 0x7f800000 | (mant << 13));
            return *(float*)&sign;
        }
    }

    // Normalized number
    exp = exp + (127 - 15);
    mant = mant << 13;

    // Assemble s, e and m.
    int32t i_val = (sign << 31) | (exp << 23) | mant;
    float f_val = *(float*) & i_val;

    return f_val;
}

///////////////////////////////////////////////////////////////////////////////
float m_half_to_float_amd(uint16t h_val)
{
    uint32t sign = (uint32t)(h_val >> 15);
    uint32t mant = (uint32t)(h_val & ((1 << 10) - 1));
    int32t  exp  = (int32t)((h_val >> 10) & 0x1F) - 15;
    uint32t i_val;

    if (exp == 16)
    {
        // Convert 16-bit FP inf/NaN To 32-bit inf/NaN Value
        exp = 128;
    }
    else if (exp == -15)
    {
        // Convert 16-bit FP zero/denorm To 32-bit zero/denorm Value
        exp = -127;
    }
    i_val = (sign << 31) | ((exp + 127) << 23) | (mant << 13);

    return *(float *) & i_val;
}

///////////////////////////////////////////////////////////////////////////////
float m_half_to_float_table(uint16t h_val)
{
    //0xAABBCCDD, 0xAABBCCDD, 0xAABBCCDD, 0xAABBCCDD, 0xAABBCCDD, 0xAABBCCDD,
    return 0.0f;
}


//////////////////////////////////////////////////////////////////////////////
inline uint16t m_float_to_half(float f_val)
{
    return m_float_to_half_amd(f_val);
}

//////////////////////////////////////////////////////////////////////////////
inline float m_half_to_float(uint16t h_val)
{
    return m_half_to_float_amd(h_val);
}


//////////////////////////////////////////////////////////////////////////////
// Conversion functions
thalf mftoh(float f_val)
{
    thalf out;
    uint16t * ptr = (uint16t *) out;
    *ptr = m_float_to_half(f_val);
}

///////////////////////////////////////////////////////////////////////////////
float mhtof(const thalf & h_val)
{
    uint16t * ptr = (uint16t *) h_val;

    return m_half_to_float(*ptr);
}

///////////////////////////////////////////////////////////////////////////////
uint16t mftoh_raw(float f_val)
{
    return m_float_to_half(f_val);
}

///////////////////////////////////////////////////////////////////////////////
float  mhtof_raw(uint16t h_val)
{
    return m_half_to_float(h_val);
}

///////////////////////////////////////////////////////////////////////////////
// Making half number
thalf mmakeh(int32t sign, int32t exp, int32t mant)
{
    uint16t val = (sign << 15) | ((exp & 31) << 10) | (mant & ((1 << 10) - 1));

    return val;
}
///////////////////////////////////////////////////////////////////////////////
// Splitting half number
void msplith(const thalf & h, int32t & sign, int32t & exp, int32t & mant)
{
    sign = (h >> 15) != 0 ? 1 : 0;
    exp =  (h >> 10) & 31;
    mant = h & ((1 << 10) - 1);
}

///////////////////////////////////////////////////////////////////////////////
// Testing value functions
///////////////////////////////////////////////////////////////////////////////
// Test if not a number
bool misnan(thalf & h)
{
    int32t s,e,m;
    msplith(h,s,e,m);

    return (e == 255) && (m != 0);
}
///////////////////////////////////////////////////////////////////////////////
// If quiet NaN
bool mqnan(thalf & h)
{
    int32t s,e,m;
    msplith(h,s,e,m);

    return (e == 255) && (m != 0) && !(m >> 9);
}
///////////////////////////////////////////////////////////////////////////////
// If signaling NaN
bool msnan(thalf & h)
{
    int32t s,e,m;
    msplith(h,s,e,m);

    return (e == 255) && (m != 0) && (m >> 9);
}

///////////////////////////////////////////////////////////////////////////////
// Is infinity
bool misinf(thalf & h)
{
    int32t s,e,m;
    msplith(h,s,e,m);

    return (e == 31) && (m == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Is positive infinity
bool mpinf(thalf & h)
{
    int32t s,e,m;
    msplith(h,s,e,m);

    return (s == 0) && (e == 31) && (m == 0);
}
///////////////////////////////////////////////////////////////////////////////
// Is negative infinity
bool mninf(thalf & h)
{
    int32t s,e,m;
    msplith(h,s,e,m);

    return (s == 1) && (e == 31) && (m == 0);
}

///////////////////////////////////////////////////////////////////////////////
// Math 16-bit half precision floating point type class
///////////////////////////////////////////////////////////////////////////////
// Class constructors
///////////////////////////////////////////////////////////////////////////////
thalf::thalf():value(0){}
///////////////////////////////////////////////////////////////////////////////
thalf::thalf(const thalf & half_val):value(half_val.value){}
///////////////////////////////////////////////////////////////////////////////
thalf::thalf(float f_val)
{
    value = m_float_to_half(f_val);
}
///////////////////////////////////////////////////////////////////////////////
thalf::thalf(double d_val)
{
    value = m_float_to_half((float)d_val);
}
///////////////////////////////////////////////////////////////////////////////
thalf::thalf(int32t i_val)
{
    value = m_float_to_half((float)i_val);
}
///////////////////////////////////////////////////////////////////////////////
thalf::thalf(int16t i_val)
{
    value = m_float_to_half((float)i_val);
}

///////////////////////////////////////////////////////////////////////////////
// Cast to pointer uint16t *
thalf::operator uint16t* () const
{
    return (uint16t*) this;
}
///////////////////////////////////////////////////////////////////////////////
thalf::operator const uint16t* () const
{
    return (const uint16t*) this;
}

///////////////////////////////////////////////////////////////////////////////
// Cast to promitive built-in types
thalf::operator float() const
{
    return m_half_to_float(value);
}
///////////////////////////////////////////////////////////////////////////////
thalf::operator double() const
{
    return (double)m_half_to_float(value);
}
///////////////////////////////////////////////////////////////////////////////
thalf::operator int32t() const
{
    return (int32t)m_half_to_float(value);
}
///////////////////////////////////////////////////////////////////////////////
thalf::operator int16t() const
{
    return (int16t)m_half_to_float(value);
}
///////////////////////////////////////////////////////////////////////////////
thalf::operator bool() const
{
    return value != 0;
}


///////////////////////////////////////////////////////////////////////////////
// Assignment operators
thalf & thalf::operator =  (const thalf & half_val)
{
    value = half_val.value;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
thalf & thalf::operator += (const thalf & half_val)
{
    float lvalue = (float) *this;
    float rvalue = (float) half_val;
    float result = lvalue + rvalue;
    value = m_float_to_half(result);

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
thalf & thalf::operator -= (const thalf & half_val)
{
    float lvalue = (float) *this;
    float rvalue = (float) half_val;
    float result = lvalue - rvalue;
    value = m_float_to_half(result);

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
thalf & thalf::operator *= (const thalf & half_val)
{
    float lvalue = (float) *this;
    float rvalue = (float) half_val;
    float result = lvalue * rvalue;
    value = m_float_to_half(result);

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
thalf & thalf::operator /= (const thalf & half_val)
{
    float lvalue = (float) *this;
    float rvalue = (float) half_val;
    float result = lvalue / rvalue;
    value = m_float_to_half(result);

    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// Unary operators
thalf thalf::operator + () const
{
    thalf out;
    out.value = value;

    return out;
}
///////////////////////////////////////////////////////////////////////////////
thalf thalf::operator - () const
{
    thalf out;
    out.value = value;

    return out;
}

///////////////////////////////////////////////////////////////////////////////
// Binary operators
thalf thalf::operator + (const thalf & half_val) const
{
    thalf out;
    float lvalue = (float) *this;
    float rvalue = (float) half_val;
    float result = lvalue + rvalue;
    out.value = m_float_to_half(result);

    return out;
}
///////////////////////////////////////////////////////////////////////////////
thalf thalf::operator - (const thalf & half_val) const
{
    thalf out;
    float lvalue = (float) *this;
    float rvalue = (float) half_val;
    float result = lvalue - rvalue;
    out.value = m_float_to_half(result);

    return out;
}
///////////////////////////////////////////////////////////////////////////////
thalf thalf::operator * (const thalf & half_val) const
{
    thalf out;
    float lvalue = (float) *this;
    float rvalue = (float) half_val;
    float result = lvalue * rvalue;
    out.value = m_float_to_half(result);

    return out;
}
///////////////////////////////////////////////////////////////////////////////
thalf thalf::operator / (const thalf & half_val) const
{
    thalf out;
    float lvalue = (float) *this;
    float rvalue = (float) half_val;
    float result = lvalue / rvalue;
    out.value = m_float_to_half(result);

    return out;
}
///////////////////////////////////////////////////////////////////////////////
bool thalf::operator == (const thalf & half_val) const
{
    return value == half_val.value;
}
///////////////////////////////////////////////////////////////////////////////
bool thalf::operator != (const thalf & half_val) const
{
    return value != half_val.value;
}
///////////////////////////////////////////////////////////////////////////////
bool thalf::operator <  (const thalf & half_val) const
{
    float lvalue = (float) *this;
    float rvalue = (float) half_val;

    return lvalue < rvalue;
}
///////////////////////////////////////////////////////////////////////////////
bool thalf::operator >  (const thalf & half_val) const
{
    float lvalue = (float) *this;
    float rvalue = (float) half_val;

    return lvalue > rvalue;
}

///////////////////////////////////////////////////////////////////////////////
bool thalf::operator <= (const thalf & half_val) const
{
    if (value == half_val.value)
        return true;

    return ((float) *this) < ((float) half_val);
}
///////////////////////////////////////////////////////////////////////////////
bool thalf::operator >= (const thalf & half_val) const
{
    if (value == half_val.value)
        return true;

    return ((float) *this) > ((float) half_val);
}

///////////////////////////////////////////////////////////////////////////////
