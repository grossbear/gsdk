////////////////////////////////////////////////////////////////////////////////
//  tfixed64.h
//
//  Math 64-bit fixed data type class with 16-bit fracture
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef __TFIXED64_H__
#define __TFIXED64_H__

#ifndef __PLATFORMTYPES_H__
//typedef int int32t;
typedef long long int64t;
#endif
////////////////////////////////////////////////////////////////////////////////
class tfixed64
{
protected:
    int64t value;
    static const int64t bits = 16; //16-bit fraction part
    static const int64t one = 1 << bits;
    static const int64t half = bits ? (1 << (bits-1)) : 0;

public:
    // Constructor declarations
    tfixed64();
    tfixed64(const tfixed64 & fx_val);
    tfixed64(const float f_val);
    tfixed64(const double d_val);
    tfixed64(const int64t i_val);


    // Cast to pointer int64t *
    operator int64t * () const;
    operator const int64t * () const;

    // Cast to primitive build-in types
    operator float() const;
    operator double() const;
    operator int64t() const;
    operator bool() const;

    // Assignment operators
    tfixed64 & operator =  (const tfixed64 & fx_val);
    tfixed64 & operator += (const tfixed64 & fx_val);
    tfixed64 & operator -= (const tfixed64 & fx_val);
    tfixed64 & operator *= (const tfixed64 & fx_val);
    tfixed64 & operator /= (const tfixed64 & fx_val);

    tfixed64 & operator >>= (int64t i);
    tfixed64 & operator <<= (int64t i);

    tfixed64 & operator &= (int64t mask);
    tfixed64 & operator |= (int64t mask);
    tfixed64 & operator ^= (int64t mask);

    // Unary operators
    tfixed64 operator + () const;
    tfixed64 operator - () const;

    // Binary operators
    tfixed64 operator + (const tfixed64 & fx_val) const;
    tfixed64 operator - (const tfixed64 & fx_val) const;
    tfixed64 operator * (const tfixed64 & fx_val) const;
    tfixed64 operator / (const tfixed64 & fx_val) const;

    tfixed64 operator >> (int64t i) const;
    tfixed64 operator << (int64t i) const;

    tfixed64 operator & (int64t mask) const;
    tfixed64 operator | (int64t mask) const;
    tfixed64 operator ^ (int64t mask) const;

    bool operator == (const tfixed64 & fx_val) const;
    bool operator != (const tfixed64 & fx_val) const;

    bool operator <  (const tfixed64 & fx_val) const;
    bool operator >  (const tfixed64 & fx_val) const;

    bool operator <= (const tfixed64 & fx_val) const;
    bool operator >= (const tfixed64 & fx_val) const;

    bool operator && (const tfixed64 & fx_val) const;
    bool operator || (const tfixed64 & fx_val) const;
};
////////////////////////////////////////////////////////////////////////////////



#endif //__TFIXED64_H__
