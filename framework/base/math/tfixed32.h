////////////////////////////////////////////////////////////////////////////////
//  tfixed32.h
//
//  Math 32-bit fixed data type class with 16-bit fraction part
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef __TFIXED32_H__
#define __TFIXED32_H__

#ifndef __PLATFORMTYPES_H__
typedef int int32t;
typedef long long int64t;
#endif
////////////////////////////////////////////////////////////////////////////////

class tfixed32
{
    protected:
        int32t value;
        static const int32t bits = 16; //16-bit fraction part
        static const int32t one = 1 << bits;
        static const int32t half = bits ? (1 << (bits-1)) : 0;

    public:
        // Constructor declarations
        tfixed32();
        tfixed32(const tfixed32 & fx_val);
        tfixed32(const float f_val);
        tfixed32(const double d_val);
        tfixed32(const int32t i_val);

        // Cast to pointer int32t *
        operator int32t * () const;
        operator const int32t * () const;

        // Cast to primitive build-in types
        operator float() const;
        operator double() const;
        operator int32t() const;
        operator bool() const;

        // Assignment operators
        tfixed32 & operator =  (const tfixed32 & fx_val);
        tfixed32 & operator += (const tfixed32 & fx_val);
        tfixed32 & operator -= (const tfixed32 & fx_val);
        tfixed32 & operator *= (const tfixed32 & fx_val);
        tfixed32 & operator /= (const tfixed32 & fx_val);

        tfixed32 & operator >>= (int32t i);
        tfixed32 & operator <<= (int32t i);

        tfixed32 & operator &= (int32t mask);
        tfixed32 & operator |= (int32t mask);
        tfixed32 & operator ^= (int32t mask);


        // Unary operators
        tfixed32 operator + () const;
        tfixed32 operator - () const;

        // Binary operators
        tfixed32 operator + (const tfixed32 & fx_val) const;
        tfixed32 operator - (const tfixed32 & fx_val) const;
        tfixed32 operator * (const tfixed32 & fx_val) const;
        tfixed32 operator / (const tfixed32 & fx_val) const;

        tfixed32 operator >> (int32t i) const;
        tfixed32 operator << (int32t i) const;

        tfixed32 operator & (int32t mask) const;
        tfixed32 operator | (int32t mask) const;
        tfixed32 operator ^ (int32t mask) const;

        bool operator == (const tfixed32 & fx_val) const;
        bool operator != (const tfixed32 & fx_val) const;

        bool operator <  (const tfixed32 & fx_val) const;
        bool operator >  (const tfixed32 & fx_val) const;

        bool operator <= (const tfixed32 & fx_val) const;
        bool operator >= (const tfixed32 & fx_val) const;

        bool operator && (const tfixed32 & fx_val) const;
        bool operator || (const tfixed32 & fx_val) const;
};
////////////////////////////////////////////////////////////////////////////////


#endif //__TFIXED32_H__
