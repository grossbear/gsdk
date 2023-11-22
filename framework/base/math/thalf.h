///////////////////////////////////////////////////////////////////////////////
//  thalf.h
//
//  Half precision floating point class declaration.
//  Conversion and utils functions declalaration available.
//
//
///////////////////////////////////////////////////////////////////////////////

#ifndef __HALF_H__
#define __HALF_H__

#include "platform.h"

///////////////////////////////////////////////////////////////////////////////
// Math 16-bit half precision floating point type class
class thalf
{
protected:
    uint16t value;

public:
    // Constructors declarations
    thalf();
    thalf(const thalf & half_val);
    thalf(float f_val);
    thalf(double d_val);
    thalf(int32t i_val);
    thalf(int16t i_val);

    // Cast to pointer uint16t *
    operator uint16t * () const;
    operator const uint16t * () const;

    // Cast to promitive built-in types
    operator float() const;
    operator double() const;
    operator int32t() const;
    operator int16t() const;
    operator bool() const;

    // Assignment operators
    thalf & operator =  (const thalf & half_val);
    thalf & operator += (const thalf & half_val);
    thalf & operator -= (const thalf & half_val);
    thalf & operator *= (const thalf & half_val);
    thalf & operator /= (const thalf & half_val);

    // Unary operators
    thalf operator + () const;
    thalf operator - () const;

    // Binary operators
    thalf operator + (const thalf & half_val) const;
    thalf operator - (const thalf & half_val) const;
    thalf operator * (const thalf & half_val) const;
    thalf operator / (const thalf & half_val) const;

    bool operator == (const thalf & half_val) const;
    bool operator != (const thalf & half_val) const;

    bool operator <  (const thalf & half_val) const;
    bool operator >  (const thalf & half_val) const;

    bool operator <= (const thalf & half_val) const;
    bool operator >= (const thalf & half_val) const;
};
///////////////////////////////////////////////////////////////////////////////

#endif //__HALF_H__
