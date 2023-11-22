///////////////////////////////////////////////////////////////////////////////////////
//  intreal.h
//
//  Union contains integer and real variable
//
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _INTREAL_H_
#define _INTREAL_H_

///////////////////////////////////////////////////////////////////////////////////////
// INTORFLOAT Union for easy access to bits of a float.
union INTFLOAT
{
    int32t  i;          // As Integer
    float   f;          // As Float
    struct              // As Bit Fields
    {
        uint32t    sign:1;
        uint32t    exponent:8;
        uint32t    mantissa:23;
    }
    bits;
};

///////////////////////////////////////////////////////////////////////////////////////
// INTORDOUBLE Union for easy access to bits of a double.
union INTDOUBLE
{
    int64t  i;          // As 64-bit Integer
    double  d;          // As Double
    struct              // As Bit Fields
    {
        uint64t    sign:1;
        uint64t    exponent:11;
        uint64t    mantissa:52;
    }
    bits;
};
///////////////////////////////////////////////////////////////////////////////////////
// INTHALF Union for easy access to bits of a half.
union INTHALF
{
    uint16t h;          // as half
    struct              // as bits field
    {
        uint16t    sign:1;
        uint16t    exponent:5;
        uint16t    mantissa:10;
    }
    bits;
};
///////////////////////////////////////////////////////////////////////////////////////

#endif //_INTREAL_H_
