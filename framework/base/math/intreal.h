///////////////////////////////////////////////////////////////////////////////////////
//  intreal.h
//
//  Union contains integer and real variable
//
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _INTREAL_H_
#define _INTREAL_H_

///////////////////////////////////////////////////////////////////////////////////////
// INTORFLOAT Union For Easy Access To Bits Of A Float.
union INTFLOAT
{
    int32t  i;          // As Integer
    float   f;          // As Float
    struct              // As Bit Fields
    {
        uint32t    sign:1;
        uint32t    biasedexponent:8;
        uint32t    significand;
    }
    bits;
};

///////////////////////////////////////////////////////////////////////////////////////
// INTORDOUBLE Union For Easy Access To Bits Of A Float.
union INTDOUBLE
{
    int64t 	i;          // As 64-bit Integer
    double	d;          // As Double
    struct              // As Bit Fields
    {
        uint64t    sign:1;
        uint64t    biasedexponent:11;
        uint64t    significand;
    }
    bits;
};
///////////////////////////////////////////////////////////////////////////////////////

#endif //_INTREAL_H_
