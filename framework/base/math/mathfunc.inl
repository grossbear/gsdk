//-----------------------------------------------------------------------------
//
// File: mathfunc.inl
// Content: inline functions used in algebra classes
//
//-----------------------------------------------------------------------------

#ifndef _CMFUNC_INL_
#define _CMFUNC_INL_

///////////////////////////////////////////////////////////////////////////////
template <class T>
class CMFunc
{
public:
    static inline void MSINCOS(T rad, T & sin_val, T & cos_val);
    static inline T MACOS(T cos_val);
    static inline T MCOS(T rad);

    static inline T MSQRT(T x, T y);
    static inline T MREVSQRT(T x, T y);

    static inline T MSQRT(T x, T y, T z)
    static inline T MSQRT(T x, T y, T z);

    static inline T MSQRT(T x, T y, T z, T w);
    static inline T MREVSQRT(T x, T y, T z, T w);
};

///////////////////////////////////////////////////////////////////////////////
template <class T>
void CMFunc<T>::MSINCOS(T rad, T & sin_val, T & cos_val)
{
    sin_val = sin(rad);
    cos_val = cos(rad);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T CMFunc<T>::MACOS(T cos_val)
{
    return acos(cos_val);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T CMFunc<T>::MCOS(T rad)
{
    return cos(rad);
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
T CMFunc<T>::MSQRT(T x, T y)
{
    return sqrt(x*x + y*y);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T CMFunc<T>::MREVSQRT(T x, T y)
{
     T root = sqrt(x*x + y*y);
     T revert = T(1.0)/root;

    return revert;
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
T CMFunc<T>::MSQRT(T x, T y, T z)
{
    return sqrt(x*x + y*y + z*z);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T CMFunc<T>::MREVSQRT(T x, T y, T z)
{
     T root = sqrt(x*x + y*y + z*z);
     T revert = T(1.0)/root;

    return revert;
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
T CMFunc<T>::MSQRT(T x, T y, T z, T w)
{
    return sqrt(x*x + y*y + z*z + w*w);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
T CMFunc<T>::MREVSQRT(T x, T y, T z, T w)
{
     T root = sqrt(x*x + y*y + z*z + w*w);
     T revert = T(1.0)/root;

    return revert;
}

///////////////////////////////////////////////////////////////////////////////


//float specialization

///////////////////////////////////////////////////////////////////////////////
template <>
void CMFunc<float>::MSINCOS(float rad, float & sin_val, float & cos_val)
{
    sin_val = sin(rad);
    cos_val = cos(rad);
}
///////////////////////////////////////////////////////////////////////////////
template <>
float CMFunc<float>::MACOS(float cos_val)
{
    return acos(cos_val);
}
///////////////////////////////////////////////////////////////////////////////
template <>
float CMFunc<float>::MCOS(float rad)
{
    return cos(rad);
}

///////////////////////////////////////////////////////////////////////////////
template <>
float CMFunc<float>::MSQRT(float x, float y)
{
    return sqrt(x*x + y*y);
}
///////////////////////////////////////////////////////////////////////////////
template <>
float CMFunc<float>::MREVSQRfloat(float x, float y)
{
     float root = sqrt(x*x + y*y);
     float revert = 1.0f / root;

    return revert;
}

V///////////////////////////////////////////////////////////////////////////////
template <>
float CMFunc<float>::MSQRT(float x, float y, float z)
{
    return sqrt(x*x + y*y + z*z);
}
///////////////////////////////////////////////////////////////////////////////
template <>
float CMFunc<float>::MREVSQR(float x, float y, float z)
{
     float root = sqrt(x*x + y*y + z*z);
     float revert = 1.0f / root;

    return revert;
}

///////////////////////////////////////////////////////////////////////////////
template <>
float CMFunc<float>::MSQRT(float x, float y, float z, float w)
{
    return sqrt(x*x + y*y + z*z + w*w);
}
///////////////////////////////////////////////////////////////////////////////
template <>
float CMFunc<float>::MREVSQRT(float x, float y, float z, float w)
{
     float root = sqrt(x*x + y*y + z*z + w*w);
     float revert = 1.0f / root;

    return revert;
}

///////////////////////////////////////////////////////////////////////////////


#endif //_CMFUNC_INL_
