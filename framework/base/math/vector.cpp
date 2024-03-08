//-----------------------------------------------------------------------------
// vector.cpp
//
// Methods operating on vector classes
//
//-----------------------------------------------------------------------------

#include "vector.h"
#include "mathval.inl"
#include "mathfunc.inl"

#ifdef _CMVECTOR2D_
///////////////////////////////////////////////////////////////////////////////
// Vector 2D methods

///////////////////////////////////////////////////////////////////////////////
// Rotate 2D vector
template <class T>
void CMVec2Rotate(  CMVector2D<T> & vOut,
                    const CMVector2D<T> & vIn, T rad)
{
    T cosA;
    T sinA;

    // Calculating sinus and cosinus angle values
    CMFunc<T>::MSINCOS(rad, sinA, cosA);

    T x =  vIn.x * cosA + vIn.y * sinA;
    T y = -vIn.x * sinA + vIn.x * cosA;

    vOut.x = x;
    vOut.y = y;

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Returns the length of a 2D vector
template <class T>
T CMVec2Length(const CMVector2D<T> & v)
{
    return CMVFunc<T>::MSQRT(v.x, v.y);
}

///////////////////////////////////////////////////////////////////////////////
// Returns the square of the length of a 2D vector
template <class T>
T CMVec2LengthSq(const CMVector2D<T> & v)
{
    return v.x * v.x + v.y * v.y;
}

///////////////////////////////////////////////////////////////////////////////
// Normalizing 2D vector
template <class T>
void CMVec2Normalize(CMVector2D<T> & vOut, const CMVector2D<T> & vIn)
{
    T vec[2] = {vIn.x, vIn.y};
    T revsq = CMFunc<T>::MREVSQRT(vIn.x, vIn.y);

    vec[0] *= revsq;
    vec[1] *= revsq;

    vOut.x = vec[0];
    vOut.y = vec[1];

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Calculate length between two 2D points
template <class T>
T CMVec2Distance(const CMVector2D<T> & p0, const CMVector2D<T> & p1)
{
    CMVector2D<T> vec;
    CMVec2Sub(vec, p1, p0);

    return CMFunc<T>::MSQRT(vec.x, vec.y);
}

///////////////////////////////////////////////////////////////////////////////
// Calculate square length between two 2D points
template <class T>
T CMVec2DistanceSq(const CMVector2D<T> & p0, const CMVector2D<T> & p1)
{
    CMVector2D<T> vec;
    CMVec2Sub(vec, p1, p0);

    return vec.x*vec.x + vec.y*vec.y;
}

///////////////////////////////////////////////////////////////////////////////
// Getting angle in radians between two 2D vectors
template <class T>
T CMVec2Angle(const CMVector2D<T> & v1, const CMVector2D<T> & v2)
{
    const T len1 = CMVec2Length(v1);
    const T len2 = CMVec2Length(v2);

    ASSERT(!mnear0(len1));
    ASSERT(!mnear0(len2));
    ASSERT(!mnear0(len1*len2));

    return CMFunc<T>::MACOS(CMVec2Dot(v1,v2) / (len1*len2));
}

///////////////////////////////////////////////////////////////////////////////
// Getting reflected 2D vector
template <class T>
void CMVec2Reflect( CMVector2D<T> & vR,
                    const CMVector2D<T> & vL,
                    const CMVector2D<T> & vN)
{
    //R = 2*N*(NoL) - L

    T dot = CMVec2Dot(vN, vL);
    T val2 = T(2.0);

    CMVector2D<T> res(vN);
    res *= val2;
    res *= dot;

    CMVec2Sub(vR, res, vL);
}

///////////////////////////////////////////////////////////////////////////////
// Getting refracted 2D vector
template <class T>
void CMVec2Refract( CMVector2D<T> &vR,
                    const CMVector2D<T> &vL,
                    const CMVector2D<T> &vN, T etaRatio)
{
    T val1 = T(1.0);

    T cosI = CMVec2Dot(vL,vN);
    T cosT2 = val1 - etaRatio*etaRatio*(val1 - cosI * cosI);

    CMVector2D<T> T1 = vL * (-etaRatio);
    CMVector2D<T> T2 = vN * (etaRatio * cosI - CMFunc<T>::MSQRT(mabs(cosT2)));

    CMVector2D<T> res;
    CMVec2Add(res, T1, T2);

    vR = res * (T) mgre0(cosT2);
}

///////////////////////////////////////////////////////////////////////////////
// Performs a square interpolation, using two 2D vectors
template <class T>
void CMVec2Sqrp(CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2, T factor)
{
    T f2 = factor * factor;

    T x = v1.x + ((v2.x - v1.x) * f2);
    T y = v1.y + ((v2.y - v1.y) * f2);

    vOut.x = x;
    vOut.y = y;

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a cosinus interpolation, using two 2D vectors
template <class T>
void CMVec2Cosrp(CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2, T factor)
{
    const T pi = CMConst<T>::MATH_PI;
    const T angle = factor * pi;
    const T c05 = CMValues<T>::HALF;
    const T c1  = CMValues<T>::ONE;

    T prc = (c1 - CMFunc<T>::MCOS(angle)) * c05;

    T x = v1.x + ((v2.x - v1.x) * prc);
    T y = v1.y + ((v2.y - v1.y) * prc);

    vOut.x = x;
    vOut.y = y;

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a point in barycentric coordinates, using the specified 2D vectors
template <class T>
void CMVec2BaryCentric( CMVector2D<T> & vOut,
                        const CMVector2D<T> & v1,
                        const CMVector2D<T> & v2,
                        const CMVector2D<T> & v3,
                        T f, T g)
{
    const T c1 = CMValues<T>::ONE;
    vOut.x = (c1 - f - g) * (v1.x) + f * (v2.x) + g * (v3.x);
    vOut.y = (c1 - f - g) * (v1.y) + f * (v2.y) + g * (v3.y);

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a quadratic interpolation, using the specified 2D vectors
template <class T>
void CMVec2Quadratic(CMVector2D<T> & vOut,
                    const CMVector2D<T> & v1,
                    const CMVector2D<T> & v2,
                    const CMVector2D<T> & v3, T factor)
{
    const T c1 = CMValues<T>::ONE;
    const T c2 = CMValues<T>::TWO;

    T f2 = factor * factor;
    T val = c1 - factor;
    T val2 = val * val;
    T facval = factor * val * c2;

    //pOut = (*pV1)*val*val + ((*pV2)*factor*val)*T(2) + (*pV3)*f2;+
    vOut.x = v1.x * val2 + v2.x * facval + v3.x * f2;
    vOut.y = v1.y * val2 + v2.y * facval + v3.y * f2;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a cubic interpolation, using the specified 2D vectors
template <class T>
void CMVec2Cubic(CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2,
                const CMVector2D<T> & v3,
                const CMVector2D<T> & v4, T factor)
{
    T x2 = factor * factor;
    T f = factor;

    CMVector2D<T> Vec12(v1 - v2);
    CMVector2D<T> P((v4 - v3) - Vec12);
    CMVector2D<T> Q(Vec12 - P);
    CMVector2D<T> R(v3 - v1);

    vOut = P*x2*f + Q*x2 + R*f + v2;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a catmull-rom spline interpolation, using the specified 2D vectors
// calculate value between points v2 and v3
template <class T>
void CMVec2CatmullRom(  CMVector2D<T> & vOut,
                        const CMVector2D<T> & v1,
                        const CMVector2D<T> & v2,
                        const CMVector2D<T> & v3,
                        const CMVector2D<T> & v4, T x)
{
    const T val05 = CMValues<T>::HALF;
    const T val2  = CMVAlues<T>::TWO;
    const T val3  = CMValues<T>::THREE;
    const T val4  = CMValues<T>::FOUR;
    const T val5  = CMValues<T>::FIVE;

    CMVector2D<T> C1( v2 * val2 );
    CMVector2D<T> C2( (v1 * (-val1) ) + v3);
    CMVector2D<T> C3( (v1 * val2) + (v2 * (-val5)) + (v3 * val4) + (v4 * (-val1)) );
    CMVector2D<T> C4( (v1 * (-val1)) + (v2 * val3) + (v3 * (-val3)) + (v4) );

    vOut = ( (((C4*x + C3)*x + C2)*x + C1) ) * val05;

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a hermite spline interpolation, using the specified 2D vectors
//
// hermite interpolate
// calculate value between points v2 and v3
// tension: 1 is high, 0 normal, -1 is low
// bias: 0 is even,
//         positive is towards first segment,
//         negative towards the other
template <class T>
void CMVec2Hermite( CMVector2D<T> & vOut,
                    const CMVector2D<T> & v1,
                    const CMVector2D<T> & v2,
                    const CMVector2D<T> & v3,
                    const CMVector2D<T> & v4, T x, T tension, T bias)
{
    T x2 = x  * x;
    T x3 = x2 * x;

    CMVector2D<T>   m1 = (v2 - v1)*(T(1.0) + bias)*(T(1.0) - tension)*T(0.5);
                        m1 += (v3 - v2)*(T(1.0) - bias)*(T(1.0) - tension)*T(0.5);
    CMVector2D<T>   m2 = (v3 - v2)*(T(1.0) + bias)*(T(1.0) - tension)*T(0.5);
                        m2 += (v4 - v3)*(T(1.0) - bias)*(T(1.0) - tension)*T(0.5);

    T a1 = x3*T(2.0) - x2*T(3.0) + T(1.0);
    T a2 = x3 - x2*T(2.0) + x;
    T a3 = x3 - x2;
    T a4 = x3*T(-2.0) + x2*T(3.0);

    vOut = v2*a1 + m1*a2 + m2*a3 + v3*a4;
}
///////////////////////////////////////////////////////////////////////////////
#endif //_CMVECTOR2D_

#ifdef _CMVECTOR3D_
///////////////////////////////////////////////////////////////////////////////
// Vector 3D methods

//////////////////////////////////////////////////////////////////////////////
// Rotate a 3D vector in X axis, using the specified angle in radians
template <class T>
void CMVec3RotateX(CMVector3D<T> & vOut, const CMVector3D<T> & vIn, T rad)
{
    T cosA;
    T sinA;

    // Calculating sinus and cosinus angle values
    CMFunc<T>::MSINCOS(rad, sinA, cosA);

    T y = vIn.y * cosA - vIn.z * sinA;//y*cosAngle - z*sinAngle
    T z = vIn.y * sinA + vIn.z * cosA;//y*sinAngle + z*cosAngle

    vOut.x = vIn.x;
    vOut.y = y;
    vOut.z = z;
}


///////////////////////////////////////////////////////////////////////////////
// Rotate a 3D vector in Y axis, using the specified angle in radians
template <class T>
void CMVec3RotateY(CMVector3D<T> & vOut, const CMVector3D<T> & vIn, T rad)
{
    T cosA;
    T sinA;

    // Calculating sinus & cosinus angle values
    CMFunc<T>::MSINCOS(rad, sinA, cosA);

    T x =  vIn.x * cosA + vIn.z * sinA;
    T z = -vIn.x * sinA + vIn.z * cosA;

    vOut.x = x;
    vOut.y = vIn.y;
    vOut.z = z;
}

///////////////////////////////////////////////////////////////////////////////
// Rotate a 3D vector in Z axis, using the specified angle in radians
template <class T>
void CMVec3RotateZ(CMVector3D<T> & vOut, const CMVector3D<T> & vIn, T rad)
{
    T cosA;
    T sinA;

    // Calculating sinus & cosinus angle values
    CMFunc<T>::MSINCOS(rad, sinA, cosA);

    T x = vIn.x * cosA - vIn.y * sinA;
    T y = vIn.x * sinA + vIn.y * cosA;

    vOut.x = x;
    vOut.y = y;
    vOut.z = vIn.z;
}


///////////////////////////////////////////////////////////////////////////////
// Returns the length of a 3D vector
template <class T>
T CMVec3Length(const CMVector3D<T> & v)
{
    return CMFunc<T>::MSQRT(v.x, v.y, v.z);
}

///////////////////////////////////////////////////////////////////////////////
// Returns the square of the length of a 3D vector
template <class T>
T CMVec3LengthSq(const CMVector3D<T> & v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

///////////////////////////////////////////////////////////////////////////////
// Calculate length between two 3D points
template <class T>
T CMVec3Distance(const CMVector3D<T> & p0, const CMVector3D<T> & p1)
{
    CMVector3D<T> vec;
    CMVec3Sub(vec, p1, p0);

    return MVEC3SQRT(vec.x, vec.y, vec.z);
}

///////////////////////////////////////////////////////////////////////////////
// Calculate square length between two 3D points
template <class T>
T CMVec3DistanceSq(const CMVector3D<T> & p0, const CMVector3D<T> & p1)
{
    CMVector3D<T> vec;
    CMVec3Sub(vec, p1, p0);

    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}

///////////////////////////////////////////////////////////////////////////////
// Normalize a 3D vector
template <class T>
void CMVec3Normalize(CMVector3D<T> & vOut, const CMVector3D<T> & vIn)
{
    T vec[3] = {vIn.x, vIn.y, vIn.z};
    T revsq = CMFunc<T>::MREVSQRT(vIn.x, vIn.y, vIn.z);

    vec[0] *= revsq;
    vec[1] *= revsq;
    vec[2] *= revsq;

    memcpy(vOut, vec, sizeof(T) * 3);

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Getting angle in radians between two 3D vectors
template <class T>
T CMVec3Angle(const CMVector3D<T> & v1, const CMVector3D<T> & v2)
{
    T len1 = CMVec3Length(v1);
    T len2 = CMVec3Length(v2);

    ASSERT(!mnear0(len1));
    ASSERT(!mnear0(len2));
    ASSERT(!mnear0(len1*len2));

    return CMFunc<T>::MACOS(CMVec3Dot(v1,v2) / (len1*len2));
}


///////////////////////////////////////////////////////////////////////////////
// Getting reflected 3D vector
template <class T>
void CMVec3Reflect(CMVector3D<T> & vR, const CMVector3D<T> & vL, const CMVector3D<T> & vN)
{
    //R = 2*N*(NoL) - L

    T dot = CMVec3Dot(vN, vL);
    T two = T(2);

    CMVector3D<T> res = vN;
    res *= two;
    res *= dot;

    CMVec3Sub(vR, res, vL);
}

///////////////////////////////////////////////////////////////////////////////
// Getting refracted 3D vector
template <class T>
void CMVec3Refract( CMVector3D<T> & vR,
                    const CMVector3D<T> & vL,
                    const CMVector3D<T> & vN, T etaRatio)
{
    T one = T(1);

    T cosI = CMVec3Dot(vL, vN);
    T cosT2 = one - etaRatio * etaRatio * (one - cosI*cosI);

    CMVector3D<T> T1 = vL * (-etaRatio);
    CMVector3D<T> T2 = vN * (etaRatio * cosI - CMFunc<T>::MSQRT(mabs(cosT2)));

    CMVector3D<T> res;
    CMVec3Add(res, T1, T2);

    vR = res * (T)mgre0(cosT2);
}

///////////////////////////////////////////////////////////////////////////////
// Performs a catmull-rom interpolation, using the specified 3D vectors
template <class T>
void CMVec3CatmullRom(  CMVector3D<T> & vOut,
                        const CMVector3D<T> & v1,
                        const CMVector3D<T> & v2,
                        const CMVector3D<T> & v3,
                        const CMVector3D<T> & v4,  T f)
{
    const T c05 = CMValues<T>::HALF;
    const T c1  = CMValues<T>::ONE;
    const T c2  = CMValues<T>::TWO;
    const T c3  = CMValues<T>::THREE;
    const T c4  = CMValues<T>::FOUR;
    const T c5  = CMValues<T>::FIVE;

    CMVector3D<T> C1, C2, C3, C4;

    C1 = v2 * c2;
    C2 = (v1 * (-c1)) + v3;
    C3 = (v1 * c2) + (v2 * (-c5)) + (v3 * c4) + (v4 * (-c1));
    C4 = (v1 * (-c1)) + (v2 * c3) + (v3 * (-c3)) + v4;
    //CMVector3D<T> C1( v2 * T(2));
    //CMVector3D<T> C2( (v1 * T(-1)) + v3);
    //CMVector3D<T> C3( (v1 * T(2)) + (v2 * T(-5)) + (v3 * T(4)) + (v4 * T(-1)) );
    //CMVector3D<T> C4( (v1 * T(-1)) + (v2 * T(3)) + (v3 * T(-3)) + (v4) );

    //vOut = ( (((C4*f + C3)*f + C2)*f + C1) ) * T(0.5);
    vOut = (((C4 * f + C3) * f + C2) * f + C1) * c05;
}

///////////////////////////////////////////////////////////////////////////////
// Returns a point in barycentric coordinates, using the specified 3D vectors
template <class T>
void CMVec3BaryCentric( CMVector3D<T> & vOut,
                        const CMVector3D<T> & v1,
                        const CMVector3D<T> & v2,
                        const CMVector3D<T> & v3, T f, T g)
{
    const T c1 = CMValues<T>::ONE;
    vOut.x = (c1 - f - g) * (v1.x) + f * (v2.x) + g * (v3.x);
    vOut.y = (c1 - f - g) * (v1.y) + f * (v2.y) + g * (v3.y);
    vOut.z = (c1 - f - g) * (v1.z) + f * (v2.z) + g * (v3.z);

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a hermite spline interpolation, using the specified 3D vectors
template <class T>
void CMVec3Hermite( CMVector3D<T> & vOut,
                    const CMVector3D<T> & v1,
                    const CMVector3D<T> & v2,
                    const CMVector3D<T> & v3,
                    const CMVector3D<T> & v4, T x, T tension, T bias)
{
    T x2 = x*x;
    T x3 = x2*x;
    //T one = T(1);
    //T halfone = T(0.5);

    CMVector3D<T>   m1 = (v2 - v1)*(T(1) + bias)*(T(1) - tension)*T(0.5);
                    m1 += (v3 - v2)*(T(1) - bias)*(T(1) - tension)*T(0.5);
    CMVector3D<T>   m2 = (v3 - v2)*(T(1) + bias)*(T(1) - tension)*T(0.5);
                    m2 += (v4 - v3)*(T(1) - bias)*(T(1) - tension)*T(0.5);

    T a1 = x3*T(2) - x2*T(3) + T(1);
    T a2 = x3 - x2*T(2) + x;
    T a3 = x3 - x2;
    T a4 = x3*T(-2) + x2*T(3);

    vOut = v2*a1 + m1*a2 + m2*a3 + v3*a4;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a quadratic interpolation, using the specified 3D vectors
template <class T>
void CMVec3Quadratic(CMVector3D<T> & vOut,
                    const CMVector3D<T> & v1,
                    const CMVector3D<T> & v2,
                    const CMVector3D<T> & v3, T factor)
{
    const c1 = CMValues<T>::ONE;
    const c2 = CMValues<T>::TWO;

    T f2 = factor * factor;
    T val = c1 - factor;
    T val2 = val * val;
    T facval = factor * val * c2;


    vOut.x = v1.x*val2 + v2.x*facval + v3.x*f2;
    vOut.y = v1.y*val2 + v2.y*facval + v3.y*f2;
    vOut.z = v1.z*val2 + v2.y*facval + v3.z*f2;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a cubic interpolation, using the specified 3D vectors
template <class T>
void CMVec3Cubic(CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2,
                const CMVector3D<T> & v3,
                const CMVector3D<T> & v4, T factor)
{
    T x2 = factor*factor;
    T f = factor;

    CMVector3D<T> Vec12(v1 - v2);
    CMVector3D<T> P((v4 - v3) - Vec12);
    CMVector3D<T> Q(Vec12 - P);
    CMVector3D<T> R(v3 - v1);

    vOut = P*x2*f + Q*x2 + R*f + v2;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a square interpolation, using two 3D vectors
template <class T>
void CMVec3Sqrp(CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2, T factor)
{
    T f2 = factor*factor;

    T x = v1.x + ((v2.x - v1.x)*f2);
    T y = v1.y + ((v2.y - v1.y)*f2);
    T z = v1.z + ((v2.z - v1.z)*f2);

    vOut.x = x;
    vOut.y = y;
    vOut.z = z;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a cosinus interpolation, using two 3D vectors
template <class T>
void CMVec3Cosrp(CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2, T factor)
{
    T angle = factor * CMConst<T>::MATH_PI;
    T prc = (T(1) - CMFunc<T>::MCOS(angle))*T(0.5);

    T x = v1.x + ((v2.x - v1.x) * prc);
    T y = v1.y + ((v2.y - v1.y) * prc);
    T z = v1.z + ((v2.z - v1.z) * prc);

    vOut.x = x;
    vOut.y = y;
    vOut.z = z;
}
#endif //_CMVECTOR3D_

#ifdef _CMVECTOR4D_
///////////////////////////////////////////////////////////////////////////////
// 4D vector

// DO POPRAKI - Cross Product in 4 dimentions
///////////////////////////////////////////////////////////////////////////////
// Cross product of two 4D vectors
template <class T>
void CMVec4Cross(CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2)
{
    vOut.x = v1.y*v2.z - v1.z*v2.y;
    vOut.y = v1.z*v2.x - v1.x*v2.z;
    vOut.z = v1.x*v2.y - v1.y*v2.x;
}

///////////////////////////////////////////////////////////////////////////////
// Rotate a 4D vector in X axis, using the specified angle in radians
template <class T>
void CMVec4RotateX(CMVector4D<T> & vOut, const CMVector4D<T> & vIn, T rad)
{
    T cosA;
    T sinA;

    // Calculating Sinus And Cosinus Angle Values
    CMFunc<T>::MSINCOS(rad, sinA, cosA);

    T y = vIn.y*cosA - vIn.z*sinA;//y*cosAngle - z*sinAngle
    T z = vIn.y*sinA + vIn.z*cosA;//y*sinAngle + z*cosAngle

    vOut.x = vIn.x;
    vOut.y = y;
    vOut.z = z;
    vOut.w = vIn.w;
}


///////////////////////////////////////////////////////////////////////////////
// Rotate a 4D vector in Y axis, using the specified angle in radians
template <class T>
void CMVec4RotateY(CMVector4D<T> & vOut, const CMVector4D<T> & vIn, T rad)
{
    T cosA;
    T sinA;

    // Calculating Sinus And Cosinus Angle Values
    CMFunc<T>::MSINCOS(rad, sinA, cosA);

    T x = vIn.x*cosA + vIn.z*sinA;//x*cosAngle + z*sinAngle
    T z = -vIn.x*sinA + vIn.z*cosA;//-x*sinAngle + z*cosAngle

    vOut.x = x;
    vOut.y = vIn.y;
    vOut.z = z;
    vOut.w = vIn.w;
}

///////////////////////////////////////////////////////////////////////////////
// Rotate a 4D vector in Z axis, using the specified angle in radians
template <class T>
void CMVec4RotateZ(CMVector4D<T> & vOut, const CMVector4D<T> & vIn, T rad)
{
    T cosA;
    T sinA;

    // Calculating sinus & cosinus angle values
    CMFunc<T>::MSINCOS(rad, sinA, cosA);

    T x = vIn.x*cosA - vIn.y*sinA;//x*cosAngle - y*sinAngle
    T y = vIn.x*sinA + vIn.y*cosA;//x*sinAngle + y*cosAngle

    vOut.x = x;
    vOut.y = y;
    vOut.z = vIn.z;
    vOut.w = vIn.w;
}

///////////////////////////////////////////////////////////////////////////////
// Returns the length of a 4D vector
template <class T>
T CMVec4Length(const CMVector4D<T> & v)
{
    return CMFunc<T>::MSQRT(v.x, v.y, v.z, v.w);
}

///////////////////////////////////////////////////////////////////////////////
// Returns the square of the length of a 4D vector
template <class T>
T CMVec4LengthSq(const CMVector4D<T> & v)
{
    return v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w;
}

///////////////////////////////////////////////////////////////////////////////
// Calculate length between two 4D points
template <class T>
T CMVec4Distance(CMVector4D<T> & p0, CMVector4D<T> & p1)
{
    CMVector4D<T> vec;
    CMVec4Sub(vec, p1, p0);

    return MVEC4SQRT(vec.x, vec.y, vec.z, vec.w);
}

///////////////////////////////////////////////////////////////////////////////
// Calculate square length between two 4D points
template <class T>
T CMVec4DistanceSq(CMVector4D<T> & p0, CMVector4D<T> & p1)
{
    CMVector4D<T> vec;
    CMVec4Sub(vec,p1,p0);

    return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z + vec.w*vec.w;
}

///////////////////////////////////////////////////////////////////////////////
// Normalize a 4D vector
template <class T>
void CMVec4Normalize(CMVector4D<T> & vOut, const CMVector4D<T> & vIn)
{
    T vec[4] = {vIn.x, vIn.y, vIn.z, vIn.w};
    T revsq = CMFunc<T>::MREVSQRT(vIn.x, vIn.y, vIn.z, vIn.w);

    vec[0] *= revsq;
    vec[1] *= revsq;
    vec[2] *= revsq;
    vec[3] *= revsq;

    memcpy(vOut, vec, sizeof(T)*4);
}

///////////////////////////////////////////////////////////////////////////////
// Getting angle in radians between two 4D vectors
template <class T>
T CMVec4Angle(const CMVector4D<T> & v1, const CMVector4D<T> & v2)
{
    T len1 = CMVec4Length(v1);
    T len2 = CMVec4Length(v2);

    ASSERT(!mnear0(len1));
    ASSERT(!mnear0(len2));
    ASSERT(!mnear0(len1*len2));

    return CMFunc<T>::MACOS(CMVec4Dot(v1,v2) / (len1*len2));
}

///////////////////////////////////////////////////////////////////////////////
// Getting reflected 4D vector
template <class T>
void CMVec4Reflect( CMVector4D<T> & vR,
                    const CMVector4D<T> & vL,
                    const CMVector4D<T> & vN)
{
    //R = 2*N*(NoL) - L

    T dot = CMVec4Dot(vN, vL);
    T two = T(2);

    CMVector4D<T> res = vN;
    res *= two;
    res *= dot;

    CMVec4Sub(vR, res, vL);
}

///////////////////////////////////////////////////////////////////////////////
// Getting refracted 4D vector
template <class T>
void CMVec4Refract( CMVector4D<T> & vR,
                    const CMVector4D<T> & vL,
                    const CMVector4D<T> & vN, T etaRatio)
{
    T one = T(1);

    T cosI = CMVec4Dot(vL,vN);
    T cosT2 = one - etaRatio*etaRatio * (one - cosI*cosI);

    CMVector4D<T> T1 = vL * (-etaRatio);
    CMVector4D<T> T2 = vN * (etaRatio * cosI - CMFunc<T>::MSQRT(mabs(cosT2)));

    CMVector4D<T> res;
    CMVec4Add(res, T1, T2);

    vR = res * (T)mgre0(cosT2);
}

///////////////////////////////////////////////////////////////////////////////
// Performs a catmull-rom interpolation, using the specified 4D vectors
template <class T>
void CMVec4CatmullRom(  CMVector4D<T> & vOut,
                        const CMVector4D<T> & v1,
                        const CMVector4D<T> & v2,
                        const CMVector4D<T> & v3,
                        const CMVector4D<T> & v4,  T f)
{
    CMVector4D<T> C1( v2 * T(2));
    CMVector4D<T> C2( (v1 * T(-1)) + v3);
    CMVector4D<T> C3( (v1 * T(2)) + (v2 * T(-5)) + (v3 * T(4)) + (v4 * T(-1)) );
    CMVector4D<T> C4( (v1 * T(-1)) + (v2 * T(3)) + (v3 * T(-3)) + (v4) );

    vOut = ( (((C4*f + C3)*f + C2)*f + C1) ) * T(0.5);
}

///////////////////////////////////////////////////////////////////////////////
// Returns a point in barycentric coordinates, using the specified 4D vectors
template <class T>
void CMVec4BaryCentric( CMVector4D<T> & vOut,
                        const CMVector4D<T> & v1,
                        const CMVector4D<T> & v2,
                        const CMVector4D<T> & v3, T f, T g)
{
    vOut.x = (T(1)-f-g) * (v1.x) + f * (v2.x) + g * (v3.x);
    vOut.y = (T(1)-f-g) * (v1.y) + f * (v2.y) + g * (v3.y);
    vOut.z = (T(1)-f-g) * (v1.z) + f * (v2.z) + g * (v3.z);
    vOut.w = (T(1)-f-g) * (v1.w) + f * (v2.w) + g * (v3.w);
}

///////////////////////////////////////////////////////////////////////////////
// Performs a hermite spline interpolation, using the specified 4D vectors
template <class T>
void CMVec4Hermite( CMVector4D<T> & vOut,
                    const CMVector4D<T> & v1,
                    const CMVector4D<T> & v2,
                    const CMVector4D<T> & v3,
                    const CMVector4D<T> & v4, T x, T tension, T bias)
{
    T x2 = x*x;
    T x3 = x2*x;

    CMVector4D<T>   m1 = (v2 - v1)*(T(1) + bias)*(T(1) - tension)*T(0.5);
                    m1 += (v3 - v2)*(T(1) - bias)*(T(1) - tension)*T(0.5);
    CMVector4D<T>   m2 = (v3 - v2)*(T(1) + bias)*(T(1) - tension)*T(0.5);
                    m2 += (v4 - v3)*(T(1) - bias)*(T(1) - tension)*T(0.5);

    T a1 = x3*T(2) - x2*T(3) + T(1);
    T a2 = x3 - x2*T(2) + x;
    T a3 = x3 - x2;
    T a4 = x3*T(-2) + x2*T(3);

    vOut = v2*a1 + m1*a2 + m2*a3 + v3*a4;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a quadratic interpolation, using the specified 4D vectors
template <class T>
void CMVec4Quadratic(CMVector4D<T> & vOut,
                    const CMVector4D<T> & v1,
                    const CMVector4D<T> & v2,
                    const CMVector4D<T> & v3, T factor)
{
    T f2 = factor * factor;
    T val = T(1) - factor;
    T val2 = val * val;
    T facval = factor * val * T(2);

    vOut.x = v1.x*val2 + v2.x*facval + v3.x*f2;
    vOut.y = v1.y*val2 + v2.y*facval + v3.y*f2;
    vOut.z = v1.z*val2 + v2.y*facval + v3.z*f2;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a cubic interpolation, using the specified 4D vectors
template <class T>
void CMVec4Cubic(CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2,
                const CMVector4D<T> & v3,
                const CMVector4D<T> & v4, T factor)
{
    T x2 = factor * factor;
    T f = factor;

    CMVector4D<T> Vec12(v1 - v2);
    CMVector4D<T> P((v4 - v3) - Vec12);
    CMVector4D<T> Q(Vec12 - P);
    CMVector4D<T> R(v3 - v1);

    vOut = P*x2*f + Q*x2 + R*f + v2;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a square interpolation, using two 4D vectors
template <class T>
void CMVec4Sqrp(CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2, T factor)
{
    T f2 = factor * factor;

    T x = v1.x + ((v2.x - v1.x)*f2);
    T y = v1.y + ((v2.y - v1.y)*f2);
    T z = v1.z + ((v2.z - v1.z)*f2);
    T w = v1.w + ((v2.w - v1.w)*f2);

    vOut.x = x;
    vOut.y = y;
    vOut.z = z;
    vOut.w = w;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a cosinus interpolation, using two 4D vectors
template <class T>
void CMVec4Cosrp(CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2, T factor)
{
    T angle = factor * CMConst<T>::MATH_PI;
    T prc = (T(1) - MVEC4COS(angle))*T(0.5);

    T x = v1.x + ((v2.x - v1.x)*prc);
    T y = v1.y + ((v2.y - v1.y)*prc);
    T z = v1.z + ((v2.z - v1.z)*prc);
    T w = v1.w + ((v2.w - v1.w)*prc);

    vOut.x = x;
    vOut.y = y;
    vOut.z = z;
    vOut.w = w;
}
///////////////////////////////////////////////////////////////////////////////
#endif //_CMVECTOR4D_
