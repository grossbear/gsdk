///////////////////////////////////////////////////////////////////////////////////////
//  quaternion.h
//
//  Quaternion class template declaration
//
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _QUATERNION_H_
#define _QUATERNION_H_


///////////////////////////////////////////////////////////////////////////////////////
// Quaternion Class Declaration
template <class T>
class CMQuaternion
{
public:
    union{
        struct{
            T x,y,z,w;
        };
        T q[4];
    };

public:
    // Constructors
    CMQuaternion();
    CMQuaternion(CMQuaternion &q);
    CMQuaternion(CMVector3D<T> &v, T s);
    CMQuaternion(const T* const p);
    CMQuaternion(T x, T y, T z, T w);

    // Cast To Pointer  (T *)
    operator T* () const;
    operator const T* () const;

    // Cast To CMVector3D Type
    operator  CMVector3D<T>();

    // Cast To CMVector4D Type
    operator  CMVector4D<T>();

    // Access Grants
    T operator [] (int32t );
    T operator [] (int32t ) const;

    // Assignment Operators
    CMQuaternion& operator =  ( const CMQuaternion & );
    CMQuaternion& operator += ( const CMQuaternion & );
    CMQuaternion& operator -= ( const CMQuaternion & );
    CMQuaternion& operator *= ( const CMQuaternion & );
    CMQuaternion& operator *= ( const T & );
    CMQuaternion& operator /= ( const T & );

    // Unary Operators
    CMQuaternion operator + () const;
    CMQuaternion operator - () const;

    // Binary Operators
    CMQuaternion operator + ( const CMQuaternion & ) const;
    CMQuaternion operator - ( const CMQuaternion & ) const;
    CMQuaternion operator * ( const CMQuaternion & ) const;
    CMQuaternion operator * ( const T & ) const;
    CMQuaternion operator / ( const T & ) const;

    boolt operator == ( const CMQuaternion & ) const;
    boolt operator != ( const CMQuaternion & ) const;
};

//---------------------------
// Inline
//---------------------------

////////////////////////////////////////////////////////////////////////////
// Adding Two Quaternion
template <class T>
void CMQuatAdd(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2);

////////////////////////////////////////////////////////////////////////////
// Subtracting Two Quaternions
template <class T>
void CMQuatSub(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2);

////////////////////////////////////////////////////////////////////////////
// Set Identity Quaternion
template <class T>
void CMQuatIdentity(CMQuaternion<T> &Q);

////////////////////////////////////////////////////////////////////////////
// Is Identity
template <class T>
boolt CMQuatIsIdentity(CMQuaternion<T> &Q);

////////////////////////////////////////////////////////////////////////////
// Quaternion Dot Product
template <class T>
T CMQuatDot(const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2);

////////////////////////////////////////////////////////////////////////////
// Conjugate Quaternion
template <class T>
void CMQuatConjugate(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q);

//---------------------------
// Non-Inline
//---------------------------

////////////////////////////////////////////////////////////////////////////
// Multiplying Two Quaternions
template <class T>
void CMQuatMult(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2);

////////////////////////////////////////////////////////////////////////////
// Multiplying 3D Vector And Quaternion Treating 3D Vector As Quaternion
// With Zero Const Part
template <class T>
void CMQuatMult(CMQuaternion<T> &qOut, const CMVector3D<T> &V1, const CMQuaternion<T> &Q1);

////////////////////////////////////////////////////////////////////////////
// Multiplying Quaternions And 3D Vector Treating 3D Vector As Quaternion
// With Zero Const Part
template <class T>
void CMQuatMult(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q1, const CMVector3D<T> &V1);

////////////////////////////////////////////////////////////////////////////
// Getting Rotating Angle
template <class T>
T CMQuatAngle(const CMQuaternion<T> &Q);

////////////////////////////////////////////////////////////////////////////
// Getting Normalised Vector Specified By Vector Part Of Quaternion
template <class T>
void CMQuatAxis(CMVector3D<T> &V, const CMQuaternion<T> &Q);

////////////////////////////////////////////////////////////////////////////
// Returns The Length Of Quaternion
template <class T>
T CMQuatLength(const CMQuaternion<T> &Q);

////////////////////////////////////////////////////////////////////////////
// Returns The Square Of The Length Of Quaternion
template <class T>
T CMQuatLengthSq(const CMQuaternion<T> &Q);

////////////////////////////////////////////////////////////////////////////
// Normalize A Quaternion
template <class T>
void CMQuatNormalize(CMQuaternion<T> &qOut, const CMQuaternion<T> &qIn);

////////////////////////////////////////////////////////////////////////////
// Rotate Quaternion Using Another Quaternion
template <class T>
void CMQuatRotate(CMQuaternion<T> &Q, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2);

////////////////////////////////////////////////////////////////////////////
// Rotating 3D Vector Using Quaternion
template <class T>
void CMQuatRotate(CMVector3D<T> &vOut, const CMQuaternion<T> &Q, const CMVector3D<T> &V);

////////////////////////////////////////////////////////////////////////////
// Create Quaternion For Rotating From Vector V1 To Vector V2
// Vectors Must Be Normalized
template <class T>
void CMQuatRotateArc(CMQuaternion<T> &Q, const CMVector3D<T> &V1, const CMVector3D<T> &V2);

////////////////////////////////////////////////////////////////////////////
// Bulding Quaternion From Three Euler Angles
template <class T>
void CMQuatFromAxisAngles(CMQuaternion<T> &qOut, T x, T y, T z);

////////////////////////////////////////////////////////////////////////////
// Getting Euler Angles From Quaternion
template <class T>
void CMQuatToAxisAngles(CMVector3D<T> &vOut, const CMQuaternion<T> &qIn);

////////////////////////////////////////////////////////////////////////////
// Bulding Quaternion From Matix
template <class T>
void CMMtxToQuat(CMQuaternion<T> &Q, const CMMatrix44<T> &M);
////////////////////////////////////////////////////////////////////////////
template <class T>
void CMMtxToQuat(CMQuaternion<T> &Q, const CMMatrix33<T> &M);

////////////////////////////////////////////////////////////////////////////
// Building Matrix From Quaternion
template <class T>
void CMQuatToMtx(CMMatrix44<T> &M, const CMQuaternion<T> &Q);
////////////////////////////////////////////////////////////////////////////
template <class T>
void CMQuatToMtx(CMMatrix33<T> &M, const CMQuaternion<T> &Q);

////////////////////////////////////////////////////////////////////////////
// Logarithm Of A Quaternion
template <class T>
void CMQuatLog(CMQuaternion<T> &qOut, const CMQuaternion<T> &qIn);

////////////////////////////////////////////////////////////////////////////
// Exponent Of A Quaternion
template <class T>
void CMQuatExp(CMQuaternion<T> &qOut, const CMQuaternion<T> &qIn);

////////////////////////////////////////////////////////////////////////////
// Linear Interpolation Between Two Quaternions
template <class T>
void CMQuatLerp(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2, T weight);

////////////////////////////////////////////////////////////////////////////
// Spherical Interpolation Between Two Quaternions
template <class T>
void CMQuatSLerp(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2, T weight);

////////////////////////////////////////////////////////////////////////////
// Spherical Cubic Interpolation
template <class T>
void CMQuatSQuad(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2,
                 const CMQuaternion<T> &A, const CMQuaternion<T> &B, T weight);

////////////////////////////////////////////////////////////////////////////
// Spline Interpolation
template <class T>
void CMQuatSpline(CMQuaternion<T> &Q, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Qm, const CMQuaternion<T> &Q2);

///////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//
//  File:       Quaternion.inl
//  Content:    Quaternion inline functions
//
//////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////
// Quaternion Structure Constructors

// Default Constructor
template <class T>
M_INLINE CMQuaternion<T>::CMQuaternion()
{
    x = y = z = 0;
    w = T(1);
}

// Copy Constructor
template <class T>
M_INLINE CMQuaternion<T>::CMQuaternion(CMQuaternion &q)
{
    x = q.x;
    y = q.y;
    z = q.z;
    w = q.w;
}

// Building Quaternion From 3D Vector And Scalar Value
template <class T>
M_INLINE CMQuaternion<T>::CMQuaternion(CMVector3D<T> &v, T s)
{
    x = v.x;
    y = v.y;
    z = v.z;
    w = s;
}

// Getting Pointer To Array
template <class T>
M_INLINE CMQuaternion<T>::CMQuaternion(const T* const p)
{
    ASSERT(p != NULL);

    x = p[0];
    y = p[1];
    z = p[2];
    w = p[3];
}

// Getting Values
template <class T>
M_INLINE CMQuaternion<T>::CMQuaternion(T x, T y, T z, T c)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

////////////////////////////////////////////////////////////////////////////
// Cast To Pointer  (T *)
template <class T>
M_INLINE CMQuaternion<T>::operator T* () const
{
    return (T*) this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T>::operator const T* () const
{
    return (const T*) this;
}

////////////////////////////////////////////////////////////////////////////
// Cast To CMVector3D Type
template <class T>
M_INLINE CMQuaternion<T>::operator CMVector3D<T>()
{
    return CMVector3D<T>(x,y,z);
}

////////////////////////////////////////////////////////////////////////////
// Cast To CMVector3D Type
template <class T>
M_INLINE CMQuaternion<T>::operator CMVector4D<T>()
{
    return CMVector4D<T>(x,y,z,w);
}

////////////////////////////////////////////////////////////////////////////
// Access Grants
template <class T>
M_INLINE T CMQuaternion<T>::operator [] (int32t i)
{
    ASSERT(i >= 0 && i <= 3);
    return q[i];
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE T CMQuaternion<T>::operator [] (int32t i) const
{
    ASSERT(i >= 0 && i <= 3);
    return q[i];
}

////////////////////////////////////////////////////////////////////////////
// Assignment Operators
template <class T>
M_INLINE CMQuaternion<T>& CMQuaternion<T>::operator = ( const CMQuaternion &quat )
{
    x = quat.x;
    y = quat.y;
    z = quat.z;
    w = quat.w;

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T>& CMQuaternion<T>::operator += ( const CMQuaternion &quat )
{
    x += quat.x;
    y += quat.y;
    z += quat.z;
    w += quat.w;

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T>& CMQuaternion<T>::operator -= ( const CMQuaternion &quat )
{
    x -= quat.x;
    y -= quat.y;
    z -= quat.z;
    w -= quat.w;

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T>& CMQuaternion<T>::operator *= ( const CMQuaternion &quat )
{
    CMQuatMult(*this,*this,quat);
    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T>& CMQuaternion<T>::operator *= ( const T &val )
{
    x *= val;
    y *= val;
    z *= val;
    w *= val;

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T>& CMQuaternion<T>::operator /= (const T &val )
{
    ASSERT(!mnear0(val));

    T inv = minvert(val);
    x *= inv;
    y *= inv;
    z *= inv;
    w *= inv;

    return *this;
}

////////////////////////////////////////////////////////////////////////////
// Unary Operators
template <class T>
M_INLINE CMQuaternion<T> CMQuaternion<T>::operator + () const
{
    return CMQuaternion(x,y,z,w);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T> CMQuaternion<T>::operator - () const
{
    return CMQuaternion(-x,-y,-z,-w);
}

////////////////////////////////////////////////////////////////////////////
// Binary Operators
template <class T>
M_INLINE CMQuaternion<T> CMQuaternion<T>::operator + ( const CMQuaternion &quat ) const
{
    return CMQuaternion(x + quat.x, y + quat.y, z + quat.z, w + quat.w);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T> CMQuaternion<T>::operator - ( const CMQuaternion &quat ) const
{
    return CMQuaternion(x - quat.x, y - quat.y, z - quat.z, w - quat.w);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T> CMQuaternion<T>::operator * ( const CMQuaternion &quat ) const
{
    CMQuaternion<T> q;

    CMQuatMult(q,*this,quat);

    return q;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T> CMQuaternion<T>::operator * ( const T &val ) const
{
    return CMQuaternion(x*val, y*val, z*val, w*val);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMQuaternion<T> CMQuaternion<T>::operator / ( const T &val ) const
{
    T inv = minvert(val);
    return CMQuaternion(x*inv, y*inv, z*inv, w*inv);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE boolt CMQuaternion<T>::operator == ( const CMQuaternion & quat ) const
{
    return (x == quat.x) && (y == quat.y) && (z == quat.z) && (w == quat.w);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE boolt CMQuaternion<T>::operator != ( const CMQuaternion & quat ) const
{
    return (x != quat.x) || (y != quat.y) || (z != quat.z) || (w != quat.w);
}


////////////////////////////////////////////////////////////////////////////
// Adding Two Quaternion
template <class T>
M_INLINE void CMQuatAdd(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2)
{
    qOut.x = Q1.x + Q2.x;
    qOut.y = Q1.y + Q2.y;
    qOut.z = Q1.z + Q2.z;
    qOut.w = Q1.w + Q2.w;
    /*
    ASSERT(pQ != NULL);
    ASSERT(pQ1 != NULL);
    ASSERT(pQ2 != NULL);

    pQ->x = pQ1->x + pQ2->x;
    pQ->y = pQ1->y + pQ2->y;
    pQ->z = pQ1->z + pQ2->z;
    pQ->w = pQ1->w + pQ2->w;
    */
}


////////////////////////////////////////////////////////////////////////////
// Subtracting Two Quaternions
template <class T>
M_INLINE void CMQuatSub(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2)
{
    qOut.x = Q1.x - Q2.x;
    qOut.y = Q1.y - Q2.y;
    qOut.z = Q1.z - Q2.z;
    qOut.w = Q1.w - Q2.w;

    /*
    ASSERT(pQ != NULL);
    ASSERT(pQ1 != NULL);
    ASSERT(pQ2 != NULL);

    pQ->x = pQ1->x - pQ2->x;
    pQ->y = pQ1->y - pQ2->y;
    pQ->z = pQ1->z - pQ2->z;
    pQ->w = pQ1->w - pQ2->w;
    */
}

////////////////////////////////////////////////////////////////////////////
// Set Identity Quaternion
template <class T>
M_INLINE void CMQuatIdentity(CMQuaternion<T> &Q)
{
    Q->x = Q->y = Q->z = 0;
    Q->w = T(1);
}

////////////////////////////////////////////////////////////////////////////
// Is Identity
template <class T>
M_INLINE boolt CMQuatIsIdentity(CMQuaternion<T> &Q)
{
    return (Q.x == T(0)) && (Q.y == T(0)) && (Q.z == T(0)) && (Q.w == T(1));
}

////////////////////////////////////////////////////////////////////////////
// Quaternion Dot Product
template <class T>
M_INLINE T CMQuatDot(const CMQuaternion<T> &Q1, const CMQuaternion<T> &Q2)
{
    return Q1.x*Q2.x + Q1.y*Q2.y + Q1.z*Q2.z + Q1.w*Q2.w;
}

////////////////////////////////////////////////////////////////////////////
// Conjugate Quaternion
template <class T>
M_INLINE void CMQuatConjugate(CMQuaternion<T> &qOut, const CMQuaternion<T> &Q)
{
    qOut.x = -Q.x;
    qOut.y = -Q.y;
    qOut.z = -Q.z;
    qOut.w = Q.w;
}


///////////////////////////////////////////////////////////////////////////////////////

#endif //_QUATERNION_H_
