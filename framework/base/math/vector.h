//-----------------------------------------------------------------------------
// vector.h
//
// Vectors template clasess declarations
//
//-----------------------------------------------------------------------------

#ifndef _CMVECTOR2D_
#define _CMVECTOR2D_

//-----------------------------------------------------------------------------
// Class template declaration
//-----------------------------------------------------------------------------

//---------------------------
// 2D Vector class template declaration
template <class T>
class CMVector2D
{
public:
    // Member variables
    union {
        struct{
            T x, y;
        };
        T v[2];
    };

public:
    // Constructors declarations
    CMVector2D();
    CMVector2D(const T * const p);
    CMVector2D(T x, T y);
    CMVector2D(const CMVector2D & v);

    // Cast to pointer  (T *)
    operator T * () const;
    operator const T * () const;

    // Assignment operators
    CMVector2D & operator =  ( const CMVector2D & );
    CMVector2D & operator += ( const CMVector2D & );
    CMVector2D & operator -= ( const CMVector2D & );
    CMVector2D & operator *= ( const T & );
    CMVector2D & operator /= ( const T & );

    // Unary operators
    CMVector2D operator + () const;
    CMVector2D operator - () const;

    T operator [] (int32t );

    // Binary operators
    CMVector2D operator + ( const CMVector2D & ) const;
    CMVector2D operator - ( const CMVector2D & ) const;
    CMVector2D operator * ( const T & ) const;
    CMVector2D operator / ( const T & ) const;

    bool operator == ( const CMVector2D & ) const;
    bool operator != ( const CMVector2D & ) const;
};


//---------------------------
// Inline
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Adding two 2D vectors
template <class T>
void CMVec2Add( CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Subtracting two 2D vectors
template <class T>
void CMVec2Sub( CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Dot product between two 2D vectors
template <class T>
T CMVec2Dot(const CMVector2D<T> & v1, const CMVector2D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Returns the z-component by taking the cross product of two 2D vectors.
template <class T>
T CMVec2CCW(const CMVector2D<T> & v1, const CMVector2D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Performs a linear interpolation between two 2D vectors
template <class T>
void CMVec2Lerp(CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2,
                T weight);

///////////////////////////////////////////////////////////////////////////////
// Returns a 2D vector that is made up of the largest components
// of two 2D vectors
template <class T>
void CMVec2Max( CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Returns a 2D vector that is made up of the smallest components
// of two 2D vectors
template <class T>
void CMVec2Min( CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Scales a 2D vector
template <class T>
void CMVec2Scale(CMVector2D<T> & vOut,
                const CMVector2D<T> & vIn,
                T scale);

///////////////////////////////////////////////////////////////////////////////
// Invert 2D vector scale
template <class T>
void CMVec2ScaleInv(CMVector2D<T> & vOut,
                    const CMVector2D<T> & vIn,
                    T scale);

//---------------------------
// Non-Inline
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Rotate a 2D vector, using the specified angle in radians
template <class T>
void CMVec2Rotate(  CMVector2D<T> & vOut,
                    const CMVector2D<T> & vIn,
                    T angle);

///////////////////////////////////////////////////////////////////////////////
// Returns the length of a 2D vector
template <class T>
T CMVec2Length(const CMVector2D<T> & vec);

///////////////////////////////////////////////////////////////////////////////
// Returns the square of the length of a 2D vector
template <class T>
T CMVec2LengthSq(const CMVector2D<T> & vec);

///////////////////////////////////////////////////////////////////////////////
// Normalize a 2D vector
template <class T>
void CMVec2Normalize(CMVector2D<T> & vOut, const CMVector2D<T> & vIn);

///////////////////////////////////////////////////////////////////////////////
// Calculate length between two 2D points
template <class T>
T CMVec2Distance(const CMVector2D<T> & p0, const CMVector2D<T> & p1);

///////////////////////////////////////////////////////////////////////////////
// Calculate square length between two 2D points
template <class T>
T CMVec2DistanceSq(const CMVector2D<T> & p0, const CMVector2D<T> & p1);

///////////////////////////////////////////////////////////////////////////////
// Getting angle in radians between two 2D vectors
template <class T>
T CMVec2Angle(const CMVector2D<T> & v1, const CMVector2D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Getting reflected 2D vector
template <class T>
void CMVec2Reflect( CMVector2D<T> & vR,
                    const CMVector2D<T> & vL,
                    const CMVector2D<T> & vN);

///////////////////////////////////////////////////////////////////////////////
// Getting refracted 2D vector
template <class T>
void CMVec2Refract( CMVector2D<T> & vR,
                    const CMVector2D<T> & vL,
                    const CMVector2D<T> & vN,
                    T etaRatio);

///////////////////////////////////////////////////////////////////////////////
// Performs a catmull-rom interpolation, using the specified 2D vectors
template <class T>
void CMVec2CatmullRom(  CMVector2D<T> & vOut,
                        const CMVector2D<T> & v1,
                        const CMVector2D<T> & v2,
                        const CMVector2D<T> & v3,
                        const CMVector2D<T> & v4,
                        T weightingFactor);

///////////////////////////////////////////////////////////////////////////////
// Returns a point in barycentric coordinates, using the specified 2D vectors
template <class T>
void CMVec2BaryCentric( CMVector2D<T> & vOut,
                        const CMVector2D<T> & v1,
                        const CMVector2D<T> & v2,
                        const CMVector2D<T> & v3,
                        T f, T g);

///////////////////////////////////////////////////////////////////////////////
// Performs a hermite spline interpolation, using the specified 2D vectors
template <class T>
void CMVec2Hermite( CMVector2D<T> &vOut,
                    const CMVector2D<T> &v1,
                    const CMVector2D<T> &vT1,
                    const CMVector2D<T> &v2,
                    const CMVector2D<T> &vT2,
                    T weight, T tension, T bias);

///////////////////////////////////////////////////////////////////////////////
// Performs a quadratic interpolation, using the specified 2D vectors
template <class T>
void CMVec2Quadratic(CMVector2D<T> & vOut,
                    const CMVector2D<T> & v1,
                    const CMVector2D<T> & v2,
                    const CMVector2D<T> & v3,
                    T factor);

///////////////////////////////////////////////////////////////////////////////
// Performs a cubic interpolation, using the specified 2D vectors
template <class T>
void CMVec2Cubic(CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2,
                const CMVector2D<T> & v3,
                const CMVector2D<T> & v4,
                T factor);

///////////////////////////////////////////////////////////////////////////////
// Performs a square interpolation, using two 2D vectors
template <class T>
void CMVec2Sqrp(CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2,
                T factor);

///////////////////////////////////////////////////////////////////////////////
// Performs a cosinus interpolation, using two 2D vectors
template <class T>
void CMVec2Cosrp(CMVector2D<T> & vOut,
                const CMVector2D<T> & v1,
                const CMVector2D<T> & v2,
                T factor);

//-----------------------------------------------------------------------------
// Class template definition
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// 2D Vector structure constructors

///////////////////////////////////////////////////////////////////////////////
// Default constructor
template <class T>
inline CMVector2D<T>::CMVector2D()
{
    x = T(0.0);
    y = T(0.0);
}

///////////////////////////////////////////////////////////////////////////////
// Getting pointer to array
template <class T>
inline CMVector2D<T>::CMVector2D(const T * const ptr)
{
    CMVector2D();
    if(ptr)
    {
        x = ptr[0];
        y = ptr[1];
    }
}

///////////////////////////////////////////////////////////////////////////////
// Getting parameters
template <class T>
inline CMVector2D<T>::CMVector2D(T x, T y)
{
    this->x = x;
    this->y = y;
}

///////////////////////////////////////////////////////////////////////////////
// Getting another 2D vector
template <class T>
inline CMVector2D<T>::CMVector2D(const CMVector2D & v)
{
    x = v.x;
    y = v.y;
}

///////////////////////////////////////////////////////////////////////////////
// Cast to pointer  (T *)
template <class T>
inline CMVector2D<T>::operator T * () const
{
    return (T *) this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T>::operator const T * () const
{
    return (const T *) this;
}

///////////////////////////////////////////////////////////////////////////////
// Assignment operators
template <class T>
inline CMVector2D<T> & CMVector2D<T>::operator = (const CMVector2D & v)
{
    x = v.x;
    y = v.y;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T> & CMVector2D<T>::operator += (const CMVector2D & v)
{
    x += v.x;
    y += v.y;

    return *this;
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T> & CMVector2D<T>::operator -= (const CMVector2D & v)
{
    x -= v.x;
    y -= v.y;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T> & CMVector2D<Treal>::operator *= (const T & val)
{
    x *= val;
    y *= val;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T> & CMVector2D<T>::operator /= (const T & val)
{
    ASSERT(!mnear0(val));

    x /= val;
    y /= val;

    return *this;
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
inline T CMVector2D<T>::operator [] (int32t i)
{
    ASSERT(i >= 0 && i <= 1);
    return v[i];
}

///////////////////////////////////////////////////////////////////////////////
// Unary operators
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T> CMVector2D<T>::operator + () const
{
    return *this;
}

// Negate
template <class T>
inline CMVector2D<T> CMVector2D<T>::operator - () const
{
    return CMVector2D(-x, -y);
}

///////////////////////////////////////////////////////////////////////////////
// Binary operators
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T> CMVector2D<T>::operator + (const CMVector2D & v) const
{
    return CMVector2D(x + v.x, y + v.y);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T> CMVector2D<T>::operator - (const CMVector2D & v) const
{
    return CMVector2D(x - v.x, y - v.y);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T> CMVector2D<T>::operator * (const T & val) const
{
    return CMVector2D(x * val, y * val);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector2D<T> CMVector2D<T>::operator / (const T & val) const
{
    ASSERT(!mnear0(val));

    T rev = T(1.0) / val;
    return CMVector2D(x * rev, y * rev);
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
inline bool CMVector2D<T>::operator == (const CMVector2D & v) const
{
    return ((x == v.x) && (y == v.y)) ? (1) : (0);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline bool CMVector2D<T>::operator != (const CMVector2D & v) const
{
    return ((x == v.x) && (y == v.y)) ? (0) : (1);
}


//-----------------------------------------------------------------------------
// Methods operating on CMVector2D class
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// Adding two vectors
template <class T>
inline void CMVec2Add(  CMVector2D<T> & vOut,
                        const CMVector2D<T> & v1,
                        const CMVector2D<T> & v2)
{
    vOut.x = v1.x + v2.x;
    vOut.y = v1.y + v2.y;

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Subtracting two vectors
template <class T>
inline void CMVec2Sub(  CMVector2D<T> & vOut,
                        const CMVector2D<T> & v1,
                        const CMVector2D<T> & v2)
{
    vOut.x = v1.x - v2.x;
    vOut.y = v1.y - v2.y;

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Dot product between two 2D vectors
template <class T>
inline T CMVec2Dot(const CMVector2D<T> & v1, const CMVector2D<T> & v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

///////////////////////////////////////////////////////////////////////////////
// Return Z component by taking cross product of two 2D vectors
template <class T>
inline T CMVec2CCW(const CMVector2D<T> & v1, const CMVector2D<T> & v2)
{
    return v1.x * v2.y - v1.y * v2.x;
}

///////////////////////////////////////////////////////////////////////////////
// Returns a 2D vector that is made up of the largest components
// of two 2D vectors
template <class T>
inline void CMVec2Max(  CMVector2D<T> & vOut,
                        const CMVector2D<T> & v1,
                        const CMVector2D<T> & v2)
{
    vOut.x = (v1.x > v2.x) ? (v1.x) : (v2.x);
    vOut.y = (v1.y > v2.y) ? (v1.y) : (v2.y);

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Returns a 2D vector that is made up of the smallest components
// of two 2D vectors
template <class T>
inline void CMVec2Min(  CMVector2D<T> & vOut,
                        const CMVector2D<T> & v1,
                        const CMVector2D<T> & v2)
{
    vOut.x = (v1.x < v2.x) ? (v1.x) : (v2.x);
    vOut.y = (v1.y < v2.y) ? (v1.y) : (v2.y);

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Scales a 2D vector
template <class T>
inline void CMVec2Scale(CMVector2D<T> & vOut,
                        const CMVector2D<T> & vIn, T scale)
{
    vOut.x = vIn.x * scale;
    vOut.y = vIn.x * scale;

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Invert 2D vector scale
template <class T>
inline void CMVec2ScaleInv( CMVector2D<T> & vOut,
                            const CMVector2D<T> & vIn, T scale)
{
    ASSERT(!mnear0(scale));

    T invScale = T(1.0) / scale;

    vOut.x = vIn.x * invScale;
    vOut.y = vIn.y * invScale;

    return ;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a linear interpolation between two 2D vectors
template <class T>
inline void CMVec2Lerp( CMVector2D<T> & vOut,
                        const CMVector2D<T> & v1,
                        const CMVector2D<T> & v2, T weight)
{
    T x = v1.x + ((v2.x - v1.x) * weight);
    T y = v1.y + ((v2.y - v1.y) * weight);

    vOut.x = x;
    vOut.y = y;

    return ;
}

#endif //_CMVECTOR2D_


///////////////////////////////////////////////////////////////////////////////
#ifndef _CMVECTOR3D_
#define _CMVECTOR3D_


// 3D vector class template declaration
template <class T>
class CMVector3D
{
public:
    // Member variables
    union {
        struct{
            T x, y, z;
        };
        T v[3];
    };

public:
    // Constructor declarations
    CMVector3D();
    CMVector3D(const T * const p);
    CMVector3D(T x, T y, T z);
    CMVector3D(const CMVector3D & v);
    CMVector3D(const CMVector2D & v);

    // Cast to pointer  (T *)
    operator T * () const;
    operator const T * () const;

    // Cast to CMVector2D type
    operator  CMVector2D<T>();


    // Assignment operators
    CMVector3D & operator =  ( const CMVector3D & );
    CMVector3D & operator += ( const CMVector3D & );
    CMVector3D & operator -= ( const CMVector3D & );
    CMVector3D & operator *= ( const T & );
    CMVector3D & operator /= ( const T & );

    // Unary operators
    CMVector3D operator + () const;
    CMVector3D operator - () const;

    T operator [] (int32t );

    // Binary operators
    CMVector3D operator + ( const CMVector3D & ) const;
    CMVector3D operator - ( const CMVector3D & ) const;
    CMVector3D operator * ( const T & ) const;
    CMVector3D operator / ( const Treal & ) const;

    bool operator == ( const CMVector3D & ) const;
    bool operator != ( const CMVector3D & ) const;

};

//---------------------------
// Inline+
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Adding two 3D vectors
template <class T>
void CMVec3Add( CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Subtracting two 3D vectors
template <class T>
void CMVec3Sub( CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Dot product between two 3D vectors
template <class T>
T CMVec3Dot(const CMVector3D<T> & v1, const CMVector3D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Cross product of two 3D vectors
template <class T>
void CMVec3Cross(CMVector3D<T> & vOut,
                const CMVector3D<T> & v1, const CMVector3D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Performs a linear interpolation between two 3D vectors
template <class T>
void CMVec3Lerp(CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2, T weight);

///////////////////////////////////////////////////////////////////////////////
// Returns a 3D vector that is made up of the largest components
// of two 3D vectors
template <class T>
void CMVec3Max( CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Returns a 3D vector that is made up of the smallest components
// of two 3D vectors
template <class T>
void CMVec3Min( CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Scales a 3D vector
template <class T>
void CMVec3Scale(CMVector3D<T> & vOut,
                const CMVector3D<T> & vIn, T scale);

///////////////////////////////////////////////////////////////////////////////
// Invert scale of 3D vector
template <class T>
void CMVec3ScaleInv(CMVector3D<T> & vOut,
                    const CMVector3D<T> & vIn, T scale);

///////////////////////////////////////////////////////////////////////////////
// Setting 3D vectors values
template <class T>
void CMVec3Set(CMVector3D<T> & v, T x, T y, T z);

//---------------------------
// Non-Inline+
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Rotate a 3D vector in X axis, using the specified angle in radians
template <class T>
void CMVec3RotateX( CMVector3D<T> & vOut,
                    const CMVector3D<T> & vIn, T rad);

///////////////////////////////////////////////////////////////////////////////
// Rotate a 3D vector in Y axis, using the specified angle in radians
template <class T>
void CMVec3RotateY( CMVector3D<T> & vOut,
                    const CMVector3D<T> & vIn, T rad);

///////////////////////////////////////////////////////////////////////////////
// Rotate a 3D vector in Z axis, using the specified angle in radians
template <class T>
void CMVec3RotateZ( CMVector3D<T> & vOut,
                    const CMVector3D<T> & vIn, T rad);

///////////////////////////////////////////////////////////////////////////////
// Returns the length of a 3D vector
template <class T>
T CMVec3Length(const CMVector3D<T> & v);

///////////////////////////////////////////////////////////////////////////////
// Returns the square of the length of a 3D vector
template <class T>
T CMVec3LengthSq(const CMVector3D<T> & v);

///////////////////////////////////////////////////////////////////////////////
// Calculate length between two 3D points
template <class T>
T CMVec3Distance(const CMVector3D<T> & p0, const CMVector3D<T> & p1);

///////////////////////////////////////////////////////////////////////////////
// Calculate square length between two 3D points
template <class T>
T CMVec3DistanceSq(const CMVector3D<T> & p0, const CMVector3D<T> & p1);

///////////////////////////////////////////////////////////////////////////////
// Normalize a 3D vector
template <class T>
void CMVec3Normalize(CMVector3D<T> & vOut, const CMVector3D<T> & vIn);

///////////////////////////////////////////////////////////////////////////////
// Getting angle in radians between two 3D vectors
template <class T>
T CMVec3Angle(const CMVector3D<T> & v1, const CMVector3D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Getting reflected 3D vector
template <class T>
void CMVec3Reflect( CMVector3D<T> & vR,
                    const CMVector3D<T> & vL,
                    const CMVector3D<T> & vN);

///////////////////////////////////////////////////////////////////////////////
// Getting refracted 3D vector
template <class T>
void CMVec3Refract( CMVector3D<T> & vR,
                    const CMVector3D<T> & vL,
                    const CMVector3D<T> & vN, T etaRatio);

///////////////////////////////////////////////////////////////////////////////
// Performs a catmull-rom interpolation, using the specified 3D vectors
template <class T>
void CMVec3CatmullRom(  CMVector3D<T> & vOut,
                        const CMVector3D<T> & v1,
                        const CMVector3D<T> & v2,
                        const CMVector3D<T> & v3,
                        const CMVector3D<T> & v4,  T factor);

///////////////////////////////////////////////////////////////////////////////
// Returns a point in barycentric coordinates, using the specified 3D vectors
template <class T>
void CMVec3BaryCentric( CMVector3D<T> & vOut,
                        const CMVector3D<T> & v1,
                        const CMVector3D<T> & v2,
                        const CMVector3D<T> & v3, T f, T g);

///////////////////////////////////////////////////////////////////////////////
// Performs a hermite spline interpolation, using the specified 3D vectors
template <class T>
void CMVec3Hermite( CMVector3D<T> & vOut,
                    const CMVector3D<T> & v1, const CMVector3D<T> & vT1,
                    const CMVector3D<T> & v2, const CMVector3D<T> & vT2,
                    T weight, T tension, T bias);

///////////////////////////////////////////////////////////////////////////////
// Performs a quadratic interpolation, using the specified 3D vectors
template <class T>
void CMVec3Quadratic(CMVector3D<T> & vOut,
                    const CMVector3D<T> & v1,
                    const CMVector3D<T> & v2,
                    const CMVector3D<T> & v3, T factor);

///////////////////////////////////////////////////////////////////////////////
// Performs a cubic interpolation, using the specified 3D vectors
template <class T>
void CMVec3Cubic(CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2,
                const CMVector3D<T> & v3,
                const CMVector3D<T> & v4, T factor);

///////////////////////////////////////////////////////////////////////////////
// Performs a square interpolation, using two 3D vectors
template <class T>
void CMVec3Sqrp(CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2, T factor);

///////////////////////////////////////////////////////////////////////////////
// Performs a cosinus interpolation, using two 3D vectors
template <class T>
void CMVec3Cosrp(CMVector3D<T> & vOut,
                const CMVector3D<T> & v1,
                const CMVector3D<T> & v2, T factor);

//-----------------------------------------------------------------------------
// Class template definition
//-----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////////////
// 3D vector structure constructors

// Default constructor
template <class T>
inline CMVector3D<T>::CMVector3D()
{
    x = T(0.0);
    y = T(0.0);
    z = T(0.0);
}

///////////////////////////////////////////////////////////////////////////////
// Getting pointer to array
template <class T>
inline CMVector3D<T>::CMVector3D(const T  * const p)
{
    CMVector3D();
    if(p)
    {
        x = p[0];
        y = p[1];
        z = p[2];
    }
}

///////////////////////////////////////////////////////////////////////////////
// Getting parameters
template <class T>
inline CMVector3D<T>::CMVector3D(T x, T y, T z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

///////////////////////////////////////////////////////////////////////////////
// Copy constructor
template <class T>
inline CMVector3D<T>::CMVector3D(const CMVector3D & v)
{
    x = v.x;
    y = v.y;
    z = v.z;
}

///////////////////////////////////////////////////////////////////////////////
// Getting 2D vector
template <class T>
inline CMVector3D<T>::CMVector3D(const CMVector2D<T> & v)

{
    x = v.x;
    y = v.y;
    z = T(0.0);
}

///////////////////////////////////////////////////////////////////////////////
//Cast to pointer  (T *)
template <class T>
inline CMVector3D<T>::operator T * () const
{
    return (T*) this;
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T>::operator const T * () const
{
    return (const T*) this;
}

///////////////////////////////////////////////////////////////////////////////
//Cast to CMVector2D type
template <class T>
inline CMVector3D<T>::operator CMVector2D<T>()
{
    return CMVector2D<T>(x,y);
}

///////////////////////////////////////////////////////////////////////////////
// Assignment operators
template <class T>
inline CMVector3D<T>& CMVector3D<T>::operator =  ( const CMVector3D & v)
{
    x = v.x;
    y = v.y;
    z = v.z;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T>& CMVector3D<T>::operator += ( const CMVector3D & v)
{
    x += v.x;
    y += v.y;
    z += v.z;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T>& CMVector3D<T>::operator -= ( const CMVector3D & v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T>& CMVector3D<T>::operator *= ( const T & val)
{
    x *= val;
    y *= val;
    z *= val;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T>& CMVector3D<T>::operator /= ( const T & val)
{
    ASSERT(!mnear0(val));

    x /= val;
    y /= val;
    z /= val;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline T CMVector3D<T>::operator [](int32t i)
{
    ASSERT(i >= 0 && i <= 2);
    return v[i];
}

///////////////////////////////////////////////////////////////////////////////
// Unary operators
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T> CMVector3D<T>::operator + () const
{
    return CMVector3D(x, y, z);
}

// Negate
template <class T>
inline CMVector3D<T> CMVector3D<T>::operator - () const
{
    return CMVector3D(-x, -y, -z);
}

///////////////////////////////////////////////////////////////////////////////
// Binary operators
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T> CMVector3D<T>::operator + ( const CMVector3D & v) const
{
    return CMVector3D(x + v.x, y + v.y, z + v.z);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T> CMVector3D<T>::operator - ( const CMVector3D & v) const
{
    return CMVector3D(x - v.x, y - v.y, z - v.z);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T> CMVector3D<T>::operator * ( const T & val) const
{
    return CMVector3D(x*val, y*val, z*val);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector3D<T> CMVector3D<T>::operator / ( const T & val) const
{
    T inv = T(1)/val;
    ASSERT(!mnear0(inv));
    return CMVector3D(x * inv, y * inv, z * inv);
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
inline bool CMVector3D<T>::operator == ( const CMVector3D & v) const
{
    return ((x == v.x) && (y == v.y) && (z == v.z)) ? (1) : (0);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline bool CMVector3D<T>::operator != ( const CMVector3D & v) const
{
    return ((x == v.x) && (y == v.y) && (z == v.z)) ? (0) : (1);
}

///////////////////////////////////////////////////////////////////////////////
// Methods operating on CMVector3D class
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Adding two vectors
template <class T>
inline void CMVec3Add(CMVector3D<T> &vOut, const CMVector3D<T> &v1, const CMVector3D<T> &v2)
{
    vOut.x = v1.x + v2.x;
    vOut.y = v1.y + v2.y;
    vOut.z = v1.z + v2.z;
}
///////////////////////////////////////////////////////////////////////////////
// Subtracting two vectors
template <class T>
inline void CMVec3Sub(CMVector3D<T> &vOut, const CMVector3D<T> &v1, const CMVector3D<T> &v2)
{
    vOut.x = v1.x - v2.x;
    vOut.y = v1.y - v2.y;
    vOut.z = v1.z - v2.z;
}
///////////////////////////////////////////////////////////////////////////////
// Dot product between two 3D vectors
template <class T>
inline T CMVec3Dot(const CMVector3D<T> &v1, const CMVector3D<T> &v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
///////////////////////////////////////////////////////////////////////////////
// Cross product of two vectors
template <class T>
inline void CMVec3Cross(CMVector3D<T> &vOut, const CMVector3D<T> &v1, const CMVector3D<T> &v2)
{
    vOut.x = v1.y*v2.z - v1.z*v2.y;
    vOut.y = v1.z*v2.x - v1.x*v2.z; //v1.x*v2.z - v1.z*v2.x;
    vOut.z = v1.x*v2.y - v1.y*v2.x;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a linear interpolation between two 3D vectors
template <class T>
inline void CMVec3Lerp(CMVector3D<T> &vOut, const CMVector3D<T> &v1, const CMVector3D<T> &v2, T weight)
{
    T x = v1.x + ((v2.x - v1.x)*weight);
    T y = v1.y + ((v2.y - v1.y)*weight);
    T z = v1.z + ((v2.z - v1.z)*weight);

    vOut.x = x;
    vOut.y = y;
    vOut.z = z;
}

///////////////////////////////////////////////////////////////////////////////
// Returns a 3D vector that is made up of the largest components of two 3D vectors
template <class T>
inline void CMVec3Max(CMVector3D<T> &vOut, const CMVector3D<T> &v1, const CMVector3D<T> &v2)
{
    vOut.x = (v1.x > v2.x) ? (v1.x) : (v2.x);
    vOut.y = (v1.y > v2.y) ? (v1.y) : (v2.y);
    vOut.z = (v1.z > v2.z) ? (v1.z) : (v2.z);
}

///////////////////////////////////////////////////////////////////////////////
// Returns a 3D vector that is made up of the smallest components of two 3D vectors
template <class T>
inline void CMVec3Min(CMVector3D<T> &vOut, const CMVector3D<T> &v1, const CMVector3D<T> &v2)
{
    vOut.x = (v1.x < v2.x) ? (v1.x) : (v2.x);
    vOut.y = (v1.y < v2.y) ? (v1.y) : (v2.y);
    vOut.z = (v1.z < v2.z) ? (v1.z) : (v2.z);
}

///////////////////////////////////////////////////////////////////////////////
// Scales a 3D vector
template <class T>
inline void CMVec3Scale(CMVector3D<T> &vOut, const CMVector3D<T> &vIn, T scale)
{
    vOut.x = vIn.x*scale;
    vOut.y = vIn.y*scale;
    vOut.z = vIn.z*scale;
}

///////////////////////////////////////////////////////////////////////////////
// Invert scale of 3D vector
template <class T>
inline void CMVec3ScaleInv(CMVector3D<T> &vOut, const CMVector3D<T> &vIn, T scale)
{
    T invscale = T(1)/scale;

    vOut.x = vIn.x*invscale;
    vOut.y = vIn.y*invscale;
    vOut.z = vIn.z*invscale;
}

///////////////////////////////////////////////////////////////////////////////
// Setting 3D vectors values
template <class T>
void CMVec3Set(CMVector3D<T> &v, T x, T y, T z)
{
    v.x = x;
    v.y = y;
    v.z = z;
}

///////////////////////////////////////////////////////////////////////////////

#endif //_CMVECTOR3D_


///////////////////////////////////////////////////////////////////////////////
#ifndef _CMVECTOR4D_
#define _CMVECTOR4D_

// 4D vector class template declaration
template <class T>
class CMVector4D
{
public:
    // Member variables
    union {
        struct{
            T x, y, z, w;
        };
        T v[4];
    };

public:
    // Constructor declarations
    CMVector4D();
    CMVector4D(const T  * const p);
    CMVector4D(T x, T y, T z, T w);
    CMVector4D(const CMVector4D & v);
    CMVector4D(const CMVector3D & v);
    CMVector4D(const CMVector2D & v);

    // Cast to pointer  (T *)
    operator T * () const;
    operator const T * () const;

    // Cast to CMVector3D type
    operator  CMVector3D<T>();

    // Cast to CMVector2D type
    operator  CMVector2D<T>();

    // Assignment operators
    CMVector4D& operator =  ( const CMVector4D & );
    CMVector4D& operator += ( const CMVector4D & );
    CMVector4D& operator -= ( const CMVector4D & );
    CMVector4D& operator *= ( const T & );
    CMVector4D& operator /= ( const T & );

    // Unary operators
    CMVector4D operator + () const;
    CMVector4D operator - () const;

    T operator [] (int32t );

    // Binary operators
    CMVector4D operator + ( const CMVector4D & ) const;
    CMVector4D operator - ( const CMVector4D & ) const;
    CMVector4D operator * ( const T & ) const;
    CMVector4D operator / ( const T & ) const;

    bool operator == ( const CMVector4D & ) const;
    bool operator != ( const CMVector4D & ) const;

};

//---------------------------
// Inline+
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Adding two 4D vectors
template <class T>
void CMVec4Add( CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Subtracting two 4D vectors
template <class T>
void CMVec4Sub( CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Dot product between two 4D vectors
template <class T>
T CMVec4Dot(const CMVector3D<T> & v1, const CMVector3D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Performs a linear interpolation between two 4D vectors
template <class T>
void CMVec4Lerp(CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2, T weight);

///////////////////////////////////////////////////////////////////////////////
// Returns a 4D vector that is made up of the largest components
// of two 4D vectors
template <class T>
void CMVec4Max( CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Returns a 4D vector that is made up of the smallest components
// of two 4D vectors
template <class T>
void CMVec4Min( CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Scales a 4D vector
template <class T>
void CMVec4Scale(CMVector4D<T> & vOut, const CMVector4D<T> & vIn, T scale);

///////////////////////////////////////////////////////////////////////////////
// Invert scale of 4D vector
template <class T>
void CMVec4ScaleInv(CMVector4D<T> & vOut, const CMVector4D<T> & vIn, T scale);

//---------------------------
// Non-Inline+
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Cross product of two 4D vectors
template <class T>
void CMVec4Cross(CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Rotate a 4D vector in X axis, using the specified angle in radians
template <class T>
void CMVec4RotateX(CMVector4D<T> & vOut, const CMVector4D<T> & vIn, T rad);

///////////////////////////////////////////////////////////////////////////////
// Rotate a 4D vector in Y axis, using the specified angle in radians
template <class T>
void CMVec4RotateY(CMVector4D<T> & vOut, const CMVector4D<T> & vIn, T rad);

///////////////////////////////////////////////////////////////////////////////
// Rotate a 4D vector in Z axis, using the specified angle in radians
template <class T>
void CMVec4RotateZ(CMVector4D<T> & vOut, const CMVector4D<T> & vIn, T rad);

///////////////////////////////////////////////////////////////////////////////
// Returns the length of a 4D vector
template <class T>
T CMVec4Length(const CMVector4D<T> & v);

///////////////////////////////////////////////////////////////////////////////
// Returns the square of the length of a 4D vector
template <class T>
T CMVec4LengthSq(const CMVector4D<T> & v);

///////////////////////////////////////////////////////////////////////////////
// Calculate length between two 4D points
template <class T>
T CMVec4Distance(CMVector4D<T> & p0, CMVector4D<T> & p1);

///////////////////////////////////////////////////////////////////////////////
// Calculate square length between two 4D points
template <class T>
T CMVec4DistanceSq(CMVector4D<T> & p0, CMVector4D<T> & p1);

///////////////////////////////////////////////////////////////////////////////
// Normalize a 4D vector
template <class T>
void CMVec4Normalize(CMVector4D<T> & vOut, const CMVector4D<T> & vIn);

///////////////////////////////////////////////////////////////////////////////
// Getting angle in radians between two 4D vectors
template <class T>
T CMVec4Angle(const CMVector4D<T> & v1, const CMVector4D<T> & v2);

///////////////////////////////////////////////////////////////////////////////
// Getting reflected 4D vector
template <class T>
void CMVec4Reflect( CMVector4D<T> & vR,
                    const CMVector4D<T> & vL, const CMVector4D<T> & vN);

///////////////////////////////////////////////////////////////////////////////
// Getting refracted 4D Vector
template <class T>
void CMVec4Refract( CMVector4D<T> & vR,
                    const CMVector4D<T> & vL,
                    const CMVector4D<T> & vN, T etaRatio);

///////////////////////////////////////////////////////////////////////////////
// Performs a catmull-rom interpolation, using the specified 4D vectors
template <class T>
void CMVec4CatmullRom(  CMVector4D<T> & vOut,
                        const CMVector4D<T> & v1,
                        const CMVector4D<T> & v2,
                        const CMVector4D<T> & v3,
                        const CMVector4D<T> & v4, T f);

///////////////////////////////////////////////////////////////////////////////
// Returns a point in barycentric coordinates, using the specified 4D vectors
template <class T>
void CMVec4BaryCentric( CMVector4D<T> & vOut,
                        const CMVector4D<T> & v1,
                        const CMVector4D<T> & v2,
                        const CMVector4D<T> & v3, T f, T g);

///////////////////////////////////////////////////////////////////////////////
// Performs a hermite spline interpolation, using the specified 3D vectors
template <class T>
void CMVec4Hermite( CMVector4D<T> & vOut,
                    const CMVector4D<T> & v1,
                    const CMVector4D<T> & vT1,
                    const CMVector4D<T> & v2,
                    const CMVector4D<T> & vT2,
                    T weight, T tension, T bias);

///////////////////////////////////////////////////////////////////////////////
// Performs a quadratic interpolation, using the specified 4D vectors
template <class T>
void CMVec4Quadratic(CMVector4D<T> & vOut,
                    const CMVector4D<T> & v1,
                    const CMVector4D<T> & v2,
                    const CMVector4D<T> & v3, T factor);

///////////////////////////////////////////////////////////////////////////////
// Performs a cubic interpolation, using the specified 4D vectors
template <class T>
void CMVec4Cubic(CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2,
                const CMVector4D<T> & v3,
                const CMVector4D<T> & v4, T factor);

///////////////////////////////////////////////////////////////////////////////
// Performs a square interpolation, using two 4D vectors
template <class T>
void CMVec4Sqrp(CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2, T factor);

///////////////////////////////////////////////////////////////////////////////
// Performs a cosinus interpolation, using two 4D vectors
template <class T>
void CMVec4Cosrp(CMVector4D<T> & vOut,
                const CMVector4D<T> & v1,
                const CMVector4D<T> & v2, T factor);



///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// 4D vector structure constructors

// Default constructor
template <class T>
inline CMVector4D<T>::CMVector4D()
{
    x = T(0);
    y = T(0);
    z = T(0);
    w = T(0);
}

///////////////////////////////////////////////////////////////////////////////
// Getting pointer to array
template <class T>
inline CMVector4D<T>::CMVector4D(const T  * const p)
{
    CMVector4D();
    if(p)
    {
        x = p[0];
        y = p[1];
        z = p[2];
        w = p[3];
    }
}

///////////////////////////////////////////////////////////////////////////////
// Getting parameters
template <class T>
inline CMVector4D<T>::CMVector4D(T x, T y, T z, T w)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

///////////////////////////////////////////////////////////////////////////////
// Copy constructor
template <class T>
inline CMVector4D<T>::CMVector4D(const CMVector4D & v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    w = v.w;
}

///////////////////////////////////////////////////////////////////////////////
// Getting 3D vector
template <class T>
inline CMVector4D<T>::CMVector4D(const CMVector3D<T> & v)

{
    x = v.x;
    y = v.y;
    z = v.z;
    w = T(0.0);
}

///////////////////////////////////////////////////////////////////////////////
// Getting 2D vector
template <class T>
inline CMVector4D<T>::CMVector4D(const CMVector2D<T> & v)

{
    x = v.x;
    y = v.y;
    z = T(0.0);
    w = T(0.0);
}

///////////////////////////////////////////////////////////////////////////////
//Cast to pointer  (T *)
template <class T>
inline CMVector4D<T>::operator T * () const
{
    return (T *) this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T>::operator const T * () const
{
    return (const T *) this;
}

///////////////////////////////////////////////////////////////////////////////
//Cast to CMVector3D type
template <class T>
inline CMVector4D<T>::operator CMVector3D<T>()
{
    return CMVector3D<T>(x,y,z);
}

///////////////////////////////////////////////////////////////////////////////
// Cast to CMVector2D type
template <class T>
inline CMVector4D<T>::operator CMVector2D<T>()
{
    return CMVector2D<T>(x,y);
}

///////////////////////////////////////////////////////////////////////////////
// Assignment operators
template <class T>
inline CMVector4D<T> & CMVector4D<T>::operator = ( const CMVector4D & v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    w = v.w;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T> & CMVector4D<T>::operator += ( const CMVector4D & v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    w += v.w;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T> & CMVector4D<T>::operator -= ( const CMVector4D & v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    w -= v.w;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T> & CMVector4D<T>::operator *= ( const T & val)
{
    x *= val;
    y *= val;
    z *= val;
    w *= val;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T> & CMVector4D<T>::operator /= ( const T & val)
{
    ASSERT(!mnear0(val));

    x /= val;
    y /= val;
    z /= val;
    w /= val;

    return *this;
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline T CMVector4D<T>::operator [] (int32t i)
{
    ASSERT(i >= 0 && i <= 3);
    return v[i];
}

///////////////////////////////////////////////////////////////////////////////
// Unary operators
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T> CMVector4D<T>::operator + () const
{
    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// Negate
template <class T>
inline CMVector4D<T> CMVector4D<T>::operator - () const
{
    return CMVector4D(-x, -y, -z, -w);
}

///////////////////////////////////////////////////////////////////////////////
// Binary Operators
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T> CMVector4D<T>::operator + ( const CMVector4D & v) const
{
    return CMVector4D(x + v.x, y + v.y, z + v.z, w + v.w);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T> CMVector4D<T>::operator - ( const CMVector4D & v) const
{
    return CMVector4D(x - v.x, y - v.y, z - v.z, w - v.w);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T> CMVector4D<T>::operator * ( const T & val) const
{
    return CMVector4D(x*val, y*val, z*val, w*val);
}
///////////////////////////////////////////////////////////////////////////////
template <class T>
inline CMVector4D<T> CMVector4D<T>::operator / ( const T & val) const
{
    T inv = T(1.0)/val;
    ASSERT(!mnear0(inv));
    return CMVector4D(x*inv, y*inv, z*inv, w*inv);
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
inline bool CMVector4D<T>::operator == ( const CMVector4D & v) const
{
    return ((x == v.x) && (y == v.y) && (z == v.z) && (w == v.w)) ? (1) : (0);
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
inline bool CMVector4D<T>::operator != ( const CMVector4D & v) const
{
    return ((x == v.x) && (y == v.y) && (z == v.z) && (w == v.w)) ? (0) : (1);
}

///////////////////////////////////////////////////////////////////////////////
// Methods operating on CMVector4D class
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Adding two vectors
template <class T>
inline void CMVec4Add(  CMVector4D<T> & vOut,
                        const CMVector4D<T> & v1,
                        const CMVector4D<T> & v2)
{
    vOut.x = v1.x + v2.x;
    vOut.y = v1.y + v2.y;
    vOut.z = v1.z + v2.z;
    vOut.w = v1.w + v2.w;
}

///////////////////////////////////////////////////////////////////////////////
// Subtracting two vectors
template <class T>
inline void CMVec4Sub(CMVector4D<T> & vOut,
const CMVector4D<T> & v1, const CMVector4D<T> & v2)
{
    vOut.x = v1.x - v2.x;
    vOut.y = v1.y - v2.y;
    vOut.z = v1.z - v2.z;
    vOut.w = v1.w - v2.w;
}

///////////////////////////////////////////////////////////////////////////////
// Dot product between two 4D vectors
template <class T>
inline T CMVec4Dot(const CMVector4D<T> & v1, const CMVector4D<T> & v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;
}

///////////////////////////////////////////////////////////////////////////////
// Performs a linear interpolation between two 4D vectors
template <class T>
inline void CMVec4Lerp( CMVector4D<T> & vOut,
                        const CMVector4D<T> & v1,
                        const CMVector4D<T> & v2, T weight)
{
    T x = v1.x + ((v2.x - v1.x) * weight);
    T y = v1.y + ((v2.y - v1.y) * weight);
    T z = v1.z + ((v2.z - v1.z) * weight);
    T w = v1.w + ((v2.w - v1.w) * weight);

    vOut.x = x;
    vOut.y = y;
    vOut.z = z;
    vOut.w = w;
}

///////////////////////////////////////////////////////////////////////////////
// Returns a 4D vector that is made up of the largest components
// of two 4D vectors
template <class T>
inline void CMVec4Max(  CMVector4D<T> & vOut,
                        const CMVector4D<T> & v1,
                        const CMVector4D<T> & v2)
{
    vOut.x = (v1.x > v2.x) ? (v1.x) : (v2.x);
    vOut.y = (v1.y > v2.y) ? (v1.y) : (v2.y);
    vOut.z = (v1.z > v2.z) ? (v1.z) : (v2.z);
    vOut.w = (v1.w > v2.w) ? (v1.w) : (v2.w);
}

///////////////////////////////////////////////////////////////////////////////
// Returns a 4D vector that is made up of the smallest components
// of two 4D vectors
template <class T>
inline void CMVec4Min(  CMVector4D<T> & vOut,
                        const CMVector4D<T> & v1,
                        const CMVector4D<T> & v2)
{
    vOut.x = (v1.x < v2.x) ? (v1.x) : (v2.x);
    vOut.y = (v1.y < v2.y) ? (v1.y) : (v2.y);
    vOut.z = (v1.z < v2.z) ? (v1.z) : (v2.z);
    vOut.w = (v1.w < v2.w) ? (v1.w) : (v2.w);
}

///////////////////////////////////////////////////////////////////////////////
// Scales a 4D vector
template <class T>
inline void CMVec4Scale(CMVector4D<T> & vOut, const CMVector4D<T> & vIn, T scale)
{
    vOut.x = vIn.x * scale;
    vOut.y = vIn.y * scale;
    vOut.z = vIn.z * scale;
    vOut.w = vIn.w * scale;
}

///////////////////////////////////////////////////////////////////////////////
// Invert scale of 4D vector
template <class T>
inline void CMVec4ScaleInv( CMVector4D<T> & vOut,
                            const CMVector4D<T> & vIn, T scale)
{
    T invscale = T(1.0)/scale;
    ASSERT(!mnear0(invscale));

    vOut.x = vIn.x * invscale;
    vOut.y = vIn.y * invscale;
    vOut.z = vIn.z * invscale;
    vOut.w = vIn.w * invscale;
}

///////////////////////////////////////////////////////////////////////////////

#endif _CMVECTOR4D_
