///////////////////////////////////////////////////////////////////////////////////////
//  plane.h
//
//  Plane class template declaration
//
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _PLANE_H_
#define _PLANE_H_


///////////////////////////////////////////////////////////////////////////////////////
// Plane Class Declaration
template <class T>
class CMPlane
{
public:
    union{
        struct{
            T a,b,c,d;
        };
        T p[4];
    };

public:
    // Constructors
    CMPlane();
    CMPlane(const CMPlane<T> &plane);
    CMPlane(const CMVector3D<T> &normal, T offset);
    CMPlane(const CMVector4D<T> &v);
    CMPlane(const T* const ptr);
    CMPlane(T a, T b, T c, T d);
    CMPlane(const CMVector3D<T> &p1, const CMVector3D<T> &p2, const CMVector3D<T> &p3);
    CMPlane(const CMVector3D<T> &normal, const CMVector3D<T> &point);



    // Cast To Pointer  (T *)
    operator T* () const;
    operator const T* () const;

    T operator [] (int32t );

    // Cast To CMVector3D Type
    operator  CMVector3D<T>();
    operator  CMVector3D<T>&();
    operator  CMVector3D<T>*();

    // Cast To CMVector4D Type
    operator  CMVector4D<T>();
    operator  CMVector4D<T>*();

    // Assignment Operators
    CMPlane& operator =  ( const CMPlane &plane );

    // Unary Operators
    CMPlane operator + () const;
    CMPlane operator - () const;
};

// Inline Function

///////////////////////////////////////////////////////////////////////////////////////
// Invert Plane
template <class T>
void CMPlaneInvert(CMPlane<T> &vOut, const CMPlane<T> &vIn);


// Non Inline Functions

///////////////////////////////////////////////////////////////////////////////////////
// Builds Plane From Three Points
template <class T>
void CMPlaneSet(CMPlane<T> &plane, const CMVector3D<T> &p1, const CMVector3D<T> &p2,
                                    const CMVector3D<T> &p3);

///////////////////////////////////////////////////////////////////////////////////////
// Builds Plane From Normal And Offset Point
template <class T>
void CMPlaneSet(CMPlane<T> &plane, const CMVector3D<T> &normal, const CMVector3D<T> &point);

///////////////////////////////////////////////////////////////////////////////////////
// Distance From Plane To Point Located In Space
template <class T>
T CMPlaneDistance(const CMPlane<T> &plane, const CMVector3D<T> &point);

///////////////////////////////////////////////////////////////////////////////////////
// Normalizes The Plane Coefficients So That The Plane Normal Has Unit Length
template <class T>
void CMPlaneNormalize(CMPlane<T> &out, const CMPlane<T> &in);

///////////////////////////////////////////////////////////////////////////////////////
// Angle Of Intersection Of The Two Planes
template <class T>
T CMPlaneIntersectionAngle(const CMPlane<T> &plane1, const CMPlane<T> &plane2);

///////////////////////////////////////////////////////////////////////////////////////
// Finds Point Of Intersection Of Plane And Line
// Return False If Line Do Not Intersect A Plane (Is Parallel To Plane)
// dir Must Be Normalized Vector
template <class T>
bool CMPlaneLineIntersect(CMVector3D<T> &vOut, const CMPlane<T> &plane, const CMVector3D<T> &point,
                        const CMVector3D<T> &dir);

///////////////////////////////////////////////////////////////////////////////////////
// Calculates Points Of Intersection Line With Plane
// Constructs Line From Given Two Points
template <class T>
bool CMPlaneLineIntersectByPoints(CMVector3D<T> &vOut, const CMPlane<T> &plane, const CMVector3D<T> &point1,
                           const CMVector3D<T> &point2);

///////////////////////////////////////////////////////////////////////////////////////
// Finds Points Of Intersection Of Ray And Plane
// dir Must Be Normalized
template <class T>
bool CMPlaneRayIntersectPt(CMVector3D<T> &vOut, const CMPlane<T> &plane, const CMVector3D<T> &orig,
                             const CMVector3D<T> &dir, bool checkTwoSide = false);

///////////////////////////////////////////////////////////////////////////////////////
// Check If Ray Intersect Plane
template <class T>
bool CMIsRayIntersectPlane(const CMPlane<T> &plane, const CMVector3D<T> &orig,
                             const CMVector3D<T> &dir, bool checkTwoSide = false);

///////////////////////////////////////////////////////////////////////////////////////
// Finds Point Of Intersection Of Line Segment And Plane
template <class T>
bool CMPlaneSegmentIntersect(CMVector3D<T> &vOut, const CMPlane<T> &plane, const CMVector3D<T> &p1,
                                const CMVector3D<T> &p2);

///////////////////////////////////////////////////////////////////////////////////////
// Calculate 3D Point Of The Intersection Of Three Planes
// Returns False If There Were No Intersection Of Planes
template <class T>
bool CMPlaneIntersect(CMVector3D<T> &vOut, const CMPlane<T> &plane1, const CMPlane<T> &plane2,
                       const CMPlane<T> &plane3);

///////////////////////////////////////////////////////////////////////////////////////
// Orthogonal Projection 3D Vector On The Plane
template <class T>
void CMPlaneOrthoProject(CMVector3D<T> &vOut, const CMPlane<T> &plane, const CMVector3D<T> &point);

///////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//
//  File:       Plane.inl
//  Content:    Plane inline functions
//
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Plane Structure Constructors

// Default Constructor
template <class T>
M_INLINE CMPlane<T>::CMPlane()
{
    a = b = c = d = T(0);
}

// Copy Constructor
template <class T>
M_INLINE CMPlane<T>::CMPlane(const CMPlane<T> &plane)
{
    a = plane.a;
    b = plane.b;
    c = plane.c;
    d = plane.d;
}

// Building Plane From Three Points
template <class T>
M_INLINE CMPlane<T>::CMPlane(const CMVector3D<T> &p1, const CMVector3D<T> &p2, const CMVector3D<T> &p3)
{
    CMPlaneSet(*this,p1,p2,p3);
}

// Building Plane From Normal And Offsert Vector
template <class T>
M_INLINE CMPlane<T>::CMPlane(const CMVector3D<T> &normal, const CMVector3D<T> &point)
{
    CMPlaneSet(*this,normal,point);
}

// Building From Normal And Offset
template <class T>
M_INLINE CMPlane<T>::CMPlane(const CMVector3D<T> &v, T offset)
{
    a = v.x;
    b = v.y;
    c = v.z;
    d = offset;
}

// Bulding From 4D Vector By Treating 4D Vector As Plane Parameters
template <class T>
M_INLINE CMPlane<T>::CMPlane(const CMVector4D<T> &v)
{
    a = v.x;
    b = v.y;
    c = v.z;
    d = v.w;
}

// Getting Pointer To Array
template <class T>
M_INLINE CMPlane<T>::CMPlane(const T* const ptr)
{
    CMPlane();
    if(ptr)
    {
        a = ptr[0];
        b = ptr[1];
        c = ptr[2];
        d = ptr[3];
    }
}

// Getting Parameters Separatly
template <class T>
M_INLINE CMPlane<T>::CMPlane(T _a, T _b, T _c, T _d):a(_a),b(_b),c(_c),d(_d)
{

}

////////////////////////////////////////////////////////////////////////////
// Assignment Operators
template <class T>
M_INLINE CMPlane<T>& CMPlane<T>::operator = ( const CMPlane &plane )
{
    a = plane.a;
    b = plane.b;
    c = plane.c;
    d = plane.d;

    return *this;
}

////////////////////////////////////////////////////////////////////////////
//Cast To Pointer  (T *)
template <class T>
M_INLINE CMPlane<T>::operator T* () const
{
    return (T*) this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMPlane<T>::operator const T*() const
{
    return (const T*) this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE T CMPlane<T>:: operator [] (int32t i)
{
    ASSERT(i >= 0 && i <= 3);

    return p[i];
}

////////////////////////////////////////////////////////////////////////////
// Cast To CMVector3D Type
template <class T>
M_INLINE CMPlane<T>::operator  CMVector3D<T>()
{
    return (CMVector3D<T>&)p[0];//return CMVector3D<T>(x,y,z);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMPlane<T>::operator  CMVector3D<T>&()
{
    return (CMVector3D<T>&)p[0];
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMPlane<T>::operator CMVector3D<T>*()
{
    return (CMVector3D<T>*)&p[0];
}

////////////////////////////////////////////////////////////////////////////
// Cast To CMVector4D Type
template <class T>
M_INLINE CMPlane<T>::operator CMVector4D<T>()
{
    return CMVector4D<T>(&p[0]);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMPlane<T>::operator  CMVector4D<T>*()
{
    return (CMVector4D<T>*)&p[0];
}

////////////////////////////////////////////////////////////////////////////
// Unary Operators
template <class T>
M_INLINE CMPlane<T>  CMPlane<T>::operator + () const
{
    return CMPlane(a,b,c,d);
}
////////////////////////////////////////////////////////////////////////////
// Negate
template <class T>
M_INLINE CMPlane<T> CMPlane<T>::operator - () const
{
    return CMPlane<T>(-a,-b,-c,-d);
}





///////////////////////////////////////////////////////////////////////////////////////

#endif //_PLANE_H_
