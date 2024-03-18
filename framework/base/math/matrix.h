//-----------------------------------------------------------------------------
//  matrix.h
//
//  Matrix classes templates declarations
//
//-----------------------------------------------------------------------------

#ifndef _MATRIX_H_
#define _MATRIX_H_

///////////////////////////////////////////////////////////////////////////////
// 2x2 Matrix Class Declaration
template <class T>
class CMMatrix22
{
public:
    union{
        struct{
            T _11, _12;
            T _21, _22;
        };
        T m[2][2];
    };

public:
    //Constructors
    CMMatrix22();
    CMMatrix22(const T* const p);
    CMMatrix22(const CMMatrix22& mtx);
    CMMatrix22(CMVector2D<T> row1, CMVector2D<T> row2);
    CMMatrix22(T _11, T _12,
               T _21, T _22);

    //Cast Operators
    operator T* () const;
    operator const T* () const;

    // Access Grants
    CMVector2D<T> operator [] (int32t );
    CMVector2D<T> operator [] (int32t ) const;

    // Assignment Operators
    CMMatrix22& operator =  ( const CMMatrix22 & );
    CMMatrix22& operator += ( const CMMatrix22 & );
    CMMatrix22& operator -= ( const CMMatrix22 & );
    CMMatrix22& operator *= ( const CMMatrix22 & );
    CMMatrix22& operator *= ( const T & );
    CMMatrix22& operator /= ( const T & );

    // Unary Operators
    CMMatrix22 operator + () const;
    CMMatrix22 operator - () const;

    // Binary Operators
    CMMatrix22 operator + ( const CMMatrix22 & ) const;
    CMMatrix22 operator - ( const CMMatrix22 & ) const;
    CMMatrix22 operator * ( const CMMatrix22 & ) const;
    CMMatrix22 operator * ( const T & ) const;
    CMMatrix22 operator / ( const T & ) const;

    bool operator == ( const CMMatrix22 & ) const;
    bool operator != ( const CMMatrix22 & ) const;
};

//---------------------------
// Inline
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Set Identity Matrix
template <class T>
void CMMtx22Identity(CMMatrix22<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Is Identity
template <class T>
bool CMMtx22IsIdentity(const CMMatrix22<T> &M);

//---------------------------
// Non-Inline
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Adding Two Matrices
template <class T>
void CMMtx22Add(CMMatrix22<T> &Out, const CMMatrix22<T> &Mtx1, const CMMatrix22<T> &Mtx2);

///////////////////////////////////////////////////////////////////////////////
// Subtracting 2x2 Matrices
template <class T>
void CMMtx22Sub(CMMatrix22<T> &Out, const CMMatrix22<T> &Mtx1, const CMMatrix22<T> &Mtx2);

///////////////////////////////////////////////////////////////////////////////
// Inverse Matrix
template <class T>
void CMMtx22Inverse(CMMatrix22<T> &Out, const CMMatrix22<T> &vM);

///////////////////////////////////////////////////////////////////////////////
// Transpose Matrix
template <class T>
void CMMtx22Transpose(CMMatrix22<T> &Out, const CMMatrix22<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Create Scale Matrix
template <class T>
void CMMtx22Scale(CMMatrix22<T> &M, T x, T y);

///////////////////////////////////////////////////////////////////////////////
// Create Rotation Matrix
template <class T>
void CMMtx22Rotate(CMMatrix22<T> &M, T rad);

///////////////////////////////////////////////////////////////////////////////
// Transform 2D Vector By 2x2 Matrix
template <class T>
void CMVec2Transform(CMVector2D<T> &Out, const CMVector2D<T> &V, const CMMatrix22<T> &M);


///////////////////////////////////////////////////////////////////////////////
// 3x3 Matrix Class Declaration
template <class T>
class CMMatrix33
{
public:
    union{
        struct{
            T _11, _12, _13;
            T _21, _22, _23;
            T _31, _32, _33;
        };
        T m[3][3];
    };

public:
    //Constructors
    CMMatrix33();
    CMMatrix33(const T* const p);
    CMMatrix33(const CMMatrix33& mtx);
    CMMatrix33(const CMVector3D<T> &row1, const CMVector3D<T> &row2, const CMVector3D<T> &row3);
    CMMatrix33(T e11, T e12, T e13,
               T e21, T e22, T e23,
               T e31, T e32, T e33);

    // Cast Operators
    operator T* () const;
    operator const T* () const;

    // Access Grants
    CMVector3D<T>& operator [] (int32t );
    CMVector3D<T>& operator [] (int32t ) const;

    T& operator () (uint32t i, uint32t j);
    T  operator () (uint32t i, uint32t j) const;

    // Assignment Operators
    CMMatrix33& operator =  ( const CMMatrix33 & );
    CMMatrix33& operator += ( const CMMatrix33 & );
    CMMatrix33& operator -= ( const CMMatrix33 & );
    CMMatrix33& operator *= ( const CMMatrix33 & );
    CMMatrix33& operator *= ( const T & );
    CMMatrix33& operator /= ( const T & );

    // Unary Operators
    CMMatrix33 operator + () const;
    CMMatrix33 operator - () const;

    // Binary Operators
    CMMatrix33 operator + ( const CMMatrix33 & ) const;
    CMMatrix33 operator - ( const CMMatrix33 & ) const;
    CMMatrix33 operator * ( const CMMatrix33 & ) const;
    CMMatrix33 operator * ( const T & ) const;
    CMMatrix33 operator / ( const T & ) const;

    bool operator == ( const CMMatrix33 & ) const;
    bool operator != ( const CMMatrix33 & ) const;
};

//---------------------------
// Inline
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Set Identity Matrix
template <class T>
void CMMtx33Identity(CMMatrix33<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Is Identity
template <class T>
bool CMMtx33IsIdentity(const CMMatrix33<T> &M);

//---------------------------
// Non-Inline
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Adding Two Matrices
template <class T>
void CMMtx33Add(CMMatrix33<T> &Out, const CMMatrix33<T> &Mtx1, const CMMatrix33<T> &Mtx2);

///////////////////////////////////////////////////////////////////////////////
// Subtracting 3x3 Matrices
template <class T>
void CMMtx33Sub(CMMatrix33<T> &Out, const CMMatrix33<T> &Mtx1, const CMMatrix33<T> &Mtx2);

///////////////////////////////////////////////////////////////////////////////
// Multiplying Two 3x3 Matrices
template <class T>
void CMMtx33Mult(CMMatrix33<T> &Out, const CMMatrix33<T> &Mtx1, const CMMatrix33<T> &Mtx2);

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Inverse Matrix
template <class T>
bool CMMtx33Inverse(CMMatrix33<T> &Out, const CMMatrix33<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Transpose Matrix
template <class T>
void CMMtx33Transpose(CMMatrix33<T> &Out, const CMMatrix33<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Create Scale Matrix In 2D Space
template <class T>
void CMMtx33Scale2D(CMMatrix33<T> &M, T x, T y);

///////////////////////////////////////////////////////////////////////////////
// Create Inverse Scale Values Matrix In 2D Space
template <class T>
void CMMtx33ScaleInv2D(CMMatrix33<T> &M, T x, T y);

///////////////////////////////////////////////////////////////////////////////
// Create Translation Matrix In 2D Space
template <class T>
void CMMtx33Translate2D(CMMatrix33<T> &M, T x, T y);

///////////////////////////////////////////////////////////////////////////////
// Create Scale Matrix
template <class T>
void CMMtx33Scale(CMMatrix33<T> &M, T x, T y, T z);

///////////////////////////////////////////////////////////////////////////////
// Create Inverse Scale Values Matrix
template <class T>
void CMMtx33ScaleInv(CMMatrix33<T> &M, T x, T y, T z);

///////////////////////////////////////////////////////////////////////////////
// Create Rotation Matrix In 2D Space
template <class T>
void CMMtx33Rotate2D(CMMatrix33<T> &M, T rad);

///////////////////////////////////////////////////////////////////////////////
// Create Rotation Matrix Around The X Axis
template <class T>
void CMMtx33RotateX(CMMatrix33<T> &M, T rad);

///////////////////////////////////////////////////////////////////////////////
// Create Rotation Matrix Around The Y Axis
template <class T>
void CMMtx33RotateY(CMMatrix33<T> &M, T rad);

///////////////////////////////////////////////////////////////////////////////
// Create Rotation Matrix Around The Z Axis
template <class T>
void CMMtx33RotateZ(CMMatrix33<T> &M, T rad);

///////////////////////////////////////////////////////////////////////////////
// Creates A Matrix That Rotates Around An Arbitrary Axis
template <class T>
void CMMtx33Rotate(CMMatrix33<T> &M, const CMVector3D<T> &vNormVec, T rad);

///////////////////////////////////////////////////////////////////////////////
// Transform 3D Vector By 3x3 Matrix
template <class T>
void CMVec3Mult(CMVector3D<T> &Out, const CMVector3D<T> &V, const CMMatrix33<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Transform 2D Vector By 2x2 Matrix
template <class T>
void CMVec2Transform(CMVector2D<T> &Out, const CMVector2D<T> &v, const CMMatrix33<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Decomposing Transformation 3x3 Matrix
template <class T>
void CMMtx33Decompose(const CMMatrix33<T> &M, CMVector3D<T> &vScale, CMVector3D<T> &vRotation);


///////////////////////////////////////////////////////////////////////////////
// 4x4 Matrix Class Declaration
template <class T>
class CMMatrix44
{
public:
    union{
        struct{
            T        _11, _12, _13, _14;
            T        _21, _22, _23, _24;
            T        _31, _32, _33, _34;
            T        _41, _42, _43, _44;
        };
        T m[4][4];
    };

public:
    //Constructors
    CMMatrix44();
    CMMatrix44(const T* const p);
    CMMatrix44(const CMMatrix44& mtx);
    CMMatrix44(const CMMatrix33<T>& mtx);
    CMMatrix44(const CMVector4D<T>& row1,const CMVector4D<T>& row2,const CMVector4D<T>& row3,const CMVector4D<T>& row4);
    CMMatrix44(const CMVector3D<T>& row1,const CMVector3D<T>& row2,const CMVector3D<T>& row3,const CMVector3D<T>& row4);
    CMMatrix44(const CMVector3D<T>& row1,const CMVector3D<T>& row2,const CMVector3D<T>& row3);
    CMMatrix44(T _11, T _12, T _13, T _14,
               T _21, T _22, T _23, T _24,
               T _31, T _32, T _33, T _34,
               T _41, T _42, T _43, T _44 );

    CMMatrix44(T _11, T _12, T _13,
               T _21, T _22, T _23,
               T _31, T _32, T _33 );

    //Cast Operators
    operator T* () const;
    operator const T* () const;

    // Access Grants
    CMVector4D<T>& operator [] (int32t );
    CMVector4D<T>& operator [] (int32t ) const;

    T& operator () (int32t i, int32t j);
    T  operator () (int32t i, int32t j) const;


    // Assignment Operators
    CMMatrix44& operator =  ( const CMMatrix44 & );
    CMMatrix44& operator += ( const CMMatrix44 & );
    CMMatrix44& operator -= ( const CMMatrix44 & );
    CMMatrix44& operator *= ( const CMMatrix44 & );
    CMMatrix44& operator *= ( const T & );
    CMMatrix44& operator /= ( const T & );

    // Unary Operators
    CMMatrix44 operator + () const;
    CMMatrix44 operator - () const;

    // Binary Operators
    CMMatrix44 operator + ( const CMMatrix44 & ) const;
    CMMatrix44 operator - ( const CMMatrix44 & ) const;
    CMMatrix44 operator * ( const CMMatrix44 & ) const;
    CMMatrix44 operator * ( const T & ) const;
    CMMatrix44 operator / ( const T & ) const;

    bool operator == ( const CMMatrix44 & ) const;
    bool operator != ( const CMMatrix44 & ) const;
};

//---------------------------
// Inline
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Set Identity Matrix
template <class T>
void CMMtx44Identity(CMMatrix44<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Is Identity
template <class T>
bool CMMtx44IsIdentity(const CMMatrix44<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Setting Matrix Row
template <class T>
void CMMtx44Set(int row, CMMatrix44<T> &M, const CMVector3D<T> &vec, T w = T(1));

///////////////////////////////////////////////////////////////////////////////
// Setting Matrix Row
template <class T>
void CMMtx44Set(int row, CMMatrix44<T> &M, const CMVector4D<T> &vec);

//---------------------------
// Non-Inline
//---------------------------

///////////////////////////////////////////////////////////////////////////////
// Adding Two Matrices
template <class T>
void CMMtx44Add(CMMatrix44<T> &Out, const CMMatrix44<T> &Mtx1, const CMMatrix44<T> &Mtx2);

///////////////////////////////////////////////////////////////////////////////
// Subtracting 4x4 Matrices
template <class T>
void CMMtx44Sub(CMMatrix44<T> &Out, const CMMatrix44<T> &Mtx1, const CMMatrix44<T> &Mtx2);

///////////////////////////////////////////////////////////////////////////////
// Multiplying Two 4x4 Matrices
template <class T>
void CMMtx44Mult(CMMatrix44<T> &Out, const CMMatrix44<T> &Mtx1, const CMMatrix44<T> &Mtx2);

///////////////////////////////////////////////////////////////////////////////
// Inverse Matrix
template <class T>
bool CMMtx44Inverse(CMMatrix44<T> &Out, const CMMatrix44<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Affine Inverse Matrix
template <class T>
bool CMMtx44InverseAffine(CMMatrix44<T> &Out, const CMMatrix44<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Transpose Matrix
template <class T>
void CMMtx44Transpose(CMMatrix44<T> &Out, const CMMatrix44<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Create Translate Matrix
template <class T>
void CMMtx44Translate(CMMatrix44<T> &M, T x, T y, T z);

///////////////////////////////////////////////////////////////////////////////
// Create Scale Matrix
template <class T>
void CMMtx44Scale(CMMatrix44<T> &M, T x, T y, T z);

///////////////////////////////////////////////////////////////////////////////
// Create Rotation Matrix Around The X Axis
template <class T>
void CMMtx44RotateX(CMMatrix44<T> &M, T rad);

///////////////////////////////////////////////////////////////////////////////
// Create Rotation Matrix Around The Y Axis
template <class T>
void CMMtx44RotateY(CMMatrix44<T> &M, T rad);

///////////////////////////////////////////////////////////////////////////////
// Create Rotation Matrix Around The Z Axis
template <class T>
void CMMtx44RotateZ(CMMatrix44<T> &M, T rad);

///////////////////////////////////////////////////////////////////////////////
// Creates A Matrix That Rotates Around An Arbitrary Axis
template <class T>
void CMMtx44Rotate(CMMatrix44<T> &M, const CMVector3D<T> &vNormVec, T rad);

///////////////////////////////////////////////////////////////////////////////
// Transform 4D Vector By 4x4 Matrix
template <class T>
void CMVec4Mult(CMVector4D<T> &Out, const CMVector4D<T> &v, const CMMatrix44<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Transform 4D Vector By 4x4 Matrix Multiplying Only By First Three Rows
template <class T>
void CMVec4Mult3(CMVector4D<T> &Out, const CMVector4D<T> &v, const CMMatrix44<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Transform 4D Vector By 3x3 Matrix
template <class T>
void CMVec4Mult(CMVector4D<T> &Out, const CMVector4D<T> &v, const CMMatrix33<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Transform 3D Vector By 4x4 Matrix
template <class T>
void CMVec3Mult(CMVector3D<T> &Out, const CMVector3D<T> &v, const CMMatrix44<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Transform 3D Vector By 4x4 Matrix Multiplying Only By First Three Rows
template <class T>
void CMVec3Mult3(CMVector3D<T> &Out, const CMVector3D<T> &v, const CMMatrix44<T> &M);

///////////////////////////////////////////////////////////////////////////////
// Builds A Transformation Matrix. NULL Arguments Are Treated As Default Transformations.
template <class T>
void CMMtx44Transformation(CMMatrix44<T> &M, const CMVector3D<T> &vTranslation, const CMVector3D<T> &vRotation,
                     const CMVector3D<T> &vScale);

///////////////////////////////////////////////////////////////////////////////
// Decomposing Transformation 4x4 Matrix
template <class T>
void CMMtx44Decompose(const CMMatrix44<T> &M, CMVector3D<T> &vScale, CMVector3D<T> &vRotation,
                      CMVector3D<T> &vTranslation);

///////////////////////////////////////////////////////////////////////////////
// Create A LookAt Matrix (Right-Handed)
template <class T>
void CMMtx44LookAtRH(CMMatrix44<T> &M, const CMVector3D<T> &vEye, const CMVector3D<T> &vAt,
                     const CMVector3D<T> &vUp);

///////////////////////////////////////////////////////////////////////////////
// Create A LookAt Matrix (Left-Handed)
template <class T>
void CMMtx44LookAtLH(CMMatrix44<T> &M, const CMVector3D<T> &vEye, const CMVector3D<T> &vAt,
                     const CMVector3D<T> &vUp);

///////////////////////////////////////////////////////////////////////////////
// Create A Perspective Projection Matrix (Right-Handed)
template <class T>
void CMMtx44PerspectiveRH(CMMatrix44<T> &M, T w, T h, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Create A Perspective Projection Matrix (Left-Handed)
template <class T>
void CMMtx44PerspectiveLH(CMMatrix44<T> &M, T w, T h, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Create A Perspective Projection Matrix (Right-Handed)
template <class T>
void CMMtx44PerspectiveFOVRH(CMMatrix44<T> &M, T fov, T aspect, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Create A Perspective Projection Matrix (Left-Handed)
template <class T>
void CMMtx44PerspectiveFOVLH(CMMatrix44<T> &M, T fov, T aspect, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Create A Perspective Projection Matrix (Right-Handed)
template <class T>
void CMMtx44PerspectiveOffCenterRH(CMMatrix44<T> &M, T l, T r, T b, T t, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Create A Perspective Projection Matrix (Left-Handed)
template <class T>
void CMMtx44PerspectiveOffCenterLH(CMMatrix44<T> &M, T l, T r, T b, T t, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Create An Ortho Projection Matrix (Right-Handed)
template <class T>
void CMMtx44OrthoRH(CMMatrix44<T> &M, T w, T h, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Create An Ortho Projection Matrix (Left-Handed)
template <class T>
void CMMtx44OrthoLH(CMMatrix44<T> &M, T w, T h, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Create An Ortho Projection Matrix (Right-Handed)
template <class T>
void CMMtx44OrthoOffCenterRH(CMMatrix44<T> &M, T l, T r, T b, T t, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Create An Ortho Projection Matrix (Left-Handed)
template <class T>
void CMMtx44OrthoOffCenterLH(CMMatrix44<T> &M, T l, T r, T b, T t, T zn, T zf);

///////////////////////////////////////////////////////////////////////////////
// Creates Matrix That Reflects Coordinate System About The Plane
// Plane Must Be Normalized
template <class T>
void CMMtx44Reflect(CMMatrix44<T> &M, const CMPlane<T> &plane);

///////////////////////////////////////////////////////////////////////////////
// Creates Matrix That Projects Geometry Into A Plane
// Plane Must Be Normalized
template <class T>
void CMMtx44Project(CMMatrix44<T> &M, const CMVector4D<T> &light, const CMPlane<T> &plane);
///////////////////////////////////////////////////////////////////////////////
template <class T>
void CMMtx44Project(CMMatrix44<T> &M, const CMVector3D<T> &light, const CMPlane<T> &plane);

///////////////////////////////////////////////////////////////////////////////







///////////////////////////////////////////////////////////////////////////////////////
// 3x3 Matrix Class Inline Functions
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// 3x3 Matrix Class Constructors

// Default Constructor
template <class T>
M_INLINE CMMatrix33<T>::CMMatrix33()
{
    memset(&_11,0,sizeof(T)*9);
}

// Getting Pointer To Array
template <class T>
M_INLINE CMMatrix33<T>::CMMatrix33(const T * const p)
{
    ASSERT(p != NULL);

    memcpy(&_11 ,p, sizeof(T)*9);
}

// Copy Constructor
template <class T>
M_INLINE CMMatrix33<T>::CMMatrix33(const CMMatrix33<T> &mtx)
{
    memcpy(&_11, &mtx._11, sizeof(T)*9);
}

// Constructing Matrix By Getting 3D Vectors
template <class T>
M_INLINE CMMatrix33<T>::CMMatrix33(const CMVector3D<T> &row1, const CMVector3D<T> &row2,
                                   const CMVector3D<T> &row3)
{
    memcpy(&_11, &row1, sizeof(T)*3);
    memcpy(&_21, &row2, sizeof(T)*3);
    memcpy(&_31, &row3, sizeof(T)*3);
}

// Building Matrix By Getting All Paremeters
template <class T>
M_INLINE CMMatrix33<T>::CMMatrix33(T _11, T _12, T _13, T _21, T _22, T _23, T _31, T _32, T _33)
{
    this->_11 = _11;    this->_12 = _12;    this->_13 = _13;
    this->_21 = _21;    this->_22 = _22;    this->_23 = _23;
    this->_31 = _31;    this->_32 = _32;    this->_33 = _33;
}

// Cast To Pointer  (T *)
template <class T>
M_INLINE CMMatrix33<T>::operator T* () const
{
    return (T*) &m[0][0];
}
///////////////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T>::operator const T* () const
{
    return (const T*) &m[0][0];
}

// Access Grants
template <class T>
M_INLINE CMVector3D<T>& CMMatrix33<T>::operator [] (int32t i)
{
    ASSERT(i >= 0 && i <= 2);

    return (CMVector3D<T> &) m[i][0];
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMVector3D<T>& CMMatrix33<T>::operator [] (int32t i) const
{
    ASSERT(i >= 0 && i <= 2);

    return (CMVector3D<T> &) m[i][0];
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE T& CMMatrix33<T>::operator () (uint32t i, uint32t j)
{
    ASSERT(i>=0 && i<=2);
    ASSERT(j>=0 && j<=2);

    return m[i][j];
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE T CMMatrix33<T>::operator () (uint32t i, uint32t j) const
{
    ASSERT(i>=0 && i<=2);
    ASSERT(j>=0 && j<=2);

    return m[i][j];
}


////////////////////////////////////////////////////////////////////////////
// Assignment Operators
template <class T>
M_INLINE CMMatrix33<T>& CMMatrix33<T>::operator = (const CMMatrix33<T> &mtx)
{
    memcpy(&_11, &mtx._11, sizeof(T)*9);

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T>& CMMatrix33<T>::operator += (const CMMatrix33<T> &mtx)
{
    _11 += mtx._11; _12 += mtx._12; _13 += mtx._13;
    _21 += mtx._21; _22 += mtx._22; _23 += mtx._23;
    _31 += mtx._31; _32 += mtx._32; _33 += mtx._33;

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T>& CMMatrix33<T>::operator -= (const CMMatrix33<T> &mtx)
{
    _11 -= mtx._11; _12 -= mtx._12; _13 -= mtx._13;
    _21 -= mtx._21; _22 -= mtx._22; _23 -= mtx._23;
    _31 -= mtx._31; _32 -= mtx._32; _33 -= mtx._33;

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T>& CMMatrix33<T>::operator *= (const T & val)
{
    _11 *= val; _12 *= val; _13 *= val;
    _21 *= val; _22 *= val; _23 *= val;
    _31 *= val; _32 *= val; _33 *= val;

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T>& CMMatrix33<T>::operator *= (const CMMatrix33<T>& mtx)
{
    CMMtx33Mult(*this,*this,mtx);

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T>& CMMatrix33<T>::operator /= (const T & val)
{
    ASSERT(!mnear0(val));

    T inv = T(1)/val;

    _11 *= inv; _12 *= inv; _13 *= inv;
    _21 *= inv; _22 *= inv; _23 *= inv;
    _31 *= inv; _32 *= inv; _33 *= inv;

    return *this;
}



////////////////////////////////////////////////////////////////////////////
// Unary Operators
template <class T>
M_INLINE CMMatrix33<T> CMMatrix33<T>::operator + () const
{
    return *this;
}

// Negate
template <class T>
M_INLINE CMMatrix33<T> CMMatrix33<T>::operator - () const
{
    return CMMatrix33(  -_11, -_12, -_13,
                        -_21, -_22, -_23,
                        -_31, -_32, -_33);
}

// Binary Operators
template <class T>
M_INLINE CMMatrix33<T> CMMatrix33<T>::operator + ( const CMMatrix33 &mtx ) const
{
    return CMMatrix33(  _11 + mtx._11, _12 + mtx._12, _13 + mtx._13,
                        _21 + mtx._21, _22 + mtx._22, _23 + mtx._23,
                        _31 + mtx._31, _32 + mtx._32, _33 + mtx._33);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T> CMMatrix33<T>::operator - ( const CMMatrix33 &mtx) const
{
    return CMMatrix33(  _11 - mtx._11, _12 - mtx._12, _13 - mtx._13,
                        _21 - mtx._21, _21 - mtx._22, _23 - mtx._23,
                        _31 - mtx._31, _32 - mtx._32, _33 - mtx._33);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T> CMMatrix33<T>::operator * ( const T & val ) const
{
    return CMMatrix33(  _11*val, _12*val, _13*val,
                        _21*val, _22*val, _23*val,
                        _31*val, _32*val, _33*val);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T> CMMatrix33<T>::operator * ( const CMMatrix33<T> &mtx ) const
{
    CMMatrix33 mtxT;
    CMMtx33Mult(mtxT,*this,mtx);
    return mtxT;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix33<T> CMMatrix33<T>::operator / ( const T & val ) const
{
    ASSERT(!mnear0(val));
    //T inv = T(1)/val;
    return CMMatrix33(  _11/val, _12/val, _13/val,
                        _21/val, _22/val, _23/val,
                        _31/val, _32/val, _33/val);
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE bool CMMatrix33<T>::operator == ( const CMMatrix33<T> & mtx ) const
{
    return 0 == memcmp(this, &mtx, sizeof(T)*9);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE bool CMMatrix33<T>::operator != ( const CMMatrix33<T> & mtx ) const
{
    return 0 != memcmp(this, &mtx, sizeof(T)*9);
}

///////////////////////////////////////////////////////////////////////////////////////
// Set Identity Matrix
template <class T>
void CMMtx33Identity(CMMatrix33<T> &mtx)
{
    memset(&mtx._11,0,sizeof(T)*9);

    mtx._11 = mtx._22 = mtx._33 = T(1);
}

///////////////////////////////////////////////////////////////////////////////////////
// Is Identity
template <class T>
bool CMMtx33IsIdentity(const CMMatrix33<T> &mtx)
{
    return  mtx.m[0][0] == T(1) && mtx.m[0][1] == T(0) && mtx.m[0][2] == T(0) &&
            mtx.m[1][0] == T(0) && mtx.m[1][1] == T(1) && mtx.m[1][2] == T(0) &&
            mtx.m[2][0] == T(0) && mtx.m[2][1] == T(0) && mtx.m[2][2] == T(1);
}

///////////////////////////////////////////////////////////////////////////////////////
// 4x4 Matrix Class Inline Functions
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// 4x4 Matrix Class Constructors

// Default Constructor
template <class T>
M_INLINE CMMatrix44<T>::CMMatrix44()
{
    memset(&this->_11,0,sizeof(T)*16);
}

// Getting Pointer To Array
template <class T>
M_INLINE CMMatrix44<T>::CMMatrix44(const T *const p)
{
    ASSERT(p != NULL);

    memcpy(&_11, p, sizeof(T)*16);
}

// Copy Constructor
template <class T>
M_INLINE CMMatrix44<T>::CMMatrix44(const CMMatrix44<T> &mtx)
{
    memcpy(&_11, &mtx._11, sizeof(T)*16);
}

// Constructing Matrix From 3x3 Matrix
template <class T>
M_INLINE CMMatrix44<T>::CMMatrix44(const CMMatrix33<T> &mtx)
{
    CMMatrix44();
    _44 = T(1);

    memcpy(&_11, &mtx._11, sizeof(T)*3);
    memcpy(&_21, &mtx._21, sizeof(T)*3);
    memcpy(&_31, &mtx._31, sizeof(T)*3);
}

// Constructing Matrix From 4D Vectors
template <class T>
M_INLINE CMMatrix44<T>::CMMatrix44(const CMVector4D<T>& row1,const CMVector4D<T>& row2,const CMVector4D<T>& row3,
                                   const CMVector4D<T>& row4)
{
    memcpy(&_11, &row1.x, sizeof(T)*4);
    memcpy(&_21, &row2.x, sizeof(T)*4);
    memcpy(&_31, &row3.x, sizeof(T)*4);
    memcpy(&_41, &row4.x, sizeof(T)*4);
}

// Constructing Matrix From Four 3D Vectors
template <class T>
M_INLINE CMMatrix44<T>::CMMatrix44(const CMVector3D<T>& row1, const CMVector3D<T>& row2, const CMVector3D<T>& row3,
                                   const CMVector3D<T>& row4)
{
    CMMatrix44();
    _44 = T(1);

    memcpy(&_11, &row1.x, sizeof(T)*3);
    memcpy(&_21, &row2.x, sizeof(T)*3);
    memcpy(&_31, &row3.x, sizeof(T)*3);
    memcpy(&_41, &row4.x, sizeof(T)*3);
}

// Constructing Matrix From 3D Vectors
template <class T>
M_INLINE CMMatrix44<T>::CMMatrix44(const CMVector3D<T>& row1, const CMVector3D<T>& row2, const CMVector3D<T>& row3)
{
    CMMatrix44();
    _44 = T(1);

    memcpy(&_11, &row1.x, sizeof(T)*3);
    memcpy(&_21, &row2.x, sizeof(T)*3);
    memcpy(&_31, &row3.x, sizeof(T)*3);
}

// Building Matrix By Getting 9 Params Of 3x3 Matrix
template <class T>
M_INLINE CMMatrix44<T>::CMMatrix44(T _11, T _12, T _13, T _21, T _22, T _23, T _31, T _32, T _33)
{
    CMMatrix44();
    m[3][3] = T(1);

    m[0][0] = _11;  m[0][1] = _12;  m[0][2] = _13;
    m[1][0] = _21;  m[1][1] = _22;  m[1][2] = _23;
    m[2][0] = _31;  m[2][1] = _32;  m[2][2] = _33;
}

// Building Matrix By Getting 16 Params Of 4x4 Matrix
template <class T>
M_INLINE CMMatrix44<T>::CMMatrix44(T _11, T _12, T _13, T _14, T _21, T _22, T _23, T _24,
                                   T _31, T _32, T _33, T _34, T _41, T _42, T _43, T _44)
{
    m[0][0] = _11;  m[0][1] = _12;  m[0][2] = _13;  m[0][3] = _14;
    m[1][0] = _21;  m[1][1] = _22;  m[1][2] = _23;  m[1][3] = _24;
    m[2][0] = _31;  m[2][1] = _32;  m[2][2] = _33;  m[2][3] = _34;
    m[3][0] = _41;  m[3][1] = _42;  m[3][2] = _43;  m[3][3] = _44;
}

// Cast To Pointer  (T *)
template <class T>
M_INLINE CMMatrix44<T>::operator T* () const
{
    return (T*) &m[0][0];
}
///////////////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T>::operator const T* () const
{
    return (const T*) &m[0][0];
}

// Acces Grant
template <class T>
M_INLINE CMVector4D<T>& CMMatrix44<T>::operator [] (int32t i)
{
    ASSERT(i >= 0 && i <= 3);

    return (CMVector4D<T> &) m[i][0];
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMVector4D<T>& CMMatrix44<T>::operator [] (int32t i) const
{
    ASSERT(i >= 0 && i <= 3);

    return (CMVector4D<T> &) m[i][0];
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE T& CMMatrix44<T>::operator () (int32t i, int32t j)
{
    ASSERT(i >= 0 && i <= 3);
    ASSERT(j >= 0 && j <= 3);

    return m[i][j];
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE T CMMatrix44<T>::operator () (int32t i, int32t j) const
{
    ASSERT(i >= 0 && j <= 3);
    ASSERT(j >= 0 && j <= 3);

    return m[i][j];
}

////////////////////////////////////////////////////////////////////////////
// Assignment Operators
template <class T>
M_INLINE CMMatrix44<T>& CMMatrix44<T>::operator = (const CMMatrix44<T> &mtx)
{
    memcpy(&_11, &mtx._11, sizeof(T)*16);

    return *this;
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T>& CMMatrix44<T>::operator += (const CMMatrix44<T> &mtx)
{
    _11 += mtx._11; _12 += mtx._12; _13 += mtx._13; _14 += mtx._14;
    _21 += mtx._21; _22 += mtx._22; _23 += mtx._23; _24 += mtx._24;
    _31 += mtx._31; _32 += mtx._32; _33 += mtx._33; _34 += mtx._34;
    _41 += mtx._41; _42 += mtx._42; _43 += mtx._43; _44 += mtx._44;

    return *this;
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T>& CMMatrix44<T>::operator -= (const CMMatrix44<T> &mtx)
{
    _11 -= mtx._11; _12 -= mtx._12; _13 -= mtx._13; _14 -= mtx._14;
    _21 -= mtx._21; _22 -= mtx._22; _23 -= mtx._23; _24 -= mtx._24;
    _31 -= mtx._31; _32 -= mtx._32; _33 -= mtx._33; _34 -= mtx._34;
    _41 -= mtx._41; _42 -= mtx._42; _43 -= mtx._43; _44 -= mtx._44;

    return *this;
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T>& CMMatrix44<T>::operator *= (const T & val)
{
    _11 *= val; _12 *= val; _13 *= val; _14 *= val;
    _21 *= val; _22 *= val; _23 *= val; _24 *= val;
    _31 *= val; _32 *= val; _33 *= val; _34 *= val;
    _41 *= val; _42 *= val; _43 *= val; _44 *= val;

    return *this;
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T>& CMMatrix44<T>::operator *= (const CMMatrix44<T>& mtx)
{
    CMMtx44Mult(*this,*this,mtx);

    return *this;
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T>& CMMatrix44<T>::operator /= (const T & val)
{
    ASSERT(!mnear0(val));

    T inv = T(1)/val;

    _11 *= inv; _12 *= inv; _13 *= inv; _14 *= inv;
    _21 *= inv; _22 *= inv; _23 *= inv; _24 *= inv;
    _31 *= inv; _32 *= inv; _33 *= inv; _34 *= inv;
    _41 *= inv; _42 *= inv; _43 *= inv; _44 *= inv;

    return *this;
}

////////////////////////////////////////////////////////////////////////////
// Unary Operators
template <class T>
M_INLINE CMMatrix44<T> CMMatrix44<T>::operator + () const
{
    return *this;
}

// Negate
template <class T>
M_INLINE CMMatrix44<T> CMMatrix44<T>::operator - () const\
{
    return CMMatrix44(  -_11, -_12, -_13, -_14,
                        -_21, -_22, -_23, -_24,
                        -_31, -_32, -_33, -_34,
                        -_41, -_42, -_43, -_44  );
}

////////////////////////////////////////////////////////////////////////////
// Binary Operators
template <class T>
M_INLINE CMMatrix44<T> CMMatrix44<T>::operator + (const CMMatrix44 &mtx) const
{
    return CMMatrix44(  _11 + mtx._11,  _12 + mtx._12,  _13 + mtx._13,  _14 + mtx._14,
                        _21 + mtx._21,  _22 + mtx._22,  _23 + mtx._23,  _24 + mtx._24,
                        _31 + mtx._31,  _32 + mtx._32,  _33 + mtx._33,  _34 + mtx._34,
                        _41 + mtx._41,  _42 + mtx._42,  _43 + mtx._43,  _44 + mtx._44   );
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T> CMMatrix44<T>::operator - (const CMMatrix44 &mtx) const
{
    return CMMatrix44(  _11 - mtx._11,  _12 - mtx._12,  _13 - mtx._13,  _14 - mtx._14,
                        _21 - mtx._21,  _22 - mtx._22,  _23 - mtx._23,  _24 - mtx._24,
                        _31 - mtx._31,  _32 - mtx._32,  _33 - mtx._33,  _34 - mtx._34,
                        _41 - mtx._41,  _42 - mtx._42,  _43 - mtx._43,  _44 - mtx._44   );
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T> CMMatrix44<T>::operator * (const T & val) const
{
    return CMMatrix44(  _11*val, _12*val, _13*val, _14*val,
                        _21*val, _22*val, _23*val, _24*val,
                        _31*val, _32*val, _33*val, _34*val,
                        _41*val, _42*val, _43*val, _44*val  );
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T> CMMatrix44<T>::operator * (const CMMatrix44<T> &mtx) const
{
    CMMatrix44 mtxT;
    CMMtx44Mult(mtxT,*this,mtx);
    return mtxT;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix44<T> CMMatrix44<T>::operator / (const T & val) const
{
    ASSERT(!mnear0(val));
    T inv = T(1)/val;
    return CMMatrix44(  _11*inv, _12*inv, _13*inv, _14*inv,
                        _21*inv, _22*inv, _23*inv, _24*inv,
                        _31*inv, _32*inv, _33*inv, _34*inv,
                        _41*inv, _42*inv, _43*inv, _44*inv  );
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE bool CMMatrix44<T>::operator == (const CMMatrix44<T> &mtx) const
{
    return memcmp(this, &mtx, sizeof(T)*16) == 0;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE bool CMMatrix44<T>::operator != (const CMMatrix44<T> &mtx) const
{
    return memcmp(this, &mtx, sizeof(T)*16) != 0;
}

///////////////////////////////////////////////////////////////////////////////////////
// Set Identity Matrix
template <class T>
M_INLINE void CMMtx44Identity(CMMatrix44<T> &M)
{
    M.m[0][1] = M.m[0][2] = M.m[0][3] =
    M.m[1][0] = M.m[1][2] = M.m[1][3] =
    M.m[2][0] = M.m[2][1] = M.m[2][3] =
    M.m[3][0] = M.m[3][1] = M.m[3][2] = T(0);

    M.m[0][0] = M.m[1][1] = M.m[2][2] = M.m[3][3] = T(1);
}

///////////////////////////////////////////////////////////////////////////////////////
// Is Identity
template <class T>
M_INLINE bool CMMtx44IsIdentity(const CMMatrix44<T> &M)
{
    return  M._11 == T(1) && M._12 == T(0) && M._13 == T(0) && M._14 == T(0) &&
            M._21 == T(0) && M._22 == T(1) && M._23 == T(0) && M._24 == T(0) &&
            M._31 == T(0) && M._32 == T(0) && M._33 == T(1) && M._34 == T(0) &&
            M._41 == T(0) && M._42 == T(0) && M._43 == T(0) && M._44 == T(1);
}

///////////////////////////////////////////////////////////////////////////////////////
// Setting Matrix Row
template <class T>
M_INLINE void CMMtx44Set(int row, CMMatrix44<T> &M, const CMVector3D<T> &vec, T w)
{
    ASSERT(row >= 0 && row <= 3);

    M[row][0] = vec.x;
    M[row][1] = vec.y;
    M[row][2] = vec.z;
    M[row][3] = w;
}

///////////////////////////////////////////////////////////////////////////////////////
// Setting Matrix Row
template <class T>
M_INLINE void CMMtx44Set(int row, CMMatrix44<T> &M, const CMVector4D<T> &vec)
{
    ASSERT(row >= 0 && row <= 3);

    M[row][0] = vec.x;
    M[row][1] = vec.y;
    M[row][2] = vec.z;
    M[row][3] = vec.w;
}

///////////////////////////////////////////////////////////////////////////////////////
// 3x4 Matrix Class Inline Functions
//////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// 3x4 Matrix Class Constructors

// Default Constructor
template <class T>
M_INLINE CMMatrix34<T>::CMMatrix34(){}

// Bulding Matrix From 3x3 Matrix And 3D Vector
template <class T>
M_INLINE CMMatrix34<T>::CMMatrix34(const CMMatrix33<T> &mtx33, const CMVector3D<T> &vec)
{
    memcpy(&this->M._11,&mtx33._11, sizeof(T)*9);
    memcpy(&this->V.v[0], &vec.v[0], sizeof(T)*3);
}

// Copy Constructor
template <class T>
M_INLINE CMMatrix34<T>::CMMatrix34(const CMMatrix34<T>& mtx34)
{
    memcpy(&this->M._11,&mtx34.M._11, sizeof(T)*9);
    memcpy(&this->V.v[0], &mtx34.V.v[0], sizeof(T)*3);
}

// Building Matrix By Getting 9 Matrix Values And 3 Vector Values
template <class T>
M_INLINE CMMatrix34<T>::CMMatrix34(T m11, T m12, T m13, T m21, T m22, T m23, T m31, T m32, T m33,
                                   T vx, T vy, T vz)
{
    M._11 = m11;    M._12 = m12;    M._13 = m13;
    M._21 = m21;    M._22 = m22;    M._23 = m23;
    M._31 = m31;    M._32 = m32;    M._33 = m33;

    V.x = vx;   V.y = vy;   V.z = vz;
}

// Access Grants
template <class T>
M_INLINE CMVector3D<T>& CMMatrix34<T>::operator [] (int32t i)
{
    ASSERT(i >= 0 && i <= 3);

    if(i < 3)
        return M[i];
    else
        return V;
}

template <class T>
M_INLINE CMVector3D<T>& CMMatrix34<T>::operator [] (int32t i) const
{
    ASSERT(i >= 0 && i <= 3);

    if(i < 3)
        return M[i];
    else
        return (CMVector3D<T> &)V;
}

////////////////////////////////////////////////////////////////////////////
// Assignment Operators
template <class T>
M_INLINE CMMatrix34<T>& CMMatrix34<T>::operator = (const CMMatrix34<T> &mtx34)
{
    M = mtx34.M;
    V = mtx34.V;

    return *this;
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T>& CMMatrix34<T>::operator += ( const CMMatrix34<T> &mtx34 )
{
    M += mtx34.M;
    V += mtx34.V;

    return *this;
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T>& CMMatrix34<T>::operator -= ( const CMMatrix34<T> &mtx34 )
{
    M -= mtx34.M;
    V -= mtx34.V;

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T>& CMMatrix34<T>::operator *= ( const CMMatrix34<T> &mtx34 )
{
    CMMtx34Mult(*this,*this,mtx34);

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T>& CMMatrix34<T>::operator *= ( const T &val )
{
    M *= val;
    V *= val;

    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T>& CMMatrix34<T>::operator /= ( const T &val )
{
    ASSERT(!mnear0(val));

    T inv = T(1)/val;

    M *= inv;
    V *= inv;

    return *this;
}

////////////////////////////////////////////////////////////////////////////
// Unary Operators
template <class T>
M_INLINE CMMatrix34<T> CMMatrix34<T>::operator + () const
{
    return *this;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T> CMMatrix34<T>::operator - () const
{
    return CMMatrix34(  -M._11, -M._12, -M._13,
                        -M._21, -M._22, -M._23,
                        -M._31, -M._32, -M._33,
                        -V.x, -V.y, -V.z);
}

////////////////////////////////////////////////////////////////////////////
// Binary Operators
template <class T>
M_INLINE CMMatrix34<T> CMMatrix34<T>::operator + ( const CMMatrix34<T> &mtx34 ) const
{
    return CMMatrix34(  M._11 + mtx34.M._11, M._12 + mtx34.M._12, M._13 + mtx34.M._13,
                        M._21 + mtx34.M._21, M._22 + mtx34.M._22, M._23 + mtx34.M._23,
                        M._31 + mtx34.M._21, M._32 + mtx34.M._22, M._33 + mtx34.M._23,
                        V.x + mtx34.V.x,  V.y + mtx34.V.y,  V.z + mtx34.V.z);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T> CMMatrix34<T>::operator - ( const CMMatrix34<T> &mtx34 ) const
{
    return CMMatrix34(  M._11 - mtx34.M._11, M._12 - mtx34.M._12, M._13 - mtx34.M._13,
                        M._21 - mtx34.M._21, M._22 - mtx34.M._22, M._23 - mtx34.M._23,
                        M._31 - mtx34.M._21, M._32 - mtx34.M._22, M._33 - mtx34.M._23,
                        V.x - mtx34.V.x,  V.y - mtx34.V.y,  V.z - mtx34.V.z);
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T> CMMatrix34<T>::operator * ( const CMMatrix34<T> &mtx34 ) const
{
    CMMatrix34 mtx;
    CMMtx34Mult(mtx,*this,mtx34);
    return mtx;
}

////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T> CMMatrix34<T>::operator * ( const T &val ) const
{
    CMMatrix34 mtx;
    mtx.M._11 = M._11*val;  mtx.M._12 = M._12*val;  mtx.M._13 = M._13*val;
    mtx.M._21 = M._21*val;  mtx.M._22 = M._22*val;  mtx.M._23 = M._23*val;
    mtx.M._31 = M._31*val;  mtx.M._32 = M._32*val;  mtx.M._33 = M._33*val;
    mtx.V.x = V.x*val;  mtx.V.y = V.y*val;  mtx.V.z = V.z*val;

    return mtx;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE CMMatrix34<T> CMMatrix34<T>::operator / ( const T &val ) const
{
    ASSERT(!mnear0(val));
    T inv = T(1)/val;
    CMMatrix34 mtx;

    mtx.M._11 = M._11*inv;  mtx.M._12 = M._12*inv;  mtx.M._13 = M._13*inv;
    mtx.M._21 = M._21*inv;  mtx.M._22 = M._22*inv;  mtx.M._23 = M._23*inv;
    mtx.M._31 = M._31*inv;  mtx.M._32 = M._32*inv;  mtx.M._33 = M._33*inv;
    mtx.V.x = V.x*inv;  mtx.V.y = V.y*inv;  mtx.V.z = V.z*inv;

    return mtx;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE bool CMMatrix34<T>::operator == ( const CMMatrix34<T> &mtx34 ) const
{
    return M == mtx34.M && V == mtx34.V;
}
////////////////////////////////////////////////////////////////////////////
template <class T>
M_INLINE bool CMMatrix34<T>::operator != ( const CMMatrix34<T> &mtx34 ) const
{
    return M != mtx34.M || V != mtx34.V;
}

///////////////////////////////////////////////////////////////////////////////////////
// Set Identity Matrix
template <class T>
M_INLINE void CMMtx34Identity(CMMatrix34<T> &M)
{
    CMMtx33Identity(&M.M);
}

///////////////////////////////////////////////////////////////////////////////////////
// Is Identity
template <class T>
M_INLINE bool CMMtx34IsIdentity(const CMMatrix34<T> &M)
{
    return CMMtx33IsIdentity(M.M);
}



#endif //_MATRIX_H_
