////////////////////////////////////////////////////////////////////////////////
// platform.h
//
// Defines common platform types
////////////////////////////////////////////////////////////////////////////////

#ifndef __PLATFORMTYPES_H__
#define __PLATFORMTYPES_H__

// C++ version
#define __CPLUSPLUS_VER_01_     1
#define __CPLUSPLUS_VER_98_     19971L
#define __CPLUSPLUS_VER_11_     201103L
#define __CPLUSPLUS_VER_14_     201402L
#define __CPLUSPLUS_VER_17_     201703L
#define __CPLUSPLUS_VER_20_     202002L

// C version
#define __C_VER_90_ 199000L
#define __C_VER_99_ 199901L
#define __C_VER_11_ 201112L
/*
_WIN32
__MINGW64__
__APPLE__
__linux__
*/


//Compiler type for C++
#if defined( __cplusplus)
#define __C_VER_ 0

// Using MS Visual Studio
#if defined(_MSC_VER)

#if (_MSC_VER >= 1700)
    #define __CPLUSPLUS_VER_    __CPLUSPLUS_VER_11_
#else
    #define __CPLUSPLUS_VER_    __CPLUSPLUS_VER_98_
#endif //_MSC_VER version

// Using GNU complier
#elif defined(__GNUC__)

// Set C++ version for GNU Compiler
#define __CPLUSPLUS_VER_    __cplusplus

#endif //__GNUC__ version

// Compiler type for C
#else //C lang
#define __CPLUSPLUS_VER_ 0

#if defined(__STDC_VERSION__)

#if (__STDC_VERSION__ >= __C_VER_11_)
#define __C_VER_    __C_VER_11_
#else
#define __C_VER_    __C_VER_99_
#endif // C version

#else
#define __C_VER_    __C_VER_90_
#endif //__STDC_VERSION__

#endif// Compiler version

////////////////////////////////////////////////////////////////////////////////
// Define basic data types
#if (__CPLUSPLUS_VER_ >= __CPLUSPLUS_VER_11_) || (__C_VER_ >= __C_VER_99_)

#include <stdint.h>
#if !defined(__cplusplus)
#include <inttypes.h>
#endif

typedef int8_t              tbyte;      // signed byte
typedef uint8_t             ubyte;      // unsigned char
typedef int8_t              int8t;      // signed 8-bit integer (char)
typedef uint8_t             uint8t;     // unsigned 8-bit integer (unsigned char)
typedef int16_t             int16t;     // signed 16-bit integer
typedef uint16_t            uint16t;    // unsigned 16-bit integer
typedef int32_t             int32t;     // signed 32-bit integer
typedef uint32_t            uint32t;    // signed 32-bit integer
typedef int64_t             int64t;     // signed 64-bit integer
typedef uint64_t            uint64t;    // unsigned 64-bit integer

// __int128_t
#else

typedef char                tbyte;      // signed byte
typedef unsigned char       ubyte;      // unsigned char
typedef char                int8t;      // signed 8-bit integer (char)
typedef unsigned char       uint8t;     // unsigned 8-bit integer (unsigned char)
typedef short               int16t;     // signed 16-bit integer
typedef unsigned short      uint16t;    // unsigned 16-bit integer
typedef int                 int32t;     // signed 32-bit integer
typedef unsigned int        uint32t;    // signed 32-bit integer

#ifdef _WIN32
typedef __int64             int64t;     // signed 64-bit integer
typedef unsigned __int64    uint64t;    // unsigned 64-bit integer
#else  //Linux
typedef long long           int64t;     // signed 64-bit integer
typedef unsigned long long  uint64t;    // unsigned 64-bit integer
#endif //os type

#endif // C and C++ version for data types declaration
////////////////////////////////////////////////////////////////////////////////

/*
MSVC++ 14.0 _MSC_VER == 1900 (Visual Studio 2015)
MSVC++ 12.0 _MSC_VER == 1800 (Visual Studio 2013)
MSVC++ 11.0 _MSC_VER == 1700 (Visual Studio 2012)
MSVC++ 10.0 _MSC_VER == 1600 (Visual Studio 2010)
MSVC++ 9.0  _MSC_VER == 1500 (Visual Studio 2008)
MSVC++ 8.0  _MSC_VER == 1400 (Visual Studio 2005)
MSVC++ 7.1  _MSC_VER == 1310 (Visual Studio 2003)
MSVC++ 7.0  _MSC_VER == 1300
MSVC++ 6.0  _MSC_VER == 1200
MSVC++ 5.0  _MSC_VER == 110
*/


/*
__GNUC__
__GNUC_MINOR__
__GNUC_PATCHLEVEL__
__GNUG__

__clang__
__clang_major__
__clang_minor__
__clang_patchlevel__
*/


//NULL pointer
#if (__CPLUSPLUS_VER_ == __CPLUSPLUS_VER_11__)
    #define NULLPTR     nullptr
#else
    #define NULLPTR     NULL
#endif //NULL pointer

enum VARDATATYPE
{
    VTYPE_BYTE,
    VTYPE_UBYTE,
    VTYPE_INT16,
    VTYPE_UINT16,
    VTYPE_INT32,
    VTYPE_UINT32,
    VTYPE_INT64,
    VTYPE_UINT64,
    VTYPE_HALF,
    VTYPE_FLOAT,
    VTYPE_DOUBLE,
};
static const unsigned int vartypesizestab[] = {1,1,2,2,4,4,8,8,2,4,8};


//ASSERTION
#ifdef _WIN32
#include <assert.h>
#define ASSERT assert
#else
#include <assert.h>
#define ASSERT assert
#endif //ASSERTION



// STATIC_ASSERT
#if defined(__CPLUSPLUS_VER_) && (__CPLUSPLUS_VER_ >= __CPLUSPLUS_VER_11_)
    #include <type_traits>
    #define STATIC_ASSERT(COND,MSG) static_assert(COND,MSG)
#elif defined(_C_VER_) && (_C_VER_ >= _C_VER_99_)
    #include <assert.h>
    #define STATIC_ASSERT(COND,MSG) static_assert(COND,MSG)
#else
    #define STATIC_ASSERT(COND,MSG)
#endif

STATIC_ASSERT (sizeof(tbyte)    == 1, "type tbyte does not match size");
STATIC_ASSERT (sizeof(ubyte)    == 1, "type ubyte does not match size");
STATIC_ASSERT (sizeof(int8t)    == 1, "type int8t does not match size");
STATIC_ASSERT (sizeof(uint8t)   == 1, "type uin8t does not match size");
STATIC_ASSERT (sizeof(int16t)   == 2, "type int16t does not match size");
STATIC_ASSERT (sizeof(uint16t)  == 2, "type uint16t does not match size");
STATIC_ASSERT (sizeof(int32t)   == 4, "type int32t does not match size");
STATIC_ASSERT (sizeof(uint32t)  == 4, "type uint32t does not match size");
STATIC_ASSERT (sizeof(int64t)   == 8, "type int64t does not match size");
STATIC_ASSERT (sizeof(uint64t)  == 8, "type uint64t does not match size");

#endif //__PLATFORMTYPES_H__

////////////////////////////////////////////////////////////////////////////////
