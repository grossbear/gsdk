////////////////////////////////////////////////////////////////////////////////
// mtypes.h
//
// Defines some simple types in library
////////////////////////////////////////////////////////////////////////////////

#ifndef __MTYPES_H__
#define __MTYPES_H__

#include "tfixed32.h"
#include "tfixed64.h"
#include "thalf.h"


////////////////////////////////////////////////////////////////////////////////
// Array struct that consists of 2 elements
template <typename T>
struct marray2
{
#if defined(__CPLUSPLUS_VER_) && (__CPLUSPLUS_VER_ >= __CPLUSPLUS_VER_11_)
    union
    {
        struct
        {
            T elem1;
            T elem2;
        };
        T array[2];
    };
#else
    T elem1;
    T elem2;

#endif
//    marray2(): elem1(T(0)), elem2(T(0)) {}
//    marray2(T e1, T e2): elem1(e1), elem2(e2) {}
//    marray2(T elems[]): elem1(elems[0]), elem2(elems[1]) {}
};

////////////////////////////////////////////////////////////////////////////////
// Array struct that consists of 3 elements
template <typename T>
struct marray3
{
#if defined(__CPLUSPLUS_VER_) && (__CPLUSPLUS_VER_ >= __CPLUSPLUS_VER_11_)
    union
    {
        struct
        {
            T elem1;
            T elem2;
            T elem3;
        };
        T array[3];
    };
#else
    T elem1;
    T elem2;
    T elem3;

//    marray3(): elem1(T(0)), elem2(T(0)), elem3(T(0)) {}
//    marray3(T e1, T e2, T e3): elem1(e1), elem2(e2), elem3(e3) {}
//    marray3(T elems[]): elem1(elems[0]), elem2(elems[1]), elem3(elems[3]) {}
#endif
};


////////////////////////////////////////////////////////////////////////////////
// Array struct that consists of 4 elements
template <typename T>
struct marray4
{
#if defined(__CPLUSPLUS_VER_) && (__CPLUSPLUS_VER_ >= __CPLUSPLUS_VER_11_)
    union
    {
        struct
        {
            T elem1;
            T elem2;
            T elem3;
            T elem4;
        };
        T array[4];
    };
#else
    T elem1;
    T elem2;
    T elem3;
    T elem4;

//    marray4(): elem1(T(0)), elem2(T(0)), elem3(T(0)), elem4(T(0)) {}
//    marray4(T e1, T e2, T e3, T e4): elem1(e1), elem2(e2), elem3(e3), elem4(e4) {}
//    marray4(T elems[]): elem1(elems[0]), elem2(elems[1]), elem3(elems[3]), elem4(elems[4]) {}
#endif
};

#endif //__MTYPES_H__
