////////////////////////////////////////////////////////////////////////////////
// elementary.cpp
//
// Math elementary functions definition
////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "mathdefs.h"
#include "mathconsts.h"
#include "primary.h"
#include "elementary.h"
#include "intreal.h"


////////////////////////////////////////////////////////////////////////////////
// Square root calc functions
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Square Root Table Values
////////////////////////////////////////////////////////////////////////////////
INTFLOAT table_sqrt[] =
{
#include "tables/sqrt_table256.h"
};

////////////////////////////////////////////////////////////////////////////////
// Table square root function
////////////////////////////////////////////////////////////////////////////////
float m_tsqrt(float f)
{
    ASSERT(!mless0(f));

    INTFLOAT fi;
    fi.f = f;
    int32t n = fi.i;
    int32t e = (n >> 1) & 0x3f800000;
    n = (n >> 16) & 0xff;
    fi.i = e + table_sqrt[n].i;

    return fi.f;
}

////////////////////////////////////////////////////////////////////////////////
// Fast reverse square root function
////////////////////////////////////////////////////////////////////////////////
float m_rfsqrt(float f)
{
    ASSERT(!mless0(f));

    float xhalf = 0.5f * f;
    int32t i = *(int32t*)&f;    // evil floating point bit level hacking
    i = 0x5f3759d5 - (i >> 1);
    f = *(float*)&i;
    f = f*(1.5f - xhalf*f*f);   // 1st iteration
    f = f*(1.5f - xhalf*f*f);   // 2nd iteration, this can be removed

    return f;
}


////////////////////////////////////////////////////////////////////////////////
// Approximation square root function
////////////////////////////////////////////////////////////////////////////////
float m_asqrt(float x)
{
    int32t iter = 12;

    float sum = 1.0f;
    float val = 1.0f;

    float nfact1 = 1.0f;
    float factVal1 = 1.0f;
    float factVal2 = 1.0f;
    float n = 0.0f;
    float cnstVal2 = 1.0f;

    for(int32t i = 1; i < iter; i++)
    {
        float tmp1, tmp2;
        val = val*x;
        tmp1 = (i%2)?(-1.0f):(1.0f);

        factVal1 *= nfact1;
        nfact1 += 1.0f;
        factVal1 *= nfact1;
        nfact1 += 1.0f;

        n += 1.0f;
        tmp2 = 1.0f - 2.0f*n;
        factVal2 *= n;
        cnstVal2 *= 4.0f;

        sum = sum + (tmp1*factVal1)/(tmp2*factVal2*factVal2*cnstVal2)*val;
    }

    return sum;
}


////////////////////////////////////////////////////////////////////////////////
// 32-bit fixed point square root function
////////////////////////////////////////////////////////////////////////////////
uint32t m_isqrt(uint32t x)
{
//
// this follows http://www.finesse.demon.co.uk/steven/sqrt.html
// 4 cycle/bit C routine
// Wilco Dijkstra also provided the following C code
// which produces optimised ARM code which takes 4 cycles per bit:

// fixed-point square root
#define ITER1(N) \
    tr = root + (1 << (N)); \
    if (n >= tr << (N))   \
    {   n -= tr << (N);   \
        root |= 2 << (N); \
    }

    unsigned int root;
    unsigned int tr;
    unsigned int n = x;
    root = 0;

    ITER1 (15);    ITER1 (14);    ITER1 (13);    ITER1 (12);
    ITER1 (11);    ITER1 (10);    ITER1 ( 9);    ITER1 ( 8);
    ITER1 ( 7);    ITER1 ( 6);    ITER1 ( 5);    ITER1 ( 4);
    ITER1 ( 3);    ITER1 ( 2);    ITER1 ( 1);    ITER1 ( 0);

    return root >> 1;
#undef ITER1
}


////////////////////////////////////////////////////////////////////////////////
// 64-bit fixed point square root function
////////////////////////////////////////////////////////////////////////////////
uint64t m_isqrt(uint64t x)
{
#define ITER1(N) \
    tr = root + (1LL << (N)); \
    if (n >= tr << (N))   \
    {   n -= tr << (N);   \
        root |= 2LL << (N); \
    }
    unsigned long long root;
    unsigned long long tr;
    unsigned long long n = x;
    root = 0LL;

    ITER1 (31LL);    ITER1 (30LL);    ITER1 (29LL);    ITER1 (28LL);
    ITER1 (27LL);    ITER1 (26LL);    ITER1 (25LL);    ITER1 (24LL);
    ITER1 (23LL);    ITER1 (22LL);    ITER1 (21LL);    ITER1 (20LL);
    ITER1 (19LL);    ITER1 (18LL);    ITER1 (17LL);    ITER1 (16LL);

    ITER1 (15LL);    ITER1 (14LL);    ITER1 (13LL);    ITER1 (12LL);
    ITER1 (11LL);    ITER1 (10LL);    ITER1 ( 9LL);    ITER1 ( 8LL);
    ITER1 ( 7LL);    ITER1 ( 6LL);    ITER1 ( 5LL);    ITER1 ( 4LL);
    ITER1 ( 3LL);    ITER1 ( 2LL);    ITER1 ( 1LL);    ITER1 ( 0LL);

    return rooti >> 1LL;
#undef ITER1
}



////////////////////////////////////////////////////////////////////////////////
// Trigonometry functions
////////////////////////////////////////////////////////////////////////////////

#define     TABLESIZE_SIN       256
#define     FIXED_TABLE_SIZE    1024
#define     TWOPISCALE          ((float)TABLESIZE_SIN * ((float)CONST_1_2PI))

////////////////////////////////////////////////////////////////////////////////
// Precalculated float sinus and cosinus values table
const INTFLOAT table_sinus [TABLESIZE_SIN] =
{
    #include "tables/sin_table256.h"
};
////////////////////////////////////////////////////////////////////////////////
// Float table sinus function
float m_tsinf(float angle)
{
    ASSERT(angle >= -10000.0f && angle <= 10000.0f);

    INTFLOAT fi;
    fi.f = angle * TWOPISCALE + bias.f;
    int32t i = fi.i & (TABLESIZE_SIN - 1);

    return table_sinus[i].f;
}

////////////////////////////////////////////////////////////////////////////////
// Float table cosinus function
float m_tcosf(float angle)
{
    ASSERT(angle >= -10000.0f && angle <= 10000.0f);

    INTFLOAT  fi;
    fi.f = angle * TWOPISCALE + bias.f;
    int32t i = (fi.i + (TABLESIZE_SIN / 4)) & (TABLESIZE_SIN - 1);

    return table_sinus[i].f;
}

////////////////////////////////////////////////////////////////////////////////
// Sinus and cosinus precalculated table function
void m_tsincosf(float angle, float & sin_val, float & cos_val)
{
    ASSERT(angle >= -10000.0f && angle <= 10000.0f);

    INTFLOAT fi;
    fi.f = angle * TWOPISCALE + bias.f;

    int i = fi.i & (TABLESIZE_SIN - 1);
    sin_val = sinus_tab[i].f;

    i = (fi.i + (TABLESIZE_SIN / 4)) & (TABLESIZE_SIN - 1);
    cos_val = sinus_tab[i].f;
}



////////////////////////////////////////////////////////////////////////////////
const uint32t fxcos_tab [FIXED_TABLE_SIZE+1] =+
{
    #include "tables/fxcos_table1024.h"
};

////////////////////////////////////////////////////////////////////////////////
const uint32t fxsin_tab [FIXED_TABLE_SIZE+1] =
{
    #include "tables/fxsin_table1024.h"
};

////////////////////////////////////////////////////////////////////////////////
const uint32t fxtan_tab [FIXED_TABLE_SIZE+1] =
{
    #include "tables/fxtan_table1024.h"
};

////////////////////////////////////////////////////////////////////////////////
const uint32t fxacos_tab [FIXED_TABLE_SIZE+1] =
{
    #include "tables/fxacos_table1024.h"
};

////////////////////////////////////////////////////////////////////////////////
const ubyte acosdeg_tab [101] =+
{
    #include "tables/acosdeg_table.h"
};

////////////////////////////////////////////////////////////////////////////////
// 32-bit fixed point table sinus function
tfixed32 m_tsinx(const tfixed32 &x)
{
    static const tfixed32 fx2pi = CMathConst< tfixed32 >::MATH_2PI;

    tfixed32 fxval = (*(int32t*)x < 0) ? (fx2pi - x) : (x);
    fxval /= fx2pi;
    fxval &= 0x0000FFFF; //fraction part
    tfixed32 fxent;
    *(int32t*)fxent = FIXED_TABLE_SIZE - 1;
    fxval *= fxent;

    int32t tabidx = *(int32t*)fxval;
    int32t intval = fxsin_tab[tabidx];

    tfixed32 outval;
    *(int32t*)outval = intval;

    return outval;
}

////////////////////////////////////////////////////////////////////////////////
// 32-bit fixed point table cosinus function
tfixed32 m_tcosx(const tfixed32 &x)
{
    static const tfixed32 fx2pi = CMathConst< tfixed32 >::MATH_2PI;

    tfixed32 fxval = (*(int32t*)x < 0) ? (-x) : (x);
    fxval /= fx2pi;
    fxval &= 0x0000FFFF;
    tfixed32 fxent;
    *(int32t*)fxent = FIXED_TABLE_SIZE - 1;
    fxval *= fxent;

    int32t tabidx = *(int32t*)fxval;
    int32t intval = fxcos_tab[tabidx];

    tfixed32 outval;
    *(int32t*)outval = intval;

    return outval;
}

////////////////////////////////////////////////////////////////////////////////
// Sinus and cosinus precalculated table function
void m_tsincosx(const tfixed32 & angle, tfixed32 & sin_val, tfixed32 & cos_val)
{
    static const tfixed32 fx2pi = CMathConst< tfixed32 >::MATH_2PI;

    sin_val = (*(int32t*)angle < 0) ? (fx2pi - angle) : (angle);
    sin_val /= fx2pi;
    sin_val &= 0x0000FFFF;

    cos_val = (*(int32t*)angle < 0) ? (-angle) : (angle);
    cos_val /= fx2pi;
    cos_val &= 0x0000FFFF;

    tfixed32 fxent;
    *(int32t*)fxent = FIXED_TABLE_SIZE - 1;

    int32t tabidx = 0;
    int32t intval = 0;

    sin_val *= fxent;
    tabidx = *(int32t*)sin_val;
    intval = fxsin_tab[tabidx];
    *(int32t*)sin_val = intval;

    cos_val *= fxent;
    tabidx = *(int32t*)cos_val;
    intval = fxcos_tab[tabidx];
    *(int32t*)cos_val = intval;
}

///////////////////////////////////////////////////////////////////////////////////////
// Arcus Cosinus Table Approxiamtion In Degrees
int32t m_tacosdeg(float f)
{
    int32t index = mftoi(f * 100);
    if(index > 100) index = 100;
    ubyte degval = acosdeg_tab[index];
    return (uint32t)degval;
}

/*
///////////////////////////////////////////////////////////////////////////////////////
int32t m_tacosdeg(const tfixed32<16> &x)
{
    static const tfixed32<16> max_idx(100);
    tfixed32<16> fx_index = max_idx * x;
    if( fx_index > max_idx) fx_index = max_idx;
    int32t index = (int32t)fx_index;
    ubyte degval = acosdeg_tab[index];
    return (int32t)degval;
}
*/


/*
////////////////////////////////////////////////////////////////////////////////
// Power functions
////////////////////////////////////////////////////////////////////////////////
float mpow(float x, float y)
{
    return (float)pow((double)x, (double)y);
}

////////////////////////////////////////////////////////////////////////////////
double mpow(double x, double y)
{
    return pow(x, y);
}
*/


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
