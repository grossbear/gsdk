//-----------------------------------------------------------------------------
//
// File: mathval.inl
// Content: inline class of const values
//
//-----------------------------------------------------------------------------

#ifndef _CMVALUES_INL_
#define _CMVALUES_INL_

///////////////////////////////////////////////////////////////////////////////
template <class TReal>
class CMValues
{
public:
/*
    static constexpr T QUATER   = T(0.25);
    static constexpr T HALF     = T(0.5);
    static constexpr T ONE      = T(1.0);
    static constexpr T TWO      = T(2.0);
    static constexpr T THREE    = T(3.0);
    static constexpr T FOUR     = T(4.0);
    static constexpr T FIVE     = T(5.0);
*/

    static const TReal QUATER;
    static const TReal HALF;
    static const TReal ONE;
    static const TReal TWO;
    static const TReal THREE;
    static const TReal FOUR;
    static const TReal FIVE;

};

///////////////////////////////////////////////////////////////////////////////
template <class TReal>
const TReal CMValues<TReal>::QUATER = TReal(0.25);

///////////////////////////////////////////////////////////////////////////////
template <class TReal>
const TReal CMValues<TReal>::HALF   = TReal(0.5);

///////////////////////////////////////////////////////////////////////////////
template <class TReal>
const TReal CMValues<TReal>::ONE    = TReal(1.0);

///////////////////////////////////////////////////////////////////////////////
template <class TReal>
const TReal CMValues<TReal>::TWO    = TReal(2.0);

///////////////////////////////////////////////////////////////////////////////
template <class TReal>
const TReal CMValues<TReal>::THREE  = TReal(3.0);

///////////////////////////////////////////////////////////////////////////////
template <class TReal>
const TReal CMValues<TReal>::FOUR   = TReal(4.0);

///////////////////////////////////////////////////////////////////////////////
template <class TReal>
const TReal CMValues<TReal>::FIVE   = TReal(5.0);

///////////////////////////////////////////////////////////////////////////////

#endif //CMVALUES_INL_
