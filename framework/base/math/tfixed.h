////////////////////////////////////////////////////////////////////////////////
//  tfixed.h
//
//  Math 32-bit fixed data type class
//
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __TFIXED_H__
#define __TFIXED_H__

////////////////////////////////////////////////////////////////////////////////
template <int bits>
class tfixed
{
    protected:
        int value;
        static const int one = 1 << bits;
        static const int half = bits ? (1 << (bits-1)) : 0;

    public:
        // Constructors Declarations
        tfixed();
        tfixed(const tfixed<bits>& fx_val);
        tfixed(const float f_val);
        tfixed(const double d_val);
        tfixed(const int i_val);


        // Cast To Pointer int *
        operator int* () const;
        operator const int* () const;

        // Cast To Primitive Build-in Types
        operator float() const;
        operator double() const;
        operator int() const;
        operator bool() const;

        // Assignment Operators
        tfixed<bits>& operator = (const tfixed<bits>& fx_val);
        tfixed<bits>& operator += (const tfixed<bits>& fx_val);
        tfixed<bits>& operator -= (const tfixed<bits>& fx_val);
        tfixed<bits>& operator *= (const tfixed<bits>& fx_val);
        tfixed<bits>& operator /= (const tfixed<bits>& fx_val);

        tfixed<bits>& operator >>= (int i);
        tfixed<bits>& operator <<= (int i);

        tfixed<bits>& operator &= (int mask);
        tfixed<bits>& operator |= (int mask);
        tfixed<bits>& operator ^= (int mask);


        // Unary Operators
        tfixed<bits> operator + () const;
        tfixed<bits> operator - () const;

        // Binary Operators
        tfixed<bits> operator + (const tfixed<bits> &fx_val) const;
        tfixed<bits> operator - (const tfixed<bits> &fx_val) const;
        tfixed<bits> operator * (const tfixed<bits> &fx_val) const;
        tfixed<bits> operator / (const tfixed<bits> &fx_val) const;

        tfixed<bits> operator >> (int i) const;
        tfixed<bits> operator << (int i) const;

        tfixed<bits> operator & (int mask) const;
        tfixed<bits> operator | (int mask) const;
        tfixed<bits> operator ^ (int mask) const;

        bool operator == (const tfixed<bits> &fx_val) const;
        bool operator != (const tfixed<bits> &fx_val) const;

        bool operator <  (const tfixed<bits> &fx_val) const;
        bool operator >  (const tfixed<bits> &fx_val) const;

        bool operator <= (const tfixed<bits> &fx_val) const;
        bool operator >= (const tfixed<bits> &fx_val) const;

        bool operator && (const tfixed<bits> &fx_val) const;
        bool operator || (const tfixed<bits> &fx_val) const;
};
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::tfixed():value(0){}
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::tfixed(const tfixed<bits>& fx_val):value(fx_val.value){}
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::tfixed(const float f_val)
{
    value = (int)(f_val * (float)one);
}
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::tfixed(const double d_val)
{
    value = (int)(d_val * (double)one);
}
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::tfixed(const int i_val)
{
    value = i_val << bits;
}

////////////////////////////////////////////////////////////////////////////////
//Cast To Pointer int *
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::operator int* () const
{
    return (int*) this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::operator const int* () const
{
    return (const int*) this;
}

////////////////////////////////////////////////////////////////////////////////
//Cast To Primitive Build-in Types
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::operator float() const
{
    return (float)value / (float)one;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::operator double() const
{
    return (double)value / (double)one;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::operator int() const
{
    return value >> bits;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>::operator bool() const
{
    return value != 0;
}

////////////////////////////////////////////////////////////////////////////////
// Assignment Operators
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator = (const tfixed<bits>& fx_val)
{
    value = fx_val.value;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator += (const tfixed<bits>& fx_val)
{
    value += fx_val.value;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator -= (const tfixed<bits>& fx_val)
{
    value -= fx_val.value;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator *= (const tfixed<bits>& fx_val)
{
    value = (int)(((long long)value * (long long)fx_val.value + half) >> bits);
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator /= (const tfixed<bits>& fx_val)
{
    value = (int)(((long long)value << bits) / (long long)fx_val.value);
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator >>= (int i)
{
    value >>= i;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator <<= (int i)
{
    value <<= i;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator &= (int mask)
{
    value &= mask;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator |= (int mask)
{
    value |= mask;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits>& tfixed<bits>::operator ^= (int mask)
{
    value ^= mask;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
// Unary Operators
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator + () const
{
    tfixed<bits> out;
    out.value = value;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator - () const
{
    tfixed<bits> out;
    out.value = -value;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
// Binary Operators
////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator + (const tfixed<bits> &fx_val) const
{
    tfixed<bits> out;
    out.value = value + fx_val.value;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator - (const tfixed<bits> &fx_val) const
{
    tfixed<bits> out;
    out.value = value - fx_val.value;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator * (const tfixed<bits> &fx_val) const
{
    tfixed<bits> out;
    out.value = (int)(((long long)value * (long long)fx_val.value + half) >> bits);
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator / (const tfixed<bits> &fx_val) const
{
    tfixed<bits> out;
    out.value = (int)(((long long)value << bits) / (long long)fx_val.value);
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator >> (int i) const
{
    tfixed<bits> out;
    out.value = value >> i;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator << (int i) const
{
    tfixed<bits> out;
    out.value = value << i;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator & (int mask) const
{
    tfixed<bits> out;
    out.value = value & mask;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator | (int mask) const
{
    tfixed<bits> out;
    out.value = value | mask;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
tfixed<bits> tfixed<bits>::operator ^ (int mask) const
{
    tfixed<bits> out;
    out.value = value ^ mask;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
bool tfixed<bits>::operator == (const tfixed<bits> &fx_val) const
{
    return value == fx_val.value;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
bool tfixed<bits>::operator != (const tfixed<bits> &fx_val) const
{
    return value != fx_val.value;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
bool tfixed<bits>::operator <  (const tfixed<bits> &fx_val) const
{
    return value < fx_val.value;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
bool tfixed<bits>::operator >  (const tfixed<bits> &fx_val) const
{
    return value > fx_val.value;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
bool tfixed<bits>::operator <= (const tfixed<bits> &fx_val) const
{
    return value <= fx_val.value;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
bool tfixed<bits>::operator >= (const tfixed<bits> &fx_val) const
{
    return value >= fx_val.value;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
bool tfixed<bits>::operator && (const tfixed<bits> &fx_val) const
{
    return value && fx_val.value;
}

////////////////////////////////////////////////////////////////////////////////
template <int bits>
bool tfixed<bits>::operator || (const tfixed<bits> &fx_val) const
{
    return value || fx_val.value;
}
////////////////////////////////////////////////////////////////////////////////

#endif //__TFIXED_H__
