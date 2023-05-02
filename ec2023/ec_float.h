//! Rohde & Schwarz Engineering Competition 2023
//!
//! (c) Rohde & Schwarz GmbH & Co. KG, Munich
//!
//! This code is needed to run your solution locally. Any changes will be lost once you upload your solution!

#pragma once

#include "ec_measurement.h"
#include "ec_float.h"

#include <vector>

namespace ec
{
  class Float final
  {
  public:
    Float() = default;
    Float(float value)noexcept;
    Float(int value) noexcept;
    Float(size_t value) noexcept;
    Float(const Float& rhs)noexcept;
    Float(Float&& rhs)noexcept;
    Float& operator=(Float&& rhs)noexcept;

    Float& operator=(const Float& rhs)noexcept;
    Float& operator=(float rhs)noexcept;

    ~Float() = default;

    friend bool operator==( const Float& lhs, const Float& rhs)noexcept;
    friend bool operator==( const Float& lhs, const Float& rhs)noexcept;
    friend bool operator==( float lhs, const Float& rhs)noexcept;
    friend bool operator< ( const Float& lhs, const Float& rhs)noexcept;
    friend bool operator<=( const Float& lhs, const Float& rhs)noexcept;
    friend bool operator> ( const Float& lhs, const Float& rhs)noexcept;
    friend bool operator>=( const Float& lhs, const Float& rhs)noexcept;
    friend bool operator!=( const Float& lhs, const Float& rhs)noexcept;

    Float operator+(const Float& rhs)const noexcept;
    Float& operator+=(const Float& rhs) noexcept;
    friend Float operator+(float rhs, const Float& lhs)noexcept;

    Float operator-(const Float& rhs)const noexcept;
    Float& operator-=(const Float& rhs) noexcept;
    friend Float operator-(float rhs, const Float& lhs)noexcept;

    Float operator*(const Float& rhs)const noexcept;
    Float& operator*=(const Float& rhs) noexcept;
    friend Float operator*(float rhs, const Float& lhs)noexcept;

    Float operator/(const Float& rhs)const noexcept;
    Float& operator/=(const Float& rhs) noexcept;
    friend Float operator/(float rhs, const Float& lhs)noexcept;

    Float& operator++() noexcept;    // prefix increment
    Float operator++(int) noexcept;  // postfix increment
    Float& operator--() noexcept;    // prefix decrement
    Float operator--(int) noexcept;  // postfix decrement

 

    friend class VecHw;
    friend class StreamHw;

    friend Float ec_acos(Float input);
    friend Float ec_cos(Float input);
    friend Float ec_sin(Float input);
    friend Float ec_asin(Float input);
    friend Float ec_tan(Float input);
    friend Float ec_atan(Float input);
    friend Float ec_sqrt(Float input);
    friend Float ec_abs(Float input);
    friend Float ec_max(Float left, Float right);
    friend Float ec_exp(Float input);
    friend Float ec_round(Float input);
    friend Float ec_log(Float input);
    friend Float ec_log10(Float input);
    friend Float ec_pow(Float base, Float exp);
    

  private:

   
    float m_value = 0.0f;
    
    float toFloat()const  noexcept;
  };


  inline Float::Float(float value) noexcept
    : m_value(value)
  {
    Measurement::inc(Measurement::MeasType::assignment);
  }



  inline Float::Float(int value) noexcept
    : m_value(static_cast<float>(value))
  {
    Measurement::inc(Measurement::MeasType::assignment);
  }


  inline Float::Float(size_t value) noexcept
    : m_value(static_cast<float>(value))
  {
    Measurement::inc(Measurement::MeasType::assignment);
  }



  inline Float::Float(const Float& rhs) noexcept
    : m_value(rhs.m_value)
  {
    Measurement::inc(Measurement::MeasType::assignment);
  }


  inline Float& Float::operator=(const Float& rhs) noexcept
  {
    Measurement::inc(Measurement::MeasType::assignment);
    m_value = rhs.m_value;
    return *this;
  }


  inline Float::Float(Float&& rhs) noexcept
    : m_value(rhs.m_value)
  {
    Measurement::inc(Measurement::MeasType::assignment);
  }


  inline Float& Float::operator=(Float&& rhs) noexcept
  {
    Measurement::inc(Measurement::MeasType::assignment);
    m_value = rhs.m_value;
    return *this;
  }


  inline bool operator==(const Float& lhs, const Float& rhs)noexcept
  {
    Measurement::inc(Measurement::MeasType::compare);
    return lhs.m_value == rhs.m_value;
  }


  inline bool operator==(float lhs, const Float& rhs)noexcept
  {
    Measurement::inc(Measurement::MeasType::compare);
    return lhs == rhs.m_value;
  }


  inline Float& Float::operator=(float rhs) noexcept
  {
    Measurement::inc(Measurement::MeasType::assignment);
    m_value = rhs;
    return *this;
  }


  inline bool operator<(const Float& lhs, const Float& rhs)noexcept
  {
    Measurement::inc(Measurement::MeasType::compare);
    return lhs.m_value < rhs.m_value;
  }


  inline bool operator<=(const Float& lhs, const Float& rhs)noexcept
  {
    Measurement::inc(Measurement::MeasType::compare);
    return lhs.m_value <= rhs.m_value;
  }


  inline bool operator>(const Float& lhs, const Float& rhs)noexcept
  {
    Measurement::inc(Measurement::MeasType::compare);
    return lhs.m_value > rhs.m_value;
  }


  inline bool operator>=(const Float& lhs, const Float& rhs)noexcept
  {
    Measurement::inc(Measurement::MeasType::compare);
    return lhs.m_value >= rhs.m_value;
  }


  inline bool operator!=(const Float& lhs, const Float& rhs)noexcept
  {
    Measurement::inc(Measurement::MeasType::compare);
    return lhs.m_value != rhs.m_value;
  }


  inline Float Float::operator+(const Float& rhs)const noexcept
  {
    Measurement::inc(Measurement::MeasType::addition);
    return { m_value + rhs.m_value };
  }


  inline Float& Float::operator+=(const Float& rhs) noexcept
  {
    Measurement::inc(Measurement::MeasType::addition);
    Measurement::inc(Measurement::MeasType::assignment);
    m_value += rhs.m_value;
    return *this;
  }


  inline Float operator+(float rhs, const Float& lhs)noexcept
  {
    return Float(rhs) + lhs;
  }


  inline Float Float::operator-(const Float& rhs)const noexcept
  {
    Measurement::inc(Measurement::MeasType::subtraction);
    return { m_value - rhs.m_value };
  }


  inline Float& Float::operator-=(const Float& rhs) noexcept
  {
    Measurement::inc(Measurement::MeasType::subtraction);
    Measurement::inc(Measurement::MeasType::assignment);
    m_value -= rhs.m_value;
    return *this;
  }


  inline Float operator-(float rhs, const Float& lhs)noexcept
  {
    return Float(rhs) - lhs;
  }


  inline Float Float::operator*(const Float& rhs)const noexcept
  {
    Measurement::inc(Measurement::MeasType::multiplication);
    return { m_value * rhs.m_value };
  }


  inline Float& Float::operator*=(const Float& rhs) noexcept
  {
    Measurement::inc(Measurement::MeasType::multiplication);
    Measurement::inc(Measurement::MeasType::assignment);
    m_value *= rhs.m_value;
    return *this;
  }


  inline Float operator*(float rhs, const Float& lhs)noexcept
  {
    return Float(rhs) * lhs;
  }


  inline Float Float::operator/(const Float& rhs)const noexcept
  {
    Measurement::inc(Measurement::MeasType::division);
    return { m_value / rhs.m_value };
  }


  inline Float& Float::operator/=(const Float& rhs) noexcept
  {
    Measurement::inc(Measurement::MeasType::division);
    Measurement::inc(Measurement::MeasType::assignment);
    m_value /= rhs.m_value;
    return *this;
  }


  inline Float operator/(float rhs, const Float& lhs)noexcept
  {
    return Float(rhs) / lhs;
  }


  inline Float& Float::operator++() noexcept
  {
    Measurement::inc(Measurement::MeasType::addition);
    ++m_value;
    return *this;
  }


  inline Float& Float::operator--() noexcept
  {
    Measurement::inc(Measurement::MeasType::subtraction);
    --m_value;
    return *this;
  }


  inline Float Float::operator++(int) noexcept
  {
    Float old = *this;
    operator++();
    return old;
  }


  inline Float Float::operator--(int) noexcept
  {
    Float old = *this;
    operator--();
    return old;
  }


  inline float Float::toFloat() const noexcept
  {
    return m_value;
  }
} // namespace ec
