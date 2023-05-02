//! Rohde & Schwarz Engineering Competition 2023
//!
//! (c) Rohde & Schwarz GmbH & Co. KG, Munich
//!
//! This code is needed to run your solution locally. Any changes will be lost once you upload your solution!

#pragma once

#include "ec_measurement.h"
#include "ec_float.h"

#include <cmath>

namespace ec
{

  inline Float ec_cos(Float input)
  {
    Measurement::inc(Measurement::MeasType::sin);
    return { std::cos(input.toFloat()) };
  }


  inline Float ec_acos(ec::Float input)
  {
    Measurement::inc(Measurement::MeasType::asin);
    return { std::acos(input.toFloat()) };
  }


  inline Float ec_sin(Float input)
  {
    Measurement::inc(Measurement::MeasType::sin);
    return { std::sin(input.toFloat()) };
  }


  inline Float ec_asin(ec::Float input)
  {
    Measurement::inc(Measurement::MeasType::asin);
    return { std::asin(input.toFloat()) };
  }

  
  inline Float ec_tan(Float input)
  {
    Measurement::inc(Measurement::MeasType::tan);
    return { std::tan(input.toFloat()) };
  }


  inline Float ec_atan(Float input)
  {
    Measurement::inc(Measurement::MeasType::atan);
    return { std::atan(input.toFloat()) };
  }


  inline Float ec_sqrt(Float input)
  {
    Measurement::inc(Measurement::MeasType::sqrt);
    return { std::sqrt(input.toFloat()) };
  }


  inline Float ec_abs(Float input)
  {
    Measurement::inc(Measurement::MeasType::abs);
    return { std::fabs(input.toFloat()) };
  }

  inline Float ec_max(Float left, Float right)
  {
    Measurement::inc(Measurement::MeasType::max);
    return { std::fmax(left.toFloat(), right.toFloat()) };
  }


  inline Float ec_round(Float input)
  {
    Measurement::inc(Measurement::MeasType::round);
    return { std::round(input.toFloat()) };
  }


  inline Float ec_exp(Float input)
  {
    Measurement::inc(Measurement::MeasType::exp);
    return { std::exp(input.toFloat()) };
  }


  inline Float ec_log(Float input)
  {
    Measurement::inc(Measurement::MeasType::log);
    return { std::log(input.toFloat()) };
  }


  inline Float ec_log10(Float input)
  {
    Measurement::inc(Measurement::MeasType::log10);
    return { std::log10(input.toFloat()) };
  }


  inline Float ec_pow(Float base, Float exp)
  {
    Measurement::inc(Measurement::MeasType::pow);
    return { std::pow(base.toFloat(), exp.toFloat()) };
  }

} // namespace ec
