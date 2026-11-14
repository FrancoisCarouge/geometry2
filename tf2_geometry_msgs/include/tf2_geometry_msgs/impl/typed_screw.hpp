// Copyright 2026 Open Source Robotics Foundation, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TF2_GEOMETRY_MSGS__IMPL__TYPED_SCREW_HPP_
#define TF2_GEOMETRY_MSGS__IMPL__TYPED_SCREW_HPP_

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

#include "fcarouge/sparse_linear_algebra.hpp"
#include "fcarouge/typed_linear_algebra.hpp"
#include "mp-units/framework/quantity.h"
#include "mp-units/systems/si.h"
#include "tf2/impl/typed_linear_algebra.hpp"

namespace tf2
{
namespace impl
{
namespace typed_screw
{
template<auto Reference>
using Quantity = mp_units::quantity<Reference, double>;

/// A 3-vector of one quantity type.
template<typename Type>
using Vector = fcarouge::typed_column_vector<
  fcarouge::sparse_matrix<double, fcarouge::sparse::dense<3, 1>>, Type, Type, Type>;

/// A rotation of 3-vectors of one quantity type: dimensionless elements.
template<typename Type>
using Rotation = fcarouge::typed_matrix<
  fcarouge::sparse_matrix<double, fcarouge::sparse::dense<3, 3>>, std::tuple<Type, Type, Type>,
  std::tuple<decltype(1. / Type{}), decltype(1. / Type{}), decltype(1. / Type{})>>;

template<typename Type>
Vector<Type> typed(const std::array<double, 3> & value)
{
  return Vector<Type>{value[0] * Type::reference, value[1] * Type::reference,
    value[2] * Type::reference};
}

template<typename Type>
std::array<double, 3> untyped(const Vector<Type> & value)
{
  return {value.template at<0>().numerical_value_in(Type::unit),
    value.template at<1>().numerical_value_in(Type::unit),
    value.template at<2>().numerical_value_in(Type::unit)};
}

template<typename Type>
Rotation<Type> rotation(const std::array<std::array<double, 3>, 3> & basis)
{
  Rotation<Type> value;
  [&]<std::size_t... Index>(std::index_sequence<Index...>) {
    (value.template at<Index / 3, Index % 3>(basis[Index / 3][Index % 3] * mp_units::one), ...);
  }(std::make_index_sequence<9>{});
  return value;
}

/// The screw (first, second) in the frame of the rotation `basis` and the translation `origin`
/// [m]: (R first, R second + origin x R first). The second half carries a length more than the
/// first: an angular velocity [1/s] and a linear velocity [m/s], or a force [N] and a torque
/// [N m]. The angle of an angular velocity is a ratio here: mp-units keeps the radian's kind
/// apart, and m rad/s does not add to m/s.
template<auto First>
std::pair<std::array<double, 3>, std::array<double, 3>> transform(
  const std::array<std::array<double, 3>, 3> & basis, const std::array<double, 3> & origin,
  const std::array<double, 3> & first, const std::array<double, 3> & second)
{
  using A = Quantity<First>;
  using B = Quantity<First * mp_units::si::metre>;
  using Length = Quantity<mp_units::si::metre>;

  const Vector<A> a{rotation<A>(basis) * typed<A>(first)};
  const Vector<B> b{rotation<B>(basis) * typed<B>(second)};
  const Vector<Length> p{typed<Length>(origin)};

  // the lever arm, origin x R first, component by component as tf2::Vector3::cross
  const Vector<B> lever{
    p.template at<1>() * a.template at<2>() - p.template at<2>() * a.template at<1>(),
    p.template at<2>() * a.template at<0>() - p.template at<0>() * a.template at<2>(),
    p.template at<0>() * a.template at<1>() - p.template at<1>() * a.template at<0>()};

  return {untyped(a), untyped(Vector<B>{b + lever})};
}
}  // namespace typed_screw
}  // namespace impl
}  // namespace tf2

#endif  // TF2_GEOMETRY_MSGS__IMPL__TYPED_SCREW_HPP_
