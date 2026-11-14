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

#ifndef TF2__TYPED__LINEAR_MATH_HPP_
#define TF2__TYPED__LINEAR_MATH_HPP_

/// A typed linear math beside tf2's untyped LinearMath types: 3-vectors of quantities, the rotation
/// of a vector, the cross product, and a rigid transform of positions, over the typed and sparse
/// linear algebra with mp-units. The untyped types carry any unit; these name it. The operations
/// keep the order of tf2's, so that a computation moved to them keeps its results.
/// A package including it finds mp-units, fcarouge-typed-linear-algebra, and
/// fcarouge-sparse-linear-algebra.

#include <cstddef>
#include <tuple>
#include <utility>

#include "fcarouge/sparse_linear_algebra.hpp"
#include "fcarouge/typed_linear_algebra.hpp"
#include "mp-units/framework/quantity.h"
#include "mp-units/systems/si.h"
#include "tf2/LinearMath/Matrix3x3.hpp"
#include "tf2/LinearMath/Transform.hpp"
#include "tf2/LinearMath/Vector3.hpp"
#include "tf2/impl/typed_linear_algebra.hpp"

namespace tf2
{
namespace typed
{
template<auto Reference>
using Quantity = mp_units::quantity<Reference, double>;
using Length = Quantity<mp_units::si::metre>;
using Duration = Quantity<mp_units::si::second>;
using Speed = Quantity<mp_units::si::metre / mp_units::si::second>;
/// An angular speed, the angle a ratio: mp-units keeps the radian's kind apart, and m rad/s does
/// not add to m/s.
using AngularSpeed = Quantity<mp_units::one / mp_units::si::second>;

/// A 3-vector of one quantity type.
template<typename Type>
using Vector = fcarouge::typed_column_vector<
  fcarouge::sparse_matrix<double, fcarouge::sparse::dense<3, 1>>, Type, Type, Type>;
using Position = Vector<Length>;

/// A rotation of 3-vectors of one quantity type: dimensionless elements.
template<typename Type>
using Rotation = fcarouge::typed_matrix<
  fcarouge::sparse_matrix<double, fcarouge::sparse::dense<3, 3>>, std::tuple<Type, Type, Type>,
  std::tuple<decltype(1. / Type{}), decltype(1. / Type{}), decltype(1. / Type{})>>;

/// The typed vector of the untyped `value`, read in the unit of `Type`.
template<typename Type>
Vector<Type> typed(const tf2::Vector3 & value)
{
  return Vector<Type>{value.x() * Type::reference, value.y() * Type::reference,
    value.z() * Type::reference};
}

/// The untyped vector of the typed `value`, in its unit.
template<typename Type>
tf2::Vector3 untyped(const Vector<Type> & value)
{
  return tf2::Vector3{value.template at<0>().numerical_value_in(Type::unit),
    value.template at<1>().numerical_value_in(Type::unit),
    value.template at<2>().numerical_value_in(Type::unit)};
}

/// The rotation of vectors of `Type` by the basis `basis`.
template<typename Type>
Rotation<Type> rotation(const tf2::Matrix3x3 & basis)
{
  Rotation<Type> value;
  [&]<std::size_t... Index>(std::index_sequence<Index...>) {
    (value.template at<Index / 3, Index % 3>(basis[Index / 3][Index % 3] * mp_units::one), ...);
  }(std::make_index_sequence<9>{});
  return value;
}

/// The vector `value` rotated by `basis`.
template<typename Type>
Vector<Type> rotate(const tf2::Matrix3x3 & basis, const Vector<Type> & value)
{
  return Vector<Type>{rotation<Type>(basis) * value};
}

/// The cross product, component by component as tf2::Vector3::cross.
template<typename Lhs, typename Rhs>
auto cross(const Vector<Lhs> & lhs, const Vector<Rhs> & rhs)
{
  using Product = decltype(Lhs {} *Rhs {});
  return Vector<Product>{
    lhs.template at<1>() * rhs.template at<2>() - lhs.template at<2>() * rhs.template at<1>(),
    lhs.template at<2>() * rhs.template at<0>() - lhs.template at<0>() * rhs.template at<2>(),
    lhs.template at<0>() * rhs.template at<1>() - lhs.template at<1>() * rhs.template at<0>()};
}

/// The position `point` transformed by the rigid transform `transform`: its basis rotates, its
/// origin in metres translates, as tf2::Transform's product.
inline Position apply(const tf2::Transform & transform, const Position & point)
{
  return Position{rotate(transform.getBasis(), point) + typed<Length>(transform.getOrigin())};
}
}  // namespace typed
}  // namespace tf2

#endif  // TF2__TYPED__LINEAR_MATH_HPP_
