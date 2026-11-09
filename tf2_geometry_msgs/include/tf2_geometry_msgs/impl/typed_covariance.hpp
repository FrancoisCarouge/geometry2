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

#ifndef TF2_GEOMETRY_MSGS__IMPL__TYPED_COVARIANCE_HPP_
#define TF2_GEOMETRY_MSGS__IMPL__TYPED_COVARIANCE_HPP_

#include <concepts>

#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

#include "fcarouge/sparse_linear_algebra.hpp"
#include "fcarouge/typed_linear_algebra.hpp"
#include "mp-units/framework/quantity.h"
#include "mp-units/systems/si.h"

namespace fcarouge
{
/// The adapter of the typed linear algebra to mp-units elements: a quantity element is stored as
/// its numerical value in its own unit, and a stored value is read back as a quantity of that unit.
template<typename To, mp_units::Quantity From>
struct element_caster<To, From>
{
  [[nodiscard]] static constexpr auto operator()(From value) -> To
  {
    static_assert(std::same_as<typename From::rep, std::remove_cvref_t<To>>);
    return value.numerical_value_in(From::unit);
  }
};

template<mp_units::Quantity To, typename From>
struct element_caster<To, From>
{
  [[nodiscard]] static constexpr auto operator()(From value) -> To
  {
    static_assert(std::same_as<typename To::rep, std::remove_cvref_t<From>>);
    return value * To::reference;
  }
};
}  // namespace fcarouge

namespace tf2
{
namespace impl
{
namespace typed_covariance
{
template<auto Reference>
using Quantity = mp_units::quantity<Reference, double>;
using Length = Quantity<mp_units::si::metre>;
using Angle = Quantity<mp_units::si::radian>;
using InverseLength = Quantity<mp_units::one / mp_units::si::metre>;
using InverseAngle = Quantity<mp_units::one / mp_units::si::radian>;

/// The elements of a pose: a position and a rotation vector.
using Pose = std::tuple<Length, Length, Length, Angle, Angle, Angle>;
using InversePose =
  std::tuple<InverseLength, InverseLength, InverseLength, InverseAngle, InverseAngle, InverseAngle>;

/// The covariance of a pose: square metres, metre radians, and square radians.
using DenseStorage = fcarouge::sparse_matrix<double, fcarouge::sparse::dense<6, 6>>;
using Covariance = fcarouge::typed_matrix<DenseStorage, Pose, Pose>;

/// The change of frame of a pose by a rotation, the derivative of the new pose by the old one:
/// the rotation of the position and the rotation of the orientation, two blocks of a pattern of
/// 18 elements.
inline constexpr fcarouge::sparse::pattern rotation_blocks{fcarouge::sparse::shape<6, 6>,
  [](std::size_t row, std::size_t column) {return (row < 3) == (column < 3);}};
using Jacobian =
  fcarouge::typed_matrix<fcarouge::sparse_matrix<double, rotation_blocks>, Pose, InversePose>;

template<std::size_t Row, std::size_t Column>
using CovarianceElement = fcarouge::tla::element_at<Covariance, Row, Column>;

/// The covariance J Σ Jᵀ of the pose in the frame of the rotation `basis`, from the row-major
/// `covariance`.
inline std::array<double, 36> rotate(
  const std::array<double, 36> & covariance, const std::array<std::array<double, 3>, 3> & basis)
{
  Covariance sigma;
  Jacobian jacobian;
  [&]<std::size_t... Index>(std::index_sequence<Index...>) {
    (sigma.template at<Index / 6, Index % 6>(
      covariance[Index] * CovarianceElement<Index / 6, Index % 6>::reference), ...);
  }(std::make_index_sequence<36>{});
  [&]<std::size_t... Index>(std::index_sequence<Index...>) {
    (jacobian.template at<Index / 3, Index % 3>(basis[Index / 3][Index % 3] * mp_units::one), ...);
    (jacobian.template at<Index / 3 + 3, Index % 3 + 3>(
      basis[Index / 3][Index % 3] * mp_units::one), ...);
  }(std::make_index_sequence<9>{});

  const Covariance result{jacobian * sigma * fcarouge::transposed(jacobian)};

  std::array<double, 36> out;
  [&]<std::size_t... Index>(std::index_sequence<Index...>) {
    ((out[Index] = result.template at<Index / 6, Index % 6>().numerical_value_in(
      CovarianceElement<Index / 6, Index % 6>::unit)), ...);
  }(std::make_index_sequence<36>{});
  return out;
}
}  // namespace typed_covariance
}  // namespace impl
}  // namespace tf2

#endif  // TF2_GEOMETRY_MSGS__IMPL__TYPED_COVARIANCE_HPP_
