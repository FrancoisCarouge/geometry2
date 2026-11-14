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

#ifndef TF2_EIGEN__IMPL__TYPED_ADJOINT_HPP_
#define TF2_EIGEN__IMPL__TYPED_ADJOINT_HPP_

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
namespace typed_adjoint
{
template<auto Reference>
using Quantity = mp_units::quantity<Reference, double>;
using Ratio = Quantity<mp_units::one>;
using Length = Quantity<mp_units::si::metre>;
using InverseLength = Quantity<mp_units::one / mp_units::si::metre>;

/// A twist (angular, linear velocity) or a wrench (force, torque), per unit of its first half:
/// its second half carries a length more, the lever arm of the rigid transform.
using Screw = std::tuple<Ratio, Ratio, Ratio, Length, Length, Length>;
using InverseScrew = std::tuple<Ratio, Ratio, Ratio, InverseLength, InverseLength, InverseLength>;
using ScrewVector = fcarouge::typed_column_vector<
  fcarouge::sparse_matrix<double, fcarouge::sparse::dense<6, 1>>,
  Ratio, Ratio, Ratio, Length, Length, Length>;

/// The adjoint of a rigid transform: the rotation on the diagonal blocks, the lever arm [p]x R in
/// metres below, zero above, a pattern of 27 elements.
inline constexpr fcarouge::sparse::pattern lower_blocks{fcarouge::sparse::shape<6, 6>,
  [](std::size_t row, std::size_t column) {return row >= 3 || column < 3;}};
using Adjoint =
  fcarouge::typed_matrix<fcarouge::sparse_matrix<double, lower_blocks>, Screw, InverseScrew>;

/// The twist or wrench `screw` in the frame of the rotation `basis` and the translation `origin`
/// [m]: Ad(T) screw, with Ad(T) = [R 0; [p]x R R].
inline std::array<double, 6> transform(
  const std::array<double, 6> & screw, const std::array<std::array<double, 3>, 3> & basis,
  const std::array<double, 3> & origin)
{
  const Length px{origin[0] * mp_units::si::metre};
  const Length py{origin[1] * mp_units::si::metre};
  const Length pz{origin[2] * mp_units::si::metre};
  const Length zero{0. * mp_units::si::metre};
  // the skew-symmetric matrix of the translation
  using Row = std::array<Length, 3>;
  const std::array<Row, 3> skew{Row{zero, -pz, py}, Row{pz, zero, -px}, Row{-py, px, zero}};

  Adjoint adjoint;
  [&]<std::size_t... Index>(std::index_sequence<Index...>) {
    (adjoint.template at<Index / 3, Index % 3>(basis[Index / 3][Index % 3] * mp_units::one), ...);
    (adjoint.template at<Index / 3 + 3, Index % 3 + 3>(
      basis[Index / 3][Index % 3] * mp_units::one), ...);
    (adjoint.template at<Index / 3 + 3, Index % 3>(
      skew[Index / 3][0] * basis[0][Index % 3] + skew[Index / 3][1] * basis[1][Index % 3] +
      skew[Index / 3][2] * basis[2][Index % 3]), ...);
  }(std::make_index_sequence<9>{});

  ScrewVector in;
  [&]<std::size_t... Index>(std::index_sequence<Index...>) {
    (in.template at<Index>(
      screw[Index] * fcarouge::tla::element_at<ScrewVector, Index, 0>::reference), ...);
  }(std::make_index_sequence<6>{});

  const ScrewVector result{adjoint * in};

  std::array<double, 6> out;
  [&]<std::size_t... Index>(std::index_sequence<Index...>) {
    ((out[Index] = result.template at<Index>().numerical_value_in(
      fcarouge::tla::element_at<ScrewVector, Index, 0>::unit)), ...);
  }(std::make_index_sequence<6>{});
  return out;
}
}  // namespace typed_adjoint
}  // namespace impl
}  // namespace tf2

#endif  // TF2_EIGEN__IMPL__TYPED_ADJOINT_HPP_
