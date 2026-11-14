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
#include <utility>

#include "mp-units/systems/si.h"
#include "tf2/LinearMath/Matrix3x3.hpp"
#include "tf2/LinearMath/Vector3.hpp"
#include "tf2/typed/linear_math.hpp"

namespace tf2
{
namespace impl
{
namespace typed_screw
{
/// The screw (first, second) in the frame of the rotation `basis` and the translation `origin`
/// [m]: (R first, R second + origin x R first), on tf2's typed linear math. The second half
/// carries a length more than the first: an angular velocity [1/s] and a linear velocity [m/s],
/// or a force [N] and a torque [N m].
template<auto First>
std::pair<std::array<double, 3>, std::array<double, 3>> transform(
  const std::array<std::array<double, 3>, 3> & basis, const std::array<double, 3> & origin,
  const std::array<double, 3> & first, const std::array<double, 3> & second)
{
  using A = tf2::typed::Quantity<First>;
  using B = tf2::typed::Quantity<First * mp_units::si::metre>;
  using tf2::typed::Length;
  using tf2::typed::Vector;

  const tf2::Matrix3x3 rotation{basis[0][0], basis[0][1], basis[0][2], basis[1][0], basis[1][1],
    basis[1][2], basis[2][0], basis[2][1], basis[2][2]};
  const Vector<A> a{tf2::typed::rotate(rotation, tf2::typed::typed<A>(
      tf2::Vector3{first[0], first[1], first[2]}))};
  const Vector<B> b{tf2::typed::rotate(rotation, tf2::typed::typed<B>(
      tf2::Vector3{second[0], second[1], second[2]}))};
  const Vector<B> lever{tf2::typed::cross(
      tf2::typed::typed<Length>(tf2::Vector3{origin[0], origin[1], origin[2]}), a)};

  const tf2::Vector3 out_first = tf2::typed::untyped(a);
  const tf2::Vector3 out_second = tf2::typed::untyped(Vector<B>{b + lever});
  return {{out_first.x(), out_first.y(), out_first.z()},
    {out_second.x(), out_second.y(), out_second.z()}};
}
}  // namespace typed_screw
}  // namespace impl
}  // namespace tf2

#endif  // TF2_GEOMETRY_MSGS__IMPL__TYPED_SCREW_HPP_
