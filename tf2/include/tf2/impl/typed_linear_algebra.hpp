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

#ifndef TF2__IMPL__TYPED_LINEAR_ALGEBRA_HPP_
#define TF2__IMPL__TYPED_LINEAR_ALGEBRA_HPP_

/// The adapter of the typed linear algebra to mp-units elements, shared by the packages that type
/// their frame changes. Header only and not used by tf2 itself: a package including it finds
/// mp-units, fcarouge-typed-linear-algebra, and fcarouge-sparse-linear-algebra.

#include <concepts>
#include <type_traits>

#include "fcarouge/sparse_linear_algebra.hpp"
#include "fcarouge/typed_linear_algebra.hpp"
#include "mp-units/framework/quantity.h"

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

#endif  // TF2__IMPL__TYPED_LINEAR_ALGEBRA_HPP_
