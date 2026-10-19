/*
 * Copyright © 2025-2026 Akreon LLC
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include "akreon/linalg/linalg.hpp"

namespace dynamics {

/** @brief Fixed-size coordinate history and inputs for a physical model.
 * Models supply T(q0, q1), V(q), and Qnc(q0, q1, inputs, h) directly.
 */
template <typename ScalarT, std::size_t DegreesOfFreedom, std::size_t InputDim>
  requires(linalg::Field<ScalarT>)
struct DynamicalSystem {
  using Scalar = ScalarT;
  static constexpr std::size_t DoF = DegreesOfFreedom;
  static constexpr std::size_t Nu = InputDim;

  linalg::Matrix<ScalarT, DoF, 3> Q{};
  linalg::Vector<ScalarT, Nu> U{};
};

template <typename ScalarT, std::size_t DoF>
constexpr auto seed_Q(const linalg::Vector<ScalarT, DoF>& ic,
                      ScalarT perturbation = 1e-4) {
  linalg::Matrix<ScalarT, DoF, 3> Q;
  for (std::size_t i = 0; i < DoF; ++i) {
    Q(i, 0) = ic(i);
    Q(i, 1) = ic(i) + perturbation;
    Q(i, 2) = Q(i, 1) + (Q(i, 1) - Q(i, 0));
  }
  return Q;
}

} // namespace dynamics
