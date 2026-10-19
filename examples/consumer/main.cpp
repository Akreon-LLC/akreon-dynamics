/*
 * Copyright © 2026 Akreon LLC
 * SPDX-License-Identifier: Apache-2.0
 */

#include <akreon/dynamics/dynamics.hpp>

#include <cmath>
#include <iostream>

// A 2 kg mass attached to a 1 N/m spring and a 0.5 N s/m damper.
class SpringMassDamper : public dynamics::DynamicalSystem<double, 1, 0> {
public:
  using Coordinates = linalg::Matrix<double, 1, 3>;
  using Position = linalg::Vector<double, 1>;

  explicit SpringMassDamper(const Coordinates& coordinates)
      : DynamicalSystem(coordinates) {}

  // The discrete Lagrangian divides this displacement-based kinetic term by
  // h^2.
  constexpr double T(const Position& q0, const Position& q1) const {
    const double displacement = q1(0) - q0(0);
    return 0.5 * mass * displacement * displacement;
  }

  constexpr double V(const Position& q) const {
    return 0.5 * stiffness * q(0) * q(0);
  }

  constexpr Position Qnc(const Position& q0, const Position& q1,
                         const linalg::Vector<double, 0>&, double h) const {
    // Match the current experimental discrete-force convention.
    return Position{-damping * (q1(0) - q0(0)) / (2.0 * h)};
  }

  static constexpr double mass = 2.0;
  static constexpr double stiffness = 1.0;
  static constexpr double damping = 0.5;
};

int main() {
  constexpr double h = 0.01;
  constexpr std::size_t steps = 100;
  // Start at a 1 m displacement with equal initial coordinate samples.
  auto coordinates = dynamics::seed_Q(linalg::Vector<double, 1>{1.0}, 0.0);
  for (std::size_t step = 0; step < steps; ++step) {
    const auto [next, A, B] = dynamics::Solve(SpringMassDamper(coordinates), h);
    coordinates = next;
    // Solve also returns local discrete-state A/B linearizations.
    if (!std::isfinite(coordinates(0, 2)) || !std::isfinite(A(0, 0))) {
      std::cerr << "Non-finite propagation result\n";
      return 1;
    }
  }

  const double position = coordinates(0, 2);
  // A coarse consumer smoke check, not an integrator accuracy guarantee.
  if (position <= 0.0 || position >= 1.0) {
    std::cerr << "Expected the damped mass to move toward equilibrium\n";
    return 1;
  }
  std::cout << "Spring-mass-damper: " << steps << " steps at h=" << h
            << " s, displacement=" << position << " m\n";
}
