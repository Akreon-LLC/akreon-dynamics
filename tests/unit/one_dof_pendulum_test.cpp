/*
 * Copyright © 2025-2026 Akreon LLC
 * SPDX-License-Identifier: Apache-2.0
 */

#include <catch2/catch_test_macros.hpp>

#include "akreon/dynamics/dynamics.hpp"
#include <akreon/linalg/test_helpers.hpp>

using namespace dynamics;
using namespace unit_test;
using namespace linalg;

template <typename ScalarT, std::size_t DoF = 1, std::size_t Nu = 0>
class Pendulum : public DynamicalSystem<ScalarT, DoF, Nu> {
public:
  constexpr explicit Pendulum(const Matrix<ScalarT, DoF, 3>& coords,
                              const ScalarT& len, const ScalarT& m)
      : DynamicalSystem<ScalarT, DoF, Nu>(coords), length(len), mass(m) {}

  constexpr ScalarT T(const Vector<ScalarT, DoF>& q0,
                      const Vector<ScalarT, DoF>& q1) const {
    const ScalarT v = length * (q1(0) - q0(0));
    return 0.5 * mass * v * v;
  }

  constexpr ScalarT V(const Vector<ScalarT, DoF>& q) const {
    return mass * G * length * (1 - constexpr_cos(q(0)));
  }

  constexpr Vector<ScalarT, DoF> Qnc(const Vector<ScalarT, DoF>& q0,
                                     const Vector<ScalarT, DoF>& q1,
                                     const Vector<ScalarT, Nu>& U,
                                     const ScalarT& h) const {
    (void)q0, (void)q1, (void)h, (void)U;
    return Vector<ScalarT, DoF>();
  }

private:
  ScalarT length;
  ScalarT mass;
  static constexpr ScalarT G = 9.81;
};

TEST_CASE("One-DOF pendulum energy, forces and integration", "[pendulum]") {
  SECTION("T, V and Qnc vanish at rest at the potential minimum") {
    constexpr auto Q = seed_Q(Vector<double, 1>{0.0});

    constexpr auto pendulum = Pendulum<double>(Q, 1.0, 2.0);

    constexpr auto T =
        pendulum.T(Vector<double, 1>{0.0}, Vector<double, 1>{0.0});
    constexpr auto V = pendulum.V(Vector<double, 1>{0.0});
    constexpr auto Qnc =
        pendulum.Qnc(Vector<double, 1>{0.0}, Vector<double, 1>{0.0},
                     Vector<double, 0>{}, 0.1);
    STATIC_REQUIRE(T == 0.0);
    STATIC_REQUIRE(V == 0.0);
    STATIC_REQUIRE(Qnc == Vector<double, 1>{});
  }

  SECTION("Solve matches ten steps from theta = pi/4 with h = 0.01 s") {
    Vector<double, 10> expected{0.784905, 0.783618, 0.781638, 0.778967,
                                0.775607, 0.77156,  0.766829, 0.761417,
                                0.755329, 0.748568};
    const double h = 0.01;
    auto Q = seed_Q(Vector<double, 1>{M_PI / 4});
    for (std::size_t i = 0; i < 10; ++i) {
      const auto result = Solve(Pendulum<double>(Q, 1.0, 1.0), h);
      Q = std::get<0>(result);
      REQUIRE(L2<Vector<double, 1>>()(
                  Vector<double, 1>{Q(0, 2) - expected(i)}) <= 1e-4);
    }
  }
}
