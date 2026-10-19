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

template <typename ScalarT, std::size_t DoF = 2, std::size_t Nu = 0>
class Pendulum : public DynamicalSystem<ScalarT, DoF, Nu> {
public:
  constexpr explicit Pendulum(const Matrix<ScalarT, DoF, 3>& coords,
                              const ScalarT& len_one, const ScalarT& m_one,
                              const ScalarT&& len_two, const ScalarT& m_two)
      : DynamicalSystem<ScalarT, DoF, Nu>(coords), length_one(len_one),
        mass_one(m_one), length_two(len_two), mass_two(m_two) {}

  constexpr ScalarT T(const Vector<ScalarT, DoF>& q0,
                      const Vector<ScalarT, DoF>& q1) const {
    const auto theta_one_dot = q1(0) - q0(0);
    const auto v1 = length_one * theta_one_dot;
    const auto T1 = 0.5 * mass_one * v1 * v1;

    const auto theta_two_dot = q1(1) - q0(1);
    const auto v2 = length_two * (theta_one_dot + theta_two_dot);
    const auto T2 =
        0.5 * mass_two *
        (v1 * v1 + v2 * v2 +
         2 * length_one * length_two * theta_one_dot *
             (theta_one_dot + theta_two_dot) * constexpr_cos(q1(1)));
    return T1 + T2;
  }

  constexpr ScalarT V(const Vector<ScalarT, DoF>& q) const {
    const auto V1 = mass_one * G * length_one * (1 - constexpr_cos(q(0)));
    const auto V2 =
        mass_two * G *
        ((length_one + length_two) - length_one * constexpr_cos(q(0)) -
         length_two * constexpr_cos(q(0) + q(1)));
    return V1 + V2;
  }

  constexpr Vector<ScalarT, DoF> Qnc(const Vector<ScalarT, DoF>& q0,
                                     const Vector<ScalarT, DoF>& q1,
                                     const Vector<ScalarT, Nu>& U,
                                     const ScalarT& h) const {
    (void)q0, (void)q1, (void)h, (void)U;
    return Vector<ScalarT, DoF>();
  }

private:
  ScalarT length_one, mass_one, length_two, mass_two;
  static constexpr ScalarT G = 9.81;
};

TEST_CASE("seed_Q offsets each coordinate across three samples",
          "[dynamics][seed_Q]") {
  constexpr auto scalar_Q = seed_Q(Vector<double, 1>{0.0});
  STATIC_REQUIRE(scalar_Q(0, 0) == 0.0);
  STATIC_REQUIRE(scalar_Q(0, 1) == 1e-4);
  STATIC_REQUIRE(scalar_Q(0, 2) == 2e-4);

  constexpr auto Q = seed_Q(Vector<double, 2>{0.0, 0.0});
  STATIC_REQUIRE(Q(0, 0) == 0.0);
  STATIC_REQUIRE(Q(0, 1) == 1e-4);
  STATIC_REQUIRE(Q(0, 2) == 2e-4);

  STATIC_REQUIRE(Q(1, 0) == 0.0);
  STATIC_REQUIRE(Q(1, 1) == 1e-4);
  STATIC_REQUIRE(Q(1, 2) == 2e-4);
}

TEST_CASE("Two-DOF pendulum energy, forces and integration", "[pendulum]") {
  SECTION("T, V and Qnc vanish at rest at the potential minimum") {
    constexpr auto zeros = Vector<double, 2>{0.0, 0.0};
    constexpr auto Q = seed_Q(zeros);

    constexpr auto pendulum = Pendulum<double>(Q, 1.0, 2.0, 1.0, 2.0);

    constexpr auto T = pendulum.T(zeros, zeros);
    constexpr auto V = pendulum.V(zeros);
    constexpr auto Qnc = pendulum.Qnc(zeros, zeros, Vector<double, 0>{}, 0.1);
    STATIC_REQUIRE(T == 0.0);
    STATIC_REQUIRE(V == 0.0);
    STATIC_REQUIRE(Qnc == Vector<double, 2>{});
  }

  SECTION("Solve matches ten steps from angles pi/4 and pi/6 with h = 0.01 s") {
    Vector<double, 10> expected_theta_one{
        0.785144, 0.784337, 0.783077, 0.781367, 0.77921,
        0.776609, 0.77357,  0.770096, 0.766195, 0.761873};

    Vector<double, 10> expected_theta_two{0.52369,  0.52359, 0.52339,  0.52308,
                                          0.522655, 0.52211, 0.521437, 0.52062,
                                          0.51966,  0.51854};
    const double h = 0.01;
    auto Q = seed_Q(Vector<double, 2>{M_PI / 4, M_PI / 6});
    for (std::size_t i = 0; i < 10; ++i) {
      const auto result = Solve(Pendulum<double>(Q, 1.0, 1.0, 1.0, 1.0), h);
      Q = std::get<0>(result);
      REQUIRE(L2<Vector<double, 1>>()(
                  Vector<double, 1>{Q(0, 2) - expected_theta_one(i)}) <= 1e-4);
      REQUIRE(L2<Vector<double, 1>>()(
                  Vector<double, 1>{Q(1, 2) - expected_theta_two(i)}) <= 1e-4);
    }
  }
}
