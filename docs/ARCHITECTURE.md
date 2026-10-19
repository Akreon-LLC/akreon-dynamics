# Akreon Dynamics Architecture

This guide describes physical model composition and propagation data flow.
For experimental status, validation limits, and future work, see the
[support roadmap](ROADMAP.md).

## Source organization

Paths below are relative to `include/akreon/dynamics/`.

| Path | Responsibility |
| --- | --- |
| `dynamical_system.hpp` | Plain model-state struct, fixed-size coordinate/input storage, and coordinate seeding |
| `del.hpp` | Discrete Lagrangian, finite-difference derivatives, discrete forces, propagation, and linearization |
| `dynamics.hpp` | Public umbrella header |

## Experimental propagation flow

1. A model inherits the `DynamicalSystem` state struct and supplies `T`, `V`, and
   `Qnc` functions over fixed-size coordinates and inputs. Coordinate history `Q`
   and inputs `U` are public values; models use ordinary value semantics.
2. Coordinate history and inputs define the discrete Euler–Lagrange residual.
3. `del.hpp` estimates derivatives using finite differences and uses Linalg QR
   solves to compute coordinate corrections and state/input linearization matrices.
4. `Solve` returns updated coordinates and matrices `A` and `B`. The caller owns
   the simulation loop and how returned values are persisted or consumed.

## Dependency boundary

Akreon Linalg owns the matrices and numerical solves. Dynamics owns physical
model interfaces and propagation logic. Returned linearization matrices can be
consumed externally without requiring control-model types here. This repository
does not depend on Akreon Convex or Akreon Control.

## API reference

Use Doxygen's namespace, class, and file indexes for the model interface and
propagation routines documented in the headers.
