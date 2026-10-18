# Support and Roadmap {#dynamics_support_roadmap}

For model composition and propagation flow, see the [architecture guide](ARCHITECTURE.md).

## Status

**Experimental.** APIs and numerical methods may change. Regression tests cover
one- and two-DoF pendulums and a spring-mass-damper model; they do not establish
general-purpose simulation support.

The current solver has a bounded iteration loop but no dedicated convergence
result, and some underlying numerical solve flags are discarded.

## Remaining work

- Surface convergence and numerical failure status to callers.
- Validate finite-difference derivatives and returned linearizations.
- Review force discretization and accuracy across time steps and physical models.

## Future directions

- Structure-preserving dynamics for control: extend the discrete Euler–Lagrange
  formulation along the research path below.
- Automatic differentiation for model derivatives and state/input linearization.
- Dynamic mode decomposition with control (DMDc) for data-driven state-space models.
- Koopman operator approximations for lifted dynamics and prediction.
- URDF import into physical models, with state-space linearization about operating points.

These are exploratory directions, not supported features.
There are no delivery commitments. Keep the Linalg-only dependency boundary;
validation should precede broader support claims.

## Structure-preserving dynamics for control

Build on the existing discrete Euler–Lagrange formulation:

- Variational integrators from discrete action principles, including consistent
  control forces and dissipation through discrete Lagrange–d'Alembert equations.
- Study symplecticity, momentum preservation under symmetry, and long-time
  energy behavior; validate forced work/energy balance separately.
- Explore constrained mechanics and Lie-group formulations for rigid bodies.
- Expose discrete transition maps and consistent state/input derivatives for
  MPC and discrete mechanics and optimal control (DMOC). Dynamics owns the
  mechanics and integration; optimization and controller synthesis stay in Control.

References: [Marsden and West, Discrete mechanics and variational integrators](https://www.cds.caltech.edu/~marsden/bib/2001/09-MaWe2001/MaWe2001.pdf)
and [Ober-Blöbaum, Junge, and Marsden, Discrete Mechanics and Optimal Control: an Analysis](https://arxiv.org/abs/0810.1386).
