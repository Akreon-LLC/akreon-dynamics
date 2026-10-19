# akreon-dynamics
Experimental dynamical-system abstractions, discrete Euler–Lagrange integration,
state propagation, and physical-system regression models.

**Status: experimental.** See the [support roadmap](docs/ROADMAP.md) for limitations.

Header-only C++26 library. Public target: `Akreon::Dynamics`; installed package:
`AkreonDynamics`.
Public dynamics APIs use the `dynamics` namespace.
Public headers live under `include/akreon/dynamics/`. Only canonical `akreon/...` include paths are supported.

## Build and test

Use CMake 3.30+ (for the exported `cxx_std_26` requirement) and a C++26-capable compiler with constexpr math support.
From this checkout:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DFETCHCONTENT_SOURCE_DIR_AKREON_LINALG=/absolute/path/akreon-linalg
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

`host-default` and `host-functional-only` presets are also available. Pass the
same local-source overrides when configuring presets. Tests default on only for
the top-level project. Disable tests using `-DBUILD_UNIT_TESTS=OFF
-DBUILD_FUNCTIONAL_TESTS=OFF`.

## Dependencies

The library depends on `akreon-linalg` (`Akreon::Linalg`).
Each dependency resolves in this order:

1. An explicit local checkout via `FETCHCONTENT_SOURCE_DIR_AKREON_LINALG`.
2. An installed CMake package found through `CMAKE_PREFIX_PATH`.
3. FetchContent from the corresponding `Akreon-LLC/akreon-<name>` Git repository,
   pinned to the planned `v0.1.0` release.

The build command above uses local checkouts. No sibling-directory layout or
Git submodules are required. Until release tags exist, use local checkouts or
installed packages. For private Git access or another compatible release, set
`AKREON_LINALG_GIT_REPOSITORY` / `AKREON_LINALG_GIT_TAG`.

Unit tests use installed Catch2 3 or fetch v3.8.1. For offline testing, pass
`-DFETCHCONTENT_SOURCE_DIR_CATCH2=/absolute/path/to/Catch2`.

## Install and consume

The [consumer example](examples/consumer) shows how a separate CMake project
uses the installed library. Ensure `akreon-linalg` is installed
into the same prefix (or include their prefixes in `CMAKE_PREFIX_PATH`). After
building above, install this package, then configure, build, and test the example:

```sh
cmake --install build --prefix /absolute/path/akreon-install
cmake -S examples/consumer -B build/consumer \
  -DCMAKE_PREFIX_PATH=/absolute/path/akreon-install
cmake --build build/consumer
ctest --test-dir build/consumer --output-on-failure
```

The example's [CMakeLists.txt](examples/consumer/CMakeLists.txt) finds the package
and links its `consumer` executable:

```cmake
find_package(AkreonDynamics 0.1 CONFIG REQUIRED)
add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE Akreon::Dynamics)
```

Its [main.cpp](examples/consumer/main.cpp) includes `<akreon/dynamics/dynamics.hpp>`
and defines a spring-mass-damper model, advances it with `dynamics::Solve`, and
checks the propagated displacement.
Follow the same pattern in your own project, setting `CMAKE_PREFIX_PATH` to the
installation prefix. `Akreon::Linalg` propagates through
`Akreon::Dynamics`, along with include directories and the C++ standard requirement.
Installed packages resolve their dependencies without fetching source.

## Documentation

The [architecture guide](docs/ARCHITECTURE.md) is the Doxygen landing page.
Contributor workflow is in [CONTRIBUTING.md](CONTRIBUTING.md); test coverage is
in [tests/README.md](tests/README.md).
API and algorithm documentation is generated from the header comments. To build
both together, install Doxygen and Graphviz, then run from this repository:

```sh
doxygen Doxyfile
```

Open `docs/doxygen/html/index.html`. The same generated HTML directory can be
published to GitHub Pages; no separate copy of the architecture guide is needed.

## Citation

If you use this software in research, please cite it using the metadata in
[CITATION.cff](CITATION.cff).

## License and notices

Licensed under [Apache 2.0](LICENSE). See [NOTICE](NOTICE) for attribution
information.
