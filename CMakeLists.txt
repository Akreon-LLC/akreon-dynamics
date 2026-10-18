cmake_minimum_required(VERSION 3.30)
project(AkreonDynamics VERSION 0.1.0 LANGUAGES CXX)
option(BUILD_UNIT_TESTS "Build Catch2 unit tests" ${PROJECT_IS_TOP_LEVEL})
option(BUILD_FUNCTIONAL_TESTS "Build functional tests" ${PROJECT_IS_TOP_LEVEL})
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(cmake/AkreonDependencies.cmake)
akreon_resolve_dependency(linalg AkreonLinalg Akreon::Linalg)

add_library(akreon_dynamics INTERFACE)
add_library(Akreon::Dynamics ALIAS akreon_dynamics)
set_target_properties(akreon_dynamics PROPERTIES EXPORT_NAME Dynamics)
target_compile_features(akreon_dynamics INTERFACE cxx_std_26)
target_include_directories(akreon_dynamics INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(akreon_dynamics INTERFACE Akreon::Linalg)

install(TARGETS akreon_dynamics EXPORT AkreonDynamicsTargets)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES LICENSE NOTICE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME})
install(EXPORT AkreonDynamicsTargets NAMESPACE Akreon:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AkreonDynamics)
configure_package_config_file(cmake/AkreonDynamicsConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/AkreonDynamicsConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AkreonDynamics)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/AkreonDynamicsConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/AkreonDynamicsConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/AkreonDynamicsConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AkreonDynamics)

if(PROJECT_IS_TOP_LEVEL AND (BUILD_UNIT_TESTS OR BUILD_FUNCTIONAL_TESTS))
  enable_testing()
  add_subdirectory(tests)
endif()
find_program(CLANG_FORMAT_BIN NAMES clang-format-22 clang-format)
if(PROJECT_IS_TOP_LEVEL AND CLANG_FORMAT_BIN)
  file(GLOB_RECURSE format_files CONFIGURE_DEPENDS include/*.hpp tests/*.cpp tests/*.hpp test-support/*.hpp examples/*.cpp)
  add_custom_target(clang_format_check COMMAND ${CLANG_FORMAT_BIN} --dry-run --Werror ${format_files} VERBATIM)
  add_custom_target(clang_format_fix COMMAND ${CLANG_FORMAT_BIN} -i ${format_files} VERBATIM)
endif()
