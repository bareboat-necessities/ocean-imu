cmake_minimum_required(VERSION 3.16)

project(ocean_imu LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Match existing Makefile optimization profile.
add_compile_options(-O3 -funroll-loops -fno-finite-math-only)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NOT WIN32)
  add_compile_options(-march=native)
endif()

set(OCEAN_IMU_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(OCEAN_IMU_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")

# Include project headers.
include_directories("${OCEAN_IMU_SRC_DIR}")

# Eigen include fallback: vendored copy first, then system path.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen")
else()
  include_directories("/usr/include/eigen3")
endif()

set(W3D_SIM_COMMON "${OCEAN_IMU_SRC_DIR}/util/W3dSimCommon.cpp")

add_executable(ahrs-qmekf-sim "${OCEAN_IMU_TESTS_DIR}/ahrs/ahrs-qmekf-sim.cpp" ${W3D_SIM_COMMON})
add_executable(detrend-basic-test "${OCEAN_IMU_TESTS_DIR}/detrend/detrend-basic-test.cpp")
add_executable(freq-track "${OCEAN_IMU_TESTS_DIR}/freq/freq-track.cpp")
add_executable(wave-direction-test "${OCEAN_IMU_TESTS_DIR}/wave_dir/wave-direction-test.cpp")
add_executable(wave-direction-iq-test "${OCEAN_IMU_TESTS_DIR}/wave_dir/wave-direction-iq-test.cpp")
add_executable(wave-encounter-test "${OCEAN_IMU_TESTS_DIR}/wave_dir/wave-encounter-test.cpp")
add_executable(imu_calibrate-test "${OCEAN_IMU_TESTS_DIR}/imu_calibrate/imu_calibrate-test.cpp")
add_executable(kalman_ou_ii-sim "${OCEAN_IMU_TESTS_DIR}/kalman_ou_ii/kalman_ou_ii-sim.cpp" ${W3D_SIM_COMMON})
add_executable(kalman_ou_iii-sim "${OCEAN_IMU_TESTS_DIR}/kalman_ou_iii/kalman_ou_iii-sim.cpp" ${W3D_SIM_COMMON})
add_executable(kalman_ou_common-test "${OCEAN_IMU_TESTS_DIR}/kalman_ou_iii/kalman_ou_common-test.cpp")
add_executable(phi_coeff_boundary-test "${OCEAN_IMU_TESTS_DIR}/kalman_ou_iii/phi_coeff_boundary-test.cpp")
add_executable(pii_observer-adaptive "${OCEAN_IMU_TESTS_DIR}/pii_observer/pii_observer-adaptive.cpp" ${W3D_SIM_COMMON})
