cmake_minimum_required(VERSION 3.5)

project(helloworld)


IF(NOT CMAKE_BUILD_TYPE MATCHES Release)
  message(STATUS "CMAKE_BUILD_TYPE not set to Release -> impacts performance")
endif()

find_package(Ceres REQUIRED)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#---------------------------------------
add_executable(testy testy.cpp)
target_link_libraries(testy Ceres::ceres)

# --------- Unconstrained KOMO ------------
add_executable(testy_komo unconstrained/testy_komo.cpp)
target_link_libraries(testy_komo Ceres::ceres)
set_target_properties(testy_komo PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unconstrained")

add_executable(komo_v2 unconstrained/komo_v2.cpp)
target_link_libraries(komo_v2 Ceres::ceres)
set_target_properties(komo_v2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unconstrained")


add_executable(komo_v3_deriv unconstrained/komo_v3_deriv.cpp)
target_link_libraries(komo_v3_deriv Ceres::ceres)
set_target_properties(komo_v3_deriv PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unconstrained")



# ---------- Constrained KOMO --------------
add_library(spm_lib constrained/nu_loss.h constrained/nu_loss.cpp)
target_link_libraries(spm_lib Ceres::ceres)

add_executable(squared_penalty_method constrained/squared_penalty_method.cpp)
target_link_libraries(squared_penalty_method spm_lib Ceres::ceres)
set_target_properties(squared_penalty_method PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/constrained")