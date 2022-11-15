cmake_minimum_required(VERSION 3.5)

project(helloworld)


IF(NOT CMAKE_BUILD_TYPE MATCHES Release)
  message(STATUS "CMAKE_BUILD_TYPE not set to Release -> impacts performance")
endif()

find_package(Ceres REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_library(komo komo.h komo.cpp)
target_link_libraries(komo Ceres::ceres Eigen3::Eigen) 

# helloworld
# add_executable(helloworld main.cpp)
# target_link_libraries(helloworld Ceres::ceres)

add_executable(tutorial tutorial.cpp)
target_link_libraries(tutorial Ceres::ceres)


add_executable(komo_main komo_main.cpp)
target_link_libraries(komo_main komo)

add_executable(komo_v2 komo_v2.cpp)
target_link_libraries(komo_v2 Ceres::ceres)

add_executable(komo_v3_deriv komo_v3_deriv.cpp)
target_link_libraries(komo_v3_deriv Ceres::ceres)
