cmake_minimum_required(VERSION 3.5)

project(helloworld)

find_package(Ceres REQUIRED)

# helloworld
add_executable(helloworld main.cpp)
target_link_libraries(helloworld Ceres::ceres)

add_executable(komo_1x komo_1x.cpp)
target_link_libraries(komo_1x Ceres::ceres)