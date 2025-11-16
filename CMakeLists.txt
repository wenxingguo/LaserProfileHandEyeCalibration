cmake_minimum_required(VERSION 3.16)
project(LaserProfileHandEyeCalibration)
set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
add_definitions(${PCL_DEFINITIONS})

add_executable(LaserProfileHandEyeCalibration test_LaserProfileHandEyeCalibration.cpp LaserProfileHandEyeCalibration.cpp)
target_link_libraries(LaserProfileHandEyeCalibration Eigen3::Eigen ${PCL_LIBRARIES})

# 对于utf-8编码的中文注释，vs的c++flag中添加/utf-8
