cmake_minimum_required(VERSION 3.16)
project(LaserProfileHandEyeCalibration)

find_package(Eigen3 REQUIRED)

add_executable(LaserProfileHandEyeCalibration test_LaserProfileHandEyeCalibration.cpp LaserProfileHandEyeCalibration.cpp)
target_link_libraries(LaserProfileHandEyeCalibration Eigen3::Eigen)

# 对于utf-8编码的中文注释，vs的c++flag中添加/utf-8
