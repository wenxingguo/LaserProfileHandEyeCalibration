/*
 * 该文件定义了线激光轮廓仪的手眼标定方法
 * 这里规定线激光轮廓仪的坐标系 —————— x轴： 沿着点标号增大的方向，z轴：沿着高度增大的方向，y轴：与x z构成右手坐标系
 * 使用标定球
 * 在标定过程中首先保持每一个**姿态不变（rx ry rz 不变）**，平移法兰坐标系，获得一组轮廓
 * ***注意球心是否在轮廓的同一测*** 这关乎球心在线激光轮廓仪的坐标y是否带符号
 */
#ifndef LaserProfileHandEyeCalibration_H
#define LaserProfileHandEyeCalibration_H

#include <vector>

struct point2d {
    double x;
    double y;
};

struct point3d {
    double x;
    double y;
    double z;
};

/******************************************************************************
 * Function:         fitCircle
 * Description:      拟合平面上的圆
 * Where:
 *                   const std::vector<std::vector<double>>& points - 轮廓点的X Z坐标
 *                   point2d& center - 返回圆心的位置
 *                   double& radius - 返回半径
 * Return:           double - 偏差 (err);
 * Error:
 *****************************************************************************/
double fitCircle(const std::vector<point2d>& pointArray, point2d& center, double& radius);

/******************************************************************************
 * Function:         calculateSphereCenterin_LPcoor
 * Description:      计算标定球球心在线激光轮廓仪坐标系下的位置，通过y_flag确定y的正负
 * Where:
 *                   const point2d& circle_center - 激光轮廓扫描仪轮廓的圆心
 *                   const double circle_radius - 轮廓的半径
 *                   const double sphere_radius - 标定球的半径
 *                   bool y_flag - 指示y坐标的正负
 * Return:           point3d - 球心在线激光轮廓仪坐标系下的坐标（x y z）
 * Error:
 *****************************************************************************/
point3d calculateSphereCenterin_LPcoor(const point2d& circle_center, const double circle_radius, const double sphere_radius, const bool y_flag);

/******************************************************************************
 * Function:         calibrationLP_R
 * Description:      计算线激光轮廓仪到机器人flange坐标系的手眼矩阵
 *                   输入的robot_pose的rx ry rz应该相同，即机器人工具末端只做平移运动
 * Where:
 *                   const std::vector<std::vector<double>>& robot_poses - 每个拍摄位置机器人tool0在base坐标下的位姿（x y z rx ry rz）rx ry rz为旋转向量，单位为rad
 *                   const std::vector<point3d>& sphere_centers_in_Lscoor - 对应位置下标定球球心在线激光轮廓仪坐标系下的坐标
 * Return:           std::vector<double> - 返回手眼变换的旋转变换（rx ry rz）;
 * Error:
 *****************************************************************************/
std::vector<double> calibrationLP_R(const std::vector<std::vector<double>>& robot_poses, const std::vector<point3d>& sphere_centers_in_Lscoor);

/******************************************************************************
 * Function:         calibrationLP_t
 * Description:      根据之前标定的R_c（线激光轮廓仪的手眼旋转部分）计算手眼的偏移矩阵t
 *                   每一个机器人位姿（R_tool02Base不变）下都会有一组在线激光轮廓仪球心坐标，计算t是应该尽量不要让同一组的数据点一起计算，
 *                   因为存在 R_tool02Base_1 - R_tool02Base_2，如果一样则会产生0
 * Where:
 *                   const std::vector<double>& R_c_rxyz - 计算得到的手眼矩阵旋转部分
 *                   const std::vector<std::vector<double>>& robot_poses -
 *                   const std::vector<point3d>& sphere_centers_in_Lscoor - TODO
 * Return:           std::vector<double> - TODO
 * Error:
 *****************************************************************************/
std::vector<double> calibrationLP_t(const std::vector<double>& R_c_rxyz,
    const std::vector<std::vector<double>>& robot_poses,
    const std::vector<point3d>& sphere_centers_in_Lscoor);

#endif /* ifndef LaserProfileHandEyeCalibration_H */
