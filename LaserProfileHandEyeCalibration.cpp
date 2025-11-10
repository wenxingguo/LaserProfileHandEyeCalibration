#include "LaserProfileHandEyeCalibration.h"
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

double fitCircle(const std::vector<point2d>& pointArray, point2d& center, double& radius)
{
    int N = pointArray.size();
    assert(N >= 3);

    double sumX = 0.0;
    double sumY = 0.0;
    double sumX2 = 0.0;
    double sumY2 = 0.0;
    double sumX3 = 0.0;
    double sumY3 = 0.0;
    double sumXY = 0.0;
    double sumXY2 = 0.0;
    double sumX2Y = 0.0;

    for (int pId = 0; pId < N; ++pId) {
        sumX += pointArray[pId].x;
        sumY += pointArray[pId].y;

        double x2 = pointArray[pId].x * pointArray[pId].x;
        double y2 = pointArray[pId].y * pointArray[pId].y;
        sumX2 += x2;
        sumY2 += y2;

        sumX3 += x2 * pointArray[pId].x;
        sumY3 += y2 * pointArray[pId].y;
        sumXY += pointArray[pId].x * pointArray[pId].y;
        sumXY2 += pointArray[pId].x * y2;
        sumX2Y += x2 * pointArray[pId].y;
    }

    double C, D, E, G, H;
    double a, b, c;

    C = N * sumX2 - sumX * sumX;
    D = N * sumXY - sumX * sumY;
    E = N * sumX3 + N * sumXY2 - (sumX2 + sumY2) * sumX;
    G = N * sumY2 - sumY * sumY;
    H = N * sumX2Y + N * sumY3 - (sumX2 + sumY2) * sumY;

    a = (H * D - E * G) / (C * G - D * D);
    b = (H * C - E * D) / (D * D - G * C);
    c = -(a * sumX + b * sumY + sumX2 + sumY2) / N;

    center.x = -a / 2.0;
    center.y = -b / 2.0;
    radius = std::sqrt(a * a + b * b - 4 * c) / 2.0;

    double err = 0.0;
    double r2 = radius * radius;
    for (int pId = 0; pId < N; ++pId) {
        double dx = pointArray[pId].x - center.x;
        double dy = pointArray[pId].y - center.y;
        err += (dx * dx + dy * dy) - r2;
    }
    return err;
}

point3d calculateSphereCenterin_LPcoor(const point2d& circle_center, const double circle_radius, const double sphere_radius, const bool y_flag)
{

    point3d sphere_center;
    sphere_center.x = circle_center.x;
    sphere_center.z = circle_center.y; // circle_center的y实际上是线激光轮廓仪的z
    double y = std::sqrt(sphere_radius * sphere_radius - circle_radius * circle_radius);
    if (y_flag)
        sphere_center.y = y;
    else
        sphere_center.y = -y;
    return sphere_center;
}

std::vector<double> calibrationLP_R(const std::vector<std::vector<double>>& robot_poses, const std::vector<point3d>& sphere_centers_in_Lscoor)
{
    assert(robot_poses.size() == sphere_centers_in_Lscoor.size());

    Eigen::Matrix3Xd robot_xyz;
    robot_xyz.resize(3, robot_poses.size());

    Eigen::Matrix3Xd robot_rxyz;
    robot_rxyz.resize(3, robot_poses.size());

    Eigen::Matrix3Xd sphere_center_xyz;
    sphere_center_xyz.resize(3, robot_poses.size());

    // 输入数据到矩阵
    for (int i = 0; i < robot_poses.size(); ++i) {
        const std::vector<double>& robot_pose = robot_poses[i];
        const point3d& sphere_center = sphere_centers_in_Lscoor[i];
        robot_xyz.col(i) << robot_pose[0], robot_pose[1], robot_pose[2];
        robot_rxyz.col(i) << robot_pose[3], robot_pose[4], robot_pose[5];
        sphere_center_xyz.row(i) << sphere_center.x, sphere_center.y, sphere_center.z;
    }

    // 姿态不变，计算平均值
    Eigen::Vector3d robot_rxyz_vec = robot_rxyz.rowwise().mean();
    // 计算旋转矩阵
    Eigen::AngleAxisd rotation_axis(robot_rxyz_vec.norm(), robot_rxyz_vec.normalized());
    Eigen::Matrix3d R_tool02Base = rotation_axis.matrix();

    // 计算差分
    Eigen::Matrix3Xd d_pc = sphere_center_xyz.rightCols(sphere_center_xyz.cols() - 1)
        - sphere_center_xyz.leftCols(sphere_center_xyz.cols() - 1); // 球心坐标的差分
    Eigen::Matrix3Xd d_tr = robot_xyz.rightCols(robot_xyz.cols() - 1)
        - robot_xyz.leftCols(robot_xyz.cols() - 1); // 机器人偏移的差分

    /*****************************************************************************************************/
    // 本质上求解   XA = B  限制X为正交矩阵
    // 求解放方程 R_tool02Base * R_c * d_pc = d_tr
    //             3 x 3         3x3   3xn     nx3
    // R_c * d_pc = R_tool02Base^T * d_tr
    // d_pc^T * R_c^T = d_tr^T * R_tool02Base
    // R_c^T = pinv(d_pc^T) * d_tr^T * R_tool02Base
    //  3x3    3xn            nx3      3x3
    // 求R_c （手眼的旋转部分）

    // svd分解计算伪逆
    // A   =  U  *  S  * V^T
    // nx3 = nxn * nx3 * 3x3
    // A^+ = V  * S^+ * U^T
    // 3xn  3x3 * 3xn * nxn

    Eigen::JacobiSVD<Eigen::MatrixX3d> svd(d_pc.transpose(),
        Eigen::ComputeFullU | Eigen::ComputeFullV);

    // 获取奇异值向量
    const Eigen::VectorXd& singularValues = svd.singularValues();
    double epsilon = 1e-5;
    // 计算S矩阵的伪逆
    Eigen::VectorXd invSingularValues(singularValues.size());
    for (int i = 0; i < singularValues.size(); ++i) {
        // 判断奇异值是否大于阈值，若大于则取其倒数，否则设为0
        if (singularValues(i) > epsilon) {
            invSingularValues(i) = 1.0 / singularValues(i);
        } else {
            invSingularValues(i) = 0.0;
        }
    }
    Eigen::Matrix3Xd S_pinv = Eigen::Matrix3Xd::Zero(3, d_pc.cols());
    for (int i = 0; i < invSingularValues.size(); ++i) {
        S_pinv(i, i) = invSingularValues(i);
    }

    // 计算原矩阵的伪逆: A^+ = V * S^+ * U^T
    Eigen::Matrix3Xd d_pc_T_pinv = svd.matrixV() * S_pinv * svd.matrixU().transpose();
    Eigen::Matrix3d R_c = (d_pc_T_pinv * d_tr.transpose() * R_tool02Base).transpose();
    /*****************************************************************************************************/
    /*
     * 可以将上述方法替换为使用克罗内科积的方法
     *
     */
    /*****************************************************************************************************/

    // 需要保证R_c为正交矩阵
    // 使用SVD分解
    // R_c = U * S * V^T
    // 3x3  3x3 3x3  3x3
    // 其中S接近于单位矩阵，这是我们想要的
    // 可以判断S与单位矩阵的差异来确定标定的效果

    Eigen::JacobiSVD<Eigen::Matrix3d> svd3x3(R_c, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d singularValues_rc = svd3x3.singularValues();
    // 判断singularValues是否与1相差很大
    double error = (Eigen::Vector3d::Ones() - singularValues_rc).squaredNorm(); // 差值的平方和
    if (error > 1e-3) {
        // 发出警告
    }

    // 矫正R_c = U * V^T 确保其为正交矩阵
    R_c = svd3x3.matrixU() * svd3x3.matrixV().transpose();

    // 转化为rx ry rz
    Eigen::AngleAxisd R_c_axis;
    R_c_axis.fromRotationMatrix(R_c);
    Eigen::Vector3d rotation_vector = R_c_axis.axis() * R_c_axis.angle();
    return { rotation_vector(0), rotation_vector(1), rotation_vector(2) };
}

Eigen::Matrix3d rotation_vector2mat(const Eigen::Vector3d& robot_rxyz_vec)
{
    Eigen::AngleAxisd rotation_axis(robot_rxyz_vec.norm(), robot_rxyz_vec.normalized());
    return rotation_axis.matrix();
}

std::vector<double> calibrationLP_t(const std::vector<double>& R_c_rxyz,
    const std::vector<std::vector<double>>& robot_poses,
    const std::vector<point3d>& sphere_centers_in_Lscoor)
{
    assert(robot_poses.size() == sphere_centers_in_Lscoor.size());
    // 每一个robot_pose的位姿都不一样
    Eigen::MatrixX3d R_tool02Base_s;
    R_tool02Base_s.resize(3 * robot_poses.size(), 3); // 3n x 3

    Eigen::VectorXd pc;
    pc.resize(3 * robot_poses.size());

    Eigen::VectorXd tr;
    tr.resize(3 * robot_poses.size());

    // 准备数据
    // 手眼矩阵的旋转部分
    Eigen::Vector3d rxyz_vec_Rc(R_c_rxyz[0], R_c_rxyz[1], R_c_rxyz[2]);
    Eigen::Matrix3d Rc_rmat = rotation_vector2mat(rxyz_vec_Rc);

    for (int i = 0; i < robot_poses.size(); ++i) {
        const std::vector<double>& robot_pose = robot_poses[i];
        // 机器人末端位置
        tr.block(3 * i, 0, 3, 1) << robot_pose[0], robot_pose[1], robot_pose[2];
        // 标定球球心坐标
        Eigen::Vector3d pc_vec(sphere_centers_in_Lscoor[i].x,
            sphere_centers_in_Lscoor[i].y,
            sphere_centers_in_Lscoor[i].z);
        // 机器人末端工具姿态
        Eigen::Vector3d R_r_rxyz(robot_pose[3], robot_pose[4], robot_pose[5]);
        Eigen::Matrix3d R_r_i = rotation_vector2mat(R_r_rxyz);
        // 计算组合
        pc.block(3 * i, 0, 3, 1) = R_r_i * Rc_rmat * pc_vec;
        R_tool02Base_s.block(3 * i, 0, 3, 3) = R_r_i;
    }

    // 计算差分
    Eigen::MatrixX3d d_R_tool02Base_s = R_tool02Base_s.topRows(R_tool02Base_s.rows() - 3) - R_tool02Base_s.bottomRows(R_tool02Base_s.rows() - 3); // A 3n-3 x 3
    Eigen::VectorXd d_tr = tr.bottomRows(tr.rows() - 3) - tr.topRows(tr.rows() - 3);
    Eigen::VectorXd d_pc = pc.bottomRows(pc.rows() - 3) - pc.topRows(pc.rows() - 3);
    Eigen::VectorXd B = d_tr + d_pc; // B 3n-3 x1

    // 解方程 A*x = B => x = A^+ * B
    // 计算A的伪逆
    Eigen::JacobiSVD<Eigen::MatrixX3d> svd(d_R_tool02Base_s, Eigen::ComputeFullU | Eigen::ComputeFullV);

    const Eigen::VectorXd& singularValues = svd.singularValues();
    double epsilon = 1e-5;
    // 计算A矩阵的伪逆
    Eigen::VectorXd invSingularValues(singularValues.size());
    for (int i = 0; i < singularValues.size(); ++i) {
        // 判断奇异值是否大于阈值，若大于则取其倒数，否则设为0
        if (singularValues(i) > epsilon) {
            invSingularValues(i) = 1.0 / singularValues(i);
        } else {
            invSingularValues(i) = 0.0;
        }
    }
    Eigen::Matrix3Xd S_pinv = Eigen::Matrix3Xd::Zero(3, d_R_tool02Base_s.rows());
    for (int i = 0; i < invSingularValues.size(); ++i) {
        S_pinv(i, i) = invSingularValues(i);
    }

    // 计算原矩阵的伪逆: A^+ = V * S^+ * U^T
    Eigen::Matrix3Xd d_R_tool02Base_s_pinv = svd.matrixV() * S_pinv * svd.matrixU().transpose();

    Eigen::Vector3d t_c = d_R_tool02Base_s_pinv * B;
    return { t_c(0), t_c(1), t_c(2) };
}
