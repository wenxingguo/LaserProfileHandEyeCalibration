#include "LaserProfileHandEyeCalibration.h"
#include <Eigen/Dense>
#include <iostream>

int main(int argc, char* argv[])
{
    /*point2d p1 = { 0, 1 + 1 };
    point2d p2 = { 1, 0 + 1 };
    point2d p3 = { 0, -1 + 1 };
    std::vector<point2d> points;
    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);

    point2d center;
    double radius;

    double err = fitCircle(points, center, radius);
    std::cout << err << " " << center.x << " " << center.y << " " << radius << std::endl;*/

    Eigen::Matrix3Xd X;
    X.resize(3, 5);
    X << 1, 3, 4, 5, 6,
        2, 3, 4, 67, 89,
        1, 2, 5, 6, 100;
    // Eigen::JacobiSVD<Eigen::Matrix3Xd> svd(X, )
    Eigen::JacobiSVD<Eigen::Matrix3Xd> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV);

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

    // Eigen::MatrixXd dig = invSingularValues.asDiagonal();

    Eigen::MatrixX3d S_pinv = Eigen::MatrixX3d::Zero(X.cols(), 3);
    for (int i = 0; i < invSingularValues.size(); ++i) {
        S_pinv(i, i) = invSingularValues(i);
    }

    // 计算原矩阵的伪逆: A^+ = V * S^+ * U^T
    Eigen::MatrixX3d X_pinv = svd.matrixV() * S_pinv * svd.matrixU().transpose();

    std::cout << X * X_pinv << std::endl;

    std::cout << "U:" << svd.matrixU().rows() << "X" << svd.matrixU().cols() << std::endl;
    std::cout << "V:" << svd.matrixV().rows() << "X" << svd.matrixV().cols() << std::endl;
    // std::cout << "dig:" << dig.rows() << "X" << dig.cols() << std::endl;

    return 0;
}
