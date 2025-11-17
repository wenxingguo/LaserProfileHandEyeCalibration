#include "LaserProfileHandEyeCalibration.h"
#include "Viewer.hpp"
#include "json.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <thread>

static double Sphere_R = 20.2086; // 不一定是对的
// std::numeric_limits<double>::quiet_NaN();

struct profile_data {
	std::vector<double> robot_pose;

	std::vector<point3d> profile;
	// 分别为轮廓圆心在线激光坐标系的位置
	// c_y 为球心在线激光坐标系下的位置，不在x-z平面
	double c_x = std::numeric_limits<double>::quiet_NaN();
	double c_y = std::numeric_limits<double>::quiet_NaN();
	double c_z = std::numeric_limits<double>::quiet_NaN();
	// 球心在世界坐标系下的位置
	double w_x = std::numeric_limits<double>::quiet_NaN();
	double w_y = std::numeric_limits<double>::quiet_NaN();
	double w_z = std::numeric_limits<double>::quiet_NaN();

	// 轮廓的半径
	double p_r = std::numeric_limits<double>::quiet_NaN();
	bool in_robot_coor = false;
	bool filtered = false;
	bool fitted = false;
	std::vector<bool> removed_point;

	// double c_y_signe = std::numeric_limits<double>::quiet_NaN(); // 保存该条轮廓y的符号，+- 1

	void filter(const double percent = 0.40) {
		// 从中间位置开始拟合一个初始的半径，然后对点进行筛选
		// 假设中间位置的点测量比较准确
		// 该函数只能在相机坐标系下进行
		assert(filtered == false);
		assert(in_robot_coor == false);
		int mid = profile.size() / 2;
		std::vector<point2d> profile_xz;
		for (int index = mid - profile.size() * percent / 2; index < mid + profile.size() * percent / 2; ++index) {
			point2d p;
			p.x = profile[index].x;
			p.y = profile[index].z;
			profile_xz.push_back(p);
		}
		double r;
		point2d c;
		fitCircle(profile_xz, c, r);
		//std::vector<bool> removed_point;
		removed_point.assign(profile.size(), false);
		double threshold = r * 0.001;
		for (int i = 0; i < profile.size(); ++i) {
			point3d& p3 = profile[i];
			double d = std::sqrt((p3.x - c.x) * (p3.x - c.x) + (p3.z - c.y) * (p3.z - c.y)) - r;
			//std::cout << d << std::endl;
			if (std::abs(d) > threshold) {
				removed_point[i] = true;
				p3.r = 255;
				p3.g = 0;
				p3.b = 0;
			}
		}
		c_x = c.x;
		c_z = c.y;
		c_y = 0;
		p_r = r;

		filtered = true;
	}

	void fit() {
		assert(in_robot_coor == false);
		std::vector<point2d> profile_xz;
		for (int index = 0; index < profile.size(); ++index) {
			if (removed_point[index] == false) {
				point2d p;
				p.x = profile[index].x;
				p.y = profile[index].z;
				profile_xz.push_back(p);
			}
		}
		double r;
		point2d c;
		fitCircle(profile_xz, c, r);
		std::cout << "c_x:" << c_x << "->" << c.x << " ";
		std::cout << "c_z:" << c_z << "->" << c.y << " ";
		std::cout << "c_y:" << c_y << "->" << 0 << " ";
		std::cout << "c_r:" << p_r << "->" << r << std::endl;
		c_x = c.x;
		c_z = c.y;
		c_y = 0;
		p_r = r;
		fitted = true;
	}
	void transform2robot_coor()
	{
		assert(in_robot_coor == false);
		Eigen::Vector3d r_vec(robot_pose[3], robot_pose[4], robot_pose[5]);
		Eigen::Vector3d t_vec(robot_pose[0], robot_pose[1], robot_pose[2]);
		Eigen::AngleAxisd r_axis(r_vec.norm(), r_vec.normalized());
		Eigen::Matrix3d r_m = r_axis.matrix();
		Eigen::Vector3d tmp_v;
		for (int i = 0; i < profile.size(); ++i) {
			tmp_v << profile[i].x, profile[i].y, profile[i].z;
			tmp_v = r_m * tmp_v + t_vec;
			profile[i].x = tmp_v(0);
			profile[i].y = tmp_v(1);
			profile[i].z = tmp_v(2);
		}
		tmp_v << c_x, c_y, c_z;
		tmp_v = r_m * tmp_v + t_vec;
		w_x = tmp_v(0);
		w_y = tmp_v(1);
		w_z = tmp_v(2);
		in_robot_coor = true;
	}
};


Sphere  fitLinearLeastSquares(const std::vector<point3d>& points) {
	int n = points.size();
	if (n < 4) {
		throw std::invalid_argument("至少需要4个点来拟合球");
	}

	// 构建矩阵 A 和向量 b
	Eigen::MatrixXd A(n, 4);
	Eigen::VectorXd b(n);

	for (int i = 0; i < n; ++i) {
		double x = points[i].x;
		double y = points[i].y;
		double z = points[i].z;

		A(i, 0) = 2 * x;
		A(i, 1) = 2 * y;
		A(i, 2) = 2 * z;
		A(i, 3) = -1;

		b(i) = x * x + y * y + z * z;
	}

	// 求解最小二乘问题: A * p = b
	Eigen::Vector4d p = A.colPivHouseholderQr().solve(b);

	// 提取球心坐标
	double x0 = p(0);
	double y0 = p(1);
	double z0 = p(2);
	double C = p(3);

	// 计算半径
	double radius = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0 - C);

	return Sphere{ x0, y0, z0, radius };


}


Sphere fitLinearLeastSquares(const std::vector<profile_data>& pose_data) {
	std::vector<point3d> points;
	for (const profile_data& p_data : pose_data) {
		assert(p_data.in_robot_coor == true);
		for (int i = 0; i < p_data.profile.size(); ++i) {
			if (p_data.removed_point[i] == false) {
				points.push_back(p_data.profile[i]);
			}
		}
		//points.insert(points.end(), p_data.profile.begin(), p_data.profile.end());
	}
	return fitLinearLeastSquares(points);
}

void compute_func(std::vector<profile_data>& pose1_data, Viewer& pc_viewer)
{
	//bool out = true;
	double prev_p_r = 0;
	for (profile_data& p_data : pose1_data) {
		p_data.filter();
		p_data.fit();
		if (p_data.p_r > prev_p_r) {

		}
		else {

		}
		//p_data.transform2robot_coor();
		pc_viewer.addPoints(p_data.profile);
	}
	// 计算每条轮廓的c_y 带有符号

	
	// 拟合球体
	// Sphere s = fitLinearLeastSquares(pose1_data);
	// Sphere_R = s.radius;
	// std::cout << "Sphere R" << s.radius << std::endl;

	std::cout << "end " << std::endl;
}

int main(int argc, char* argv[])
{
	nlohmann::json J;

	std::ifstream f("../pose1.json");
	if (f.is_open()) {
		J = nlohmann::json::parse(f);
	}
	else {
		std::cerr << "open file error" << std::endl;
		std::exit(-1);
	}
	// std::cout << "ok" << std::endl;
	std::vector<profile_data> pose1_data;
	int i = 0;
	for (const nlohmann::json& pose_frofile : J) {
		if (!pose_frofile["profile"].is_null()) {
			const nlohmann::json& pose = pose_frofile["robot pose"];
			assert(!pose.is_null());
			assert(pose.is_array());
			const nlohmann::json& profile = pose_frofile["profile"];
			if (profile.size() < 300) {
				std::cout << "profile point size less than 5, pass" << std::endl;

			}
			else {
				profile_data P_data;
				P_data.robot_pose = pose.get<std::vector<double>>();
				//for (const nlohmann::json& point : profile) {
				for (int i = 0; i < profile.size(); ++i) {
					const nlohmann::json& point = profile.at(i);
					std::vector<double> point_xyz = point.get<std::vector<double>>();
					point3d p_xyz;
					p_xyz.x = point_xyz[0];
					p_xyz.y = point_xyz[1];
					p_xyz.z = point_xyz[2];
					P_data.profile.push_back(p_xyz);
				}
				pose1_data.push_back(P_data);
			}
		}
	}

	Viewer pc_viewer;
	std::thread get_data_thread(compute_func, std::ref(pose1_data), std::ref(pc_viewer));
	pc_viewer.run();
	get_data_thread.join();

	return 0;
}
