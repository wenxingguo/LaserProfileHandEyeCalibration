#include "LaserProfileHandEyeCalibration.h"
#include "Viewer.hpp"
#include "json.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <thread>

struct profile_data {
    std::vector<double> robot_pose;

    std::vector<point3d> profile;
    double c_x;
    double c_y;
    double c_z;
    double r;
    bool in_robot_coor;

    void transform2robot_coor()
    {
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
        c_x = tmp_v(0);
        c_y = tmp_v(1);
        c_z = tmp_v(2);
        in_robot_coor = true;
    }
    // std::vector<double> profile_x;
    //  std::vector<double> profile_y; // 恒为 0
    // std::vector<double> profile_z;
};

void add_data_func(std::vector<profile_data>& pose1_data, Viewer& pc_viewer)
{
    bool out = true;
    for (profile_data& p_data : pose1_data) {
        std::vector<point2d> profile_xz;
        for (int index = 0; index < p_data.profile.size(); ++index) {
            point2d p;
            p.x = p_data.profile[index].x;
            p.y = p_data.profile[index].z;
            profile_xz.push_back(p);
        }
        point2d c;
        double r;
        fitCircle(profile_xz, c, r);
        p_data.c_x = c.x;
        p_data.c_z = c.y;
        p_data.c_y = 0;
        p_data.transform2robot_coor();

        if (out) {
            pc_viewer.addPoints(p_data.profile, { 255, 0, 0 });
            out = false;
        } else {
            pc_viewer.addPoints(p_data.profile, { 0, 0, 255 });
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    std::cout << "end " << std::endl;
}

int main(int argc, char* argv[])
{
    nlohmann::json J;

    std::ifstream f("../pose1.json");
    if (f.is_open()) {
        J = nlohmann::json::parse(f);
    } else {
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
            if (profile.size() < 5) {
                std::cout << "profile point size less than 5, pass" << std::endl;

            } else {
                profile_data P_data;
                P_data.robot_pose = pose.get<std::vector<double>>();
                for (const nlohmann::json& point : profile) {
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
    std::thread get_data_thread(add_data_func, std::ref(pose1_data), std::ref(pc_viewer));
    pc_viewer.run();
    get_data_thread.join();
    /*std::thread viewer_thread(&Viewer::run, &pc_viewer);
    for (const profile_data& p_data : pose1_data) {
        std::vector<point2d> profile_xz;
        for (int index = 0; index < p_data.profile.size(); ++index) {
            point2d p;
            p.x = p_data.profile[index].x;
            p.y = p_data.profile[index].z;
            profile_xz.push_back(p);
            // if (out) {
            //     std::cout << p.x << ", " << p.y << std::endl;
            // }
        }

        pc_viewer.addPoints(p_data.profile, { 0, 0, 255 });
        if (out) {
            point2d c;
            double r;
            pc_viewer.addPoints(p_data.profile, { 255, 0, 0 });
            fitCircle(profile_xz, c, r);
            std::cout << "profile point szie: " << profile_xz.size() << " " << c.x << " " << c.y << " " << r << std::endl;
            out = false;
            // pcl::visualization::PCLVisualizer::Ptr viewer_ptr = pc_viewer.getViewer();
            // pcl::ModelCoefficients m;
            // m.values.push_back(c.x);
            // m.values.push_back(c.y);
            // m.values.push_back(r);
            // viewer_ptr->addCircle(m);
        }
        std::this_thread::sleep_for(std::chrono::seconds(2));
    }*/
    // pc_viewer.run();

    return 0;
}
