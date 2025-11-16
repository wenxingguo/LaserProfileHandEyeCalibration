#include "LaserProfileHandEyeCalibration.h"
#include "json.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

struct profile_data {
    std::vector<double> robot_pose;
    std::vector<double> profile_x;
    // std::vector<double> profile_y; // 恒为 0
    std::vector<double> profile_z;
};

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
                    P_data.profile_x.push_back(point_xyz[0]);
                    // P_data.profile_y.push_back(point_xyz[1]);
                    P_data.profile_z.push_back(point_xyz[2]);
                }
                pose1_data.push_back(P_data);
            }
        }
    }
    bool out = true;
    for (const profile_data& p_data : pose1_data) {
        assert(p_data.profile_x.size() == p_data.profile_z.size());
        // assert(p_data.profile_x.size() == p_data.profile_z.size());

        std::vector<point2d> profile_xz;
        for (int index = 0; index < p_data.profile_x.size(); ++index) {
            point2d p;
            p.x = p_data.profile_x[index];
            p.y = p_data.profile_z[index];
            profile_xz.push_back(p);
            if (out) {
                std::cout << p.x << ", " << p.y << std::endl;
            }
        }
        out = false;
        point2d c;
        double r;
        fitCircle(profile_xz, c, r);
        // std::cout << "profile point szie: " << profile_xz.size() << " " << c.x << " " << c.y << " " << r << std::endl;
    }

    return 0;
}
