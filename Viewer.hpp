#ifndef VIEWER_HPP
#define VIEWER_HPP
#include "geometry.h"
#include <atomic>
#include <mutex>
#include <pcl/visualization/pcl_visualizer.h>

class Viewer {
public:
    Viewer(const std::string& window_title_ = "pointcloud viewer",
        const std::string& pc_id_ = "showed-pointcloud")
        : update(false)
        , stop_flag(false)
        , pc_id(pc_id_)
    {
        viewer_ptr.reset(new pcl::visualization::PCLVisualizer(window_title_));
        temp_pc.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
    }
    ~Viewer()
    {
    }

    void run()
    {
        if (!viewer_ptr->contains(pc_id)) {
            viewer_ptr->addPointCloud(temp_pc, pc_id);
        }
        while (!stop_flag && !viewer_ptr->wasStopped()) {
            if (update) {
                std::lock_guard<std::mutex> lock_guard(data_mutex);
                viewer_ptr->updatePointCloud(temp_pc, pc_id);
                // temp_pc->clear();
            }
            viewer_ptr->spinOnce();
        }
        stop_flag = false;
    }

    void stop()
    {
        stop_flag = true;
    }

    void addPoints(const std::vector<point3d>& points)
        //const std::vector<unsigned char>& rgb = { 255, 255, 255 })
    {
        std::lock_guard<std::mutex> lock_guard(data_mutex);
        for (const point3d& point : points) {
            pcl::PointXYZRGB point_xyzrgb;
            point_xyzrgb.x = point.x;
            point_xyzrgb.y = point.y;
            point_xyzrgb.z = point.z;
            point_xyzrgb.r = point.r;
            point_xyzrgb.g = point.g;
            point_xyzrgb.b = point.b;
            temp_pc->push_back(point_xyzrgb);
        }
        if (temp_pc->size() > max_point) {
            temp_pc->erase(temp_pc->begin(), temp_pc->begin() + max_point / 3);
        }
        update = true;
    }
    pcl::visualization::PCLVisualizer::Ptr getViewer() const
    {
        return viewer_ptr;
    }

private:
    std::atomic_bool update;
    std::atomic_bool stop_flag;
    std::mutex data_mutex;
    pcl::visualization::PCLVisualizer::Ptr viewer_ptr;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr temp_pc;
    std::string pc_id;
    int max_point = 10000000;
};

#endif /* ifndef VIEWER_HPP */
