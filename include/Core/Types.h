#pragma once

#include "../libs/eigen/Eigen/Dense"
#include <cstdint>
#include <vector>
#include <functional>

namespace V4PCS {

    using Vector3 = Eigen::Vector3d;

    // Hash pour Eigen::Vector3i (utilisable avec std::unordered_map)
    struct Vector3iHash {
        std::size_t operator()(const Eigen::Vector3i& k) const noexcept {
            std::size_t h1 = std::hash<int>()(k.x());
            std::size_t h2 = std::hash<int>()(k.y());
            std::size_t h3 = std::hash<int>()(k.z());
            return ((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1);
        }
    };

    struct Point {
        Vector3 position;
        std::uint8_t r = 255;
        std::uint8_t g = 255;
        std::uint8_t b = 255;
    };
    using PointCloud = std::vector<Point>;

    struct Voxel {
        int id;
        Eigen::Vector3i grid_index;
        PointCloud points;
        Vector3 centroid;
    };

    struct Plane {
        Vector3 normal;
        Vector3 centroid;
        Eigen::Vector3i grid_index;
        int id_voxel;
    };

    struct PlanarPatch {
        int id;
        Vector3 center;
        Vector3 normal;
        std::vector<Plane*> component_planes;
        double area = 0.0;

        std::uint8_t r, g, b;
    };

} // namespace V4PCS
