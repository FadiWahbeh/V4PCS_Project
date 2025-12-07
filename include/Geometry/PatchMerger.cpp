#include "Geometry/PatchMerger.h"
#include <queue>
#include <cmath>
#include <iostream>
#include <random>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace V4PCS {
    namespace Geometry {

        PatchMerger::PatchMerger(double angle_thresh_deg, double distance_thresh)
            : distance_thresh_(distance_thresh) {
            angle_thresh_cos_ = std::cos(angle_thresh_deg * M_PI / 180.0);
        }

        bool PatchMerger::isSmoothConnection(const Plane& p1, const Plane& p2) const {
            double dot = std::abs(p1.normal.dot(p2.normal));
            if (dot < angle_thresh_cos_) return false;

            double dist1 = std::abs(p1.normal.dot(p1.centroid - p2.centroid));
            double dist2 = std::abs(p2.normal.dot(p2.centroid - p1.centroid));

            if (dist1 > distance_thresh_ || dist2 > distance_thresh_) return false;
            return true;
        }

        std::vector<PlanarPatch> PatchMerger::merge(
            const std::vector<Plane>& voxel_planes,
            const std::unordered_map<Eigen::Vector3i, int, Vector3iHash>& grid_map)
        {
            std::vector<PlanarPatch> patches;
            std::vector<bool> visited(voxel_planes.size(), false);

            // Directions de voisinage (26-connexité)
            std::vector<Eigen::Vector3i> neighbors_offsets;
            for (int x = -1; x <= 1; ++x)
                for (int y = -1; y <= 1; ++y)
                    for (int z = -1; z <= 1; ++z)
                        if (x != 0 || y != 0 || z != 0)
                            neighbors_offsets.emplace_back(x, y, z);

            std::mt19937 rng(42);
            std::uniform_int_distribution<int> color_dist(50, 255);

            for (size_t i = 0; i < voxel_planes.size(); ++i) {
                if (visited[i]) continue;

                PlanarPatch current_patch;
                current_patch.id = static_cast<int>(patches.size());
                current_patch.r = static_cast<std::uint8_t>(color_dist(rng));
                current_patch.g = static_cast<std::uint8_t>(color_dist(rng));
                current_patch.b = static_cast<std::uint8_t>(color_dist(rng));

                std::queue<int> q;
                q.push(static_cast<int>(i));
                visited[i] = true;

                Vector3 sum_centers = Vector3::Zero();
                Vector3 sum_normals = Vector3::Zero();

                while (!q.empty()) {
                    int current_idx = q.front();
                    q.pop();

                    const Plane& current_plane = voxel_planes[current_idx];
                    current_patch.component_planes.push_back(
                        const_cast<Plane*>(&current_plane)
                    );

                    sum_centers += current_plane.centroid;

                    if (sum_normals.isZero(1e-8)) {
                        sum_normals += current_plane.normal;
                    }
                    else {
                        if (current_plane.normal.dot(sum_normals) < 0.0)
                            sum_normals -= current_plane.normal;
                        else
                            sum_normals += current_plane.normal;
                    }

                    // retrouver l'indice de grille associé
                    Eigen::Vector3i current_grid_idx;
                    bool found = false;
                    for (const auto& kv : grid_map) {
                        if (kv.second == current_idx) {
                            current_grid_idx = kv.first;
                            found = true;
                            break;
                        }
                    }
                    if (!found) continue;

                    for (const auto& offset : neighbors_offsets) {
                        Eigen::Vector3i neighbor_grid = current_grid_idx + offset;

                        auto it = grid_map.find(neighbor_grid);
                        if (it != grid_map.end()) {
                            int neighbor_idx = it->second;
                            if (neighbor_idx >= 0 &&
                                neighbor_idx < static_cast<int>(voxel_planes.size()) &&
                                !visited[neighbor_idx])
                            {
                                if (isSmoothConnection(current_plane,
                                    voxel_planes[neighbor_idx])) {
                                    visited[neighbor_idx] = true;
                                    q.push(neighbor_idx);
                                }
                            }
                        }
                    }
                }

                if (!current_patch.component_planes.empty()) {
                    double n = static_cast<double>(current_patch.component_planes.size());
                    current_patch.center = sum_centers / n;
                    if (!sum_normals.isZero(1e-8))
                        current_patch.normal = sum_normals.normalized();
                    else
                        current_patch.normal = current_patch.component_planes.front()->normal;

                    patches.push_back(current_patch);
                }
            }

            return patches;
        }

    } // namespace Geometry
} // namespace V4PCS
