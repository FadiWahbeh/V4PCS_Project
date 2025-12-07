#pragma once
#include "Core/Types.h"

namespace V4PCS {
	class PlaneExtractor {
	public:
		PlaneExtractor(double voxel_size, double linearity_thresh, double curvature_thresh);
		std::vector<Plane> extractPlanes(const std::vector<Voxel>& voxels);
	
	private:
		double voxel_size_;
		double linearity_thresh_;
		double curvature_thresh_;
		bool analyzeVoxel(const PointCloud& voxel_points, Plane& out_plane);
	};
}