#pragma once

#include "Core/Types.h"
#include <String>
#include <fstream>
#include <iostream>
#include <sstream>

namespace V4PCS {
		namespace IO {

            class FileManager {
            public:
                // entré sortie PLY
                static PointCloud loadPLY(const std::string& path);
                static void savePLY(const std::string& path, const PointCloud& cloud);
                static PointCloud loadPointCloudAuto(const std::string& path);

				// lire les fichiers ASCII génériques x y z
                static PointCloud loadASCII_XYZ(const std::string& path);
            };
		}
}