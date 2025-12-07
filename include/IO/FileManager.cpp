// IO/FileManager.cpp
#include "IO/FileManager.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <cctype>


namespace fs = std::filesystem;

namespace V4PCS {
    namespace IO {

        // Helper local pour créer un point à partir de x,y,z
        static inline Point makePoint(double x, double y, double z) {
            Point p;
            p.position = Vector3(x, y, z);
            return p;
        }

        // PLY ASCII
        PointCloud FileManager::loadPLY(const std::string& path) {
            PointCloud cloud;

            std::ifstream file(path);
            if (!file.is_open()) {
                std::cerr << "Erreur: Impossible d'ouvrir " << path << std::endl;
                return cloud;
            }

            std::string line;
            bool header_end = false;
            while (std::getline(file, line)) {
                if (line.find("end_header") != std::string::npos) {
                    header_end = true;
                    continue;
                }
                if (!header_end) continue;
                if (line.empty()) continue;

                std::stringstream ss(line);
                double x, y, z;
                if (ss >> x >> y >> z) {
                    cloud.push_back(makePoint(x, y, z));
                }
            }

            std::cout << "[FileManager] PLY : chargé " << cloud.size()
                << " points depuis " << path << std::endl;
            return cloud;
        }

        void FileManager::savePLY(const std::string& path, const PointCloud& cloud) {
            std::ofstream file(path);
            if (!file.is_open()) {
                std::cerr << "Erreur: Impossible d'ouvrir " << path << " en écriture" << std::endl;
                return;
            }

            file << "ply\n";
            file << "format ascii 1.0\n";
            file << "element vertex " << cloud.size() << "\n";
            file << "property float x\n";
            file << "property float y\n";
            file << "property float z\n";
            file << "property uchar red\n";
            file << "property uchar green\n";
            file << "property uchar blue\n";
            file << "end_header\n";

            for (const auto& p : cloud) {
                file << p.position.x() << " "
                    << p.position.y() << " "
                    << p.position.z() << " "
                    << static_cast<int>(p.r) << " "
                    << static_cast<int>(p.g) << " "
                    << static_cast<int>(p.b) << "\n";
            }

            file.close();
            std::cout << "Sauvegarde PLY dans: " << path << std::endl;
        }


        // ASCII générique au moins "x y z" par ligne
        PointCloud FileManager::loadASCII_XYZ(const std::string& path) {
            PointCloud cloud;
            std::ifstream f(path);
            if (!f.is_open()) {
                std::cerr << "Erreur: Impossible d'ouvrir " << path << std::endl;
                return cloud;
            }

            std::cout << "[FileManager] Début lecture ASCII de " << path << std::endl;

            std::string line;
            cloud.reserve(1000000); 

            std::size_t lineCount = 0;
            std::size_t validPoints = 0;

            while (std::getline(f, line)) {
                ++lineCount;
                if (line.empty()) continue;
                if (line[0] == '#') continue;

                std::stringstream ss(line);
                double x, y, z;

                // On lit UNIQUEMENT x y z et on ignore le reste
                if (!(ss >> x >> y >> z)) {
                    continue;
                }

                cloud.push_back(makePoint(x, y, z));
                ++validPoints;

                // Affiche de la progression tous les 1 million de lignes
                if (lineCount % 1000000 == 0) {
                    std::cout << "  lignes lues : " << lineCount
                              << " | points valides : " << validPoints << std::endl;
                }
            }

            std::cout << "[FileManager] Fin lecture ASCII : "
                      << validPoints << " points valides sur ~" << lineCount << " lignes."
                      << std::endl;

            return cloud;
		}

        // -----------------------------------------------------------------------------
        // PCD avec PCL (TUM MLS, etc.)
        // On lit pcl::PointXYZ (ou PointXYZI si tu veux l'intensité), puis on convertit.
        // -----------------------------------------------------------------------------
        
        /*
        PointCloud FileManager::loadPCD(const std::string& path) {
            PointCloud cloud;
            
            pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_cloud(new pcl::PointCloud<pcl::PointXYZ>);
            if (pcl::io::loadPCDFile(path, *pcl_cloud) < 0) {
                std::cerr << "Erreur: Impossible de lire le fichier PCD " << path << std::endl;
                return cloud;
            }

            cloud.reserve(pcl_cloud->size());
            for (const auto& p : pcl_cloud->points) {
                cloud.push_back(makePoint(p.x, p.y, p.z));
            }

            std::cout << "[FileManager] PCD : chargé " << cloud.size()
                << " points depuis " << path << std::endl;
            return cloud;
        }
        */
        // -----------------------------------------------------------------------------
        // Choix automatique du loader en fonction de l'extension
        // -----------------------------------------------------------------------------
        PointCloud FileManager::loadPointCloudAuto(const std::string& path) {
            PointCloud cloud;

            fs::path p(path);
            if (!fs::exists(p)) {
                std::cerr << "Erreur: Fichier inexistant : " << path << std::endl;
                return cloud;
            }

            std::string ext = p.extension().string();
            for (auto& c : ext) c = static_cast<char>(std::tolower(c));

            if (ext == ".ply") {
                return loadPLY(path);
            }
            /*
            else if (ext == ".pcd") {
                return loadPCD(path);
            }*/
            else if (ext == ".txt" || ext == ".xyz") {
                return loadASCII_XYZ(path);
            }
            else {
                std::cerr << "[FileManager] Extension \"" << ext
                    << "\" non reconnue, tentative de lecture ASCII générique.\n";
                return loadASCII_XYZ(path);
            }
        }

    } // namespace IO
} // namespace V4PCS
