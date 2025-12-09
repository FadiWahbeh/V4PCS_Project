//
// Created by Enzo Gallet on 09/12/2025.
//

#ifndef FICHIER_MANAGER_H
#define FICHIER_MANAGER_H

#include <string>
#include <vector>
#include <filesystem>



class fichier_manager {
    public:
    explicit fichier_manager(const std::filesystem::path &chemin);

    std::vector<std::filesystem::path> lister_fichiers(const std::string& extension = "") const;


private:
    std::filesystem::path m_chemin;
};



#endif //FICHIER_MANAGER_H
