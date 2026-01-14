//
// Created by Enzo Gallet on 09/12/2025.
//

#include "fichier_manager.h"

fichier_manager::fichier_manager(const std::filesystem::path &chemin)
    : m_chemin(chemin)
{
    if (!std::filesystem::exists(chemin)) {
        throw std::runtime_error("Le répertoire suivante n'existe pas : " + chemin.string());
    }
}

std::vector<std::filesystem::path> fichier_manager::lister_fichiers(const std::string &extension) const {
    std::vector<std::filesystem::path> résultat;

    for (const auto& entry : std::filesystem::directory_iterator(m_chemin)) {
        if (!entry.is_regular_file())
            continue;

        if (extension.empty() || entry.path().extension() == extension)
            résultat.push_back(entry.path());
    }
    return résultat;
}
