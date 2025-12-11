//
// Created by Enzo Gallet on 09/12/2025.
//

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>


#include "fichier_manager.h"
#include "grille_voxel.h"
#include "patch_builder.h"
#include "voxel_planaire_extraction.h"

int main() {
    auto input = fichier_manager(DATA_DIR_INPUT);
    auto res = input.lister_fichiers(".txt");

    for (const auto &temp : res) {
        std::cout << "Fichier : " << temp << std::endl;
        std::cout << "Construction de la grille de voxel" << std::endl;
        auto grille = grille_voxel(temp);
        double vs = grille.get_voxel_taille();
        std::cout << "Taille de voxel utilisée : " << vs << " m" << std::endl;
        std::cout << "Extraction des voxel planaires" << std::endl;
        double linearity_thresh  = 0.8;
        double curvature_thresh  = 0.1;
        voxel_planaire_extraction extraction = voxel_planaire_extraction(linearity_thresh, curvature_thresh);
        extraction.extraire(grille);

        auto grille_simple = grille.get_grille();
        auto grille_planaire = extraction.get_grille_planaire();


        std::cout << "Construction des patches" << std::endl;
        double angle_thresh_deg = 10.0;      // seuil angulaire
        double distance_thresh  = vs * 1.5;  // seuil de distance basé sur taille voxel
        auto builder = patch_builder(angle_thresh_deg, distance_thresh);
        builder.construire_patches(grille_planaire,grille_simple);
        std::cout << "Nombre de patches obtenus : " << builder.get_patches_taille() << std::endl;
        builder.exporter_patches_ply(temp);


    }


    return 0;
}
