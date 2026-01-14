cmake_minimum_required(VERSION 3.30)
project(V4PCS_Project)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



include(FetchContent)

# Déclarer Eigen
FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        5.0.0
)

# Récupérer
FetchContent_MakeAvailable(eigen)

add_executable(V4PCS_Project
        src/main.cpp
        src/fichier_manager.cpp
        src/fichier_manager.h
        src/grille_voxel.cpp
        src/grille_voxel.h
        src/voxel_planaire_extraction.cpp
        src/voxel_planaire_extraction.h
        src/patch_builder.cpp
        src/patch_builder.h
        src/patch_pair_generator.cpp
        src/patch_pair_generator.h
        src/v4pcs_algo.cpp
        src/v4pcs_algo.h
)



# Eigen fournit une cible IMPORTED appelée Eigen3::Eigen
target_link_libraries(${PROJECT_NAME}
        PRIVATE
        Eigen3::Eigen
)



target_compile_definitions(${PROJECT_NAME}
        PRIVATE DATA_DIR_INPUT="${CMAKE_SOURCE_DIR}/data/input"
        PRIVATE DATA_DIR_OUTPUT="${CMAKE_SOURCE_DIR}/data/output"
)
