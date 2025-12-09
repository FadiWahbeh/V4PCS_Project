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

add_executable(${PROJECT_NAME}
        src/main.cpp
)

# Eigen fournit une cible IMPORTED appelée Eigen3::Eigen
target_link_libraries(${PROJECT_NAME}
        PRIVATE
        Eigen3::Eigen
)
