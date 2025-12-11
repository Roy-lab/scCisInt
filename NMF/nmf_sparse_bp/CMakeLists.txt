cmake_minimum_required(VERSION 3.12)
project(nmf_sparse_bp CXX C)

set(CMAKE_CXX_STANDARD 11)

# --- Find GSL automatically ---
find_package(GSL REQUIRED)

# --- Compiler flags ---
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

# --- Random SVD C files ---
set(RANDOM_SVD_C
        modules/random_svd/matrix_vector_functions_gsl.c
        modules/random_svd/low_rank_svd_algorithms_gsl.c
)
add_library(random_svd_objs OBJECT ${RANDOM_SVD_C})
set_target_properties(random_svd_objs PROPERTIES LINKER_LANGUAGE C)
# Add include paths for GSL
target_include_directories(random_svd_objs
        PRIVATE
        modules/random_svd
        ${GSL_INCLUDE_DIRS}   # <-- THIS IS THE FIX
)

# --- All module C++ files ---
file(GLOB MODULE_CPP modules/*.cpp)

# --- Main executable ---
add_executable(run_nmf_bp
        run_nmf.cpp
        ${MODULE_CPP}
        $<TARGET_OBJECTS:random_svd_objs>
)

# --- Include directories ---
target_include_directories(run_nmf_bp PRIVATE
        modules
        modules/random_svd
)

# --- Link GSL ---
target_link_libraries(run_nmf_bp PRIVATE
        GSL::gsl
        GSL::gslcblas
)
