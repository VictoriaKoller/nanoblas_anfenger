# Demos for nanoblas

set(NANOBLAS_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src")

# Demo: demo_vector
add_executable(demo_vector demo_vector.cpp)
target_include_directories(demo_vector PRIVATE "${NANOBLAS_SRC_DIR}")
target_compile_features(demo_vector PRIVATE cxx_std_20)

# Demo: demo_matrix
add_executable(demo_matrix demo_matrix.cpp)
target_include_directories(demo_matrix PRIVATE "${NANOBLAS_SRC_DIR}")
target_link_libraries(demo_matrix PRIVATE LAPACK::LAPACK)
target_compile_features(demo_matrix PRIVATE cxx_std_20)

# Demo: demo_lapack
add_executable(demo_lapack demo_lapack.cpp)
target_include_directories(demo_lapack PRIVATE "${NANOBLAS_SRC_DIR}")
target_link_libraries(demo_lapack PRIVATE LAPACK::LAPACK)

target_compile_features(demo_lapack PRIVATE cxx_std_20)

# Demo: demo_parallel
add_executable(demo_parallel demo_parallel.cpp)
target_include_directories(demo_parallel PRIVATE "${NANOBLAS_SRC_DIR}")
target_include_directories(demo_parallel PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../HP_Anfenger/src")
target_link_libraries(demo_parallel PRIVATE LAPACK::LAPACK ASC_HPC)
target_compile_features(demo_parallel PRIVATE cxx_std_20)

# Windows: copy openblas DLL for demo_lapack after build
if(WIN32)
    add_custom_command(TARGET demo_lapack POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_BINARY_DIR}/openblas/bin/libopenblas.dll"
            $<TARGET_FILE_DIR:demo_lapack>
    )
endif()

# Install demo executables (optional)
install(TARGETS demo_vector demo_matrix demo_lapack demo_parallel
    RUNTIME DESTINATION nanoblas/demo
)