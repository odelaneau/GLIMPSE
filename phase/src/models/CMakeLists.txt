## Create model module
add_library(phase_models imputation_hmm.cpp
                         phasing_hmm.cpp
)

## Enable AVX 
target_compile_options(phase_models PRIVATE -mavx -mfma -mavx2 -Wno-ignored-attributes)

## Link to common library
target_link_libraries(phase_models PRIVATE Utils Containers)

## Link header file
target_include_directories(phase_models PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/.." )
