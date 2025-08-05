## Create model module
add_library(phase_models imputation_hmm.cpp
                         phasing_hmm.cpp
)

## Include header files
target_include_directories(phase_models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})