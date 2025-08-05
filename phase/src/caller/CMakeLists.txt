## Create phase caller library
add_library(phase_caller caller_algorithm.cpp
                         caller_finalise.cpp
                         caller_initialise.cpp
                         caller_management.cpp
                         caller_parameters.cpp
)

## Link header file
target_include_directories(phase_caller PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

