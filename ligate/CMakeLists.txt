## Create Ligater module
add_library(Ligater src/ligater/ligater_algorithm.cpp
                    src/ligater/ligater_finalise.cpp
                    src/ligater/ligater_initialise.cpp
                    src/ligater/ligater_management.cpp
                    src/ligater/ligater_parameters.cpp
)

## Link dependency
target_link_libraries(Ligater PUBLIC Utils Version_tools ${HTS_LIB})

## Link to header files
target_include_directories(Ligater PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)