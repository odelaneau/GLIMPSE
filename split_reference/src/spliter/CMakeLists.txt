# Create spliter library
add_library(SplitCaller spliter_algorithm.cpp
                        spliter_finalise.cpp
                        spliter_initialise.cpp
                        spliter_management.cpp
                        spliter_parameters.cpp
)

# Link to library
target_link_libraries(SplitCaller PRIVATE Utils Containers Version_tools)

# Link to header files
target_include_directories(SplitCaller  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})