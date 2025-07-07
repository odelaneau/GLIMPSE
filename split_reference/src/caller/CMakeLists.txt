# Create caller library
add_library(SplitCaller caller_algorithm.cpp
                        caller_finalise.cpp
                        caller_initialise.cpp
                        caller_management.cpp
                        caller_parameters.cpp
)

# Link to library
target_link_libraries(SplitCaller PUBLIC Utils Containers Version_tools)

# Link to header files
target_include_directories(SplitCaller  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})