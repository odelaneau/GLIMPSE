# Create IO library
add_library(Objects  genotype.cpp
                     variant.cpp
)

# Link to library
target_link_libraries(Objects Utils)

# Link to header files
target_include_directories(Objects  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})