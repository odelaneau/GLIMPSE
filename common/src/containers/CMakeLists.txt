# Create caller library
add_library(Containers bitmatrix.cpp
                      ref_haplotype_set.cpp
                      variant_map.cpp
)

# Link to library
target_link_libraries(Containers PUBLIC Utils Objects IO)

# Link to header files
target_include_directories(Containers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})