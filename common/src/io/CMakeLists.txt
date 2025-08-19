# Create IO library
add_library(IO  gmap_reader.cpp
                ref_genotype_reader.cpp
)

# Link to library
target_link_libraries(IO Utils Containers)

# Link to header files
target_include_directories(IO PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})