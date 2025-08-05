## Create phase IO module
add_library(phase_io genotype_bam_caller.cpp
                     genotype_reader.cpp
                     genotype_writer.cpp
                     gmap_reader.cpp
)

## Link to common library
target_link_libraries(phase_io PRIVATE IO)

## Include header files
target_include_directories(phase_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})