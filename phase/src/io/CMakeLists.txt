## Create phase IO module
add_library(phase_io genotype_bam_caller.cpp
                     genotype_reader.cpp
                     genotype_writer.cpp
)

## Link to common library
target_link_libraries(phase_io PRIVATE IO)