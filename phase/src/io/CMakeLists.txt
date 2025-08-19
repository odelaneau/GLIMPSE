## Create phase IO module
add_library(phase_io genotype_bam_caller.cpp
                     genotype_reader.cpp
                     genotype_writer.cpp
)

## Link to common library
target_link_libraries(phase_io PRIVATE IO Version_tools)

## Ignore warning 
target_compile_options(phase_io PRIVATE -Wno-ignored-attributes)

## Link header file
target_include_directories(phase_io PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/.." )