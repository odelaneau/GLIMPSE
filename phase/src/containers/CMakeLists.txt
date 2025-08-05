## Create phase containers header
add_library(phase_containers conditioning_set.cpp
                             genotype_set.cpp
                             glimpse_mpileup.cpp
                             haplotype_set.cpp
)

## Link to dependencies
target_link_libraries(phase_containers PRIVATE Containers Utils Objects)