## Create phase object
add_library(Phase INTERFACE)

## Link to dependencies
add_subdirectory(src/caller)
add_subdirectory(src/containers)
add_subdirectory(src/io)
add_subdirectory(src/models)

target_link_libraries(Phase INTERFACE   phase_caller
                                        phase_containers
                                        phase_io
                                        phase_models
                                        ${HTS_LIB}
)