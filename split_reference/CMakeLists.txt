# Create Split_ref library
add_library(SPLIT_REFERENCE INTERFACE)

# Add caller
add_subdirectory(src/caller)

# Add Libraries
target_link_libraries(SPLIT_REFERENCE INTERFACE Containers
                                                IO
                                                Objects
                                                SplitCaller
                                                ${HTS_LIB}
)
