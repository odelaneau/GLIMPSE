# Create Split_ref library
add_library(SPLIT_REFERENCE INTERFACE)

# Add caller
add_subdirectory(src/caller)

# Add Libraries
target_link_libraries(SPLIT_REFERENCE INTERFACE Containers
                                      INTERFACE IO
                                      INTERFACE Objects
                                      INTERFACE SplitCaller
                                      INTERFACE Boost::serialization
                                      INTERFACE Boost::program_options
                                      INTERFACE Boost::iostreams
                                      INTERFACE Boost::regex
                                      INTERFACE ZLIB::ZLIB
                                      INTERFACE ${HTS_LIB}
)
