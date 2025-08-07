# Create chunker library
add_library(Chunker src/chunker/chunker_management.cpp
                    src/chunker/chunker_algorithm.cpp
                    src/chunker/chunker_io.cpp
                    src/chunker/chunker_map.cpp
                    src/chunker/chunker_parameters.cpp
)

# Link to packages
target_link_libraries(Chunker PUBLIC Utils 
                                     IO 
                                     Version_tools
                                     Boost::program_options
                                     Boost::iostreams
                                     Boost::regex
                                     ZLIB::ZLIB
                                     ${HTS_LIB}
)

# Link to header files
target_include_directories(Chunker PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/chunker)