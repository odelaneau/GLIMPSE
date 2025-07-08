# Create concordance library
add_library(Concordance INTERFACE)

# Add directories
add_subdirectory(src/checker)
add_subdirectory(src/containers)


# Link to libraries
target_link_libraries(Concordance INTERFACE ConcorChecker 
                                  INTERFACE ConcorContainer
                                  INTERFACE Boost::program_options
                                  INTERFACE Boost::iostreams
                                  INTERFACE Boost::regex
                                  INTERFACE z
                                  INTERFACE ${HTS_LIB}
)

target_include_directories(Concordance INTERFACE ${HTS_INCLUDE_DIR})


