# Create concordance library
add_library(Concordance INTERFACE)

# Link packages
find_package(Boost REQUIRED COMPONENTS program_options iostreams regex)
set(HTS_INCLUDE_DIR /usr/include)
set(HTS_LIB /usr/lib/x86_64-linux-gnu/libhts.so)

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


