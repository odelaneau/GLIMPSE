cmake_minimum_required(VERSION 3.14)
project(GLIMPSE2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set Policy
cmake_policy(SET CMP0167 NEW)

# Define Git information
execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE __COMMIT_ID__
    OUTPUT_STRIP_TRAILING_WHITESPACE
) 

execute_process(
    COMMAND git log -1 --format=%cd --date=short
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE __COMMIT_DATE__
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
    ${CMAKE_SOURCE_DIR}/versions/versions.h.in
    ${CMAKE_BINARY_DIR}/versions/versions.h
)

# Set HTS path
set(HTS_INCLUDE_DIR /usr/include)
set(HTS_LIB /usr/lib/x86_64-linux-gnu/libhts.so)

# Find libraries
find_package(Boost REQUIRED COMPONENTS serialization program_options iostreams regex)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Create Version packages
add_library(Version_tools INTERFACE)
target_include_directories(Version_tools INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/versions)

# Add common packages
add_subdirectory(common/src/utils)
add_subdirectory(common/src/containers)
add_subdirectory(common/src/io)
add_subdirectory(common/src/objects)

# Add modules
add_subdirectory(concordance)
add_subdirectory(chunk)
add_subdirectory(split_reference)
add_subdirectory(phase)
add_subdirectory(ligate)

# Switch Test modules or run main
option(TEST "DO GTEST" OFF)
if (TEST)
    # Workaround macro
    set(TOOLBOX_DEFINE "#ifdef _DECLARE_TOOLBOX_HERE")

    # Run tests
    add_subdirectory(tests)
else()
    # Workaround macro
    set(TOOLBOX_DEFINE "#ifdef BUILD_MAIN")

    # Create main exe
    add_executable(glimpse2 main.cpp)
    
    # Link to library
    target_link_libraries(glimpse2  PRIVATE Chunker
                                            Spliter
                                            Phase
                                            Ligater
                                            Concordance
                                            Threads::Threads
                                            Utils
    )

    # Link to compile options
    target_compile_options(glimpse2 PRIVATE -Wno-ignored-attributes) 
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/supports/otools.h.in
    ${PROJECT_BINARY_DIR}/supports/otools.h
)

## Create Doxygen (option)
option(DOXY "CREATE DOXYGEN DOCUMENT" OFF)
if(DOXY)
    find_package(Doxygen REQUIRED)

    ## Format css
    include(FetchContent)
    FetchContent_Declare(
        doxygen-awesome-css
        URL https://github.com/jothepro/doxygen-awesome-css/archive/refs/heads/main.zip
    )
    FetchContent_MakeAvailable(doxygen-awesome-css)
    FetchContent_GetProperties(doxygen-awesome-css SOURCE_DIR AWESOME_CSS_DIR)

    # Paths for input/output
    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # Replace variables in template
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    # Create custom target
    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating Doxygen documentation"
        DEPENDS 
            ${PROJECT_BINARY_DIR}/versions/versions.h 
            ${PROJECT_BINARY_DIR}/supports/otools.h
        VERBATIM
    )

endif()