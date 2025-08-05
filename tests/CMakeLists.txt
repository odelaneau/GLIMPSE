include(FetchContent)

# Invoke Gtest
fetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)
enable_testing()

# Set up test_Concordance
add_executable(test_Concordance test_Concordance.cpp TestFile.cpp)
target_link_libraries(test_Concordance GTest::gtest_main Concordance)
add_test(NAME test_Concordance COMMAND test_Concordance)

# Set up test_Chunk
add_executable(test_Chunk test_Chunk.cpp TestFile.cpp)
target_link_libraries(test_Chunk GTest::gtest_main Chunker)
target_compile_definitions(test_Chunk PRIVATE MAP_FOLDER="${CMAKE_SOURCE_DIR}/maps/genetic_maps.b38")
add_test(NAME test_Chunk COMMAND test_Chunk)

# Set up test_Split_reference
add_executable(test_Split_reference test_Split_reference.cpp TestFile.cpp)
target_link_libraries(test_Split_reference GTest::gtest_main SPLIT_REFERENCE)
target_compile_definitions(test_Split_reference PRIVATE MAP_FOLDER="${CMAKE_SOURCE_DIR}/maps/genetic_maps.b38")
add_test(NAME test_Split_reference COMMAND test_Split_reference)

# Set up test_Phase
add_executable(test_Phase test_Phase.cpp TestFile.cpp)
target_link_libraries(test_Phase GTest::gtest_main Phase)
dd_test(NAME test_Phase COMMAND test_Phase)