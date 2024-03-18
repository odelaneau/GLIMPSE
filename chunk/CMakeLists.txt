include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB CHUNK_SOURCE src/*.cpp src/chunker/*.cpp)

add_executable(GLIMPSE_chunk ${CHUNK_SOURCE})

target_link_libraries(GLIMPSE_chunk
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
  CURL::libcurl
  OpenSSL::Crypto
  ${HTSlib_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  )

install(TARGETS GLIMPSE_chunk DESTINATION bin COMPONENT runtime)
