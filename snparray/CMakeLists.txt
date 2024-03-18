include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB SNPARRAY_SOURCE src/*.cpp src/snparray/*.cpp src/utils/*.cpp)

add_executable(GLIMPSE_snparray ${SNPARRAY_SOURCE})

target_link_libraries(GLIMPSE_snparray
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
  CURL::libcurl
  OpenSSL::Crypto
  ${HTSlib_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  )

install(TARGETS GLIMPSE_snparray DESTINATION bin COMPONENT runtime)
