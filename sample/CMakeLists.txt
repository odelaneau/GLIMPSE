include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB SAMPLE_SOURCE src/*.cpp src/sampler/*.cpp src/utils/*.cpp)

add_executable(GLIMPSE_sample ${SAMPLE_SOURCE})

target_link_libraries(GLIMPSE_sample
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
  CURL::libcurl
  OpenSSL::Crypto
  ${HTSlib_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  )

install(TARGETS GLIMPSE_sample DESTINATION bin COMPONENT runtime)
