include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB STATS_SOURCE src/*.cpp src/stater/*.cpp src/utils/*.cpp)

add_executable(GLIMPSE_stats ${STATS_SOURCE})

target_link_libraries(GLIMPSE_stats
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
  CURL::libcurl
  OpenSSL::Crypto
  ${HTSlib_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  )

install(TARGETS GLIMPSE_stats DESTINATION bin COMPONENT runtime)
