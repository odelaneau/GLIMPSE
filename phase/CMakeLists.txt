include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB PHASE_SOURCE
  src/*.cpp
  src/caller/*.cpp
  src/containers/*.cpp
  src/io/*.cpp
  src/models/*.cpp
  src/objects/*.cpp
  src/utils/*.cpp
  )

add_executable(GLIMPSE_phase ${PHASE_SOURCE})

target_link_libraries(GLIMPSE_phase
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
  CURL::libcurl
  OpenSSL::Crypto
  ${HTSlib_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  )

install(TARGETS GLIMPSE_phase DESTINATION bin COMPONENT runtime)
