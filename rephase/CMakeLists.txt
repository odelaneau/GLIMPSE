include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB REPHASE_SOURCE
  main/*.cpp
  src/containers/*.cpp
  src/io/*.cpp
  src/models/*.cpp
  src/objects/*.cpp
  src/rephaser/*.cpp
  src/utils/*.cpp
  )

add_executable(GLIMPSE_rephase ${REPHASE_SOURCE})

target_link_libraries(GLIMPSE_rephase
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
  CURL::libcurl
  OpenSSL::Crypto
  ${HTSlib_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  )

install(TARGETS GLIMPSE_rephase DESTINATION bin COMPONENT runtime)
