include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB LIGATE_SOURCE src/*.cpp src/ligater/*.cpp src/utils/*.cpp)

add_executable(GLIMPSE_ligate ${LIGATE_SOURCE})

target_link_libraries(GLIMPSE_ligate
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
  CURL::libcurl
  OpenSSL::Crypto
  ${HTSlib_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  )

install(TARGETS GLIMPSE_ligate DESTINATION bin COMPONENT runtime)
