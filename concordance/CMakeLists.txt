include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB CONCORDANCE_SOURCE src/*.cpp src/checker/*.cpp src/containers/*.cpp)

add_executable(GLIMPSE_concordance ${CONCORDANCE_SOURCE})

target_link_libraries(GLIMPSE_concordance
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
  CURL::libcurl
  OpenSSL::Crypto
  ${HTSlib_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  )

install(TARGETS GLIMPSE_concordance DESTINATION bin COMPONENT runtime)
