cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(PSH10)

# GEANT4
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CMAKE_INSTALL_PREFIX}/include)

add_library(PSH10 "SHARED" ${sources} ${headers}) 

target_link_libraries(PSH10 ${Geant4_LIBRARIES})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/PSH10Config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/PSH10Config.cmake" @ONLY IMMEDIATE)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/PSH10Config.cmake" DESTINATION "lib/cmake/PSH10")
install(FILES ${headers} DESTINATION include)
install(TARGETS PSH10 DESTINATION lib)
