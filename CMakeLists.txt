cmake_minimum_required(VERSION 3.10.0)

set (CMAKE_CXX_STANDARD 17)
project(img2spec VERSION 0.1.0 LANGUAGES C CXX)


file (GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME} PRIVATE ${SRC_FILES})

target_include_directories(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")