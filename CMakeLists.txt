cmake_minimum_required(VERSION 3.6.0)

project(exampleDoPri5 LANGUAGES C CXX VERSION 0.4.2)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")

find_package(VecCore REQUIRED COMPONENTS Vc UMESIMD)
find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cc ${sources} ${headers})

target_link_libraries(${PROJECT_NAME} PRIVATE ${Geant4_LIBRARIES} VecCore::VecCore)